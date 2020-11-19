package dbfv

import (
	"encoding/binary"
	"errors"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// RTGProtocol is the structure storing the parameters for the collective rotation-keys generation.
type RTGProtocol struct {
	context *dbfvContext

	tmpSwitchKey    [][2]*ring.Poly
	tmpPoly         [2]*ring.Poly
	gaussianSampler *ring.GaussianSampler
}

// RTGShare is the structure storing the shares of the RTG protocol
type RTGShare struct {
	Type  bfv.Rotation
	K     uint64
	Value []*ring.Poly
}

// MarshalBinary encode the target element on a slice of byte.
func (share *RTGShare) MarshalBinary() ([]byte, error) {
	lenRing := share.Value[0].GetDataLen(true)
	data := make([]byte, 3*8+lenRing*uint64(len(share.Value)))
	binary.BigEndian.PutUint64(data[0:8], share.K)
	binary.BigEndian.PutUint64(data[8:16], uint64(share.Type))
	binary.BigEndian.PutUint64(data[16:24], lenRing)
	ptr := uint64(24)
	for _, val := range share.Value {
		cnt, err := val.WriteTo(data[ptr : ptr+lenRing])
		if err != nil {
			return []byte{}, err
		}
		ptr += cnt
	}

	return data, nil
}

// UnmarshalBinary decodes a slice of bytes on the target element.
func (share *RTGShare) UnmarshalBinary(data []byte) error {
	if len(data) <= 24 {
		return errors.New("Unsufficient data length")
	}
	share.K = binary.BigEndian.Uint64(data[0:8])
	share.Type = bfv.Rotation(binary.BigEndian.Uint64(data[8:16]))
	lenRing := binary.BigEndian.Uint64(data[16:24])
	valLength := uint64(len(data)-3*8) / lenRing

	share.Value = make([]*ring.Poly, valLength)
	ptr := uint64(24)
	for i := range share.Value {
		share.Value[i] = new(ring.Poly)
		err := share.Value[i].UnmarshalBinary(data[ptr : ptr+lenRing])
		if err != nil {
			return err
		}
		ptr += lenRing

	}

	return nil
}

// AllocateShare allocates the shares of the RTG protocol.
func (rtg *RTGProtocol) AllocateShare() (rtgShare RTGShare) {
	rtgShare.Value = make([]*ring.Poly, rtg.context.params.Beta())
	for i := uint64(0); i < rtg.context.params.Beta(); i++ {
		rtgShare.Value[i] = rtg.context.ringQP.NewPoly()
	}
	return
}

// NewRotKGProtocol creates a new rotkg object and will be used to generate collective rotation-keys from a shared secret-key among j parties.
func NewRotKGProtocol(params *bfv.Parameters) (rtg *RTGProtocol) {

	context := newDbfvContext(params)

	rtg = new(RTGProtocol)
	rtg.context = context

	rtg.tmpSwitchKey = make([][2]*ring.Poly, rtg.context.params.Beta())
	for i := range rtg.tmpSwitchKey {
		rtg.tmpSwitchKey[i][0] = context.ringQP.NewPoly()
		rtg.tmpSwitchKey[i][1] = context.ringQP.NewPoly()
	}

	rtg.tmpPoly = [2]*ring.Poly{context.ringQP.NewPoly(), context.ringQP.NewPoly()}

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	rtg.gaussianSampler = ring.NewGaussianSampler(prng, context.ringQP, params.Sigma(), uint64(6*params.Sigma()))

	return rtg
}

// GenShare is the first and unique round of the rotkg protocol. Each party, using its secret share of the collective secret-key
// and a collective random polynomial, a public share of the rotation-key by computing :
//
// [a*s_i + (pi(s_i) - s_i) + e]
//
// and broadcasts it to the other j-1 parties. The protocol must be repeated for each desired rotation.
func (rtg *RTGProtocol) GenShare(rotType bfv.Rotation, k uint64, sk *ring.Poly, crp []*ring.Poly, shareOut *RTGShare) {

	ringQP := rtg.context.ringQP

	shareOut.Type = rotType
	shareOut.K = k
	switch rotType {
	case bfv.RotationRight:
		rtg.genShare(sk, bfv.GaloisGen, k&(ringQP.N/2-1), crp, shareOut.Value)
		return
	case bfv.RotationLeft:
		rtg.genShare(sk, bfv.GaloisGen, (ringQP.N/2-k)&(ringQP.N/2-1), crp, shareOut.Value)
		return
	case bfv.RotationRow:
		rtg.genShare(sk, ringQP.NthRoot-1, 1, crp, shareOut.Value)
		return
	}
}

func (rtg *RTGProtocol) genShare(sk *ring.Poly, gen, k uint64, crp []*ring.Poly, evakey []*ring.Poly) {

	ringQP := rtg.context.ringQP

	ring.PermuteNTT(sk, gen, k, ringQP.N, ringQP.NthRoot, rtg.tmpPoly[1])

	ringQP.MulScalarBigint(sk, rtg.context.ringP.ModulusBigint, rtg.tmpPoly[0])

	var index uint64

	for i := uint64(0); i < rtg.context.params.Beta(); i++ {

		// e
		rtg.gaussianSampler.Read(evakey[i])
		ringQP.MForm(evakey[i], evakey[i])
		ringQP.NTT(evakey[i], evakey[i])

		// a is the CRP

		// e + sk_in * (qiBarre*qiStar) * 2^w
		// (qiBarre*qiStar)%qi = 1, else 0
		for j := uint64(0); j < rtg.context.params.Alpha(); j++ {

			index = i*rtg.context.params.Alpha() + j

			qi := ringQP.Modulus[index]
			tmp0 := rtg.tmpPoly[0].Coeffs[index]
			tmp1 := evakey[i].Coeffs[index]

			for w := uint64(0); w < ringQP.N; w++ {
				tmp1[w] = ring.CRed(tmp1[w]+tmp0[w], qi)
			}

			// Handles the case where nb pj does not divides nb qi

			if index >= rtg.context.params.QiCount() {
				break
			}
		}

		// sk_in * (qiBarre*qiStar) * 2^w - a*sk + e
		ringQP.MulCoeffsMontgomeryAndSub(crp[i], rtg.tmpPoly[1], evakey[i])
	}

	rtg.tmpPoly[0].Zero()
	rtg.tmpPoly[1].Zero()

	return
}

// Aggregate is the second part of the unique round of the rotkg protocol. Uppon receiving the j-1 public shares,
// each party computes  :
//
// [sum(a*a_j + (pi(a_j) - a_j) + e_j), a]
func (rtg *RTGProtocol) Aggregate(share1, share2, shareOut RTGShare) {
	ringQP := rtg.context.ringQP

	if share1.Type != share2.Type || share1.K != share2.K {
		panic("cannot aggregate shares of different types")
	}

	shareOut.Type = share1.Type
	shareOut.K = share1.K
	for i := uint64(0); i < rtg.context.params.Beta(); i++ {
		ringQP.Add(share1.Value[i], share2.Value[i], shareOut.Value[i])
	}
}

// Finalize populates the input RotationKeys struture with the Switching key computed from the protocol.
func (rtg *RTGProtocol) Finalize(share RTGShare, crp []*ring.Poly, rotKey *bfv.RotationKeys) {

	ringQP := rtg.context.ringQP

	k := share.K & ((rtg.context.n >> 1) - 1)

	for i := uint64(0); i < rtg.context.params.Beta(); i++ {
		ringQP.Copy(share.Value[i], rtg.tmpSwitchKey[i][0])
		ringQP.Copy(crp[i], rtg.tmpSwitchKey[i][1])
	}

	rotKey.SetRotKey(share.Type, k, rtg.tmpSwitchKey)
}
