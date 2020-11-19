package drckks

import (
	"github.com/ldsec/lattigo/v2/rckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// RTGProtocol is the structure storing the parameters for the collective rotation-keys generation.
type RTGProtocol struct {
	drckksContext *drckksContext

	galElRotRow uint64
	galElRotCol map[rckks.Rotation][]uint64

	tmpSwitchKey    [][2]*ring.Poly
	tmpPoly         [2]*ring.Poly
	gaussianSampler *ring.GaussianSampler
}

// RTGShare is a struct storing the share of the RTG protocol.
type RTGShare struct {
	Type  rckks.Rotation
	K     uint64
	Value []*ring.Poly
}

// AllocateShare allocates the share the the RTG protocol.
func (rtg *RTGProtocol) AllocateShare() (rtgShare RTGShare) {
	rtgShare.Value = make([]*ring.Poly, rtg.drckksContext.beta)
	for i := uint64(0); i < rtg.drckksContext.beta; i++ {
		rtgShare.Value[i] = rtg.drckksContext.ringQP.NewPoly()
	}
	return
}

// NewRotKGProtocol creates a new rotkg object and will be used to generate collective rotation-keys from a shared secret-key among j parties.
func NewRotKGProtocol(params *rckks.Parameters) (rtg *RTGProtocol) {

	rtg = new(RTGProtocol)

	drckksContext := newDrckksContext(params)

	rtg.drckksContext = drckksContext

	rtg.tmpSwitchKey = make([][2]*ring.Poly, rtg.drckksContext.beta)
	for i := range rtg.tmpSwitchKey {
		rtg.tmpSwitchKey[i][0] = drckksContext.ringQP.NewPoly()
		rtg.tmpSwitchKey[i][1] = drckksContext.ringQP.NewPoly()
	}

	rtg.tmpPoly = [2]*ring.Poly{drckksContext.ringQP.NewPoly(), drckksContext.ringQP.NewPoly()}

	N := drckksContext.n

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	rtg.gaussianSampler = ring.NewGaussianSampler(prng, drckksContext.ringQP, params.Sigma(), uint64(6*params.Sigma()))

	rtg.galElRotRow = (N << 1) - 1

	return rtg
}

// GenShare is the first and unique round of the rotkg protocol. Each party, using its secret share of the collective secret-key
// and a collective random polynomial, a public share of the rotation-key by computing :
//
// [a*s_i + pi(s_i) + e]
//
// and broadcasts it to the other j-1 parties. The protocol must be repeated for each desired rotation.
func (rtg *RTGProtocol) GenShare(rotType rckks.Rotation, k uint64, sk *ring.Poly, crp []*ring.Poly, shareOut *RTGShare) {
	shareOut.Type = rotType
	shareOut.K = k
	switch rotType {
	case rckks.RotationRight:
		rtg.genShare(sk, k, crp, shareOut.Value)
		return
	case rckks.RotationLeft:
		rtg.genShare(sk, rtg.drckksContext.ringQP.N-k, crp, shareOut.Value)
		return
	}
}

// genswitchkey is a generic method to generate the public-share of the collective rotation-key.
func (rtg *RTGProtocol) genShare(sk *ring.Poly, k uint64, crp []*ring.Poly, evakey []*ring.Poly) {

	ringQP := rtg.drckksContext.ringQP

	ring.PermuteNTT(sk, rckks.GaloisGen, k, ringQP.N, ringQP.NthRoot, rtg.tmpPoly[1])

	ringQP.MulScalarBigint(sk, rtg.drckksContext.ringP.ModulusBigint, rtg.tmpPoly[0])

	var index uint64

	for i := uint64(0); i < rtg.drckksContext.beta; i++ {

		// e
		rtg.gaussianSampler.Read(evakey[i])
		ringQP.MForm(evakey[i], evakey[i])
		rckks.NTTRCKKS(ringQP, evakey[i], evakey[i])

		// a is the CRP

		// e + sk_in * (qiBarre*qiStar) * 2^w
		// (qiBarre*qiStar)%qi = 1, else 0
		for j := uint64(0); j < rtg.drckksContext.alpha; j++ {

			index = i*rtg.drckksContext.alpha + j

			qi := ringQP.Modulus[index]
			tmp0 := rtg.tmpPoly[0].Coeffs[index]
			tmp1 := evakey[i].Coeffs[index]

			for w := uint64(0); w < ringQP.N; w++ {
				tmp1[w] = ring.CRed(tmp1[w]+tmp0[w], qi)
			}

			// Handles the case where nb pj does not divides nb qi
			if index >= rtg.drckksContext.params.QiCount() {
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
// [sum(a*a_j + pi(s_j) + e_j), a]
func (rtg *RTGProtocol) Aggregate(share1, share2, shareOut RTGShare) {
	ringQP := rtg.drckksContext.ringQP

	if share1.Type != share2.Type || share1.K != share2.K {
		panic("cannot aggregate shares of different types")
	}

	shareOut.Type = share1.Type
	shareOut.K = share1.K
	for i := uint64(0); i < rtg.drckksContext.beta; i++ {
		ringQP.Add(share1.Value[i], share2.Value[i], shareOut.Value[i])
	}
}

// Finalize finalizes the RTG protocol and populates the input RotationKey with the computed collective SwitchingKey.
func (rtg *RTGProtocol) Finalize(params *rckks.Parameters, share RTGShare, crp []*ring.Poly, rotKey *rckks.RotationKeys) {

	ringQP := rtg.drckksContext.ringQP

	k := share.K & ((rtg.drckksContext.n >> 1) - 1)

	for i := uint64(0); i < rtg.drckksContext.beta; i++ {
		ringQP.Copy(share.Value[i], rtg.tmpSwitchKey[i][0])
		ringQP.Copy(crp[i], rtg.tmpSwitchKey[i][1])
	}

	rotKey.SetRotKey(params, rtg.tmpSwitchKey, share.Type, k)
}
