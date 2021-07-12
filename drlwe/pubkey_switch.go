package drlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// PublicKeySwitchingProtocol is an interface describing the local steps of a generic RLWE PCKS protocol.
type PublicKeySwitchingProtocol interface {
	AllocateShares(level int) *PCKSShare
	GenShare(skInput *rlwe.SecretKey, pkOutput *rlwe.PublicKey, ct rlwe.Ciphertext, shareOut *PCKSShare)
	AggregateShares(share1, share2, shareOut *PCKSShare)
	KeySwitch(combined *PCKSShare, ct rlwe.Ciphertext, ctOut rlwe.Ciphertext)
}

// PCKSShare represents a party's share in the PCKS protocol.
type PCKSShare struct {
	Value [2]*ring.Poly
}

// PCKSProtocol is the structure storing the parameters for the collective public key-switching.
type PCKSProtocol struct {
	ringQ         *ring.Ring
	ringP         *ring.Ring
	ringQP        *ring.Ring
	sigmaSmudging float64

	tmpQ       *ring.Poly
	tmpP       *ring.Poly
	share0tmpQ *ring.Poly
	share1tmpQ *ring.Poly
	share0tmpP *ring.Poly
	share1tmpP *ring.Poly

	baseconverter             *ring.FastBasisExtender
	gaussianSampler           *ring.GaussianSampler
	ternarySamplerMontgomeryQ *ring.TernarySampler
}

// NewPCKSProtocol creates a new PCKSProtocol object and will be used to re-encrypt a ciphertext ctx encrypted under a secret-shared key among j parties under a new
// collective public-key.
func NewPCKSProtocol(params rlwe.Parameters, sigmaSmudging float64) *PCKSProtocol {
	pcks := new(PCKSProtocol)
	pcks.ringQ = params.RingQ()
	pcks.ringP = params.RingP()
	pcks.sigmaSmudging = sigmaSmudging

	pcks.tmpQ = pcks.ringQ.NewPoly()
	pcks.tmpP = pcks.ringP.NewPoly()
	pcks.share0tmpQ = pcks.ringQ.NewPoly()
	pcks.share1tmpQ = pcks.ringQ.NewPoly()
	pcks.share0tmpP = pcks.ringP.NewPoly()
	pcks.share1tmpP = pcks.ringP.NewPoly()

	pcks.baseconverter = ring.NewFastBasisExtender(pcks.ringQ, pcks.ringP)
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	pcks.gaussianSampler = ring.NewGaussianSampler(prng, pcks.ringQ, sigmaSmudging, int(6*sigmaSmudging))
	pcks.ternarySamplerMontgomeryQ = ring.NewTernarySampler(prng, pcks.ringQ, 0.5, false)

	return pcks
}

// AllocateShares allocates the shares of the PCKS protocol
func (pcks *PCKSProtocol) AllocateShares(level int) (s *PCKSShare) {
	s = new(PCKSShare)
	s.Value[0] = pcks.ringQ.NewPolyLvl(level)
	s.Value[1] = pcks.ringQ.NewPolyLvl(level)
	return
}

// GenShare is the first part of the unique round of the PCKSProtocol protocol. Each party computes the following :
//
// [s_i * ctx[0] + (u_i * pk[0] + e_0i)/P, (u_i * pk[1] + e_1i)/P]
//
// and broadcasts the result to the other j-1 parties.
func (pcks *PCKSProtocol) GenShare(sk *rlwe.SecretKey, pk *rlwe.PublicKey, ct rlwe.Ciphertext, shareOut *PCKSShare) {

	el := ct.RLWEElement()

	ringQ := pcks.ringQ
	ringP := pcks.ringP

	level := el.Level()

	pk0Q := new(ring.Poly)
	pk0P := new(ring.Poly)
	pk1Q := new(ring.Poly)
	pk1P := new(ring.Poly)

	// Splits pk[0] and pk[1] with their respective modulus Qlvl and P
	pk0Q.Coeffs = pk.Value[0].Coeffs[:level+1]
	pk1Q.Coeffs = pk.Value[1].Coeffs[:level+1]
	pk0P.Coeffs = pk.Value[0].Coeffs[len(ringQ.Modulus):]
	pk1P.Coeffs = pk.Value[1].Coeffs[len(ringQ.Modulus):]

	// samples MForm(u_i) in Q and P separately
	pcks.ternarySamplerMontgomeryQ.ReadLvl(level, pcks.tmpQ)
	extendBasisSmallNormAndCenter(ringQ.Modulus[0], ringP.Modulus, pcks.tmpQ.Coeffs[0], pcks.tmpP.Coeffs)
	ringQ.MFormLvl(level, pcks.tmpQ, pcks.tmpQ)
	ringP.MForm(pcks.tmpP, pcks.tmpP)
	ringQ.NTTLvl(level, pcks.tmpQ, pcks.tmpQ)
	ringP.NTT(pcks.tmpP, pcks.tmpP)

	// h_0 = NTT(u_i * pk_0)
	ringQ.MulCoeffsMontgomeryLvl(level, pcks.tmpQ, pk0Q, pcks.share0tmpQ)
	ringP.MulCoeffsMontgomery(pcks.tmpP, pk0P, pcks.share0tmpP)
	ringQ.InvNTTLvl(level, pcks.share0tmpQ, pcks.share0tmpQ)
	ringP.InvNTT(pcks.share0tmpP, pcks.share0tmpP)

	// h_1 = NTT(u_i * pk_1)
	ringQ.MulCoeffsMontgomeryLvl(level, pcks.tmpQ, pk1Q, pcks.share1tmpQ)
	ringP.MulCoeffsMontgomery(pcks.tmpP, pk1P, pcks.share1tmpP)
	ringQ.InvNTTLvl(level, pcks.share1tmpQ, pcks.share1tmpQ)
	ringP.InvNTT(pcks.share1tmpP, pcks.share1tmpP)

	// h_0 = u_i * pk_0 + e0
	pcks.gaussianSampler.ReadLvl(level, pcks.tmpQ)
	extendBasisSmallNormAndCenter(ringQ.Modulus[0], ringP.Modulus, pcks.tmpQ.Coeffs[0], pcks.tmpP.Coeffs)
	ringQ.AddLvl(level, pcks.share0tmpQ, pcks.tmpQ, pcks.share0tmpQ)
	ringP.Add(pcks.share0tmpP, pcks.tmpP, pcks.share0tmpP)

	// h_1 = u_i * pk_1 + e1
	pcks.gaussianSampler.ReadLvl(level, pcks.tmpQ)
	extendBasisSmallNormAndCenter(ringQ.Modulus[0], ringP.Modulus, pcks.tmpQ.Coeffs[0], pcks.tmpP.Coeffs)
	ringQ.AddLvl(level, pcks.share1tmpQ, pcks.tmpQ, pcks.share1tmpQ)
	ringP.Add(pcks.share1tmpP, pcks.tmpP, pcks.share1tmpP)

	// h_0 = (u_i * pk_0 + e0)/P
	pcks.baseconverter.ModDownSplitPQ(level, pcks.share0tmpQ, pcks.share0tmpP, shareOut.Value[0])

	// h_1 = (u_i * pk_1 + e1)/P
	pcks.baseconverter.ModDownSplitPQ(level, pcks.share1tmpQ, pcks.share1tmpP, shareOut.Value[1])

	// h_0 = s_i*c_1 + (u_i * pk_0 + e0)/P
	if el.Value[0].IsNTT {
		ringQ.NTTLvl(level, shareOut.Value[0], shareOut.Value[0])
		ringQ.NTTLvl(level, shareOut.Value[1], shareOut.Value[1])
		ringQ.MulCoeffsMontgomeryAndAddLvl(level, el.Value[1], sk.Value, shareOut.Value[0])
	} else {
		// tmp = s_i*c_1
		ringQ.NTTLazyLvl(level, el.Value[1], pcks.tmpQ)
		ringQ.MulCoeffsMontgomeryConstantLvl(level, pcks.tmpQ, sk.Value, pcks.tmpQ)
		ringQ.InvNTTLvl(level, pcks.tmpQ, pcks.tmpQ)

		// h_0 = s_i*c_1 + (u_i * pk_0 + e0)/P
		ringQ.AddLvl(level, shareOut.Value[0], pcks.tmpQ, shareOut.Value[0])
	}
}

// AggregateShares is the second part of the first and unique round of the PCKSProtocol protocol. Each party uppon receiving the j-1 elements from the
// other parties computes :
//
// [ctx[0] + sum(s_i * ctx[0] + u_i * pk[0] + e_0i), sum(u_i * pk[1] + e_1i)]
func (pcks *PCKSProtocol) AggregateShares(share1, share2, shareOut *PCKSShare) {
	level1, level2 := len(share1.Value[0].Coeffs)-1, len(share2.Value[0].Coeffs)-1
	if level1 != level2 {
		panic("cannot aggreate two shares at different levels.")
	}
	pcks.ringQ.AddLvl(level1, share1.Value[0], share2.Value[0], shareOut.Value[0])
	pcks.ringQ.AddLvl(level1, share1.Value[1], share2.Value[1], shareOut.Value[1])
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (pcks *PCKSProtocol) KeySwitch(combined *PCKSShare, ct, ctOut rlwe.Ciphertext) {
	el, elOut := ct.RLWEElement(), ctOut.RLWEElement()
	pcks.ringQ.AddLvl(el.Level(), el.Value[0], combined.Value[0], elOut.Value[0])
	ring.CopyValuesLvl(el.Level(), combined.Value[1], elOut.Value[1])
}

// MarshalBinary encodes a PCKS share on a slice of bytes.
func (share *PCKSShare) MarshalBinary() ([]byte, error) {
	lenR1 := share.Value[0].GetDataLen(true)
	lenR2 := share.Value[1].GetDataLen(true)

	data := make([]byte, lenR1+lenR2)
	_, err := share.Value[0].WriteTo(data[0:lenR1])
	if err != nil {
		return []byte{}, err
	}

	_, err = share.Value[1].WriteTo(data[lenR1 : lenR1+lenR2])
	if err != nil {
		return []byte{}, err
	}

	return data, nil
}

// UnmarshalBinary decodes marshaled PCKS share on the target PCKS share.
func (share *PCKSShare) UnmarshalBinary(data []byte) error {

	if share.Value[0] == nil {
		share.Value[0] = new(ring.Poly)
	}

	if share.Value[1] == nil {
		share.Value[1] = new(ring.Poly)
	}

	err := share.Value[0].UnmarshalBinary(data[0 : len(data)/2])
	if err != nil {
		return err
	}

	err = share.Value[1].UnmarshalBinary(data[len(data)/2:])
	if err != nil {
		return err
	}

	return nil
}
