package drlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// KeySwitchingProtocol is an interface describing the local steps of a generic RLWE CKS protocol
type KeySwitchingProtocol interface {
	AllocateShare(level int) *CKSShare
	GenShare(skInput, skOutput *rlwe.SecretKey, ct *rlwe.Ciphertext, shareOut *CKSShare)
	AggregateShares(share1, share2, shareOut *CKSShare)
	KeySwitch(combined *CKSShare, ct, ctOut *rlwe.Ciphertext)
}

// CKSProtocol is the structure storing the parameters and and precomputations for the collective key-switching protocol.
type CKSProtocol struct {
	ringQ           *ring.Ring
	ringP           *ring.Ring
	gaussianSampler *ring.GaussianSampler
	baseconverter   *ring.FastBasisExtender

	tmpP     *ring.Poly
	tmpQ     *ring.Poly
	tmpDelta *ring.Poly
}

// CKSShare is a type for the CKS protocol shares.
type CKSShare struct {
	Value *ring.Poly
}

// MarshalBinary encodes a CKS share on a slice of bytes.
func (ckss *CKSShare) MarshalBinary() ([]byte, error) {
	return ckss.Value.MarshalBinary()
}

// UnmarshalBinary decodes marshaled CKS share on the target CKS share.
func (ckss *CKSShare) UnmarshalBinary(data []byte) error {
	if ckss.Value == nil {
		ckss.Value = new(ring.Poly)
	}
	return ckss.Value.UnmarshalBinary(data)
}

// NewCKSProtocol creates a new CKSProtocol that will be used to operate a collective key-switching on a ciphertext encrypted under a collective public-key, whose
// secret-shares are distributed among j parties, re-encrypting the ciphertext under another public-key, whose secret-shares are also known to the
// parties.
func NewCKSProtocol(params rlwe.Parameters, sigmaSmudging float64) *CKSProtocol {
	cks := new(CKSProtocol)
	cks.ringQ = params.RingQ()
	cks.ringP = params.RingP() // TODO this assumes that P larger than 1

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	cks.gaussianSampler = ring.NewGaussianSampler(prng, cks.ringQ, sigmaSmudging, int(6*sigmaSmudging))
	cks.baseconverter = ring.NewFastBasisExtender(cks.ringQ, cks.ringP)

	cks.tmpQ = cks.ringQ.NewPoly()
	cks.tmpP = cks.ringP.NewPoly()
	cks.tmpDelta = cks.ringQ.NewPoly()

	return cks
}

// AllocateShare allocates the shares of the CKSProtocol
func (cks *CKSProtocol) AllocateShare(level int) *CKSShare {
	return &CKSShare{cks.ringQ.NewPolyLvl(level)}
}

// GenShare computes a party's share in the CKS protocol.
// ct.Value[0] can be nil, computations are only done using ct.Value[1]
// NTT flag for ct.Value[1] is expected to be set correctly
func (cks *CKSProtocol) GenShare(skInput, skOutput *rlwe.SecretKey, ct *rlwe.Ciphertext, shareOut *CKSShare) {

	el := ct.RLWEElement()

	ringQ := cks.ringQ
	ringP := cks.ringP

	level := utils.MinInt(len(ringQ.Modulus)-1, el.Value[1].Level())

	ringQ.SubLvl(level, skInput.Value, skOutput.Value, cks.tmpDelta)

	ct1 := el.Value[1]
	if !el.Value[1].IsNTT {
		ringQ.NTTLazyLvl(level, el.Value[1], cks.tmpQ)
		ct1 = cks.tmpQ
	}

	// a * (skIn - skOut) mod Q
	ringQ.MulCoeffsMontgomeryConstantLvl(level, ct1, cks.tmpDelta, shareOut.Value)

	// P * a * (skIn - skOut) mod QP (mod P = 0)
	ringQ.MulScalarBigintLvl(level, shareOut.Value, cks.ringP.ModulusBigint, shareOut.Value)

	if !el.Value[1].IsNTT {
		// InvNTT(P * a * (skIn - skOut)) mod QP (mod P = 0)
		ringQ.InvNTTLazyLvl(level, shareOut.Value, shareOut.Value)

		// Samples e in Q
		cks.gaussianSampler.ReadLvl(level, cks.tmpQ)

		// Extend e to P (assumed to have norm < qi)
		extendBasisSmallNormAndCenter(ringQ.Modulus[0], ringP.Modulus, cks.tmpQ.Coeffs[0], cks.tmpP.Coeffs)

		// InvNTT(P * a * (skIn - skOut) + e) mod QP (mod P = e)
		ringQ.AddNoModLvl(level, shareOut.Value, cks.tmpQ, shareOut.Value)

		// InvNTT(P * a * (skIn - skOut) + e) * (1/P) mod QP (mod P = e)
		cks.baseconverter.ModDownSplitPQ(level, shareOut.Value, cks.tmpP, shareOut.Value)

	} else {
		// Sample e in Q
		cks.gaussianSampler.ReadLvl(level, cks.tmpQ)

		// Extend e to P (assumed to have norm < qi)
		extendBasisSmallNormAndCenter(ringQ.Modulus[0], ringP.Modulus, cks.tmpQ.Coeffs[0], cks.tmpP.Coeffs)

		// Takes the error to the NTT domain
		ringQ.NTTLvl(level, cks.tmpQ, cks.tmpQ)
		ringP.NTT(cks.tmpP, cks.tmpP)

		// P * a * (skIn - skOut) + e mod Q (mod P = 0, so P = e)
		ringQ.AddLvl(level, shareOut.Value, cks.tmpQ, shareOut.Value)

		// (P * a * (skIn - skOut) + e) * (1/P) mod QP (mod P = e)
		cks.baseconverter.ModDownSplitNTTPQ(level, shareOut.Value, cks.tmpP, shareOut.Value)
	}

	shareOut.Value.Coeffs = shareOut.Value.Coeffs[:level+1]
}

// AggregateShares is the second part of the unique round of the CKSProtocol protocol. Upon receiving the j-1 elements each party computes :
//
// [ctx[0] + sum((skInput_i - skOutput_i) * ctx[0] + e_i), ctx[1]]
func (cks *CKSProtocol) AggregateShares(share1, share2, shareOut *CKSShare) {
	cks.ringQ.AddLvl(share1.Value.Level(), share1.Value, share2.Value, shareOut.Value)
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (cks *CKSProtocol) KeySwitch(combined *CKSShare, ct, ctOut *rlwe.Ciphertext) {
	el, elOut := ct.RLWEElement(), ctOut.RLWEElement()
	cks.ringQ.AddLvl(el.Level(), el.Value[0], combined.Value, elOut.Value[0])
	ring.CopyValuesLvl(el.Level(), el.Value[1], elOut.Value[1])
}
