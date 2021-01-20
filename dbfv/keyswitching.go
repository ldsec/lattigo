package dbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// CKSProtocol is a structure storing the parameters for the collective key-switching protocol.
type CKSProtocol struct {
	context *dbfvContext

	sigmaSmudging float64

	tmpNtt   *ring.Poly
	tmpDelta *ring.Poly
	hP       *ring.Poly

	baseconverter   *ring.FastBasisExtender
	gaussianSampler *ring.GaussianSampler
}

// CKSShare is a type for the CKS protocol shares.
type CKSShare struct {
	*ring.Poly
}

// UnmarshalBinary decodes a previouls marshaled share on the target share.
func (share *CKSShare) UnmarshalBinary(data []byte) error {
	share.Poly = new(ring.Poly)
	err := share.Poly.UnmarshalBinary(data)
	return err

}

// NewCKSProtocol creates a new CKSProtocol that will be used to operate a collective key-switching on a ciphertext encrypted under a collective public-key, whose
// secret-shares are distributed among j parties, re-encrypting the ciphertext under another public-key, whose secret-shares are also known to the
// parties.
func NewCKSProtocol(params *bfv.Parameters, sigmaSmudging float64) *CKSProtocol {

	context := newDbfvContext(params)

	cks := new(CKSProtocol)

	cks.context = context

	cks.sigmaSmudging = sigmaSmudging

	cks.tmpNtt = cks.context.ringQP.NewPoly()
	cks.tmpDelta = cks.context.ringQ.NewPoly()
	cks.hP = cks.context.ringP.NewPoly()

	cks.baseconverter = ring.NewFastBasisExtender(cks.context.ringQ, cks.context.ringP)
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	cks.gaussianSampler = ring.NewGaussianSampler(prng)

	return cks
}

// AllocateShare allocates the shares of the CKSProtocol
func (cks *CKSProtocol) AllocateShare() CKSShare {

	return CKSShare{cks.context.ringQ.NewPoly()}

}

// GenShare is the first and unique round of the CKSProtocol protocol. Each party holding a ciphertext ctx encrypted under a collective publick-key musth
// compute the following :
//
// [(skInput_i - skOutput_i) * ctx[0] + e_i]
//
// Each party then broadcast the result of this computation to the other j-1 parties.
func (cks *CKSProtocol) GenShare(skInput, skOutput *ring.Poly, ct *bfv.Ciphertext, shareOut CKSShare) {

	cks.context.ringQ.Sub(skInput, skOutput, cks.tmpDelta)

	cks.genShareDelta(cks.tmpDelta, ct, shareOut)
}

func (cks *CKSProtocol) genShareDelta(skDelta *ring.Poly, ct *bfv.Ciphertext, shareOut CKSShare) {

	level := uint64(len(ct.Value()[1].Coeffs) - 1)

	ringQ := cks.context.ringQ
	ringQP := cks.context.ringQP

	ringQ.NTTLazy(ct.Value()[1], cks.tmpNtt)
	ringQ.MulCoeffsMontgomeryConstant(cks.tmpNtt, skDelta, shareOut.Poly)
	ringQ.MulScalarBigint(shareOut.Poly, cks.context.ringP.ModulusBigint, shareOut.Poly)

	ringQ.InvNTTLazy(shareOut.Poly, shareOut.Poly)

	cks.gaussianSampler.ReadLvl(uint64(len(ringQP.Modulus)-1), cks.tmpNtt, ringQP, cks.sigmaSmudging, uint64(6*cks.sigmaSmudging))
	ringQ.AddNoMod(shareOut.Poly, cks.tmpNtt, shareOut.Poly)

	for x, i := 0, uint64(len(ringQ.Modulus)); i < uint64(len(cks.context.ringQP.Modulus)); x, i = x+1, i+1 {
		tmphP := cks.hP.Coeffs[x]
		tmpNTT := cks.tmpNtt.Coeffs[i]
		for j := uint64(0); j < ringQ.N; j++ {
			tmphP[j] += tmpNTT[j]
		}
	}

	cks.baseconverter.ModDownSplitPQ(level, shareOut.Poly, cks.hP, shareOut.Poly)

	cks.tmpNtt.Zero()
	cks.hP.Zero()
}

// AggregateShares is the second part of the unique round of the CKSProtocol protocol. Upon receiving the j-1 elements each party computes :
//
// [ctx[0] + sum((skInput_i - skOutput_i) * ctx[0] + e_i), ctx[1]]
func (cks *CKSProtocol) AggregateShares(share1, share2, shareOut CKSShare) {
	cks.context.ringQ.Add(share1.Poly, share2.Poly, shareOut.Poly)
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (cks *CKSProtocol) KeySwitch(combined CKSShare, ct *bfv.Ciphertext, ctOut *bfv.Ciphertext) {
	cks.context.ringQ.Add(ct.Value()[0], combined.Poly, ctOut.Value()[0])
	cks.context.ringQ.Copy(ct.Value()[1], ctOut.Value()[1])
}
