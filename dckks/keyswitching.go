package dckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// CKSProtocol is a structure storing the parameters for the collective key-switching protocol.
type CKSProtocol struct {
	dckksContext *dckksContext

	sigmaSmudging float64

	tmpQ     *ring.Poly
	tmpP     *ring.Poly
	tmpDelta *ring.Poly
	hP       *ring.Poly

	baseconverter   *ring.FastBasisExtender
	gaussianSampler *ring.GaussianSampler
}

// CKSShare is a struct holding a share of the CKS protocol.
type CKSShare *ring.Poly

// NewCKSProtocol creates a new CKSProtocol that will be used to operate a collective key-switching on a ciphertext encrypted under a collective public-key, whose
// secret-shares are distributed among j parties, re-encrypting the ciphertext under another public-key, whose secret-shares are also known to the
// parties.
func NewCKSProtocol(params *ckks.Parameters, sigmaSmudging float64) (cks *CKSProtocol) {

	cks = new(CKSProtocol)

	dckksContext := newDckksContext(params)

	cks.dckksContext = dckksContext

	cks.tmpQ = dckksContext.ringQ.NewPoly()
	cks.tmpP = dckksContext.ringP.NewPoly()
	cks.tmpDelta = dckksContext.ringQ.NewPoly()
	cks.hP = dckksContext.ringP.NewPoly()

	cks.sigmaSmudging = sigmaSmudging

	cks.baseconverter = ring.NewFastBasisExtender(dckksContext.ringQ, dckksContext.ringP)
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	cks.gaussianSampler = ring.NewGaussianSampler(prng)

	return cks
}

// AllocateShare allocates the share of the CKS protocol.
func (cks *CKSProtocol) AllocateShare() CKSShare {
	return cks.dckksContext.ringQ.NewPoly()
}

// GenShare is the first and unique round of the CKSProtocol protocol. Each party holding a ciphertext ctx encrypted under a collective publick-key must
// compute the following :
//
// [(skInput_i - skOutput_i) * ctx[0] + e_i]
//
// Each party then broadcasts the result of this computation to the other j-1 parties.
func (cks *CKSProtocol) GenShare(skInput, skOutput *ring.Poly, ct *ckks.Ciphertext, shareOut CKSShare) {

	cks.dckksContext.ringQ.SubNoMod(skInput, skOutput, cks.tmpDelta)

	cks.genShareDelta(cks.tmpDelta, ct, shareOut)
}

func (cks *CKSProtocol) genShareDelta(skDelta *ring.Poly, ct *ckks.Ciphertext, shareOut CKSShare) {

	ringQ := cks.dckksContext.ringQ
	ringP := cks.dckksContext.ringP

	ringQ.MulCoeffsMontgomeryConstantLvl(ct.Level(), ct.Value()[1], skDelta, shareOut)

	ringQ.MulScalarBigintLvl(ct.Level(), shareOut, ringP.ModulusBigint, shareOut)

	cks.gaussianSampler.ReadLvl(ct.Level(), cks.tmpQ, cks.dckksContext.ringQP, cks.sigmaSmudging, uint64(6*cks.sigmaSmudging))
	extendBasisSmallNormAndCenter(ringQ, ringP, cks.tmpQ, cks.tmpP)

	ringQ.NTTLvl(ct.Level(), cks.tmpQ, cks.tmpQ)
	ringP.NTT(cks.tmpP, cks.tmpP)

	ringQ.AddLvl(ct.Level(), shareOut, cks.tmpQ, shareOut)

	cks.baseconverter.ModDownSplitNTTPQ(ct.Level(), shareOut, cks.tmpP, shareOut)

	cks.tmpQ.Zero()
	cks.tmpP.Zero()
}

// AggregateShares is the second part of the unique round of the CKSProtocol protocol. Upon receiving the j-1 elements each party computes :
//
// [ctx[0] + sum((skInput_i - skOutput_i) * ctx[0] + e_i), ctx[1]]
func (cks *CKSProtocol) AggregateShares(share1, share2, shareOut CKSShare) {
	cks.dckksContext.ringQ.AddLvl(uint64(len(share1.Coeffs)-1), share1, share2, shareOut)
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (cks *CKSProtocol) KeySwitch(combined CKSShare, ct *ckks.Ciphertext, ctOut *ckks.Ciphertext) {
	ctOut.SetScale(ct.Scale())
	cks.dckksContext.ringQ.AddLvl(ct.Level(), ct.Value()[0], combined, ctOut.Value()[0])
	cks.dckksContext.ringQ.CopyLvl(ct.Level(), ct.Value()[1], ctOut.Value()[1])
}
