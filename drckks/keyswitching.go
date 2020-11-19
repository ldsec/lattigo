package drckks

import (
	"github.com/ldsec/lattigo/v2/rckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// CKSProtocol is a structure storing the parameters for the collective key-switching protocol.
type CKSProtocol struct {
	drckksContext *drckksContext

	sigmaSmudging float64

	tmp      *ring.Poly
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
func NewCKSProtocol(params *rckks.Parameters, sigmaSmudging float64) (cks *CKSProtocol) {

	cks = new(CKSProtocol)

	drckksContext := newDrckksContext(params)

	cks.drckksContext = drckksContext

	cks.tmp = drckksContext.ringQP.NewPoly()
	cks.tmpDelta = drckksContext.ringQ.NewPoly()
	cks.hP = drckksContext.ringP.NewPoly()

	cks.baseconverter = ring.NewFastBasisExtender(drckksContext.ringQ, drckksContext.ringP)
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	cks.gaussianSampler = ring.NewGaussianSampler(prng, drckksContext.ringQP, params.Sigma(), uint64(6*params.Sigma()))

	return cks
}

// AllocateShare allocates the share of the CKS protocol.
func (cks *CKSProtocol) AllocateShare() CKSShare {
	return cks.drckksContext.ringQ.NewPoly()
}

// GenShare is the first and unique round of the CKSProtocol protocol. Each party holding a ciphertext ctx encrypted under a collective publick-key must
// compute the following :
//
// [(skInput_i - skOutput_i) * ctx[0] + e_i]
//
// Each party then broadcasts the result of this computation to the other j-1 parties.
func (cks *CKSProtocol) GenShare(skInput, skOutput *ring.Poly, ct *rckks.Ciphertext, shareOut CKSShare) {

	cks.drckksContext.ringQ.Sub(skInput, skOutput, cks.tmpDelta)

	cks.genShareDelta(cks.tmpDelta, ct, shareOut)
}

func (cks *CKSProtocol) genShareDelta(skDelta *ring.Poly, ct *rckks.Ciphertext, shareOut CKSShare) {

	ringQ := cks.drckksContext.ringQ
	ringP := cks.drckksContext.ringP

	ringQ.MulCoeffsMontgomeryLvl(ct.Level(), ct.Value()[1], skDelta, shareOut)

	ringQ.MulScalarBigintLvl(ct.Level(), shareOut, ringP.ModulusBigint, shareOut)

	// TODO : improve by only computing the NTT for the required primes
	cks.gaussianSampler.Read(cks.tmp)
	rckks.NTTRCKKS(cks.drckksContext.ringQP, cks.tmp, cks.tmp)

	ringQ.AddLvl(ct.Level(), shareOut, cks.tmp, shareOut)

	for x, i := 0, uint64(len(ringQ.Modulus)); i < uint64(len(cks.drckksContext.ringQP.Modulus)); x, i = x+1, i+1 {
		tmp0 := cks.tmp.Coeffs[i]
		tmp1 := cks.hP.Coeffs[x]
		for j := uint64(0); j < ringQ.N; j++ {
			tmp1[j] += tmp0[j]
		}
	}

	rckks.ModDownSplitNTTPQRCKKS(cks.baseconverter, ct.Level(), shareOut, cks.hP, shareOut)

	cks.hP.Zero()
	cks.tmp.Zero()
}

// AggregateShares is the second part of the unique round of the CKSProtocol protocol. Upon receiving the j-1 elements each party computes :
//
// [ctx[0] + sum((skInput_i - skOutput_i) * ctx[0] + e_i), ctx[1]]
func (cks *CKSProtocol) AggregateShares(share1, share2, shareOut CKSShare) {
	cks.drckksContext.ringQ.AddLvl(uint64(len(share1.Coeffs)-1), share1, share2, shareOut)
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (cks *CKSProtocol) KeySwitch(combined CKSShare, ct *rckks.Ciphertext, ctOut *rckks.Ciphertext) {
	ctOut.SetScale(ct.Scale())
	cks.drckksContext.ringQ.AddLvl(ct.Level(), ct.Value()[0], combined, ctOut.Value()[0])
	cks.drckksContext.ringQ.CopyLvl(ct.Level(), ct.Value()[1], ctOut.Value()[1])
}
