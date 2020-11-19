package drckks

import (
	"github.com/ldsec/lattigo/v2/rckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// PCKSProtocol is the structure storing the parameters for the collective public key-switching.
type PCKSProtocol struct {
	drckksContext *drckksContext

	sigmaSmudging float64

	tmpQ *ring.Poly
	tmpP *ring.Poly

	share0tmpQ *ring.Poly
	share1tmpQ *ring.Poly
	share0tmpP *ring.Poly
	share1tmpP *ring.Poly

	baseconverter   *ring.FastBasisExtender
	gaussianSampler *ring.GaussianSampler
	ternarySampler  *ring.TernarySampler
}

// PCKSShare is a struct storing the share of the PCKS protocol.
type PCKSShare [2]*ring.Poly

// NewPCKSProtocol creates a new PCKSProtocol object and will be used to re-encrypt a ciphertext ctx encrypted under a secret-shared key mong j parties under a new
// collective public-key.
func NewPCKSProtocol(params *rckks.Parameters, sigmaSmudging float64) *PCKSProtocol {

	pcks := new(PCKSProtocol)

	drckksContext := newDrckksContext(params)

	pcks.drckksContext = drckksContext

	pcks.tmpQ = drckksContext.ringQ.NewPoly()
	pcks.tmpP = drckksContext.ringP.NewPoly()
	pcks.share0tmpQ = drckksContext.ringQ.NewPoly()
	pcks.share1tmpQ = drckksContext.ringQ.NewPoly()
	pcks.share0tmpP = drckksContext.ringP.NewPoly()
	pcks.share1tmpP = drckksContext.ringP.NewPoly()

	pcks.baseconverter = ring.NewFastBasisExtender(drckksContext.ringQ, drckksContext.ringP)
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	pcks.gaussianSampler = ring.NewGaussianSampler(prng, drckksContext.ringQ, params.Sigma(), uint64(6*params.Sigma()))
	pcks.ternarySampler = ring.NewTernarySampler(prng, drckksContext.ringQ, 0.5, false)

	return pcks
}

// AllocateShares allocates the share of the PCKS protocol.
func (pcks *PCKSProtocol) AllocateShares(level uint64) (s PCKSShare) {
	s[0] = pcks.drckksContext.ringQ.NewPolyLvl(level)
	s[1] = pcks.drckksContext.ringQ.NewPolyLvl(level)
	return
}

// GenShare is the first part of the unique round of the PCKSProtocol protocol. Each party computes the following :
//
// [s_i * ctx[0] + u_i * pk[0] + e_0i, u_i * pk[1] + e_1i]
//
// and broadcasts the result to the other j-1 parties.
func (pcks *PCKSProtocol) GenShare(sk *ring.Poly, pk *rckks.PublicKey, ct *rckks.Ciphertext, shareOut PCKSShare) {

	lvl := ct.Level()

	ringQ := pcks.drckksContext.ringQ
	ringP := pcks.drckksContext.ringP

	pcks.ternarySampler.ReadLvl(lvl, pcks.tmpQ)
	extendBasisSmallNormAndCenter(ringQ, ringP, pcks.tmpQ, pcks.tmpP)

	rckks.NTTRCKKSLvl(ringQ, lvl, pcks.tmpQ, pcks.tmpQ)
	rckks.NTTRCKKS(ringP, pcks.tmpP, pcks.tmpP)

	ringQ.MFormLvl(lvl, pcks.tmpQ, pcks.tmpQ)
	ringP.MForm(pcks.tmpP, pcks.tmpP)

	pk0P := new(ring.Poly)
	pk1P := new(ring.Poly)
	pk0P.Coeffs = pk.Get()[0].Coeffs[len(ringQ.Modulus):]
	pk1P.Coeffs = pk.Get()[1].Coeffs[len(ringQ.Modulus):]

	// h_0 = u_i * pk_0
	ringQ.MulCoeffsMontgomeryLvl(lvl, pcks.tmpQ, pk.Get()[0], pcks.share0tmpQ)
	ringQ.MulCoeffsMontgomeryLvl(lvl, pcks.tmpQ, pk.Get()[1], pcks.share1tmpQ)
	// h_1 = u_i * pk_1
	ringP.MulCoeffsMontgomery(pcks.tmpP, pk0P, pcks.share0tmpP)
	ringP.MulCoeffsMontgomery(pcks.tmpP, pk1P, pcks.share1tmpP)

	// h_0 = u_i * pk_0 + e0
	pcks.gaussianSampler.ReadLvl(lvl, pcks.tmpQ)
	extendBasisSmallNormAndCenter(ringQ, ringP, pcks.tmpQ, pcks.tmpP)
	rckks.NTTRCKKSLvl(ringQ, lvl, pcks.tmpQ, pcks.tmpQ)
	rckks.NTTRCKKS(ringP, pcks.tmpP, pcks.tmpP)
	ringQ.AddLvl(lvl, pcks.share0tmpQ, pcks.tmpQ, pcks.share0tmpQ)
	ringP.Add(pcks.share0tmpP, pcks.tmpP, pcks.share0tmpP)

	// h_1 = u_i * pk_1 + e1
	pcks.gaussianSampler.ReadLvl(lvl, pcks.tmpQ)
	extendBasisSmallNormAndCenter(ringQ, ringP, pcks.tmpQ, pcks.tmpP)
	rckks.NTTRCKKSLvl(ringQ, lvl, pcks.tmpQ, pcks.tmpQ)
	rckks.NTTRCKKS(ringP, pcks.tmpP, pcks.tmpP)
	ringQ.AddLvl(lvl, pcks.share1tmpQ, pcks.tmpQ, pcks.share1tmpQ)
	ringP.Add(pcks.share1tmpP, pcks.tmpP, pcks.share1tmpP)

	// h_0 = (u_i * pk_0 + e0)/P
	rckks.ModDownSplitNTTPQRCKKS(pcks.baseconverter, lvl, pcks.share0tmpQ, pcks.share0tmpP, shareOut[0])

	// h_1 = (u_i * pk_1 + e1)/P
	// Cound be moved to the keyswitch part of the protocol, but the second element of the shares will be larger.
	rckks.ModDownSplitNTTPQRCKKS(pcks.baseconverter, lvl, pcks.share1tmpQ, pcks.share1tmpP, shareOut[1])

	// h_0 = s_i*c_1 + (u_i * pk_0 + e0)/P
	ringQ.MulCoeffsMontgomeryAndAddLvl(lvl, ct.Value()[1], sk, shareOut[0])

	pcks.tmpQ.Zero()
	pcks.tmpP.Zero()
}

// AggregateShares is the second part of the first and unique round of the PCKSProtocol protocol. Each party uppon receiving the j-1 elements from the
// other parties computes :
//
// [ctx[0] + sum(s_i * ctx[0] + u_i * pk[0] + e_0i), sum(u_i * pk[1] + e_1i)]
func (pcks *PCKSProtocol) AggregateShares(share1, share2, shareOut PCKSShare) {

	level := uint64(len(share1[0].Coeffs)) - 1
	pcks.drckksContext.ringQ.AddLvl(level, share1[0], share2[0], shareOut[0])
	pcks.drckksContext.ringQ.AddLvl(level, share1[1], share2[1], shareOut[1])
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (pcks *PCKSProtocol) KeySwitch(combined PCKSShare, ct, ctOut *rckks.Ciphertext) {

	ctOut.SetScale(ct.Scale())

	pcks.drckksContext.ringQ.AddLvl(ct.Level(), ct.Value()[0], combined[0], ctOut.Value()[0])
	pcks.drckksContext.ringQ.CopyLvl(ct.Level(), combined[1], ctOut.Value()[1])
}
