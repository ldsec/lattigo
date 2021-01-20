package dckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// PCKSProtocol is the structure storing the parameters for the collective public key-switching.
type PCKSProtocol struct {
	dckksContext *dckksContext

	sigmaSmudging float64

	tmp *ring.Poly

	share0tmp *ring.Poly
	share1tmp *ring.Poly

	baseconverter            *ring.FastBasisExtender
	gaussianSampler          *ring.GaussianSampler
	ternarySamplerMontgomery *ring.TernarySampler
}

// PCKSShare is a struct storing the share of the PCKS protocol.
type PCKSShare [2]*ring.Poly

// NewPCKSProtocol creates a new PCKSProtocol object and will be used to re-encrypt a ciphertext ctx encrypted under a secret-shared key mong j parties under a new
// collective public-key.
func NewPCKSProtocol(params *ckks.Parameters, sigmaSmudging float64) *PCKSProtocol {

	pcks := new(PCKSProtocol)

	dckksContext := newDckksContext(params)

	pcks.dckksContext = dckksContext

	pcks.tmp = dckksContext.ringQP.NewPoly()
	pcks.share0tmp = dckksContext.ringQP.NewPoly()
	pcks.share1tmp = dckksContext.ringQP.NewPoly()

	pcks.sigmaSmudging = sigmaSmudging

	pcks.baseconverter = ring.NewFastBasisExtender(dckksContext.ringQ, dckksContext.ringP)
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	pcks.gaussianSampler = ring.NewGaussianSampler(prng)
	pcks.ternarySamplerMontgomery = ring.NewTernarySampler(prng, dckksContext.ringQP, 0.5, true)

	return pcks
}

// AllocateShares allocates the share of the PCKS protocol.
func (pcks *PCKSProtocol) AllocateShares(level uint64) (s PCKSShare) {
	s[0] = pcks.dckksContext.ringQ.NewPolyLvl(level)
	s[1] = pcks.dckksContext.ringQ.NewPolyLvl(level)
	return
}

// GenShare is the first part of the unique round of the PCKSProtocol protocol. Each party computes the following :
//
// [s_i * ctx[0] + u_i * pk[0] + e_0i, u_i * pk[1] + e_1i]
//
// and broadcasts the result to the other j-1 parties.
func (pcks *PCKSProtocol) GenShare(sk *ring.Poly, pk *ckks.PublicKey, ct *ckks.Ciphertext, shareOut PCKSShare) {

	// Planned improvement : adapt share size to ct.Level() to improve efficiency.

	ringQ := pcks.dckksContext.ringQ
	ringQP := pcks.dckksContext.ringQP

	pcks.ternarySamplerMontgomery.Read(pcks.tmp)
	ringQP.NTTLazy(pcks.tmp, pcks.tmp)

	// h_0 = u_i * pk_0
	ringQP.MulCoeffsMontgomeryConstant(pcks.tmp, pk.Get()[0], pcks.share0tmp)
	// h_1 = u_i * pk_1
	ringQP.MulCoeffsMontgomeryConstant(pcks.tmp, pk.Get()[1], pcks.share1tmp)

	// h_0 = u_i * pk_0 + e0
	pcks.gaussianSampler.Read(pcks.tmp, ringQP, pcks.sigmaSmudging, uint64(6*pcks.sigmaSmudging))
	ringQP.NTT(pcks.tmp, pcks.tmp)
	ringQP.Add(pcks.share0tmp, pcks.tmp, pcks.share0tmp)

	// h_1 = u_i * pk_1 + e1
	pcks.gaussianSampler.Read(pcks.tmp, ringQP, pcks.sigmaSmudging, uint64(6*pcks.sigmaSmudging))
	ringQP.NTT(pcks.tmp, pcks.tmp)
	ringQP.Add(pcks.share1tmp, pcks.tmp, pcks.share1tmp)

	// h_0 = (u_i * pk_0 + e0)/P
	pcks.baseconverter.ModDownNTTPQ(ct.Level(), pcks.share0tmp, shareOut[0])

	// h_1 = (u_i * pk_1 + e1)/P
	// Cound be moved to the keyswitch part of the protocol, but the second element of the shares will be larger.
	pcks.baseconverter.ModDownNTTPQ(ct.Level(), pcks.share1tmp, shareOut[1])

	// h_0 = s_i*c_1 + (u_i * pk_0 + e0)/P
	ringQ.MulCoeffsMontgomeryAndAddLvl(ct.Level(), ct.Value()[1], sk, shareOut[0])

	pcks.tmp.Zero()
}

// AggregateShares is the second part of the first and unique round of the PCKSProtocol protocol. Each party uppon receiving the j-1 elements from the
// other parties computes :
//
// [ctx[0] + sum(s_i * ctx[0] + u_i * pk[0] + e_0i), sum(u_i * pk[1] + e_1i)]
func (pcks *PCKSProtocol) AggregateShares(share1, share2, shareOut PCKSShare) {

	level := uint64(len(share1[0].Coeffs)) - 1
	pcks.dckksContext.ringQ.AddLvl(level, share1[0], share2[0], shareOut[0])
	pcks.dckksContext.ringQ.AddLvl(level, share1[1], share2[1], shareOut[1])
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (pcks *PCKSProtocol) KeySwitch(combined PCKSShare, ct, ctOut *ckks.Ciphertext) {

	ctOut.SetScale(ct.Scale())

	pcks.dckksContext.ringQ.AddLvl(ct.Level(), ct.Value()[0], combined[0], ctOut.Value()[0])
	pcks.dckksContext.ringQ.CopyLvl(ct.Level(), combined[1], ctOut.Value()[1])
}
