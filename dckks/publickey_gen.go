//Package dckks implements a distributed (or threshold) version of the CKKS scheme that enables secure multiparty computation solutions with secret-shared secret keys.
package dckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// CKGProtocol is the structure storing the parameters and state for a party in the collective key generation protocol.
type CKGProtocol struct {
	dckksContext    *dckksContext
	gaussianSampler *ring.GaussianSampler
	sigma           float64
}

// CKGShare is a struct storing the CKG protocol's share.
type CKGShare *ring.Poly

// NewCKGProtocol creates a new CKGProtocol instance
func NewCKGProtocol(params *ckks.Parameters) *CKGProtocol {

	ckg := new(CKGProtocol)
	ckg.dckksContext = newDckksContext(params)
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ckg.gaussianSampler = ring.NewGaussianSampler(prng)
	ckg.sigma = params.Sigma()
	return ckg
}

// AllocateShares allocates the share of the CKG protocol.
func (ckg *CKGProtocol) AllocateShares() CKGShare {
	return ckg.dckksContext.ringQP.NewPoly()
}

// GenShare generates the party's public key share from its secret key as:
//
// crs*s_i + e_i
//
// for the receiver protocol. Has no effect is the share was already generated.
func (ckg *CKGProtocol) GenShare(sk *ring.Poly, crs *ring.Poly, shareOut CKGShare) {
	ringQP := ckg.dckksContext.ringQP

	ckg.gaussianSampler.Read(shareOut, ringQP, ckg.sigma, uint64(6*ckg.sigma))
	ringQP.NTT(shareOut, shareOut)
	ringQP.MulCoeffsMontgomeryAndSub(sk, crs, shareOut)
}

// AggregateShares aggregates a new share to the aggregate key
func (ckg *CKGProtocol) AggregateShares(share1, share2, shareOut CKGShare) {
	ckg.dckksContext.ringQP.Add(share1, share2, shareOut)
}

// GenPublicKey return the current aggregation of the received shares as a bfv.PublicKey.
func (ckg *CKGProtocol) GenPublicKey(roundShare CKGShare, crs *ring.Poly, pubkey *ckks.PublicKey) {
	pubkey.Set([2]*ring.Poly{roundShare, crs})
}
