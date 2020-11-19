//Package drckks implements a distributed (or threshold) version of the CKKS scheme that enables secure multiparty computation solutions with secret-shared secret keys.
package drckks

import (
	"github.com/ldsec/lattigo/v2/rckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// CKGProtocol is the structure storing the parameters and state for a party in the collective key generation protocol.
type CKGProtocol struct {
	drckksContext   *drckksContext
	gaussianSampler *ring.GaussianSampler
}

// CKGShare is a struct storing the CKG protocol's share.
type CKGShare *ring.Poly

// NewCKGProtocol creates a new CKGProtocol instance
func NewCKGProtocol(params *rckks.Parameters) *CKGProtocol {

	ckg := new(CKGProtocol)
	ckg.drckksContext = newDrckksContext(params)
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ckg.gaussianSampler = ring.NewGaussianSampler(prng, ckg.drckksContext.ringQP, params.Sigma(), uint64(6*params.Sigma()))
	return ckg
}

// AllocateShares allocates the share of the CKG protocol.
func (ckg *CKGProtocol) AllocateShares() CKGShare {
	return ckg.drckksContext.ringQP.NewPoly()
}

// GenShare generates the party's public key share from its secret key as:
//
// crs*s_i + e_i
//
// for the receiver protocol. Has no effect is the share was already generated.
func (ckg *CKGProtocol) GenShare(sk *ring.Poly, crs *ring.Poly, shareOut CKGShare) {
	ringQP := ckg.drckksContext.ringQP

	ckg.gaussianSampler.Read(shareOut)
	rckks.NTTRCKKS(ringQP, shareOut, shareOut)
	ringQP.MulCoeffsMontgomeryAndSub(sk, crs, shareOut)
}

// AggregateShares aggregates a new share to the aggregate key
func (ckg *CKGProtocol) AggregateShares(share1, share2, shareOut CKGShare) {
	ckg.drckksContext.ringQP.Add(share1, share2, shareOut)
}

// GenPublicKey return the current aggregation of the received shares as a bfv.PublicKey.
func (ckg *CKGProtocol) GenPublicKey(roundShare CKGShare, crs *ring.Poly, pubkey *rckks.PublicKey) {
	pubkey.Set([2]*ring.Poly{roundShare, crs})
}
