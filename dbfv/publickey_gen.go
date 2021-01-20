//Package dbfv implements a distributed (or threshold) version of the BFV scheme that enables secure multiparty computation solutions with secret-shared secret keys.
package dbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// CKGProtocol is the structure storing the parameters and state for a party in the collective key generation protocol.
type CKGProtocol struct {
	ringQP          *ring.Ring
	sigma           float64
	gaussianSampler *ring.GaussianSampler
}

// CKGShare is a struct holding a CKG share.
type CKGShare struct {
	*ring.Poly
}

// UnmarshalBinary decode a marshaled CKG share on the target CKG share.
func (share *CKGShare) UnmarshalBinary(data []byte) error {
	if share.Poly == nil {
		share.Poly = new(ring.Poly)
	}
	err := share.Poly.UnmarshalBinary(data)
	return err

}

// NewCKGProtocol creates a new CKGProtocol instance
func NewCKGProtocol(params *bfv.Parameters) *CKGProtocol {
	context := newDbfvContext(params)
	ckg := new(CKGProtocol)
	ckg.ringQP = context.ringQP
	ckg.sigma = params.Sigma()
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ckg.gaussianSampler = ring.NewGaussianSampler(prng)
	return ckg
}

// AllocateShares allocates the CKG shares.
func (ckg *CKGProtocol) AllocateShares() CKGShare {
	return CKGShare{ckg.ringQP.NewPoly()}
}

// GenShare generates the party's public key share from its secret key as:
//
// crs*s_i + e_i
//
// for the receiver protocol. Has no effect is the share was already generated.
func (ckg *CKGProtocol) GenShare(sk *ring.Poly, crs *ring.Poly, shareOut CKGShare) {
	ckg.gaussianSampler.Read(shareOut.Poly, ckg.ringQP, ckg.sigma, uint64(6*ckg.sigma))
	ckg.ringQP.NTT(shareOut.Poly, shareOut.Poly)
	ckg.ringQP.MulCoeffsMontgomeryAndSub(sk, crs, shareOut.Poly)
}

// AggregateShares aggregates a new share to the aggregate key
func (ckg *CKGProtocol) AggregateShares(share1, share2, shareOut CKGShare) {
	ckg.ringQP.Add(share1.Poly, share2.Poly, shareOut.Poly)
}

// GenPublicKey return the current aggregation of the received shares as a bfv.PublicKey.
func (ckg *CKGProtocol) GenPublicKey(roundShare CKGShare, crs *ring.Poly, pubkey *bfv.PublicKey) {
	pubkey.Set([2]*ring.Poly{roundShare.Poly, crs})
}
