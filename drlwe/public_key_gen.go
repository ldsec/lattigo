//Package drlwe implements a distributed (or threshold) version of the CKKS scheme that enables secure multiparty computation solutions with secret-shared secret keys.
package drlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// CollectivePublicKeyGenerator is an interface describing the local steps of a generic RLWE CKG protocol.
type CollectivePublicKeyGenerator interface {
	AllocateShares() *CKGShare
	GenShare(sk *rlwe.SecretKey, crs *ring.Poly, shareOut *CKGShare)
	AggregateShares(share1, share2, shareOut *CKGShare)
	GenPublicKey(aggregatedShare *CKGShare, crs *ring.Poly, pubkey *rlwe.PublicKey)
}

// CKGProtocol is the structure storing the parameters and and precomputations for the collective key generation protocol.
type CKGProtocol struct {
	params rlwe.Parameters

	ringQP          *ring.Ring
	gaussianSampler *ring.GaussianSampler
}

// CKGShare is a struct storing the CKG protocol's share.
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
func NewCKGProtocol(params rlwe.Parameters) *CKGProtocol { // TODO drlwe.Params

	ckg := new(CKGProtocol)
	ckg.params = params
	ckg.ringQP = params.RingQP()

	var err error
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ckg.gaussianSampler = ring.NewGaussianSampler(prng, ckg.ringQP, params.Sigma(), int(6*params.Sigma()))
	return ckg
}

// AllocateShares allocates the share of the CKG protocol.
func (ckg *CKGProtocol) AllocateShares() *CKGShare {
	return &CKGShare{ckg.ringQP.NewPoly()}
}

// GenShare generates the party's public key share from its secret key as:
//
// crs*s_i + e_i
//
// for the receiver protocol. Has no effect is the share was already generated.
func (ckg *CKGProtocol) GenShare(sk *rlwe.SecretKey, crs *ring.Poly, shareOut *CKGShare) {
	ringQP := ckg.ringQP
	ckg.gaussianSampler.Read(shareOut.Poly)
	ringQP.NTT(shareOut.Poly, shareOut.Poly)
	ringQP.MulCoeffsMontgomeryAndSub(sk.Value, crs, shareOut.Poly)
}

// AggregateShares aggregates a new share to the aggregate key
func (ckg *CKGProtocol) AggregateShares(share1, share2, shareOut *CKGShare) {
	ckg.ringQP.Add(share1.Poly, share2.Poly, shareOut.Poly)
}

// GenPublicKey return the current aggregation of the received shares as a bfv.PublicKey.
func (ckg *CKGProtocol) GenPublicKey(roundShare *CKGShare, crs *ring.Poly, pubkey *rlwe.PublicKey) {
	pubkey.Value[0].Copy(roundShare.Poly)
	pubkey.Value[1].Copy(crs)
}
