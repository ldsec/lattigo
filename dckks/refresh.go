package dckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// RefreshProtocol is a struct storing the relevant parameters for the Refresh protocol.
type RefreshProtocol struct {
	MaskedTransformProtocol
}

// RefreshShare is a struct storing a party's share in the Refresh protocol.
type RefreshShare struct {
	MaskedTransformShare
}

// NewRefreshProtocol creates a new Refresh protocol instance.
func NewRefreshProtocol(params ckks.Parameters, precision int, sigmaSmudging float64) (rfp *RefreshProtocol) {
	rfp = new(RefreshProtocol)
	rfp.MaskedTransformProtocol = *NewMaskedTransformProtocol(params, precision, sigmaSmudging)
	return
}

// AllocateShare allocates the shares of the PermuteProtocol
func (rfp *RefreshProtocol) AllocateShare(minLevel, maxLevel int) *RefreshShare {
	return &RefreshShare{*rfp.MaskedTransformProtocol.AllocateShare(minLevel, maxLevel)}
}

// GenShares generates a share for the Refresh protocol.
// This protocol requires additional inputs which are :
// logBound : the bit length of the masks
// logSlots : the bit length of the number of slots
//
// The method "GetMinimumLevelForBootstrapping" should be used to get the minimum level at which the refresh can be called while still ensure 128-bits of security, as well as the
// value for logBound.
func (rfp *RefreshProtocol) GenShares(sk *rlwe.SecretKey, logBound, logSlots int, ciphertext *ckks.Ciphertext, crs *ring.Poly, shareOut *RefreshShare) {
	rfp.MaskedTransformProtocol.GenShares(sk, logBound, logSlots, ciphertext, crs, nil, &shareOut.MaskedTransformShare)
}

// Aggregate aggregates two parties' shares in the Refresh protocol.
func (rfp *RefreshProtocol) Aggregate(share1, share2, shareOut *RefreshShare) {
	rfp.MaskedTransformProtocol.Aggregate(&share1.MaskedTransformShare, &share2.MaskedTransformShare, &shareOut.MaskedTransformShare)
}

// Finalize applies Decrypt, Recode and Recrypt on the input ciphertext.
func (rfp *RefreshProtocol) Finalize(ciphertext *ckks.Ciphertext, logSlots int, crs *ring.Poly, share *RefreshShare, ciphertextOut *ckks.Ciphertext) {
	rfp.MaskedTransformProtocol.Transform(ciphertext, logSlots, nil, crs, &share.MaskedTransformShare, ciphertextOut)
}
