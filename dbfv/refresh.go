package dbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
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
func NewRefreshProtocol(params bfv.Parameters, sigmaSmudging float64) (rfp *RefreshProtocol) {
	rfp = new(RefreshProtocol)
	rfp.MaskedTransformProtocol = *NewMaskedTransformProtocol(params, sigmaSmudging)
	return
}

// AllocateShare allocates the shares of the PermuteProtocol
func (rfp *RefreshProtocol) AllocateShare() *RefreshShare {
	return &RefreshShare{*rfp.MaskedTransformProtocol.AllocateShare()}
}

// GenShares generates a share for the Refresh protocol.
func (rfp *RefreshProtocol) GenShares(sk *rlwe.SecretKey, ciphertext *bfv.Ciphertext, crs *ring.Poly, shareOut *RefreshShare) {
	rfp.MaskedTransformProtocol.GenShares(sk, ciphertext, crs, nil, &shareOut.MaskedTransformShare)
}

// Aggregate aggregates two parties' shares in the Refresh protocol.
func (rfp *RefreshProtocol) Aggregate(share1, share2, shareOut *RefreshShare) {
	rfp.MaskedTransformProtocol.Aggregate(&share1.MaskedTransformShare, &share2.MaskedTransformShare, &shareOut.MaskedTransformShare)
}

// Finalize applies Decrypt, Recode and Recrypt on the input ciphertext.
func (rfp *RefreshProtocol) Finalize(ciphertext *bfv.Ciphertext, crs *ring.Poly, share *RefreshShare, ciphertextOut *bfv.Ciphertext) {
	rfp.MaskedTransformProtocol.Transform(ciphertext, nil, crs, &share.MaskedTransformShare, ciphertextOut)
}
