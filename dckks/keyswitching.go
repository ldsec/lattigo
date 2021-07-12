package dckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// CKSProtocol is a structure storing the parameters for the collective key-switching protocol.
type CKSProtocol struct {
	drlwe.CKSProtocol
}

// NewCKSProtocol creates a new CKSProtocol that will be used to operate a collective key-switching on a ciphertext encrypted under a collective public-key, whose
// secret-shares are distributed among j parties, re-encrypting the ciphertext under another public-key, whose secret-shares are also known to the
// parties.
func NewCKSProtocol(params ckks.Parameters, sigmaSmudging float64) (cks *CKSProtocol) {
	return &CKSProtocol{*drlwe.NewCKSProtocol(params.Parameters, sigmaSmudging)}
}

// KeySwitch performs the actual keyswitching operation on a ciphertext ct and put the result in ctOut
func (cks *CKSProtocol) KeySwitch(combined *drlwe.CKSShare, ct rlwe.Ciphertext, ctOut rlwe.Ciphertext) {
	ctOut.(*ckks.Ciphertext).Scale = (ct.(*ckks.Ciphertext).Scale)
	cks.CKSProtocol.KeySwitch(combined, ct, ctOut)
}
