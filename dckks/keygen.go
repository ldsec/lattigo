//Package dckks implements a distributed (or threshold) version of the CKKS scheme that enables secure multiparty computation solutions with secret-shared secret keys.
package dckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/drlwe"
)

// CKGProtocol is the structure storing the parameters and state for a party in the collective key generation protocol.
type CKGProtocol struct {
	drlwe.CKGProtocol
}

// NewCKGProtocol creates a new CKGProtocol instance
func NewCKGProtocol(params ckks.Parameters) *CKGProtocol {

	ckg := new(CKGProtocol)
	ckg.CKGProtocol = *drlwe.NewCKGProtocol(params.Parameters)
	return ckg
}

// RKGProtocol is the structure storing the parameters and state for a party in the collective relinearization key
// generation protocol.
type RKGProtocol struct {
	drlwe.RKGProtocol
}

// NewRKGProtocol creates a new RKGProtocol object that will be used to generate a collective evaluation-key
// among j parties in the given context with the given bit-decomposition.
func NewRKGProtocol(params ckks.Parameters) *RKGProtocol {
	return &RKGProtocol{*drlwe.NewRKGProtocol(params.Parameters, 0.5)}
}

// RTGProtocol is the structure storing the parameters for the collective rotation-keys generation.
type RTGProtocol struct {
	drlwe.RTGProtocol
}

// NewRotKGProtocol creates a new rotkg object and will be used to generate collective rotation-keys from a shared secret-key among j parties.
func NewRotKGProtocol(params ckks.Parameters) (rtg *RTGProtocol) {
	return &RTGProtocol{*drlwe.NewRTGProtocol(params.Parameters)}
}
