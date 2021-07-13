package bfv

import "github.com/ldsec/lattigo/v2/rlwe"

// NewKeyGenerator creates a rlwe.KeyGenerator instance from the BFV parameters.
func NewKeyGenerator(params Parameters) rlwe.KeyGenerator {
	return rlwe.NewKeyGenerator(params.Parameters)
}

// NewSecretKey returns an allocated BFV secret key with zero values.
func NewSecretKey(params Parameters) (sk *rlwe.SecretKey) {
	return rlwe.NewSecretKey(params.Parameters)
}

// NewPublicKey returns an allocated BFV public with zero values.
func NewPublicKey(params Parameters) (pk *rlwe.PublicKey) {
	return rlwe.NewPublicKey(params.Parameters)
}

// NewSwitchingKey returns an allocated BFV public switching key with zero values.
func NewSwitchingKey(params Parameters) *rlwe.SwitchingKey {
	return rlwe.NewSwitchingKey(params.Parameters)
}

// NewRelinearizationKey returns an allocated BFV public relinearization key with zero value for each degree in [2 < maxRelinDegree].
func NewRelinearizationKey(params Parameters, maxRelinDegree int) *rlwe.RelinearizationKey {
	return rlwe.NewRelinKey(params.Parameters, maxRelinDegree)
}

// NewRotationKeySet returns an allocated set of BFV public rotation keys with zero values for each galois element
// (i.e., for each supported rotation).
func NewRotationKeySet(params Parameters, galoisElements []uint64) *rlwe.RotationKeySet {
	return rlwe.NewRotationKeySet(params.Parameters, galoisElements)
}
