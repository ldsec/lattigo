package ckks

import "github.com/ldsec/lattigo/v2/rlwe"

// NewKeyGenerator creates a rlwe.KeyGenerator instance from the CKKS parameters.
func NewKeyGenerator(params Parameters) rlwe.KeyGenerator {
	return rlwe.NewKeyGenerator(params.Parameters)
}

// BootstrappingKey is a type for a CKKS bootstrapping key, wich regroups the necessary public relinearization
// and rotation keys (i.e., an EvaluationKey).
type BootstrappingKey rlwe.EvaluationKey

// NewSecretKey returns an allocated CKKS secret key with zero values.
func NewSecretKey(params Parameters) (sk *rlwe.SecretKey) {
	return rlwe.NewSecretKey(params.Parameters)
}

// NewPublicKey returns an allocated CKKS public with zero values.
func NewPublicKey(params Parameters) (pk *rlwe.PublicKey) {
	return rlwe.NewPublicKey(params.Parameters)
}

// NewSwitchingKey returns an allocated CKKS public switching key with zero values.
func NewSwitchingKey(params Parameters) *rlwe.SwitchingKey {
	return rlwe.NewSwitchingKey(params.Parameters)
}

// NewRelinearizationKey returns an allocated CKKS public relinearization key with zero value.
func NewRelinearizationKey(params Parameters) *rlwe.RelinearizationKey {
	return rlwe.NewRelinKey(params.Parameters, 2)
}

// NewRotationKeySet returns an allocated set of CKKS public rotation keys with zero values for each galois element
// (i.e., for each supported rotation).
func NewRotationKeySet(params Parameters, galoisElements []uint64) *rlwe.RotationKeySet {
	return rlwe.NewRotationKeySet(params.Parameters, galoisElements)
}
