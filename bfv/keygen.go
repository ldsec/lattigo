package bfv

import (
	"github.com/ldsec/lattigo/v2/rlwe"
)

// KeyGenerator is an interface implementing the methods of the keyGenerator.
type KeyGenerator interface {
	GenSecretKey() (sk *rlwe.SecretKey)
	GenSecretKeyGaussian() (sk *rlwe.SecretKey)
	GenSecretKeyWithDistrib(p float64) (sk *rlwe.SecretKey)
	GenSecretKeySparse(hw int) (sk *rlwe.SecretKey)
	GenPublicKey(sk *rlwe.SecretKey) (pk *rlwe.PublicKey)
	GenKeyPair() (sk *rlwe.SecretKey, pk *rlwe.PublicKey)
	GenKeyPairSparse(hw int) (sk *rlwe.SecretKey, pk *rlwe.PublicKey)
	GenRelinearizationKey(sk *rlwe.SecretKey, maxDegree int) (evk *rlwe.RelinearizationKey)
	GenSwitchingKey(skInput, skOutput *rlwe.SecretKey) (newevakey *rlwe.SwitchingKey)
	GenSwitchingKeyForGalois(galEl uint64, sk *rlwe.SecretKey) (swk *rlwe.SwitchingKey)
	GenRotationKeys(galEls []uint64, sk *rlwe.SecretKey) (rks *rlwe.RotationKeySet)
	GenSwitchingKeyForRotationBy(k int, sk *rlwe.SecretKey) (swk *rlwe.SwitchingKey)
	GenRotationKeysForRotations(ks []int, inclueSwapRows bool, sk *rlwe.SecretKey) (rks *rlwe.RotationKeySet)
	GenSwitchingKeyForRowRotation(sk *rlwe.SecretKey) (swk *rlwe.SwitchingKey)
	GenRotationKeysForInnerSum(sk *rlwe.SecretKey) (rks *rlwe.RotationKeySet)
}

// keyGenerator is a structure that stores the elements required to create new keys,
// as well as a small memory pool for intermediate values.
type keyGenerator struct {
	rlwe.KeyGenerator
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as the evaluation,
// rotation and switching keys can be generated.
func NewKeyGenerator(params Parameters) KeyGenerator {
	return &keyGenerator{KeyGenerator: rlwe.NewKeyGenerator(params.Parameters)}
}
