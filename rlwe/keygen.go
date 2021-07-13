package rlwe

import (
	"math/big"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// KeyGenerator is an interface implementing the methods of the KeyGenerator.
type KeyGenerator interface {
	GenSecretKey() (sk *SecretKey)
	GenSecretKeyGaussian() (sk *SecretKey)
	GenSecretKeyWithDistrib(p float64) (sk *SecretKey)
	GenSecretKeySparse(hw int) (sk *SecretKey)
	GenPublicKey(sk *SecretKey) (pk *PublicKey)
	GenKeyPair() (sk *SecretKey, pk *PublicKey)
	GenKeyPairSparse(hw int) (sk *SecretKey, pk *PublicKey)
	GenRelinearizationKey(sk *SecretKey, maxDegree int) (evk *RelinearizationKey)
	GenSwitchingKey(skInput, skOutput *SecretKey) (newevakey *SwitchingKey)
	GenSwitchingKeyForGalois(galEl uint64, sk *SecretKey) (swk *SwitchingKey)
	GenRotationKeys(galEls []uint64, sk *SecretKey) (rks *RotationKeySet)
	GenSwitchingKeyForRotationBy(k int, sk *SecretKey) (swk *SwitchingKey)
	GenRotationKeysForRotations(ks []int, inclueSwapRows bool, sk *SecretKey) (rks *RotationKeySet)
	GenSwitchingKeyForRowRotation(sk *SecretKey) (swk *SwitchingKey)
	GenRotationKeysForInnerSum(sk *SecretKey) (rks *RotationKeySet)
}

// KeyGenerator is a structure that stores the elements required to create new keys,
// as well as a small memory pool for intermediate values.
type keyGenerator struct {
	params          Parameters
	ringQP          *ring.Ring
	pBigInt         *big.Int
	polypool        [2]*ring.Poly
	gaussianSampler *ring.GaussianSampler
	uniformSampler  *ring.UniformSampler
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as the evaluation,
// rotation and switching keys can be generated.
func NewKeyGenerator(params Parameters) KeyGenerator {

	ringQP := params.RingQP()

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	return &keyGenerator{
		params:          params,
		ringQP:          ringQP,
		pBigInt:         params.PBigInt(),
		polypool:        [2]*ring.Poly{ringQP.NewPoly(), ringQP.NewPoly()},
		gaussianSampler: ring.NewGaussianSampler(prng, ringQP, params.Sigma(), int(6*params.Sigma())),
		uniformSampler:  ring.NewUniformSampler(prng, ringQP),
	}
}

// GenSecretKey generates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *keyGenerator) GenSecretKey() (sk *SecretKey) {
	return keygen.GenSecretKeyWithDistrib(1.0 / 3)
}

func (keygen *keyGenerator) GenSecretKeyGaussian() (sk *SecretKey) {
	sk = new(SecretKey)

	sk.Value = keygen.gaussianSampler.ReadNew()
	keygen.ringQP.NTT(sk.Value, sk.Value)
	return sk
}

// GenSecretKeyWithDistrib generates a new SecretKey with the distribution [(p-1)/2, p, (p-1)/2].
func (keygen *keyGenerator) GenSecretKeyWithDistrib(p float64) (sk *SecretKey) {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ternarySamplerMontgomery := ring.NewTernarySampler(prng, keygen.ringQP, p, true)

	sk = new(SecretKey)
	sk.Value = ternarySamplerMontgomery.ReadNew()
	keygen.ringQP.NTT(sk.Value, sk.Value)
	return sk
}

// GenSecretKeySparse generates a new SecretKey with exactly hw non-zero coefficients.
func (keygen *keyGenerator) GenSecretKeySparse(hw int) (sk *SecretKey) {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ternarySamplerMontgomery := ring.NewTernarySamplerSparse(prng, keygen.ringQP, hw, true)

	sk = new(SecretKey)
	sk.Value = ternarySamplerMontgomery.ReadNew()
	keygen.ringQP.NTT(sk.Value, sk.Value)
	return sk
}

// GenPublicKey generates a new public key from the provided SecretKey.
func (keygen *keyGenerator) GenPublicKey(sk *SecretKey) (pk *PublicKey) {

	pk = new(PublicKey)

	ringQP := keygen.ringQP

	//pk[0] = [-as + e]
	//pk[1] = [a]

	pk.Value[0] = keygen.gaussianSampler.ReadNew()
	ringQP.NTT(pk.Value[0], pk.Value[0])
	pk.Value[1] = keygen.uniformSampler.ReadNew()

	ringQP.MulCoeffsMontgomeryAndSub(sk.Value, pk.Value[1], pk.Value[0])

	return pk
}

// GenKeyPair generates a new SecretKey with distribution [1/3, 1/3, 1/3] and a corresponding public key.
func (keygen *keyGenerator) GenKeyPair() (sk *SecretKey, pk *PublicKey) {
	sk = keygen.GenSecretKey()
	return sk, keygen.GenPublicKey(sk)
}

// GenKeyPairSparse generates a new SecretKey with exactly hw non zero coefficients [1/2, 0, 1/2].
func (keygen *keyGenerator) GenKeyPairSparse(hw int) (sk *SecretKey, pk *PublicKey) {
	sk = keygen.GenSecretKeySparse(hw)
	return sk, keygen.GenPublicKey(sk)
}

// GenRelinKey generates a new EvaluationKey that will be used to relinearize Ciphertexts during multiplication.
func (keygen *keyGenerator) GenRelinearizationKey(sk *SecretKey, maxDegree int) (evk *RelinearizationKey) {

	if keygen.ringQP == nil {
		panic("modulus P is empty")
	}

	evk = new(RelinearizationKey)
	evk.Keys = make([]*SwitchingKey, maxDegree)
	for i := range evk.Keys {
		evk.Keys[i] = NewSwitchingKey(keygen.params)
	}

	keygen.polypool[0].CopyValues(sk.Value) // TODO Remove ?

	ringQP := keygen.ringQP

	keygen.polypool[1].CopyValues(sk.Value)
	for i := 0; i < maxDegree; i++ {
		ringQP.MulCoeffsMontgomery(keygen.polypool[1], sk.Value, keygen.polypool[1])
		keygen.newSwitchingKey(keygen.polypool[1], sk.Value, evk.Keys[i])
	}

	keygen.polypool[0].Zero()
	keygen.polypool[1].Zero()

	return
}

// GenRotationKeys generates a RotationKeySet from a list of galois element corresponding to the desired rotations
// See also GenRotationKeysForRotations.
func (keygen *keyGenerator) GenRotationKeys(galEls []uint64, sk *SecretKey) (rks *RotationKeySet) {
	rks = NewRotationKeySet(keygen.params, galEls)
	for _, galEl := range galEls {
		keygen.genrotKey(sk.Value, keygen.params.InverseGaloisElement(galEl), rks.Keys[galEl])
	}
	return rks
}

func (keygen *keyGenerator) GenSwitchingKeyForRotationBy(k int, sk *SecretKey) (swk *SwitchingKey) {
	swk = NewSwitchingKey(keygen.params)
	galElInv := keygen.params.GaloisElementForColumnRotationBy(-int(k))
	keygen.genrotKey(sk.Value, galElInv, swk)
	return
}

// GenRotationKeysForRotations generates a RotationKeySet supporting left rotations by k positions for all k in ks.
// Negative k is equivalent to a right rotation by k positions
// If includeConjugate is true, the resulting set contains the conjugation key.
func (keygen *keyGenerator) GenRotationKeysForRotations(ks []int, includeConjugate bool, sk *SecretKey) (rks *RotationKeySet) {
	galEls := make([]uint64, len(ks), len(ks)+1)
	for i, k := range ks {
		galEls[i] = keygen.params.GaloisElementForColumnRotationBy(k)
	}
	if includeConjugate {
		galEls = append(galEls, keygen.params.GaloisElementForRowRotation())
	}
	return keygen.GenRotationKeys(galEls, sk)
}

// GenSwitchingKey generates a new key-switching key, that will re-encrypt a Ciphertext encrypted under the input key into the output key.
func (keygen *keyGenerator) GenSwitchingKey(skInput, skOutput *SecretKey) (newevakey *SwitchingKey) {

	if keygen.params.PCount() == 0 {
		panic("Cannot GenSwitchingKey: modulus P is empty")
	}

	ring.CopyValues(skInput.Value, keygen.polypool[0])
	newevakey = NewSwitchingKey(keygen.params)
	keygen.newSwitchingKey(keygen.polypool[0], skOutput.Value, newevakey)
	keygen.polypool[0].Zero()
	return
}

func (keygen *keyGenerator) GenSwitchingKeyForRowRotation(sk *SecretKey) (swk *SwitchingKey) {
	swk = NewSwitchingKey(keygen.params)
	keygen.genrotKey(sk.Value, keygen.params.GaloisElementForRowRotation(), swk)
	return
}

func (keygen *keyGenerator) GenSwitchingKeyForGalois(galoisEl uint64, sk *SecretKey) (swk *SwitchingKey) {
	swk = NewSwitchingKey(keygen.params)
	keygen.genrotKey(sk.Value, keygen.params.InverseGaloisElement(galoisEl), swk)
	return
}

// GenRotationKeysForInnerSum generates a RotationKeySet supporting the InnerSum operation of the Evaluator
func (keygen *keyGenerator) GenRotationKeysForInnerSum(sk *SecretKey) (rks *RotationKeySet) {
	return keygen.GenRotationKeys(keygen.params.GaloisElementsForRowInnerSum(), sk)
}

func (keygen *keyGenerator) genrotKey(sk *ring.Poly, galEl uint64, swk *SwitchingKey) {

	skIn := sk
	skOut := keygen.polypool[1]

	index := ring.PermuteNTTIndex(galEl, uint64(keygen.ringQP.N))
	ring.PermuteNTTWithIndexLvl(keygen.params.QPCount()-1, skIn, index, skOut)

	keygen.newSwitchingKey(skIn, skOut, swk)

	keygen.polypool[0].Zero()
	keygen.polypool[1].Zero()
}

func (keygen *keyGenerator) newSwitchingKey(skIn, skOut *ring.Poly, swk *SwitchingKey) {

	ringQP := keygen.ringQP

	// Computes P * skIn
	ringQP.MulScalarBigint(skIn, keygen.pBigInt, keygen.polypool[0])

	alpha := keygen.params.PCount()
	beta := keygen.params.Beta()

	var index int
	for i := 0; i < beta; i++ {

		// e

		keygen.gaussianSampler.Read(swk.Value[i][0])
		ringQP.NTTLazy(swk.Value[i][0], swk.Value[i][0])
		ringQP.MForm(swk.Value[i][0], swk.Value[i][0])

		// a (since a is uniform, we consider we already sample it in the NTT and Montgomery domain)
		keygen.uniformSampler.Read(swk.Value[i][1])

		// e + (skIn * P) * (q_star * q_tild) mod QP
		//
		// q_prod = prod(q[i*alpha+j])
		// q_star = Q/qprod
		// q_tild = q_star^-1 mod q_prod
		//
		// Therefore : (skIn * P) * (q_star * q_tild) = sk*P mod q[i*alpha+j], else 0
		for j := 0; j < alpha; j++ {

			index = i*alpha + j

			qi := ringQP.Modulus[index]
			p0tmp := keygen.polypool[0].Coeffs[index]
			p1tmp := swk.Value[i][0].Coeffs[index]

			for w := 0; w < ringQP.N; w++ {
				p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
			}

			// It handles the case where nb pj does not divide nb qi
			if index >= keygen.params.QCount() {
				break
			}
		}

		// (skIn * P) * (q_star * q_tild) - a * skOut + e mod QP
		ringQP.MulCoeffsMontgomeryAndSub(swk.Value[i][1], skOut, swk.Value[i][0])
	}
}
