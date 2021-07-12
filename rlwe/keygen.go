package rlwe

import (
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
	params           Parameters
	ringQ            *ring.Ring
	ringP            *ring.Ring
	poolQ            [2]*ring.Poly
	poolP            [2]*ring.Poly
	gaussianSamplerQ *ring.GaussianSampler
	uniformSamplerQ  *ring.UniformSampler
	uniformSamplerP  *ring.UniformSampler
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as the evaluation,
// rotation and switching keys can be generated.
func NewKeyGenerator(params Parameters) KeyGenerator {

	ringQ := params.RingQ()
	ringP := params.RingP()

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	return &keyGenerator{
		params:           params,
		ringQ:            ringQ,
		ringP:            ringP,
		poolQ:            [2]*ring.Poly{ringQ.NewPoly(), ringQ.NewPoly()},
		poolP:            [2]*ring.Poly{ringP.NewPoly(), ringP.NewPoly()},
		gaussianSamplerQ: ring.NewGaussianSampler(prng, ringQ, params.Sigma(), int(6*params.Sigma())),
		uniformSamplerQ:  ring.NewUniformSampler(prng, ringQ),
		uniformSamplerP:  ring.NewUniformSampler(prng, ringP),
	}
}

// GenSecretKey generates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *keyGenerator) GenSecretKey() (sk *SecretKey) {
	return keygen.GenSecretKeyWithDistrib(1.0 / 3)
}

func (keygen *keyGenerator) GenSecretKeyGaussian() (sk *SecretKey) {
	sk = new(SecretKey)
	sk.Value[0] = keygen.gaussianSamplerQ.ReadNew()
	ExtendBasisSmallNormAndCenter(keygen.ringQ, keygen.ringP, sk.Value[0], sk.Value[1])
	keygen.ringQ.NTT(sk.Value[0], sk.Value[0])
	keygen.ringP.NTT(sk.Value[1], sk.Value[1])
	return sk
}

// GenSecretKeyWithDistrib generates a new SecretKey with the distribution [(p-1)/2, p, (p-1)/2].
func (keygen *keyGenerator) GenSecretKeyWithDistrib(p float64) (sk *SecretKey) {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ternarySamplerMontgomery := ring.NewTernarySampler(prng, keygen.ringQ, p, false)

	sk = new(SecretKey)
	sk.Value[0] = ternarySamplerMontgomery.ReadNew()
	sk.Value[1] = keygen.ringP.NewPoly()
	ExtendBasisSmallNormAndCenter(keygen.ringQ, keygen.ringP, sk.Value[0], sk.Value[1])
	keygen.ringQ.MForm(sk.Value[0], sk.Value[0])
	keygen.ringP.MForm(sk.Value[1], sk.Value[1])
	keygen.ringQ.NTT(sk.Value[0], sk.Value[0])
	keygen.ringP.NTT(sk.Value[1], sk.Value[1])
	return sk
}

// GenSecretKeySparse generates a new SecretKey with exactly hw non-zero coefficients.
func (keygen *keyGenerator) GenSecretKeySparse(hw int) (sk *SecretKey) {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ternarySamplerMontgomery := ring.NewTernarySamplerSparse(prng, keygen.ringQ, hw, false)

	sk = new(SecretKey)
	sk.Value[0] = ternarySamplerMontgomery.ReadNew()
	sk.Value[1] = keygen.ringP.NewPoly()
	ExtendBasisSmallNormAndCenter(keygen.ringQ, keygen.ringP, sk.Value[0], sk.Value[1])
	keygen.ringQ.MForm(sk.Value[0], sk.Value[0])
	keygen.ringP.MForm(sk.Value[1], sk.Value[1])
	keygen.ringQ.NTT(sk.Value[0], sk.Value[0])
	keygen.ringP.NTT(sk.Value[1], sk.Value[1])
	return sk
}

// GenPublicKey generates a new public key from the provided SecretKey.
func (keygen *keyGenerator) GenPublicKey(sk *SecretKey) (pk *PublicKey) {

	pk = new(PublicKey)

	ringQ := keygen.ringQ
	ringP := keygen.ringP

	//pk[0][0] = [-(a*s + e)] mod Q
	//pk[0][1] = [-(a*s + e)] mod P
	//pk[1][0] = [a] mod Q
	//pk[1][1] = [a] mod P

	pk.Value[0][0] = keygen.gaussianSamplerQ.ReadNew() // mod Q
	pk.Value[0][1] = ringP.NewPoly()
	ExtendBasisSmallNormAndCenter(ringQ, ringP, pk.Value[0][0], pk.Value[0][1])  // mod P
	ringQ.NTT(pk.Value[0][0], pk.Value[0][0])                                    // mod Q
	ringP.NTT(pk.Value[0][1], pk.Value[0][1])                                    // mod P
	pk.Value[1][0] = keygen.uniformSamplerQ.ReadNew()                            // mod Q
	pk.Value[1][1] = keygen.uniformSamplerP.ReadNew()                            // mod P
	ringQ.MulCoeffsMontgomeryAndSub(sk.Value[0], pk.Value[1][0], pk.Value[0][0]) // mod Q
	ringP.MulCoeffsMontgomeryAndSub(sk.Value[1], pk.Value[1][1], pk.Value[0][1]) // mod P
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

	if keygen.ringP == nil {
		panic("modulus P is empty")
	}

	evk = new(RelinearizationKey)
	evk.Keys = make([]*SwitchingKey, maxDegree)
	for i := range evk.Keys {
		evk.Keys[i] = NewSwitchingKey(keygen.params)
	}

	keygen.poolQ[1].CopyValues(sk.Value[0])
	keygen.poolP[1].CopyValues(sk.Value[1])
	for i := 0; i < maxDegree; i++ {
		keygen.ringQ.MulCoeffsMontgomery(keygen.poolQ[1], sk.Value[0], keygen.poolQ[1])
		keygen.ringP.MulCoeffsMontgomery(keygen.poolP[1], sk.Value[1], keygen.poolP[1])
		keygen.newSwitchingKey([2]*ring.Poly{keygen.poolQ[1], keygen.poolP[1]}, sk.Value, evk.Keys[i])
	}

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
	newevakey = NewSwitchingKey(keygen.params)
	keygen.newSwitchingKey(skInput.Value, skOutput.Value, newevakey)
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

func (keygen *keyGenerator) genrotKey(skIn [2]*ring.Poly, galEl uint64, swk *SwitchingKey) {
	index := ring.PermuteNTTIndex(galEl, uint64(keygen.ringQ.N))
	ring.PermuteNTTWithIndexLvl(keygen.params.QCount()-1, skIn[0], index, keygen.poolQ[1])
	ring.PermuteNTTWithIndexLvl(keygen.params.PCount()-1, skIn[1], index, keygen.poolP[1])
	keygen.newSwitchingKey(skIn, [2]*ring.Poly{keygen.poolQ[1], keygen.poolP[1]}, swk)
}

func (keygen *keyGenerator) newSwitchingKey(skIn, skOut [2]*ring.Poly, swk *SwitchingKey) {

	ringQ := keygen.ringQ
	ringP := keygen.ringP

	// Computes P * skIn
	ringQ.MulScalarBigint(skIn[0], ringP.ModulusBigint, keygen.poolQ[0])

	alpha := keygen.params.PCount()
	beta := keygen.params.Beta()

	var index int
	for i := 0; i < beta; i++ {

		// e

		keygen.gaussianSamplerQ.Read(swk.Value[i][0][0])
		ExtendBasisSmallNormAndCenter(ringQ, ringP, swk.Value[i][0][0], swk.Value[i][0][1])
		ringQ.NTTLazy(swk.Value[i][0][0], swk.Value[i][0][0])
		ringP.NTTLazy(swk.Value[i][0][1], swk.Value[i][0][1])
		ringQ.MForm(swk.Value[i][0][0], swk.Value[i][0][0])
		ringP.MForm(swk.Value[i][0][1], swk.Value[i][0][1])

		// a (since a is uniform, we consider we already sample it in the NTT and Montgomery domain)
		keygen.uniformSamplerQ.Read(swk.Value[i][1][0])
		keygen.uniformSamplerP.Read(swk.Value[i][1][1])

		// e + (skIn * P) * (q_star * q_tild) mod QP
		//
		// q_prod = prod(q[i*alpha+j])
		// q_star = Q/qprod
		// q_tild = q_star^-1 mod q_prod
		//
		// Therefore : (skIn * P) * (q_star * q_tild) = sk*P mod q[i*alpha+j], else 0
		for j := 0; j < alpha; j++ {

			index = i*alpha + j

			// It handles the case where nb pj does not divide nb qi
			if index >= keygen.params.QCount() {
				break
			}

			qi := ringQ.Modulus[index]
			p0tmp := keygen.poolQ[0].Coeffs[index]
			p1tmp := swk.Value[i][0][0].Coeffs[index]

			for w := 0; w < ringQ.N; w++ {
				p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
			}
		}

		// (skIn * P) * (q_star * q_tild) - a * skOut + e mod QP
		ringQ.MulCoeffsMontgomeryAndSub(swk.Value[i][1][0], skOut[0], swk.Value[i][0][0])
		ringP.MulCoeffsMontgomeryAndSub(swk.Value[i][1][1], skOut[1], swk.Value[i][0][1])
	}
}

/*
// GenSwitchingKeyDimensionFrom generates a new key-switching key, that will re-encrypt a ciphertext encrypted
// under skIn of dimension n to a ciphertext encrypted under sKOut of dimension N > n.
// [-a*SkOut + w*P*skIn_{Y^{N/n}} + e, a] in X^{N}
func (keygen *keyGenerator) GenSwitchingKeyDimensionFrom(paramsFrom *Parameters, skIn, skOut *ring.Poly) (swk *SwitchingKey) {

	// From a small N to a larger N
	if !(paramsFrom.N() < keygen.params.N()) {
		panic("paramsFrom dimension must be smaller than the keygenerator parameters dimension")
	}

	if keygen.params.PCount() == 0 {
		panic("Cannot GenSwitchingKey: modulus P is empty")
	}

	ringQP := keygen.params.ringQP

	swk = NewSwitchingKey(keygen.params)

	// Maps skIn Y = X^{N/n} -> X
	keygen.polypool[0].Zero()
	paramsFrom.RingQP().InvNTT(skOut, keygen.polypool[0])
	gap := paramsFrom.N() / keygen.params.N()
	for j := 0; j < len(ringQP.Modulus); j++ {
		tmp := keygen.polypool[0].Coeffs[j]
		for i := paramsFrom.RingQP().N - 1; i >= 0; i-- {
			tmp[i*gap], tmp[i] = tmp[i], tmp[i*gap]
		}
	}
	ringQP.NTT(keygen.polypool[0], keygen.polypool[0])

	ringQP.MulScalarBigint(keygen.polypool[0], paramsFrom.RingP().ModulusBigint, keygen.polypool[0])

	alpha := keygen.params.PCount()
	beta := keygen.params.Beta()

	var index int
	for i := 0; i < beta; i++ {

		keygen.gaussianSampler.Read(swk.Value[i][0])
		ringQP.NTTLazy(swk.Value[i][0], swk.Value[i][0])
		ringQP.MForm(swk.Value[i][0], swk.Value[i][0])

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

		keygen.uniformSampler.Read(swk.Value[i][1])
		ringQP.MulCoeffsMontgomeryAndSub(swk.Value[i][1], skOut, swk.Value[i][0])
	}
	return
}

// GenSwitchingKeyDimensionFrom generates a new key-switching key, that will re-encrypt a ciphertext encrypted
// under skIn of dimension N to a ciphertext encrypted under sKOut of dimension n < N.
// [-a*skOut_{Y^{N/n}} + w*P*skIn + e_{N}, a_{N}] in X^{N}
// The ciphetext moduli of paramsTo must be shared by the key generator.
func (keygen *keyGenerator) NewSwitchingKeyDimensionTo(paramsTo *Parameters, skIn, skOut *ring.Poly) (swk *SwitchingKey) {

	if !(paramsTo.N() < keygen.params.N()) {
		panic("paramsTo dimension must be smaller than the keygenerator parameters dimension")
	}

	if keygen.params.PCount() == 0 {
		panic("Cannot GenSwitchingKey: modulus P is empty")
	}

	swk = NewSwitchingKey(keygen.params)

	ringQP := paramsTo.RingQP() // uses paramsTo ringQP to ensure the output security

	// Concatenates modulus Q of small params with modulus P of large params
	skInSmallQ := new(ring.Poly)
	skInSmallQ.Coeffs = append(skIn.Coeffs[:paramsTo.MaxLevel()+1], skIn.Coeffs[keygen.params.MaxLevel()+1:]...)

	// Computes P * skIn
	ringQP.MulScalarBigint(skIn, keygen.pBigInt, keygen.polypool[0])

	// Maps skOut Y = X^{N/n} -> X
	keygen.polypool[1].Zero()
	ringQP.InvNTT(skOut, keygen.polypool[1])
	gap := paramsTo.N() / keygen.params.N()
	for j := 0; j < len(ringQP.Modulus); j++ {
		tmp := keygen.polypool[1].Coeffs[j]
		for i := paramsTo.RingQP().N - 1; i >= 0; i-- {
			tmp[i*gap], tmp[i] = tmp[i], tmp[i*gap]
		}
	}
	ringQP.NTT(keygen.polypool[1], keygen.polypool[1])

	var index int
	alpha := keygen.params.PCount()
	beta := keygen.params.Beta()
	for i := 0; i < beta; i++ {

		keygen.gaussianSampler.Read(swk.Value[i][0])

		ringQP.NTTLazy(swk.Value[i][0], swk.Value[i][0])
		ringQP.MForm(swk.Value[i][0], swk.Value[i][0])

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

		keygen.uniformSampler.Read(swk.Value[i][1])
		ringQP.MulCoeffsMontgomeryAndSub(swk.Value[i][1], keygen.polypool[1], swk.Value[i][0])
	}
	return
}
*/
