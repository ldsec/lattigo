package rckks

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
	"math/big"
)

// KeyGenerator is an interface implementing the methods of the KeyGenerator.
type KeyGenerator interface {
	GenSecretKey() (sk *SecretKey)
	GenSecretKeyGaussian() (sk *SecretKey)
	GenSecretKeyWithDistrib(p float64) (sk *SecretKey)
	GenSecretKeySparse(hw uint64) (sk *SecretKey)
	GenPublicKey(sk *SecretKey) (pk *PublicKey)
	GenKeyPair() (sk *SecretKey, pk *PublicKey)
	GenKeyPairSparse(hw uint64) (sk *SecretKey, pk *PublicKey)
	GenRelinKey(sk *SecretKey) (evakey *EvaluationKey)
	GenSwitchingKey(skInput, skOutput *SecretKey) (newevakey *SwitchingKey)
	GenRotationKey(rotType Rotation, sk *SecretKey, k uint64, rotKey *RotationKeys)
	GenRotationKeysPow2(skOutput *SecretKey) (rotKey *RotationKeys)
}

// KeyGenerator is a structure that stores the elements required to create new keys,
// as well as a small memory pool for intermediate values.
type keyGenerator struct {
	params          *Parameters
	ringQP          *ring.Ring
	pBigInt         *big.Int
	polypool        [2]*ring.Poly
	gaussianSampler *ring.GaussianSampler
	uniformSampler  *ring.UniformSampler
}

// SecretKey is a structure that stores the SecretKey
type SecretKey struct {
	sk *ring.Poly
}

// PublicKey is a structure that stores the PublicKey
type PublicKey struct {
	pk [2]*ring.Poly
}

// Rotation is a type used to represent the rotations types.
type Rotation int

// Constants for rotation types
const (
	RotationRight = iota + 1
	RotationLeft
)

// RotationKeys is a structure that stores the switching-keys required during the homomorphic rotations.
type RotationKeys struct {
	permuteNTTLeftIndex  map[uint64][]uint64
	permuteNTTRightIndex map[uint64][]uint64

	evakeyRotColLeft  map[uint64]*SwitchingKey
	evakeyRotColRight map[uint64]*SwitchingKey
}

// EvaluationKey is a structure that stores the switching-keys required during the relinearization.
type EvaluationKey struct {
	evakey *SwitchingKey
}

// SwitchingKey is a structure that stores the switching-keys required during the key-switching.
type SwitchingKey struct {
	evakey [][2]*ring.Poly
}

// Get returns the switching key backing slice
func (swk *SwitchingKey) Get() [][2]*ring.Poly {
	return swk.evakey
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as the evaluation,
// rotation and switching keys can be generated.
func NewKeyGenerator(params *Parameters) KeyGenerator {

	var qp *ring.Ring
	var err error
	if qp, err = ring.NewRingWithNthRoot(params.N(), params.N()<<2, append(params.qi, params.pi...)); err != nil {
		panic(err)
	}

	var pBigInt *big.Int
	if len(params.pi) != 0 {
		pBigInt = ring.NewUint(1)
		for _, pi := range params.pi {
			pBigInt.Mul(pBigInt, ring.NewUint(pi))
		}
	}

	prng, err := utils.NewKeyedPRNG(nil)
	if err != nil {
		panic(err)
	}

	return &keyGenerator{
		params:          params.Copy(),
		ringQP:          qp,
		pBigInt:         pBigInt,
		polypool:        [2]*ring.Poly{qp.NewPoly(), qp.NewPoly()},
		gaussianSampler: ring.NewGaussianSampler(prng, qp, params.sigma, uint64(6*params.sigma)),
		uniformSampler:  ring.NewUniformSampler(prng, qp),
	}
}

// GenSecretKey generates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *keyGenerator) GenSecretKey() (sk *SecretKey) {
	return keygen.GenSecretKeyWithDistrib(1.0 / 3)
}

func (keygen *keyGenerator) GenSecretKeyGaussian() (sk *SecretKey) {
	sk = new(SecretKey)
	sk.sk = keygen.gaussianSampler.ReadNew()
	NTTRCKKS(keygen.ringQP, sk.sk, sk.sk)
	return sk
}

// GenSecretKeyWithDistrib generates a new SecretKey with the distribution [(p-1)/2, p, (p-1)/2].
func (keygen *keyGenerator) GenSecretKeyWithDistrib(p float64) (sk *SecretKey) {
	prng, err := utils.NewKeyedPRNG(nil)
	if err != nil {
		panic(err)
	}
	ternarySamplerMontgomery := ring.NewTernarySampler(prng, keygen.ringQP, p, true)

	sk = new(SecretKey)
	sk.sk = ternarySamplerMontgomery.ReadNew()
	NTTRCKKS(keygen.ringQP, sk.sk, sk.sk)
	return sk
}

// GenSecretKeySparse generates a new SecretKey with exactly hw non-zero coefficients.
func (keygen *keyGenerator) GenSecretKeySparse(hw uint64) (sk *SecretKey) {
	prng, err := utils.NewKeyedPRNG(nil)
	if err != nil {
		panic(err)
	}
	ternarySamplerMontgomery := ring.NewTernarySamplerSparse(prng, keygen.ringQP, hw, true)

	sk = new(SecretKey)
	sk.sk = ternarySamplerMontgomery.ReadNew()
	NTTRCKKS(keygen.ringQP, sk.sk, sk.sk)
	return sk
}

// NewSecretKey generates a new SecretKey with zero values.
func NewSecretKey(params *Parameters) *SecretKey {

	sk := new(SecretKey)
	sk.sk = params.NewPolyQP()
	return sk
}

// Get returns the value of the SecretKey.
func (sk *SecretKey) Get() *ring.Poly {
	return sk.sk
}

// Set sets the value of the SecretKey to the provided value.
func (sk *SecretKey) Set(poly *ring.Poly) {
	sk.sk = poly.CopyNew()
}

// GenPublicKey generates a new public key from the provided SecretKey.
func (keygen *keyGenerator) GenPublicKey(sk *SecretKey) (pk *PublicKey) {

	pk = new(PublicKey)

	ringQP := keygen.ringQP

	//pk[0] = [-(a*s + e)]
	//pk[1] = [a]
	pk.pk[0] = keygen.gaussianSampler.ReadNew()
	NTTRCKKS(ringQP, pk.pk[0], pk.pk[0])
	pk.pk[1] = keygen.uniformSampler.ReadNew()
	ringQP.MulCoeffsMontgomeryAndAdd(sk.sk, pk.pk[1], pk.pk[0])
	ringQP.Neg(pk.pk[0], pk.pk[0])

	return pk
}

// NewPublicKey returns a new PublicKey with zero values.
func NewPublicKey(params *Parameters) (pk *PublicKey) {

	pk = new(PublicKey)

	pk.pk[0] = params.NewPolyQP()
	pk.pk[1] = params.NewPolyQP()

	return
}

// Get returns the value of the the public key.
func (pk *PublicKey) Get() [2]*ring.Poly {
	return pk.pk
}

// Set sets the value of the public key to the provided value.
func (pk *PublicKey) Set(poly [2]*ring.Poly) {
	pk.pk[0] = poly[0].CopyNew()
	pk.pk[1] = poly[1].CopyNew()
}

// GenKeyPair generates a new SecretKey with distribution [1/3, 1/3, 1/3] and a corresponding public key.
func (keygen *keyGenerator) GenKeyPair() (sk *SecretKey, pk *PublicKey) {
	sk = keygen.GenSecretKey()
	return sk, keygen.GenPublicKey(sk)
}

// GenKeyPairSparse generates a new SecretKey with exactly hw non zero coefficients [1/2, 0, 1/2].
func (keygen *keyGenerator) GenKeyPairSparse(hw uint64) (sk *SecretKey, pk *PublicKey) {
	sk = keygen.GenSecretKeySparse(hw)
	return sk, keygen.GenPublicKey(sk)
}

// GenRelinKey generates a new EvaluationKey that will be used to relinearize Ciphertexts during multiplication.
func (keygen *keyGenerator) GenRelinKey(sk *SecretKey) (evakey *EvaluationKey) {

	if len(keygen.params.pi) == 0 {
		panic("Cannot GenRelinKey: modulus P is empty")
	}

	evakey = new(EvaluationKey)
	keygen.ringQP.MulCoeffsMontgomery(sk.Get(), sk.Get(), keygen.polypool[0])
	evakey.evakey = keygen.newSwitchingKey(keygen.polypool[0], sk.Get())
	keygen.polypool[0].Zero()

	return
}

// NewRelinKey returns a new EvaluationKey with zero values.
func NewRelinKey(params *Parameters) (evakey *EvaluationKey) {

	evakey = new(EvaluationKey)
	evakey.evakey = new(SwitchingKey)

	// delta_sk = skInput - skOutput = GaloisEnd(skOutput, rotation) - skOutput
	evakey.evakey.evakey = make([][2]*ring.Poly, params.Beta())
	for i := uint64(0); i < params.Beta(); i++ {

		evakey.evakey.evakey[i][0] = params.NewPolyQP()
		evakey.evakey.evakey[i][1] = params.NewPolyQP()
	}

	return
}

// Get returns the slice of switching keys of the evaluation-key.
func (evk *EvaluationKey) Get() *SwitchingKey {
	return evk.evakey
}

// Set sets the target Evaluation key with the input polynomials.
func (evk *EvaluationKey) Set(rlk [][2]*ring.Poly) {

	evk.evakey = new(SwitchingKey)
	evk.evakey.evakey = make([][2]*ring.Poly, len(rlk))
	for j := range rlk {
		evk.evakey.evakey[j][0] = rlk[j][0].CopyNew()
		evk.evakey.evakey[j][1] = rlk[j][1].CopyNew()
	}
}

// GenSwitchingKey generates a new key-switching key, that will re-encrypt a Ciphertext encrypted under the input key into the output key.
func (keygen *keyGenerator) GenSwitchingKey(skInput, skOutput *SecretKey) (newevakey *SwitchingKey) {

	if len(keygen.params.pi) == 0 {
		panic("Cannot GenSwitchingKey: modulus P is empty")
	}

	keygen.ringQP.Copy(skInput.Get(), keygen.polypool[0])
	newevakey = keygen.newSwitchingKey(keygen.polypool[0], skOutput.Get())
	keygen.polypool[0].Zero()
	return
}

// NewSwitchingKey returns a new SwitchingKey with zero values.
func NewSwitchingKey(params *Parameters) (evakey *SwitchingKey) {

	evakey = new(SwitchingKey)

	// delta_sk = skInput - skOutput = GaloisEnd(skOutput, rotation) - skOutput

	evakey.evakey = make([][2]*ring.Poly, params.Beta())

	for i := uint64(0); i < params.Beta(); i++ {
		evakey.evakey[i][0] = params.NewPolyQP()
		evakey.evakey[i][1] = params.NewPolyQP()
	}

	return
}

// NewRotationKeys generates a new instance of RotationKeys, with the provided rotation to the left, right and conjugation if requested.
func NewRotationKeys() (rotKey *RotationKeys) {
	rotKey = new(RotationKeys)
	return
}

// GenRot populates the input RotationKeys with a SwitchingKey for the given rotation type and amount.
func (keygen *keyGenerator) GenRotationKey(rotType Rotation, sk *SecretKey, k uint64, rotKey *RotationKeys) {

	if len(keygen.params.pi) == 0 {
		panic("Cannot GenRot: modulus P is empty")
	}

	ringQP := keygen.ringQP
	N := ringQP.N
	NthRoot := ringQP.NthRoot

	if rotKey.permuteNTTLeftIndex == nil {
		rotKey.permuteNTTLeftIndex = make(map[uint64][]uint64)
	}

	if rotKey.permuteNTTRightIndex == nil {
		rotKey.permuteNTTRightIndex = make(map[uint64][]uint64)
	}

	if _, inMap := rotKey.permuteNTTLeftIndex[k]; !inMap {
		rotKey.permuteNTTLeftIndex[k] = ring.PermuteNTTIndex(GaloisGen, k, N, NthRoot)
	}

	if _, inMap := rotKey.permuteNTTRightIndex[k]; !inMap {
		rotKey.permuteNTTRightIndex[k] = ring.PermuteNTTIndex(GaloisGen, NthRoot-k, N, NthRoot)
	}

	switch rotType {
	case RotationLeft:

		if rotKey.evakeyRotColLeft == nil {
			rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
		}

		if rotKey.evakeyRotColLeft[k] == nil && k != 0 {
			rotKey.evakeyRotColLeft[k] = keygen.genrotKey(sk.Get(), rotKey.permuteNTTRightIndex[k])
		}

	case RotationRight:

		if rotKey.evakeyRotColRight == nil {
			rotKey.evakeyRotColRight = make(map[uint64]*SwitchingKey)
		}

		if rotKey.evakeyRotColRight[k] == nil && k != 0 {
			rotKey.evakeyRotColRight[k] = keygen.genrotKey(sk.Get(), rotKey.permuteNTTLeftIndex[k])
		}
	}
}

// GenRotationKeysPow2 generates a new rotation key with all the power-of-two rotations to the left and right, as well as the conjugation.
func (keygen *keyGenerator) GenRotationKeysPow2(skOutput *SecretKey) (rotKey *RotationKeys) {

	if len(keygen.params.pi) == 0 {
		panic("Cannot GenRotationKeysPow2: modulus P is empty")
	}

	rotKey = NewRotationKeys()

	for n := uint64(1); n < keygen.params.N(); n <<= 1 {
		keygen.GenRotationKey(RotationLeft, skOutput, n, rotKey)
		keygen.GenRotationKey(RotationRight, skOutput, n, rotKey)
	}

	return
}

// SetRotKey sets the target RotationKeys' SwitchingKey for the specified rotation type and amount with the input polynomials.
func (rotKey *RotationKeys) SetRotKey(params *Parameters, evakey [][2]*ring.Poly, rotType Rotation, k uint64) {

	switch rotType {
	case RotationLeft:

		if rotKey.evakeyRotColLeft == nil {
			rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
		}

		if rotKey.permuteNTTLeftIndex == nil {
			rotKey.permuteNTTLeftIndex = make(map[uint64][]uint64)
		}

		if rotKey.evakeyRotColLeft[k] == nil && k != 0 {

			rotKey.permuteNTTLeftIndex[k] = ring.PermuteNTTIndex(GaloisGen, k, params.N(), 4*params.N())

			rotKey.evakeyRotColLeft[k] = new(SwitchingKey)
			rotKey.evakeyRotColLeft[k].evakey = make([][2]*ring.Poly, len(evakey))
			for j := range evakey {
				rotKey.evakeyRotColLeft[k].evakey[j][0] = evakey[j][0].CopyNew()
				rotKey.evakeyRotColLeft[k].evakey[j][1] = evakey[j][1].CopyNew()
			}
		}

	case RotationRight:

		if rotKey.evakeyRotColRight == nil {
			rotKey.evakeyRotColRight = make(map[uint64]*SwitchingKey)
		}

		if rotKey.permuteNTTRightIndex == nil {
			rotKey.permuteNTTRightIndex = make(map[uint64][]uint64)
		}

		if rotKey.evakeyRotColRight[k] == nil && k != 0 {

			rotKey.permuteNTTRightIndex[k] = ring.PermuteNTTIndex(GaloisGen, 4*params.N()-k, params.N(), 4*params.N())

			rotKey.evakeyRotColRight[k] = new(SwitchingKey)
			rotKey.evakeyRotColRight[k].evakey = make([][2]*ring.Poly, len(evakey))
			for j := range evakey {
				rotKey.evakeyRotColRight[k].evakey[j][0] = evakey[j][0].CopyNew()
				rotKey.evakeyRotColRight[k].evakey[j][1] = evakey[j][1].CopyNew()
			}
		}
	}
}

func (keygen *keyGenerator) genrotKey(sk *ring.Poly, index []uint64) (switchingkey *SwitchingKey) {

	skIn := sk
	skOut := keygen.polypool[1]

	ring.PermuteNTTWithIndexLvl(keygen.params.QPiCount()-1, skIn, index, skOut)

	switchingkey = keygen.newSwitchingKey(skIn, skOut)

	keygen.polypool[0].Zero()
	keygen.polypool[1].Zero()

	return
}

func (keygen *keyGenerator) newSwitchingKey(skIn, skOut *ring.Poly) (switchingkey *SwitchingKey) {

	switchingkey = new(SwitchingKey)

	ringQP := keygen.ringQP

	// Computes P * skIn
	ringQP.MulScalarBigint(skIn, keygen.pBigInt, keygen.polypool[0])

	alpha := keygen.params.Alpha()
	beta := keygen.params.Beta()

	var index uint64

	switchingkey.evakey = make([][2]*ring.Poly, beta)

	for i := uint64(0); i < beta; i++ {

		// e
		switchingkey.evakey[i][0] = keygen.gaussianSampler.ReadNew()
		NTTRCKKS(ringQP, switchingkey.evakey[i][0], switchingkey.evakey[i][0])
		ringQP.MForm(switchingkey.evakey[i][0], switchingkey.evakey[i][0])

		// a (since a is uniform, we consider we already sample it in the NTT and Montgomery domain)
		switchingkey.evakey[i][1] = keygen.uniformSampler.ReadNew()

		// e + (skIn * P) * (q_star * q_tild) mod QP
		//
		// q_prod = prod(q[i*alpha+j])
		// q_star = Q/qprod
		// q_tild = q_star^-1 mod q_prod
		//
		// Therefore : (skIn * P) * (q_star * q_tild) = sk*P mod q[i*alpha+j], else 0
		for j := uint64(0); j < alpha; j++ {

			index = i*alpha + j

			qi := ringQP.Modulus[index]
			p0tmp := keygen.polypool[0].Coeffs[index]
			p1tmp := switchingkey.evakey[i][0].Coeffs[index]

			for w := uint64(0); w < ringQP.N; w++ {
				p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
			}

			// It handles the case where nb pj does not divide nb qi
			if index >= keygen.params.QiCount() {
				break
			}
		}

		// (skIn * P) * (q_star * q_tild) - a * skOut + e mod QP
		ringQP.MulCoeffsMontgomeryAndSub(switchingkey.evakey[i][1], skOut, switchingkey.evakey[i][0])
	}

	return
}
