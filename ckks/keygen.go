package ckks

import (
	"math"
	"math/big"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
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
	GenBootstrappingKey(logSlots uint64, btpParams *BootstrappingParameters, sk *SecretKey) (btpKey *BootstrappingKey)
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
	Conjugate
)

// RotationKeys is a structure that stores the switching-keys required during the homomorphic rotations.
type RotationKeys struct {
	permuteNTTLeftIndex      map[uint64][]uint64
	permuteNTTRightIndex     map[uint64][]uint64
	permuteNTTConjugateIndex []uint64

	evakeyRotColLeft  map[uint64]*SwitchingKey
	evakeyRotColRight map[uint64]*SwitchingKey
	evakeyConjugate   *SwitchingKey
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

// BootstrappingKey is a structure that stores the switching-keys required during the bootstrapping.
type BootstrappingKey struct {
	relinkey *EvaluationKey // Relinearization key
	rotkeys  *RotationKeys  // Rotation and conjugation keys
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as the evaluation,
// rotation and switching keys can be generated.
func NewKeyGenerator(params *Parameters) KeyGenerator {

	var qp *ring.Ring
	var err error
	if qp, err = ring.NewRing(params.N(), append(params.qi, params.pi...)); err != nil {
		panic(err)
	}

	var pBigInt *big.Int
	if len(params.pi) != 0 {
		pBigInt = ring.NewUint(1)
		for _, pi := range params.pi {
			pBigInt.Mul(pBigInt, ring.NewUint(pi))
		}
	}

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	return &keyGenerator{
		params:          params.Copy(),
		ringQP:          qp,
		pBigInt:         pBigInt,
		polypool:        [2]*ring.Poly{qp.NewPoly(), qp.NewPoly()},
		gaussianSampler: ring.NewGaussianSampler(prng),
		uniformSampler:  ring.NewUniformSampler(prng, qp),
	}
}

// GenSecretKey generates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *keyGenerator) GenSecretKey() (sk *SecretKey) {
	return keygen.GenSecretKeyWithDistrib(1.0 / 3)
}

func (keygen *keyGenerator) GenSecretKeyGaussian() (sk *SecretKey) {
	sk = new(SecretKey)

	sk.sk = keygen.gaussianSampler.ReadNew(keygen.ringQP, keygen.params.sigma, uint64(6*keygen.params.sigma))
	keygen.ringQP.NTT(sk.sk, sk.sk)
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
	sk.sk = ternarySamplerMontgomery.ReadNew()
	keygen.ringQP.NTT(sk.sk, sk.sk)
	return sk
}

// GenSecretKeySparse generates a new SecretKey with exactly hw non-zero coefficients.
func (keygen *keyGenerator) GenSecretKeySparse(hw uint64) (sk *SecretKey) {
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	ternarySamplerMontgomery := ring.NewTernarySamplerSparse(prng, keygen.ringQP, hw, true)

	sk = new(SecretKey)
	sk.sk = ternarySamplerMontgomery.ReadNew()
	keygen.ringQP.NTT(sk.sk, sk.sk)
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
	pk.pk[0] = keygen.gaussianSampler.ReadNew(keygen.ringQP, keygen.params.sigma, uint64(6*keygen.params.sigma))
	ringQP.NTT(pk.pk[0], pk.pk[0])
	pk.pk[1] = keygen.uniformSampler.ReadNew()

	ringQP.MulCoeffsMontgomeryAndSub(sk.sk, pk.pk[1], pk.pk[0])

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

	if rotType != Conjugate {

		if rotKey.permuteNTTLeftIndex == nil {
			rotKey.permuteNTTLeftIndex = make(map[uint64][]uint64)
		}

		if rotKey.permuteNTTRightIndex == nil {
			rotKey.permuteNTTRightIndex = make(map[uint64][]uint64)
		}

		if _, inMap := rotKey.permuteNTTLeftIndex[k]; !inMap {
			rotKey.permuteNTTLeftIndex[k] = ring.PermuteNTTIndex(GaloisGen, k, ringQP.N)
		}

		if _, inMap := rotKey.permuteNTTRightIndex[k]; !inMap {
			rotKey.permuteNTTRightIndex[k] = ring.PermuteNTTIndex(GaloisGen, 2*ringQP.N-k, ringQP.N)
		}
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

	case Conjugate:
		rotKey.permuteNTTConjugateIndex = ring.PermuteNTTIndex(2*ringQP.N-1, 1, ringQP.N)
		rotKey.evakeyConjugate = keygen.genrotKey(sk.Get(), rotKey.permuteNTTConjugateIndex)
	}
}

// GenRotationKeysPow2 generates a new rotation key with all the power-of-two rotations to the left and right, as well as the conjugation.
func (keygen *keyGenerator) GenRotationKeysPow2(skOutput *SecretKey) (rotKey *RotationKeys) {

	if len(keygen.params.pi) == 0 {
		panic("Cannot GenRotationKeysPow2: modulus P is empty")
	}

	rotKey = NewRotationKeys()

	for n := uint64(1); n < keygen.params.N()/2; n <<= 1 {
		keygen.GenRotationKey(RotationLeft, skOutput, n, rotKey)
		keygen.GenRotationKey(RotationRight, skOutput, n, rotKey)
	}

	//keygen.GenRot(Conjugate, skOutput, 0, rotKey)

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

			rotKey.permuteNTTLeftIndex[k] = ring.PermuteNTTIndex(GaloisGen, k, params.N())

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

			rotKey.permuteNTTRightIndex[k] = ring.PermuteNTTIndex(GaloisGen, 2*params.N()-1-k, params.N())

			rotKey.evakeyRotColRight[k] = new(SwitchingKey)
			rotKey.evakeyRotColRight[k].evakey = make([][2]*ring.Poly, len(evakey))
			for j := range evakey {
				rotKey.evakeyRotColRight[k].evakey[j][0] = evakey[j][0].CopyNew()
				rotKey.evakeyRotColRight[k].evakey[j][1] = evakey[j][1].CopyNew()
			}
		}

	case Conjugate:

		if rotKey.evakeyConjugate == nil {

			rotKey.permuteNTTConjugateIndex = ring.PermuteNTTIndex(2*params.N()-1, 1, params.N())

			rotKey.evakeyConjugate = new(SwitchingKey)
			rotKey.evakeyConjugate.evakey = make([][2]*ring.Poly, len(evakey))
			for j := range evakey {
				rotKey.evakeyConjugate.evakey[j][0] = evakey[j][0].CopyNew()
				rotKey.evakeyConjugate.evakey[j][1] = evakey[j][1].CopyNew()
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
		switchingkey.evakey[i][0] = keygen.gaussianSampler.ReadNew(keygen.ringQP, keygen.params.sigma, uint64(6*keygen.params.sigma))
		ringQP.NTTLazy(switchingkey.evakey[i][0], switchingkey.evakey[i][0])
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

// GenKeys generates the bootstrapping keys
func (keygen *keyGenerator) GenBootstrappingKey(logSlots uint64, btpParams *BootstrappingParameters, sk *SecretKey) (btpKey *BootstrappingKey) {

	btpKey = &BootstrappingKey{
		relinkey: keygen.GenRelinKey(sk),
		rotkeys:  NewRotationKeys(),
	}

	rotKeyIndex := computeBootstrappingDFTRotationList(keygen.params.logN, logSlots, btpParams)

	/*
		nbKeys := uint64(len(rotKeyIndex)) + 2 //rot keys + conj key + relin key
		nbPoly := keygen.params.Beta()
		nbCoefficients := 2 * keygen.params.N() * keygen.params.QPiCount()
		bytesPerCoeff := uint64(8)
		log.Println("Switching-Keys size (GB) :", float64(nbKeys*nbPoly*nbCoefficients*bytesPerCoeff)/float64(1000000000), "(", nbKeys, "keys)")
	*/

	keygen.GenRotationKey(Conjugate, sk, 0, btpKey.rotkeys)
	for _, i := range rotKeyIndex {
		keygen.GenRotationKey(RotationLeft, sk, uint64(i), btpKey.rotkeys)
	}

	return
}

func addMatrixRotToList(pVec map[uint64]bool, rotations []uint64, N1, slots uint64, repack bool) []uint64 {

	var index uint64
	for j := range pVec {

		index = (j / N1) * N1

		if repack {
			// Sparse repacking, occurring during the first DFT matrix of the CoeffsToSlots.
			index &= (2*slots - 1)
		} else {
			// Other cases
			index &= (slots - 1)
		}

		if index != 0 && !utils.IsInSliceUint64(index, rotations) {
			rotations = append(rotations, index)
		}

		index = j & (N1 - 1)

		if index != 0 && !utils.IsInSliceUint64(index, rotations) {
			rotations = append(rotations, index)
		}
	}

	return rotations
}

func computeBootstrappingDFTRotationList(logN, logSlots uint64, btpParams *BootstrappingParameters) (rotKeyIndex []uint64) {

	// List of the rotation key values to needed for the bootstrapp
	rotKeyIndex = []uint64{}

	var slots uint64 = 1 << logSlots
	var dslots uint64 = slots
	if logSlots < logN-1 {
		dslots <<= 1
	}

	//SubSum rotation needed X -> Y^slots rotations
	for i := logSlots; i < logN-1; i++ {
		if !utils.IsInSliceUint64(1<<i, rotKeyIndex) {
			rotKeyIndex = append(rotKeyIndex, 1<<i)
		}
	}

	indexCtS := computeBootstrappingDFTIndexMap(logN, logSlots, btpParams.CtSDepth(), true)

	// Coeffs to Slots rotations
	for _, pVec := range indexCtS {

		N1 := findbestbabygiantstepsplitIndexMap(pVec, dslots, btpParams.MaxN1N2Ratio)

		rotKeyIndex = addMatrixRotToList(pVec, rotKeyIndex, N1, slots, false)
	}

	indexStC := computeBootstrappingDFTIndexMap(logN, logSlots, btpParams.StCDepth(), false)

	// Slots to Coeffs rotations
	for i, pVec := range indexStC {

		N1 := findbestbabygiantstepsplitIndexMap(pVec, dslots, btpParams.MaxN1N2Ratio)

		if logSlots < logN-1 && i == 0 {
			rotKeyIndex = addMatrixRotToList(pVec, rotKeyIndex, N1, slots, true)
		} else {
			rotKeyIndex = addMatrixRotToList(pVec, rotKeyIndex, N1, slots, false)
		}
	}

	return rotKeyIndex
}

func computeBootstrappingDFTIndexMap(logN, logSlots, maxDepth uint64, forward bool) (rotationMap []map[uint64]bool) {

	var level, depth, nextLevel uint64

	level = logSlots

	rotationMap = make([]map[uint64]bool, maxDepth)

	// We compute the chain of merge in order or reverse order depending if its DFT or InvDFT because
	// the way the levels are collapsed has an inpact on the total number of rotations and keys to be
	// stored. Ex. instead of using 255 + 64 plaintext vectors, we can use 127 + 128 plaintext vectors
	// by reversing the order of the merging.
	merge := make([]uint64, maxDepth)
	for i := uint64(0); i < maxDepth; i++ {

		depth = uint64(math.Ceil(float64(level) / float64(maxDepth-i)))

		if forward {
			merge[i] = depth
		} else {
			merge[uint64(len(merge))-i-1] = depth

		}

		level -= depth
	}

	level = logSlots
	for i := uint64(0); i < maxDepth; i++ {

		if logSlots < logN-1 && !forward && i == 0 {

			// Special initial matrix for the repacking before SlotsToCoeffs
			rotationMap[i] = genWfftRepackIndexMap(logSlots, level)

			// Merges this special initial matrix with the first layer of SlotsToCoeffs DFT
			rotationMap[i] = nextLevelfftIndexMap(rotationMap[i], logSlots, 2<<logSlots, level, forward)

			// Continues the merging with the next layers if the total depth requires it.
			nextLevel = level - 1
			for j := uint64(0); j < merge[i]-1; j++ {
				rotationMap[i] = nextLevelfftIndexMap(rotationMap[i], logSlots, 2<<logSlots, nextLevel, forward)
				nextLevel--
			}

		} else {
			// First layer of the i-th level of the DFT
			rotationMap[i] = genWfftIndexMap(logSlots, level, forward)

			// Merges the layer with the next levels of the DFT if the total depth requires it.
			nextLevel = level - 1
			for j := uint64(0); j < merge[i]-1; j++ {
				rotationMap[i] = nextLevelfftIndexMap(rotationMap[i], logSlots, 1<<logSlots, nextLevel, forward)
				nextLevel--
			}
		}

		level -= merge[i]
	}

	return
}

// Finds the best N1*N2 = N for the baby-step giant-step algorithm for matrix multiplication.
func findbestbabygiantstepsplitIndexMap(vector map[uint64]bool, maxN uint64, maxRatio float64) (minN uint64) {

	for N1 := uint64(1); N1 < maxN; N1 <<= 1 {

		index := make(map[uint64][]uint64)

		for key := range vector {

			idx1 := key / N1
			idx2 := key & (N1 - 1)

			if index[idx1] == nil {
				index[idx1] = []uint64{idx2}
			} else {
				index[idx1] = append(index[idx1], idx2)
			}
		}

		if len(index[0]) > 0 {

			hoisted := len(index[0]) - 1
			normal := len(index) - 1

			// The matrice is very sparse already
			if normal == 0 {
				return N1 / 2
			}

			if hoisted > normal {
				// Finds the next split that has a ratio hoisted/normal greater or equal to maxRatio
				for float64(hoisted)/float64(normal) < maxRatio {

					if normal/2 == 0 {
						break
					}
					N1 *= 2
					hoisted = hoisted*2 + 1
					normal = normal / 2
				}
				return N1
			}
		}
	}

	return 1
}

func genWfftIndexMap(logL, level uint64, forward bool) (vectors map[uint64]bool) {

	var rot uint64

	if forward {
		rot = 1 << (level - 1)
	} else {
		rot = 1 << (logL - level)
	}

	vectors = make(map[uint64]bool)
	vectors[0] = true
	vectors[rot] = true
	vectors[(1<<logL)-rot] = true
	return
}

func genWfftRepackIndexMap(logL, level uint64) (vectors map[uint64]bool) {
	vectors = make(map[uint64]bool)
	vectors[0] = true
	vectors[(1 << logL)] = true
	return
}

func nextLevelfftIndexMap(vec map[uint64]bool, logL, N, nextLevel uint64, forward bool) (newVec map[uint64]bool) {

	var rot uint64

	newVec = make(map[uint64]bool)

	if forward {
		rot = (1 << (nextLevel - 1)) & (N - 1)
	} else {
		rot = (1 << (logL - nextLevel)) & (N - 1)
	}

	for i := range vec {
		newVec[i] = true
		newVec[(i+rot)&(N-1)] = true
		newVec[(i+N-rot)&(N-1)] = true
	}

	return
}
