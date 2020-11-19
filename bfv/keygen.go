package bfv

import (
	"math/big"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// KeyGenerator is an interface implementing the methods of the keyGenerator.
type KeyGenerator interface {
	GenSecretKey() (sk *SecretKey)
	GenSecretkeyWithDistrib(p float64) (sk *SecretKey)
	GenPublicKey(sk *SecretKey) (pk *PublicKey)
	GenKeyPair() (sk *SecretKey, pk *PublicKey)
	GenRelinKey(sk *SecretKey, maxDegree uint64) (evk *EvaluationKey)
	GenSwitchingKey(skIn, skOut *SecretKey) (evk *SwitchingKey)
	GenRot(rotType Rotation, sk *SecretKey, k uint64, rotKey *RotationKeys)
	GenRotationKeysPow2(sk *SecretKey) (rotKey *RotationKeys)
}

// keyGenerator is a structure that stores the elements required to create new keys,
// as well as a small memory pool for intermediate values.
type keyGenerator struct {
	params           *Parameters
	ringQP           *ring.Ring
	pBigInt          *big.Int
	polypool         [2]*ring.Poly
	gaussianSampler  *ring.GaussianSampler
	uniformSampler   *ring.UniformSampler
	galElRotRow      uint64   // Rows rotation generator
	galElRotColLeft  []uint64 // Columns right rotations generators
	galElRotColRight []uint64 // Columsn left rotations generators
}

// SecretKey is a structure that stores the SecretKey.
type SecretKey struct {
	sk *ring.Poly
}

// PublicKey is a structure that stores the PublicKey.
type PublicKey struct {
	pk [2]*ring.Poly
}

// Rotation is a type used to represent the rotations types.
type Rotation int

// Constants for rotation types
const (
	RotationRight = iota + 1
	RotationLeft
	RotationRow
)

// RotationKeys is a structure that stores the switching-keys required during the homomorphic rotations.
type RotationKeys struct {
	permuteNTTLeftIndex  map[uint64][]uint64
	permuteNTTRightIndex map[uint64][]uint64
	permuteNTTRowIndex   []uint64

	evakeyRotColLeft  map[uint64]*SwitchingKey
	evakeyRotColRight map[uint64]*SwitchingKey
	evakeyRotRow      *SwitchingKey
}

// EvaluationKey is a structure that stores the switching-keys required during the relinearization.
type EvaluationKey struct {
	evakey []*SwitchingKey
}

// SwitchingKey is a structure that stores the switching-keys required during the key-switching.
type SwitchingKey struct {
	evakey [][2]*ring.Poly
}

// Get returns the switching key backing slice.
func (swk *SwitchingKey) Get() [][2]*ring.Poly {
	return swk.evakey
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as the evaluation,
// rotation and switching keys can be generated.
func NewKeyGenerator(params *Parameters) KeyGenerator {

	var ringQP *ring.Ring
	var err error
	if ringQP, err = ring.NewRing(params.N(), append(params.qi, params.pi...)); err != nil {
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
		params:           params.Copy(),
		ringQP:           ringQP,
		pBigInt:          pBigInt,
		polypool:         [2]*ring.Poly{ringQP.NewPoly(), ringQP.NewPoly()},
		gaussianSampler:  ring.NewGaussianSampler(prng, ringQP, params.Sigma(), uint64(6*params.Sigma())),
		uniformSampler:   ring.NewUniformSampler(prng, ringQP),
		galElRotColLeft:  ring.GenGaloisParams(params.N(), GaloisGen),
		galElRotColRight: ring.GenGaloisParams(params.N(), ring.ModExp(GaloisGen, 2*params.N()-1, 2*params.N())),
		galElRotRow:      2*params.N() - 1,
	}
}

// GenSecretKey creates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *keyGenerator) GenSecretKey() (sk *SecretKey) {
	return keygen.GenSecretkeyWithDistrib(1.0 / 3)
}

// GenSecretkeyWithDistrib creates a new SecretKey with the distribution [(1-p)/2, p, (1-p)/2].
func (keygen *keyGenerator) GenSecretkeyWithDistrib(p float64) (sk *SecretKey) {
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

// NewSecretKey generates a new SecretKey with zero values.
func NewSecretKey(params *Parameters) *SecretKey {

	sk := new(SecretKey)
	sk.sk = ring.NewPoly(uint64(1<<params.logN), uint64(len(params.qi)+len(params.pi)))
	return sk
}

// Get returns the polynomial of the target SecretKey.
func (sk *SecretKey) Get() *ring.Poly {
	return sk.sk
}

// Set sets the polynomial of the target secret key as the input polynomial.
func (sk *SecretKey) Set(poly *ring.Poly) {
	sk.sk = poly.CopyNew()
}

// GenPublicKey generates a new PublicKey from the provided SecretKey.
func (keygen *keyGenerator) GenPublicKey(sk *SecretKey) (pk *PublicKey) {

	pk = new(PublicKey)

	ringQP := keygen.ringQP

	//pk[0] = [-(a*s + e)]
	//pk[1] = [a]

	pk.pk[0] = keygen.gaussianSampler.ReadNew()
	ringQP.NTT(pk.pk[0], pk.pk[0])
	pk.pk[1] = keygen.uniformSampler.ReadNew()

	ringQP.MulCoeffsMontgomeryAndAdd(sk.sk, pk.pk[1], pk.pk[0])
	ringQP.Neg(pk.pk[0], pk.pk[0])

	return pk
}

// NewPublicKey returns a new PublicKey with zero values.
func NewPublicKey(params *Parameters) (pk *PublicKey) {

	pk = new(PublicKey)

	pk.pk[0] = ring.NewPoly(uint64(1<<params.logN), uint64(len(params.qi)+len(params.pi)))
	pk.pk[1] = ring.NewPoly(uint64(1<<params.logN), uint64(len(params.qi)+len(params.pi)))

	return
}

// Get returns the polynomials of the PublicKey.
func (pk *PublicKey) Get() [2]*ring.Poly {
	return pk.pk
}

// Set sets the polynomial of the PublicKey as the input polynomials.
func (pk *PublicKey) Set(p [2]*ring.Poly) {
	pk.pk[0] = p[0].CopyNew()
	pk.pk[1] = p[1].CopyNew()
}

// NewKeyPair generates a new SecretKey with distribution [1/3, 1/3, 1/3] and a corresponding PublicKey.
func (keygen *keyGenerator) GenKeyPair() (sk *SecretKey, pk *PublicKey) {
	sk = keygen.GenSecretKey()
	return sk, keygen.GenPublicKey(sk)
}

// NewRelinKey generates a new evaluation key from the provided SecretKey. It will be used to relinearize a ciphertext (encrypted under a PublicKey generated from the provided SecretKey)
// of degree > 1 to a ciphertext of degree 1. Max degree is the maximum degree of the ciphertext allowed to relinearize.
func (keygen *keyGenerator) GenRelinKey(sk *SecretKey, maxDegree uint64) (evk *EvaluationKey) {

	if keygen.ringQP == nil {
		panic("Cannot GenRelinKey: modulus P is empty")
	}

	evk = new(EvaluationKey)

	evk.evakey = make([]*SwitchingKey, maxDegree)

	keygen.polypool[0].Copy(sk.Get())

	ringQP := keygen.ringQP

	ringQP.MulCoeffsMontgomery(sk.Get(), sk.Get(), keygen.polypool[1])
	evk.evakey[0] = keygen.newSwitchingKey(keygen.polypool[1], sk.Get())

	for i := uint64(1); i < maxDegree; i++ {
		ringQP.MulCoeffsMontgomery(keygen.polypool[1], sk.Get(), keygen.polypool[1])
		evk.evakey[i] = keygen.newSwitchingKey(keygen.polypool[0], sk.Get())
	}

	keygen.polypool[0].Zero()
	keygen.polypool[1].Zero()

	return
}

// NewRelinKey creates a new EvaluationKey with zero values.
func NewRelinKey(params *Parameters, maxDegree uint64) (evakey *EvaluationKey) {

	evakey = new(EvaluationKey)

	beta := params.Beta()

	evakey.evakey = make([]*SwitchingKey, maxDegree)

	for w := uint64(0); w < maxDegree; w++ {

		evakey.evakey[w] = new(SwitchingKey)

		evakey.evakey[w].evakey = make([][2]*ring.Poly, beta)

		for i := uint64(0); i < beta; i++ {

			evakey.evakey[w].evakey[i][0] = ring.NewPoly(uint64(1<<params.logN), uint64(len(params.qi)+len(params.pi)))
			evakey.evakey[w].evakey[i][1] = ring.NewPoly(uint64(1<<params.logN), uint64(len(params.qi)+len(params.pi)))
		}
	}

	return
}

// Get returns the slice of SwitchingKeys of the target EvaluationKey.
func (evk *EvaluationKey) Get() []*SwitchingKey {
	return evk.evakey
}

// SetRelinKeys sets the polynomial of the target EvaluationKey as the input polynomials.
func (evk *EvaluationKey) SetRelinKeys(rlk [][][2]*ring.Poly) {

	evk.evakey = make([]*SwitchingKey, len(rlk))
	for i := range rlk {
		evk.evakey[i] = new(SwitchingKey)
		evk.evakey[i].evakey = make([][2]*ring.Poly, len(rlk[i]))
		for j := range rlk[i] {
			evk.evakey[i].evakey[j][0] = rlk[i][j][0].CopyNew()
			evk.evakey[i].evakey[j][1] = rlk[i][j][1].CopyNew()
		}
	}
}

// GenSwitchingKey generates a new key-switching key, that will allow to re-encrypt under the output-key a ciphertext encrypted under the input-key.
func (keygen *keyGenerator) GenSwitchingKey(skInput, skOutput *SecretKey) (newevakey *SwitchingKey) {

	if keygen.ringQP == nil {
		panic("Cannot GenRelinKey: modulus P is empty")
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
		evakey.evakey[i][0] = ring.NewPoly(uint64(1<<params.logN), uint64(len(params.qi)+len(params.pi)))
		evakey.evakey[i][1] = ring.NewPoly(uint64(1<<params.logN), uint64(len(params.qi)+len(params.pi)))
	}

	return
}

// NewRotationKeys returns a new empty RotationKeys struct.
func NewRotationKeys() (rotKey *RotationKeys) {
	rotKey = new(RotationKeys)
	return
}

// GenRot populates the target RotationKeys with a SwitchingKey for the desired rotation type and amount.
func (keygen *keyGenerator) GenRot(rotType Rotation, sk *SecretKey, k uint64, rotKey *RotationKeys) {

	ringQP := keygen.ringQP
	if ringQP == nil {
		panic("Cannot GenRot: modulus P is empty")
	}

	switch rotType {
	case RotationLeft:

		if rotKey.evakeyRotColLeft == nil {
			rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
		}

		if rotKey.evakeyRotColLeft[k] == nil && k != 0 {
			rotKey.evakeyRotColLeft[k] = keygen.genrotKey(sk.Get(), GaloisGen, ringQP.N-k)
		}

	case RotationRight:

		if rotKey.evakeyRotColRight == nil {
			rotKey.evakeyRotColRight = make(map[uint64]*SwitchingKey)
		}

		if rotKey.evakeyRotColRight[k] == nil && k != 0 {
			rotKey.evakeyRotColRight[k] = keygen.genrotKey(sk.Get(), GaloisGen, k)
		}

	case RotationRow:
		rotKey.evakeyRotRow = keygen.genrotKey(sk.Get(), ringQP.NthRoot-1, 1)
	}
}

// GenRotationKeysPow2 generates a new rotation key with all the power-of-two rotations to the left and right, as well as the conjugation.
func (keygen *keyGenerator) GenRotationKeysPow2(skOutput *SecretKey) (rotKey *RotationKeys) {

	if keygen.ringQP == nil {
		panic("Cannot GenRotationKeysPow2: modulus P is empty")
	}

	rotKey = NewRotationKeys()

	for n := uint64(1); n < 1<<(keygen.params.LogN()-1); n <<= 1 {
		keygen.GenRot(RotationLeft, skOutput, n, rotKey)
		keygen.GenRot(RotationRight, skOutput, n, rotKey)
	}

	keygen.GenRot(RotationRow, skOutput, 0, rotKey)

	return
}

// SetRotKey sets the target RotationKeys' SwitchingKey for the specified rotation type and amount with the input polynomials.
func (rotKey *RotationKeys) SetRotKey(rotType Rotation, k uint64, evakey [][2]*ring.Poly) {

	switch rotType {
	case RotationLeft:

		if rotKey.evakeyRotColLeft == nil {
			rotKey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
		}

		if rotKey.evakeyRotColLeft[k] == nil && k != 0 {

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

		if rotKey.evakeyRotColRight[k] == nil && k != 0 {

			rotKey.evakeyRotColRight[k] = new(SwitchingKey)
			rotKey.evakeyRotColRight[k].evakey = make([][2]*ring.Poly, len(evakey))
			for j := range evakey {
				rotKey.evakeyRotColRight[k].evakey[j][0] = evakey[j][0].CopyNew()
				rotKey.evakeyRotColRight[k].evakey[j][1] = evakey[j][1].CopyNew()
			}
		}

	case RotationRow:

		if rotKey.evakeyRotRow == nil {

			rotKey.evakeyRotRow = new(SwitchingKey)
			rotKey.evakeyRotRow.evakey = make([][2]*ring.Poly, len(evakey))
			for j := range evakey {
				rotKey.evakeyRotRow.evakey[j][0] = evakey[j][0].CopyNew()
				rotKey.evakeyRotRow.evakey[j][1] = evakey[j][1].CopyNew()
			}
		}
	}
}

func (keygen *keyGenerator) genrotKey(sk *ring.Poly, gen, k uint64) (switchingkey *SwitchingKey) {

	skIn := sk
	skOut := keygen.polypool[1]

	ring.PermuteNTT(skIn, gen, k, keygen.ringQP.N, keygen.ringQP.NthRoot, skOut)

	switchingkey = keygen.newSwitchingKey(skIn, skOut)

	keygen.polypool[0].Zero()
	keygen.polypool[1].Zero()

	return
}

func (keygen *keyGenerator) newSwitchingKey(skIn, skOut *ring.Poly) (switchingkey *SwitchingKey) {

	switchingkey = new(SwitchingKey)

	ringQP := keygen.ringQP

	alpha := keygen.params.Alpha()
	beta := keygen.params.Beta()

	var index uint64

	// delta_sk = skIn - skOut = GaloisEnd(skOut, rotation) - skOut

	ringQP.MulScalarBigint(skIn, keygen.pBigInt, keygen.polypool[0])

	switchingkey.evakey = make([][2]*ring.Poly, beta)

	for i := uint64(0); i < beta; i++ {

		// e
		switchingkey.evakey[i][0] = keygen.gaussianSampler.ReadNew()
		ringQP.NTTLazy(switchingkey.evakey[i][0], switchingkey.evakey[i][0])
		ringQP.MForm(switchingkey.evakey[i][0], switchingkey.evakey[i][0])
		// a
		switchingkey.evakey[i][1] = keygen.uniformSampler.ReadNew()

		// e + skIn * (qiBarre*qiStar) * 2^w
		// (qiBarre*qiStar)%qi = 1, else 0

		for j := uint64(0); j < alpha; j++ {

			index = i*alpha + j

			qi := ringQP.Modulus[index]
			p0tmp := keygen.polypool[0].Coeffs[index]
			p1tmp := switchingkey.evakey[i][0].Coeffs[index]

			for w := uint64(0); w < ringQP.N; w++ {
				p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
			}

			// Handles the case where nb pj does not divide nb qi
			if index >= keygen.params.QiCount() {
				break
			}

		}

		// skIn * (qiBarre*qiStar) * 2^w - a*sk + e
		ringQP.MulCoeffsMontgomeryAndSub(switchingkey.evakey[i][1], skOut, switchingkey.evakey[i][0])
	}

	return
}
