package dckks

import (
	"flag"
	"fmt"
	"math"
	"sort"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters). Overrides -short and requires -timeout=0.")
var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")
var minPrec float64 = 15.0
var parties uint64 = 3

func testString(opname string, parties uint64, params *ckks.Parameters) string {
	return fmt.Sprintf("%sparties=%d/logN=%d/logQ=%d/levels=%d/a=%d/b=%d",
		opname,
		parties,
		params.LogN(),
		params.LogQP(),
		params.MaxLevel()+1,
		params.Alpha(),
		params.Beta())
}

type testContext struct {
	params *ckks.Parameters

	dckksContext *dckksContext

	prng utils.PRNG

	encoder   ckks.Encoder
	evaluator ckks.Evaluator

	encryptorPk0 ckks.Encryptor
	decryptorSk0 ckks.Decryptor
	decryptorSk1 ckks.Decryptor

	pk0 *ckks.PublicKey
	pk1 *ckks.PublicKey

	sk0 *ckks.SecretKey
	sk1 *ckks.SecretKey

	sk0Shards []*ckks.SecretKey
	sk1Shards []*ckks.SecretKey
}

func TestDCKKS(t *testing.T) {

	var err error

	var defaultParams = ckks.DefaultParams[ckks.PN12QP109 : ckks.PN12QP109+4] // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15
	if testing.Short() {
		defaultParams = ckks.DefaultParams[ckks.PN12QP109 : ckks.PN12QP109+2] // the short test runs for ring degree N=2^12, 2^13
	}
	if *flagLongTest {
		defaultParams = ckks.DefaultParams // the long test suite runs for all default parameters
	}

	for _, p := range defaultParams {

		var testCtx *testContext
		if testCtx, err = genTestParams(p); err != nil {
			panic(err)
		}

		testPublicKeyGen(testCtx, t)
		testRelinKeyGen(testCtx, t)
		testKeyswitching(testCtx, t)
		testPublicKeySwitching(testCtx, t)
		testRotKeyGenConjugate(testCtx, t)
		testRotKeyGenCols(testCtx, t)
		testRefresh(testCtx, t)
		testRefreshAndPermute(testCtx, t)
	}
}

func genTestParams(defaultParams *ckks.Parameters) (testCtx *testContext, err error) {

	testCtx = new(testContext)

	testCtx.params = defaultParams.Copy()

	testCtx.dckksContext = newDckksContext(testCtx.params)

	if testCtx.prng, err = utils.NewPRNG(); err != nil {
		return nil, err
	}

	testCtx.encoder = ckks.NewEncoder(testCtx.params)
	testCtx.evaluator = ckks.NewEvaluator(testCtx.params)

	kgen := ckks.NewKeyGenerator(testCtx.params)

	// SecretKeys
	testCtx.sk0Shards = make([]*ckks.SecretKey, parties)
	testCtx.sk1Shards = make([]*ckks.SecretKey, parties)
	tmp0 := testCtx.dckksContext.ringQP.NewPoly()
	tmp1 := testCtx.dckksContext.ringQP.NewPoly()

	for j := uint64(0); j < parties; j++ {
		testCtx.sk0Shards[j] = kgen.GenSecretKey()
		testCtx.sk1Shards[j] = kgen.GenSecretKey()
		testCtx.dckksContext.ringQP.Add(tmp0, testCtx.sk0Shards[j].Get(), tmp0)
		testCtx.dckksContext.ringQP.Add(tmp1, testCtx.sk1Shards[j].Get(), tmp1)
	}

	testCtx.sk0 = new(ckks.SecretKey)
	testCtx.sk1 = new(ckks.SecretKey)

	testCtx.sk0.Set(tmp0)
	testCtx.sk1.Set(tmp1)

	// Publickeys
	testCtx.pk0 = kgen.GenPublicKey(testCtx.sk0)
	testCtx.pk1 = kgen.GenPublicKey(testCtx.sk1)

	testCtx.encryptorPk0 = ckks.NewEncryptorFromPk(testCtx.params, testCtx.pk0)
	testCtx.decryptorSk0 = ckks.NewDecryptor(testCtx.params, testCtx.sk0)
	testCtx.decryptorSk1 = ckks.NewDecryptor(testCtx.params, testCtx.sk1)

	return
}

func testPublicKeyGen(testCtx *testContext, t *testing.T) {

	decryptorSk0 := testCtx.decryptorSk0
	sk0Shards := testCtx.sk0Shards

	crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.dckksContext.ringQP)

	t.Run(testString("PublicKeyGen/", parties, testCtx.params), func(t *testing.T) {

		crp := crpGenerator.ReadNew()

		type Party struct {
			*CKGProtocol
			s  *ring.Poly
			s1 CKGShare
		}

		ckgParties := make([]*Party, parties)
		for i := uint64(0); i < parties; i++ {
			p := new(Party)
			p.CKGProtocol = NewCKGProtocol(testCtx.params)
			p.s = sk0Shards[i].Get()
			p.s1 = p.AllocateShares()
			ckgParties[i] = p
		}
		P0 := ckgParties[0]

		// Each party creates a new CKGProtocol instance
		for i, p := range ckgParties {
			p.GenShare(p.s, crp, p.s1)
			if i > 0 {
				P0.AggregateShares(p.s1, P0.s1, P0.s1)
			}
		}

		pk := &ckks.PublicKey{}
		P0.GenPublicKey(P0.s1, crp, pk)

		// Verifies that decrypt((encryptp(collectiveSk, m), collectivePk) = m
		encryptorTest := ckks.NewEncryptorFromPk(testCtx.params, pk)

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorTest, 1, t)

		verifyTestVectors(testCtx, decryptorSk0, coeffs, ciphertext, t)
	})

}

func testRelinKeyGen(testCtx *testContext, t *testing.T) {

	evaluator := testCtx.evaluator
	encryptorPk0 := testCtx.encryptorPk0
	decryptorSk0 := testCtx.decryptorSk0
	sk0Shards := testCtx.sk0Shards

	t.Run(testString("RelinKeyGen/", parties, testCtx.params), func(t *testing.T) {

		type Party struct {
			*RKGProtocol
			u      *ring.Poly
			s      *ring.Poly
			share1 RKGShare
			share2 RKGShare
		}

		rkgParties := make([]*Party, parties)

		for i := range rkgParties {
			p := new(Party)
			p.RKGProtocol = NewEkgProtocol(testCtx.params)
			p.u = p.NewEphemeralKey()
			p.s = sk0Shards[i].Get()
			p.share1, p.share2 = p.AllocateShares()
			rkgParties[i] = p
		}

		P0 := rkgParties[0]

		crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.dckksContext.ringQP)
		crp := make([]*ring.Poly, testCtx.params.Beta())

		for i := uint64(0); i < testCtx.params.Beta(); i++ {
			crp[i] = crpGenerator.ReadNew()
		}

		// ROUND 1
		for i, p := range rkgParties {
			p.GenShareRoundOne(p.u, p.s, crp, p.share1)
			if i > 0 {
				P0.AggregateShareRoundOne(p.share1, P0.share1, P0.share1)
			}
		}

		//ROUND 2
		for i, p := range rkgParties {
			p.GenShareRoundTwo(P0.share1, p.u, p.s, crp, p.share2)
			if i > 0 {
				P0.AggregateShareRoundTwo(p.share2, P0.share2, P0.share2)
			}
		}

		evk := ckks.NewRelinKey(testCtx.params)
		P0.GenRelinearizationKey(P0.share1, P0.share2, evk)

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, 1, t)

		for i := range coeffs {
			coeffs[i] *= coeffs[i]
		}

		evaluator.MulRelin(ciphertext, ciphertext, evk, ciphertext)

		evaluator.Rescale(ciphertext, testCtx.params.Scale(), ciphertext)

		require.Equal(t, ciphertext.Degree(), uint64(1))

		verifyTestVectors(testCtx, decryptorSk0, coeffs, ciphertext, t)

	})

}

func testKeyswitching(testCtx *testContext, t *testing.T) {

	encryptorPk0 := testCtx.encryptorPk0
	decryptorSk1 := testCtx.decryptorSk1
	sk0Shards := testCtx.sk0Shards
	sk1Shards := testCtx.sk1Shards

	t.Run(testString("Keyswitching/", parties, testCtx.params), func(t *testing.T) {

		type Party struct {
			*CKSProtocol
			s0    *ring.Poly
			s1    *ring.Poly
			share CKSShare
		}

		cksParties := make([]*Party, parties)
		for i := uint64(0); i < parties; i++ {
			p := new(Party)
			p.CKSProtocol = NewCKSProtocol(testCtx.params, 6.36)
			p.s0 = sk0Shards[i].Get()
			p.s1 = sk1Shards[i].Get()
			p.share = p.AllocateShare()
			cksParties[i] = p
		}
		P0 := cksParties[0]

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, 1, t)

		// Each party creates its CKSProtocol instance with tmp = si-si'
		for i, p := range cksParties {
			p.GenShare(p.s0, p.s1, ciphertext, p.share)
			if i > 0 {
				P0.AggregateShares(p.share, P0.share, P0.share)
			}
		}

		ksCiphertext := ckks.NewCiphertext(testCtx.params, 1, ciphertext.Level(), ciphertext.Scale())

		P0.KeySwitch(P0.share, ciphertext, ksCiphertext)

		verifyTestVectors(testCtx, decryptorSk1, coeffs, ksCiphertext, t)

		P0.KeySwitch(P0.share, ciphertext, ciphertext)

		verifyTestVectors(testCtx, decryptorSk1, coeffs, ksCiphertext, t)

	})
}

func testPublicKeySwitching(testCtx *testContext, t *testing.T) {

	encryptorPk0 := testCtx.encryptorPk0
	decryptorSk1 := testCtx.decryptorSk1
	sk0Shards := testCtx.sk0Shards
	pk1 := testCtx.pk1

	t.Run(testString("PublicKeySwitching/", parties, testCtx.params), func(t *testing.T) {

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, 1, t)

		testCtx.evaluator.DropLevel(ciphertext, 1)

		type Party struct {
			*PCKSProtocol
			s     *ring.Poly
			share PCKSShare
		}

		pcksParties := make([]*Party, parties)
		for i := uint64(0); i < parties; i++ {
			p := new(Party)
			p.PCKSProtocol = NewPCKSProtocol(testCtx.params, 6.36)
			p.s = sk0Shards[i].Get()
			p.share = p.AllocateShares(ciphertext.Level())
			pcksParties[i] = p
		}
		P0 := pcksParties[0]

		ciphertextSwitched := ckks.NewCiphertext(testCtx.params, 1, ciphertext.Level(), ciphertext.Scale())

		for i, p := range pcksParties {
			p.GenShare(p.s, pk1, ciphertext, p.share)
			if i > 0 {
				P0.AggregateShares(p.share, P0.share, P0.share)
			}
		}

		P0.KeySwitch(P0.share, ciphertext, ciphertextSwitched)

		verifyTestVectors(testCtx, decryptorSk1, coeffs, ciphertextSwitched, t)
	})
}

func testRotKeyGenConjugate(testCtx *testContext, t *testing.T) {

	ringQP := testCtx.dckksContext.ringQP
	evaluator := testCtx.evaluator
	encryptorPk0 := testCtx.encryptorPk0
	decryptorSk0 := testCtx.decryptorSk0
	sk0Shards := testCtx.sk0Shards

	t.Run(testString("RotKeyGenConjugate/", parties, testCtx.params), func(t *testing.T) {

		type Party struct {
			*RTGProtocol
			s     *ring.Poly
			share RTGShare
		}

		pcksParties := make([]*Party, parties)
		for i := uint64(0); i < parties; i++ {
			p := new(Party)
			p.RTGProtocol = NewRotKGProtocol(testCtx.params)
			p.s = sk0Shards[i].Get()
			p.share = p.AllocateShare()
			pcksParties[i] = p
		}
		P0 := pcksParties[0]

		crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.dckksContext.ringQP)
		crp := make([]*ring.Poly, testCtx.params.Beta())

		for i := uint64(0); i < testCtx.params.Beta(); i++ {
			crp[i] = crpGenerator.ReadNew()
		}

		for i, p := range pcksParties {
			p.GenShare(ckks.Conjugate, 0, p.s, crp, &p.share)
			if i > 0 {
				P0.Aggregate(p.share, P0.share, P0.share)
			}
		}

		rotkey := ckks.NewRotationKeys()
		P0.Finalize(testCtx.params, P0.share, crp, rotkey)

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, 1, t)

		evaluator.Conjugate(ciphertext, rotkey, ciphertext)

		coeffsWant := make([]complex128, ringQP.N>>1)

		for i := uint64(0); i < ringQP.N>>1; i++ {
			coeffsWant[i] = complex(real(coeffs[i]), -imag(coeffs[i]))
		}

		verifyTestVectors(testCtx, decryptorSk0, coeffsWant, ciphertext, t)

	})
}

func testRotKeyGenCols(testCtx *testContext, t *testing.T) {

	ringQP := testCtx.dckksContext.ringQP
	evaluator := testCtx.evaluator
	encryptorPk0 := testCtx.encryptorPk0
	decryptorSk0 := testCtx.decryptorSk0
	sk0Shards := testCtx.sk0Shards

	t.Run(testString("RotKeyGenCols/", parties, testCtx.params), func(t *testing.T) {

		type Party struct {
			*RTGProtocol
			s     *ring.Poly
			share RTGShare
		}

		pcksParties := make([]*Party, parties)
		for i := uint64(0); i < parties; i++ {
			p := new(Party)
			p.RTGProtocol = NewRotKGProtocol(testCtx.params)
			p.s = sk0Shards[i].Get()
			p.share = p.AllocateShare()
			pcksParties[i] = p
		}

		P0 := pcksParties[0]

		crpGenerator := ring.NewUniformSampler(testCtx.prng, ringQP)
		crp := make([]*ring.Poly, testCtx.params.Beta())

		for i := uint64(0); i < testCtx.params.Beta(); i++ {
			crp[i] = crpGenerator.ReadNew()
		}

		mask := (ringQP.N >> 1) - 1

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, 1, t)

		receiver := ckks.NewCiphertext(testCtx.params, ciphertext.Degree(), ciphertext.Level(), ciphertext.Scale())

		for k := uint64(1); k < ringQP.N>>1; k <<= 1 {

			for i, p := range pcksParties {
				p.GenShare(ckks.RotationLeft, k, p.s, crp, &p.share)
				if i > 0 {
					P0.Aggregate(p.share, P0.share, P0.share)
				}
			}

			rotkey := ckks.NewRotationKeys()
			P0.Finalize(testCtx.params, P0.share, crp, rotkey)

			evaluator.Rotate(ciphertext, k, rotkey, receiver)

			coeffsWant := make([]complex128, ringQP.N>>1)

			for i := uint64(0); i < ringQP.N>>1; i++ {
				coeffsWant[i] = coeffs[(i+k)&mask]
			}

			verifyTestVectors(testCtx, decryptorSk0, coeffsWant, receiver, t)
		}
	})
}

func testRefresh(testCtx *testContext, t *testing.T) {

	evaluator := testCtx.evaluator
	encryptorPk0 := testCtx.encryptorPk0
	decryptorSk0 := testCtx.decryptorSk0
	sk0Shards := testCtx.sk0Shards

	levelStart := uint64(3)

	t.Run(testString("Refresh/", parties, testCtx.params), func(t *testing.T) {

		if testCtx.params.MaxLevel() < 3 {
			t.Skip("skipping test for params max level < 3")
		}

		type Party struct {
			*RefreshProtocol
			s      *ring.Poly
			share1 RefreshShareDecrypt
			share2 RefreshShareRecrypt
		}

		RefreshParties := make([]*Party, parties)
		for i := uint64(0); i < parties; i++ {
			p := new(Party)
			p.RefreshProtocol = NewRefreshProtocol(testCtx.params)
			p.s = sk0Shards[i].Get()
			p.share1, p.share2 = p.AllocateShares(levelStart)
			RefreshParties[i] = p
		}

		P0 := RefreshParties[0]

		crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.dckksContext.ringQ)
		crp := crpGenerator.ReadNew()

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, 1.0, t)

		for ciphertext.Level() != levelStart {
			evaluator.DropLevel(ciphertext, 1)
		}

		for i, p := range RefreshParties {
			p.GenShares(p.s, levelStart, parties, ciphertext, crp, p.share1, p.share2)
			if i > 0 {
				P0.Aggregate(p.share1, P0.share1, P0.share1)
				P0.Aggregate(p.share2, P0.share2, P0.share2)
			}
		}

		// We refresh the ciphertext with the simulated error
		P0.Decrypt(ciphertext, P0.share1)      // Masked decryption
		P0.Recode(ciphertext)                  // Masked re-encoding
		P0.Recrypt(ciphertext, crp, P0.share2) // Masked re-encryption

		require.Equal(t, ciphertext.Level(), testCtx.params.MaxLevel())

		verifyTestVectors(testCtx, decryptorSk0, coeffs, ciphertext, t)

	})
}

func testRefreshAndPermute(testCtx *testContext, t *testing.T) {

	evaluator := testCtx.evaluator
	encryptorPk0 := testCtx.encryptorPk0
	decryptorSk0 := testCtx.decryptorSk0
	sk0Shards := testCtx.sk0Shards

	levelStart := uint64(3)

	t.Run(testString("RefreshAndPermute/", parties, testCtx.params), func(t *testing.T) {

		if testCtx.params.MaxLevel() < 3 {
			t.Skip("skipping test for params max level < 3")
		}

		type Party struct {
			*PermuteProtocol
			s      *ring.Poly
			share1 RefreshShareDecrypt
			share2 RefreshShareRecrypt
		}

		RefreshParties := make([]*Party, parties)
		for i := uint64(0); i < parties; i++ {
			p := new(Party)
			p.PermuteProtocol = NewPermuteProtocol(testCtx.params)
			p.s = sk0Shards[i].Get()
			p.share1, p.share2 = p.AllocateShares(levelStart)
			RefreshParties[i] = p
		}

		P0 := RefreshParties[0]

		crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.dckksContext.ringQ)
		crp := crpGenerator.ReadNew()

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, 1.0, t)

		for ciphertext.Level() != levelStart {
			evaluator.DropLevel(ciphertext, 1)
		}

		permutation := make([]uint64, testCtx.params.Slots())

		for i := range permutation {
			permutation[i] = ring.RandUniform(testCtx.prng, testCtx.params.Slots(), testCtx.params.Slots()-1)
		}

		for i, p := range RefreshParties {
			p.GenShares(p.s, levelStart, parties, ciphertext, crp, testCtx.params.Slots(), permutation, p.share1, p.share2)
			if i > 0 {
				P0.Aggregate(p.share1, P0.share1, P0.share1)
				P0.Aggregate(p.share2, P0.share2, P0.share2)
			}
		}

		// We refresh the ciphertext with the simulated error
		P0.Decrypt(ciphertext, P0.share1)                           // Masked decryption
		P0.Permute(ciphertext, permutation, testCtx.params.Slots()) // Masked re-encoding
		P0.Recrypt(ciphertext, crp, P0.share2)                      // Masked re-encryption

		coeffsPermute := make([]complex128, len(coeffs))

		for i := range coeffs {
			coeffsPermute[i] = coeffs[permutation[i]]
		}

		require.Equal(t, ciphertext.Level(), testCtx.params.MaxLevel())

		verifyTestVectors(testCtx, decryptorSk0, coeffsPermute, ciphertext, t)

	})
}

func newTestVectors(testCtx *testContext, encryptor ckks.Encryptor, a float64, t *testing.T) (values []complex128, plaintext *ckks.Plaintext, ciphertext *ckks.Ciphertext) {

	slots := testCtx.params.Slots()

	values = make([]complex128, slots)

	for i := uint64(0); i < slots; i++ {
		values[i] = utils.RandComplex128(-a, a)
	}

	values[0] = complex(0.607538, 0.555668)

	plaintext = testCtx.encoder.EncodeNew(values, testCtx.params.LogSlots())

	ciphertext = encryptor.EncryptNew(plaintext)

	return values, plaintext, ciphertext
}

func verifyTestVectors(testCtx *testContext, decryptor ckks.Decryptor, valuesWant []complex128, element interface{}, t *testing.T) {

	var plaintextTest *ckks.Plaintext
	var valuesTest []complex128

	switch element := element.(type) {
	case *ckks.Ciphertext:
		plaintextTest = decryptor.DecryptNew(element)
	case *ckks.Plaintext:
		plaintextTest = element
	}

	slots := testCtx.params.Slots()

	valuesTest = testCtx.encoder.Decode(plaintextTest, testCtx.params.LogSlots())

	var deltaReal, deltaImag float64

	var minprec, maxprec, meanprec, medianprec complex128

	diff := make([]complex128, slots)

	minprec = complex(0, 0)
	maxprec = complex(1, 1)

	meanprec = complex(0, 0)

	distribReal := make(map[uint64]uint64)
	distribImag := make(map[uint64]uint64)

	for i := range valuesWant {

		deltaReal = math.Abs(real(valuesTest[i]) - real(valuesWant[i]))
		deltaImag = math.Abs(imag(valuesTest[i]) - imag(valuesWant[i]))

		diff[i] += complex(deltaReal, 0)
		diff[i] += complex(0, deltaImag)

		meanprec += diff[i]

		if real(diff[i]) > real(minprec) {
			minprec = complex(real(diff[i]), 0)
		}

		if imag(diff[i]) > imag(minprec) {
			minprec = complex(real(minprec), imag(diff[i]))
		}

		if real(diff[i]) < real(maxprec) {
			maxprec = complex(real(diff[i]), 0)
		}

		if imag(diff[i]) < imag(maxprec) {
			maxprec = complex(real(maxprec), imag(diff[i]))
		}

		distribReal[uint64(math.Floor(math.Log2(1/real(diff[i]))))]++
		distribImag[uint64(math.Floor(math.Log2(1/imag(diff[i]))))]++
	}

	meanprec /= complex(float64(slots), 0)
	medianprec = calcmedian(diff)

	if *printPrecisionStats {
		t.Logf("Minimum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(minprec)), math.Log2(1/imag(minprec)))
		t.Logf("Maximum precision : (%.2f, %.2f) bits \n", math.Log2(1/real(maxprec)), math.Log2(1/imag(maxprec)))
		t.Logf("Mean    precision : (%.2f, %.2f) bits \n", math.Log2(1/real(meanprec)), math.Log2(1/imag(meanprec)))
		t.Logf("Median  precision : (%.2f, %.2f) bits \n", math.Log2(1/real(medianprec)), math.Log2(1/imag(medianprec)))
		t.Log()
	}

	require.GreaterOrEqual(t, math.Log2(1/real(medianprec)), minPrec)
	require.GreaterOrEqual(t, math.Log2(1/imag(medianprec)), minPrec)
}

func calcmedian(values []complex128) (median complex128) {

	tmp := make([]float64, len(values))

	for i := range values {
		tmp[i] = real(values[i])
	}

	sort.Float64s(tmp)

	for i := range values {
		values[i] = complex(tmp[i], imag(values[i]))
	}

	for i := range values {
		tmp[i] = imag(values[i])
	}

	sort.Float64s(tmp)

	for i := range values {
		values[i] = complex(real(values[i]), tmp[i])
	}

	index := len(values) / 2

	if len(values)&1 == 1 {
		return values[index]
	}

	if index+1 == len(values) {
		return values[index]
	}

	return (values[index] + values[index+1]) / 2
}
