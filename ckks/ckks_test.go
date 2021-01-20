package ckks

import (
	"flag"
	"fmt"
	"math"
	"math/cmplx"
	"runtime"
	"testing"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters + secure bootstrapping). Overrides -short and requires -timeout=0.")
var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")
var testBootstrapping = flag.Bool("test-bootstrapping", false, "run the bootstrapping tests (memory intensive)")

var minPrec float64 = 15.0

func testString(testContext *testParams, opname string) string {
	return fmt.Sprintf("%slogN=%d/LogSlots=%d/logQP=%d/levels=%d/a=%d/b=%d",
		opname,
		testContext.params.LogN(),
		testContext.params.LogSlots(),
		testContext.params.LogQP(),
		testContext.params.MaxLevel()+1,
		testContext.params.Alpha(),
		testContext.params.Beta())
}

type testParams struct {
	params      *Parameters
	ringQ       *ring.Ring
	ringP       *ring.Ring
	ringQP      *ring.Ring
	prng        utils.PRNG
	encoder     Encoder
	kgen        KeyGenerator
	sk          *SecretKey
	pk          *PublicKey
	rlk         *EvaluationKey
	encryptorPk Encryptor
	encryptorSk Encryptor
	decryptor   Decryptor
	evaluator   Evaluator
}

func TestCKKS(t *testing.T) {

	var err error

	var defaultParams = DefaultParams[PN12QP109 : PN12QP109+4] // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15

	if testing.Short() {
		defaultParams = DefaultParams[PN12QP109 : PN12QP109+2] // the short test suite runs for ring degree N=2^12, 2^13
	}

	if *flagLongTest {
		defaultParams = DefaultParams // the long test suite runs for all default parameters
	}

	for _, defaultParam := range defaultParams {
		var testContext *testParams
		if testContext, err = genTestParams(defaultParam, 0); err != nil {
			panic(err)
		}

		for _, testSet := range []func(testContext *testParams, t *testing.T){
			testParameters,
			testEncoder,
			testEncryptor,
			testEvaluatorAdd,
			testEvaluatorSub,
			testEvaluatorRescale,
			testEvaluatorAddConst,
			testEvaluatorMultByConst,
			testEvaluatorMultByConstAndAdd,
			testEvaluatorMul,
			testFunctions,
			testDecryptPublic,
			testEvaluatePoly,
			testChebyshevInterpolator,
			testSwitchKeys,
			testConjugate,
			testRotateColumns,
			testMarshaller,
		} {
			testSet(testContext, t)
			runtime.GC()
		}
	}
}

func genTestParams(defaultParam *Parameters, hw uint64) (testContext *testParams, err error) {

	testContext = new(testParams)

	testContext.params = defaultParam.Copy()

	testContext.kgen = NewKeyGenerator(testContext.params)

	if hw == 0 {
		testContext.sk, testContext.pk = testContext.kgen.GenKeyPair()
	} else {
		testContext.sk, testContext.pk = testContext.kgen.GenKeyPairSparse(hw)
	}

	if testContext.ringQ, err = ring.NewRing(testContext.params.N(), testContext.params.qi); err != nil {
		return nil, err
	}

	if testContext.ringQP, err = ring.NewRing(testContext.params.N(), append(testContext.params.qi, testContext.params.pi...)); err != nil {
		return nil, err
	}

	if testContext.params.PiCount() != 0 {
		if testContext.ringP, err = ring.NewRing(testContext.params.N(), testContext.params.pi); err != nil {
			return nil, err
		}

		testContext.rlk = testContext.kgen.GenRelinKey(testContext.sk)
	}

	if testContext.prng, err = utils.NewPRNG(); err != nil {
		return nil, err
	}

	testContext.encoder = NewEncoder(testContext.params)

	testContext.encryptorPk = NewEncryptorFromPk(testContext.params, testContext.pk)
	testContext.encryptorSk = NewEncryptorFromSk(testContext.params, testContext.sk)
	testContext.decryptor = NewDecryptor(testContext.params, testContext.sk)

	testContext.evaluator = NewEvaluator(testContext.params)

	return testContext, nil

}

func newTestVectors(testContext *testParams, encryptor Encryptor, a, b complex128, t *testing.T) (values []complex128, plaintext *Plaintext, ciphertext *Ciphertext) {

	logSlots := testContext.params.LogSlots()

	values = make([]complex128, 1<<logSlots)

	for i := uint64(0); i < 1<<logSlots; i++ {
		values[i] = complex(utils.RandFloat64(real(a), real(b)), utils.RandFloat64(imag(a), imag(b)))
	}

	values[0] = complex(0.607538, 0)

	plaintext = testContext.encoder.EncodeNTTAtLvlNew(testContext.params.MaxLevel(), values, logSlots)

	if encryptor != nil {

		switch encryptor := encryptor.(type) {
		case *pkEncryptor:
			if testContext.params.PiCount() != 0 {
				ciphertext = encryptor.EncryptNew(plaintext)
			} else {
				ciphertext = encryptor.EncryptFastNew(plaintext)
			}
		case *skEncryptor:
			ciphertext = encryptor.EncryptNew(plaintext)
		}
	}

	return values, plaintext, ciphertext
}

func verifyTestVectors(testContext *testParams, decryptor Decryptor, valuesWant []complex128, element interface{}, t *testing.T, bound float64) {

	precStats := GetPrecisionStats(testContext.params, testContext.encoder, decryptor, valuesWant, element, bound)

	if *printPrecisionStats {
		t.Log(precStats.String())
	}

	require.GreaterOrEqual(t, real(precStats.MeanPrecision), minPrec)
	require.GreaterOrEqual(t, imag(precStats.MeanPrecision), minPrec)
}

func testParameters(testContext *testParams, t *testing.T) {

	t.Run("Parameters/NewParametersFromModuli/", func(t *testing.T) {
		p, err := NewParametersFromModuli(testContext.params.LogN(), testContext.params.Moduli())
		p.SetLogSlots(testContext.params.LogSlots())
		p.SetScale(testContext.params.Scale())
		assert.NoError(t, err)
		assert.True(t, p.Equals(testContext.params))
	})

	t.Run("Parameters/NewParametersFromLogModuli/", func(t *testing.T) {
		p, err := NewParametersFromLogModuli(testContext.params.LogN(), testContext.params.LogModuli())
		p.SetLogSlots(testContext.params.LogSlots())
		p.SetScale(testContext.params.Scale())
		assert.NoError(t, err)
		assert.True(t, p.Equals(testContext.params))
	})
}

func testEncoder(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Encoder/Slots/"), func(t *testing.T) {

		values, plaintext, _ := newTestVectors(testContext, nil, complex(-1, -1), complex(1, 1), t)

		verifyTestVectors(testContext, testContext.decryptor, values, plaintext, t, 0)
	})

	t.Run(testString(testContext, "Encoder/Coeffs/"), func(t *testing.T) {

		slots := testContext.params.N()

		valuesWant := make([]float64, slots)

		for i := uint64(0); i < slots; i++ {
			valuesWant[i] = utils.RandFloat64(-1, 1)
		}

		valuesWant[0] = 0.607538

		plaintext := NewPlaintext(testContext.params, testContext.params.MaxLevel(), testContext.params.Scale())

		testContext.encoder.EncodeCoeffs(valuesWant, plaintext)

		valuesTest := testContext.encoder.DecodeCoeffs(plaintext)

		var meanprec float64

		for i := range valuesWant {
			meanprec += math.Abs(valuesTest[i] - valuesWant[i])
		}

		meanprec /= float64(slots)

		if *printPrecisionStats {
			t.Log(fmt.Sprintf("\nMean    precision : %.2f \n", math.Log2(1/meanprec)))
		}

		require.GreaterOrEqual(t, math.Log2(1/meanprec), minPrec)
	})

}

func testEncryptor(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Encryptor/EncryptFromPk/Lvl=Max/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorPk, complex(-1, -1), complex(1, 1), t)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})

	t.Run(testString(testContext, "Encryptor/EncryptFromPkFast/Lvl=Max/"), func(t *testing.T) {

		logSlots := testContext.params.LogSlots()

		values := make([]complex128, 1<<logSlots)

		for i := uint64(0); i < 1<<logSlots; i++ {
			values[i] = utils.RandComplex128(-1, 1)
		}

		values[0] = complex(0.607538, 0.555668)

		plaintext := testContext.encoder.EncodeNew(values, logSlots)

		verifyTestVectors(testContext, testContext.decryptor, values, testContext.encryptorPk.EncryptFastNew(plaintext), t, 0)
	})

	t.Run(testString(testContext, "Encryptor/EncryptFromSk/Lvl=Max/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})

	t.Run(testString(testContext, "Encryptor/EncryptFromPk/Lvl=1/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		if testContext.params.MaxLevel() < 1 {
			t.Skip("skipping test for params max level < 1")
		}

		logSlots := testContext.params.LogSlots()

		values := make([]complex128, 1<<logSlots)

		for i := uint64(0); i < 1<<logSlots; i++ {
			values[i] = utils.RandComplex128(-1, 1)
		}

		values[0] = complex(0.607538, 0.555668)

		plaintext := testContext.encoder.EncodeAtLvlNew(1, values, logSlots)

		verifyTestVectors(testContext, testContext.decryptor, values, testContext.encryptorPk.EncryptNew(plaintext), t, 0)
	})

	t.Run(testString(testContext, "Encryptor/EncryptFromPkFast/Lvl=1/"), func(t *testing.T) {

		if testContext.params.MaxLevel() < 1 {
			t.Skip("skipping test for params max level < 1")
		}

		logSlots := testContext.params.LogSlots()

		values := make([]complex128, 1<<logSlots)

		for i := uint64(0); i < 1<<logSlots; i++ {
			values[i] = utils.RandComplex128(-1, 1)
		}

		values[0] = complex(0.607538, 0.555668)

		plaintext := testContext.encoder.EncodeAtLvlNew(1, values, logSlots)

		verifyTestVectors(testContext, testContext.decryptor, values, testContext.encryptorPk.EncryptFastNew(plaintext), t, 0)
	})

	t.Run(testString(testContext, "Encryptor/EncryptFromSk/Lvl=1/"), func(t *testing.T) {

		if testContext.params.MaxLevel() < 1 {
			t.Skip("skipping test for params max level < 1")
		}

		logSlots := testContext.params.LogSlots()

		values := make([]complex128, 1<<logSlots)

		for i := uint64(0); i < 1<<logSlots; i++ {
			values[i] = utils.RandComplex128(-1, 1)
		}

		values[0] = complex(0.607538, 0.555668)

		plaintext := testContext.encoder.EncodeAtLvlNew(1, values, logSlots)

		verifyTestVectors(testContext, testContext.decryptor, values, testContext.encryptorSk.EncryptNew(plaintext), t, 0)
	})

}

func testEvaluatorAdd(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "EvaluatorAdd/CtCtInPlace/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		testContext.evaluator.Add(ciphertext1, ciphertext2, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t, 0)
	})

	t.Run(testString(testContext, "EvaluatorAdd/CtCtNew/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		ciphertext3 := testContext.evaluator.AddNew(ciphertext1, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext3, t, 0)
	})

	t.Run(testString(testContext, "EvaluatorAdd/CtPlainInPlace/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, plaintext2, _ := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		testContext.evaluator.Add(ciphertext1, plaintext2, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t, 0)

		for i := range values1 {
			values1[i] += values2[i]
		}

		testContext.evaluator.Add(plaintext2, ciphertext1, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t, 0)
	})

	t.Run(testString(testContext, "EvaluatorAdd/CtPlainInPlaceNew/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, plaintext2, _ := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		ciphertext3 := testContext.evaluator.AddNew(ciphertext1, plaintext2)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext3, t, 0)

		ciphertext3 = testContext.evaluator.AddNew(plaintext2, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext3, t, 0)
	})

}

func testEvaluatorSub(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "EvaluatorSub/CtCtInPlace/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] -= values2[i]
		}

		testContext.evaluator.Sub(ciphertext1, ciphertext2, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t, 0)
	})

	t.Run(testString(testContext, "EvaluatorSub/CtCtNew/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] -= values2[i]
		}

		ciphertext3 := testContext.evaluator.SubNew(ciphertext1, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext3, t, 0)
	})

	t.Run(testString(testContext, "EvaluatorSub/CtPlainInPlace/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, plaintext2, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		valuesTest := make([]complex128, len(values1))
		for i := range values1 {
			valuesTest[i] = values1[i] - values2[i]
		}

		testContext.evaluator.Sub(ciphertext1, plaintext2, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, valuesTest, ciphertext2, t, 0)

		for i := range values1 {
			valuesTest[i] = values2[i] - values1[i]
		}

		testContext.evaluator.Sub(plaintext2, ciphertext1, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, valuesTest, ciphertext2, t, 0)
	})

	t.Run(testString(testContext, "EvaluatorSub/CtPlainNew/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, plaintext2, _ := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		valuesTest := make([]complex128, len(values1))
		for i := range values1 {
			valuesTest[i] = values1[i] - values2[i]
		}

		ciphertext3 := testContext.evaluator.SubNew(ciphertext1, plaintext2)

		verifyTestVectors(testContext, testContext.decryptor, valuesTest, ciphertext3, t, 0)

		for i := range values1 {
			valuesTest[i] = values2[i] - values1[i]
		}

		ciphertext3 = testContext.evaluator.SubNew(plaintext2, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, valuesTest, ciphertext3, t, 0)
	})

}

func testEvaluatorRescale(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "EvaluatorRescale/Single/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		constant := testContext.ringQ.Modulus[ciphertext.Level()]

		testContext.evaluator.MultByConst(ciphertext, constant, ciphertext)

		ciphertext.MulScale(float64(constant))

		testContext.evaluator.Rescale(ciphertext, testContext.params.Scale(), ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})

	t.Run(testString(testContext, "EvaluatorRescale/Many/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		nbRescales := testContext.params.MaxLevel()
		if nbRescales > 5 {
			nbRescales = 5
		}

		for i := uint64(0); i < nbRescales; i++ {
			constant := testContext.ringQ.Modulus[ciphertext.Level()]
			testContext.evaluator.MultByConst(ciphertext, constant, ciphertext)
			ciphertext.MulScale(float64(constant))
		}

		testContext.evaluator.RescaleMany(ciphertext, nbRescales, ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})
}

func testEvaluatorAddConst(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "EvaluatorAddConst/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		constant := complex(3.1415, -1.4142)

		for i := range values {
			values[i] += constant
		}

		testContext.evaluator.AddConst(ciphertext, constant, ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})

}

func testEvaluatorMultByConst(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "EvaluatorMultByConst/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		constant := 1.0 / complex(3.1415, -1.4142)

		for i := range values {
			values[i] *= constant
		}

		testContext.evaluator.MultByConst(ciphertext, constant, ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})

}

func testEvaluatorMultByConstAndAdd(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "EvaluatorMultByConstAndAdd/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		constant := 1.0 / complex(3.1415, -1.4142)

		for i := range values1 {
			values2[i] += (constant * values1[i])
		}

		testContext.evaluator.MultByConstAndAdd(ciphertext1, constant, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext2, t, 0)
	})

}

func testEvaluatorMul(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "EvaluatorMul/ct0*pt->ct0/"), func(t *testing.T) {

		values1, plaintext1, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, plaintext1, nil, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t, 0)
	})

	t.Run(testString(testContext, "EvaluatorMul/pt*ct0->ct0/"), func(t *testing.T) {

		values1, plaintext1, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, plaintext1, nil, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t, 0)
	})

	t.Run(testString(testContext, "EvaluatorMul/ct0*pt->ct1/"), func(t *testing.T) {

		values1, plaintext1, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		ciphertext2 := testContext.evaluator.MulRelinNew(ciphertext1, plaintext1, nil)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext2, t, 0)
	})

	t.Run(testString(testContext, "EvaluatorMul/ct0*ct1->ct0/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values2[i] *= values1[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, ciphertext2, nil, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext1, t, 0)
	})

	t.Run(testString(testContext, "EvaluatorMul/ct0*ct1->ct1/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values2[i] *= values1[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, ciphertext2, nil, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext2, t, 0)
	})

	t.Run(testString(testContext, "EvaluatorMul/ct0*ct1->ct2/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values2[i] *= values1[i]
		}

		ciphertext3 := testContext.evaluator.MulRelinNew(ciphertext1, ciphertext2, nil)

		verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext3, t, 0)
	})

	t.Run(testString(testContext, "EvaluatorMul/ct0*ct0->ct0/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, ciphertext1, nil, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t, 0)
	})

	t.Run(testString(testContext, "EvaluatorMul/ct0*ct0->ct1/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		ciphertext2 := testContext.evaluator.MulRelinNew(ciphertext1, ciphertext1, testContext.rlk)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext2, t, 0)
	})

	t.Run(testString(testContext, "EvaluatorMul/Relinearize(ct0*ct1->ct0)/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values2[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, ciphertext2, nil, ciphertext1)

		testContext.evaluator.Relinearize(ciphertext1, testContext.rlk, ciphertext1)

		require.Equal(t, ciphertext1.Degree(), uint64(1))

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t, 0)
	})

	t.Run(testString(testContext, "EvaluatorMul/Relinearize(ct0*ct1->ct1)/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values2[i] *= values1[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, ciphertext2, nil, ciphertext2)

		testContext.evaluator.Relinearize(ciphertext2, testContext.rlk, ciphertext2)

		require.Equal(t, ciphertext1.Degree(), uint64(1))

		verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext2, t, 0)
	})

}

func testFunctions(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Functions/PowerOf2/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		if testContext.params.MaxLevel() < 3 {
			t.Skip("skipping test for params max level < 3")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		n := uint64(2)

		valuesWant := make([]complex128, len(values))
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values[i]
		}

		for i := uint64(0); i < n; i++ {
			for j := 0; j < len(valuesWant); j++ {
				valuesWant[j] *= valuesWant[j]
			}
		}

		testContext.evaluator.PowerOf2(ciphertext, n, testContext.rlk, ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, valuesWant, ciphertext, t, 0)
	})

	t.Run(testString(testContext, "Functions/Power/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		if testContext.params.MaxLevel() < 4 {
			t.Skip("skipping test for params max level < 4")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		n := uint64(3)

		for i := range values {
			values[i] = cmplx.Pow(values[i], complex(float64(n), 0))
		}

		testContext.evaluator.Power(ciphertext, n, testContext.rlk, ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})

	t.Run(testString(testContext, "Functions/Inverse/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		if testContext.params.MaxLevel() < 7 {
			t.Skip("skipping test for params max level < 7")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(0.1, 0), complex(1, 0), t)

		n := uint64(7)

		for i := range values {
			values[i] = 1.0 / values[i]
		}

		ciphertext = testContext.evaluator.InverseNew(ciphertext, n, testContext.rlk)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})
}

func testEvaluatePoly(testContext *testParams, t *testing.T) {

	var err error

	t.Run(testString(testContext, "EvaluatePoly/Exp/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		if testContext.params.MaxLevel() < 3 {
			t.Skip("skipping test for params max level < 3")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, 0), complex(1, 0), t)

		coeffs := []complex128{
			complex(1.0, 0),
			complex(1.0, 0),
			complex(1.0/2, 0),
			complex(1.0/6, 0),
			complex(1.0/24, 0),
			complex(1.0/120, 0),
			complex(1.0/720, 0),
			complex(1.0/5040, 0),
		}

		poly := NewPoly(coeffs)

		for i := range values {
			values[i] = cmplx.Exp(values[i])
		}

		if ciphertext, err = testContext.evaluator.EvaluatePoly(ciphertext, poly, testContext.rlk); err != nil {
			t.Error(err)
		}

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})
}

func testChebyshevInterpolator(testContext *testParams, t *testing.T) {

	var err error

	t.Run(testString(testContext, "ChebyshevInterpolator/Sin/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		if testContext.params.MaxLevel() < 5 {
			t.Skip("skipping test for params max level < 5")
		}

		eval := testContext.evaluator

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, 0), complex(1, 0), t)

		cheby := Approximate(cmplx.Sin, complex(-1.5, 0), complex(1.5, 0), 15)

		for i := range values {
			values[i] = cmplx.Sin(values[i])
		}

		eval.MultByConst(ciphertext, 2/(cheby.b-cheby.a), ciphertext)
		eval.AddConst(ciphertext, (-cheby.a-cheby.b)/(cheby.b-cheby.a), ciphertext)
		eval.Rescale(ciphertext, eval.(*evaluator).scale, ciphertext)

		if ciphertext, err = eval.EvaluateCheby(ciphertext, cheby, testContext.rlk); err != nil {
			t.Error(err)
		}

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})
}

func testDecryptPublic(testContext *testParams, t *testing.T) {

	var err error

	t.Run(testString(testContext, "DecryptPublic/Sin/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		if testContext.params.MaxLevel() < 5 {
			t.Skip("skipping test for params max level < 5")
		}

		eval := testContext.evaluator

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, 0), complex(1, 0), t)

		cheby := Approximate(cmplx.Sin, complex(-1.5, 0), complex(1.5, 0), 15)

		for i := range values {
			values[i] = cmplx.Sin(values[i])
		}

		eval.MultByConst(ciphertext, 2/(cheby.b-cheby.a), ciphertext)
		eval.AddConst(ciphertext, (-cheby.a-cheby.b)/(cheby.b-cheby.a), ciphertext)
		eval.Rescale(ciphertext, eval.(*evaluator).scale, ciphertext)

		if ciphertext, err = eval.EvaluateCheby(ciphertext, cheby, testContext.rlk); err != nil {
			t.Error(err)
		}

		plaintext := testContext.decryptor.DecryptNew(ciphertext)

		valuesHave := testContext.encoder.Decode(plaintext, testContext.params.LogSlots())

		verifyTestVectors(testContext, nil, values, valuesHave, t, 0)

		sigma := testContext.encoder.GetErrSTDTimeDom(values, valuesHave, plaintext.Scale())

		valuesHave = testContext.encoder.DecodePublic(plaintext, testContext.params.LogSlots(), sigma)

		verifyTestVectors(testContext, nil, values, valuesHave, t, 0)
	})
}

func testSwitchKeys(testContext *testParams, t *testing.T) {

	var sk2 *SecretKey
	var decryptorSk2 Decryptor
	var switchingKey *SwitchingKey

	if testContext.params.PiCount() != 0 {
		sk2 = testContext.kgen.GenSecretKey()
		decryptorSk2 = NewDecryptor(testContext.params, sk2)
		switchingKey = testContext.kgen.GenSwitchingKey(testContext.sk, sk2)
	}

	t.Run(testString(testContext, "SwitchKeys/InPlace/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		testContext.evaluator.SwitchKeys(ciphertext, switchingKey, ciphertext)

		verifyTestVectors(testContext, decryptorSk2, values, ciphertext, t, 0)
	})

	t.Run(testString(testContext, "SwitchKeys/New/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		ciphertext = testContext.evaluator.SwitchKeysNew(ciphertext, switchingKey)

		verifyTestVectors(testContext, decryptorSk2, values, ciphertext, t, 0)
	})

}

func testConjugate(testContext *testParams, t *testing.T) {

	var rotKey *RotationKeys
	if testContext.params.PiCount() != 0 {
		rotKey = NewRotationKeys()
		testContext.kgen.GenRotationKey(Conjugate, testContext.sk, 0, rotKey)
	}

	t.Run(testString(testContext, "Conjugate/InPlace/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values {
			values[i] = complex(real(values[i]), -imag(values[i]))
		}

		testContext.evaluator.Conjugate(ciphertext, rotKey, ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})

	t.Run(testString(testContext, "Conjugate/New/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values {
			values[i] = complex(real(values[i]), -imag(values[i]))
		}

		ciphertext = testContext.evaluator.ConjugateNew(ciphertext, rotKey)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t, 0)
	})

}

func testRotateColumns(testContext *testParams, t *testing.T) {

	var rotKey *RotationKeys
	if testContext.params.PiCount() != 0 {
		rotKey = testContext.kgen.GenRotationKeysPow2(testContext.sk)
	}

	t.Run(testString(testContext, "RotateColumns/InPlace/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		values2 := make([]complex128, len(values1))
		ciphertext2 := NewCiphertext(testContext.params, ciphertext1.Degree(), ciphertext1.Level(), ciphertext1.Scale())

		for n := 1; n < len(values1); n <<= 1 {

			// Applies the column rotation to the values
			for i := range values1 {
				values2[i] = values1[(i+n)%len(values1)]
			}

			testContext.evaluator.Rotate(ciphertext1, uint64(n), rotKey, ciphertext2)

			verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext2, t, 0)
		}
	})

	t.Run(testString(testContext, "RotateColumns/New/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		values2 := make([]complex128, len(values1))

		for n := 1; n < len(values1); n <<= 1 {

			// Applies the column rotation to the values
			for i := range values1 {
				values2[i] = values1[(i+n)%len(values1)]
			}

			verifyTestVectors(testContext, testContext.decryptor, values2, testContext.evaluator.RotateNew(ciphertext1, uint64(n), rotKey), t, 0)

		}

	})

	t.Run(testString(testContext, "RotateColumns/Random/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		values2 := make([]complex128, len(values1))
		ciphertext2 := NewCiphertext(testContext.params, ciphertext1.Degree(), ciphertext1.Level(), ciphertext1.Scale())

		for n := 1; n < 4; n++ {

			rand := utils.RandUint64() % uint64(len(values1))

			// Applies the column rotation to the values
			for i := range values1 {
				values2[i] = values1[(i+int(rand))%len(values1)]
			}

			testContext.evaluator.Rotate(ciphertext1, rand, rotKey, ciphertext2)

			verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext2, t, 0)
		}

	})

	t.Run(testString(testContext, "RotateColumns/Hoisted/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		values2 := make([]complex128, len(values1))
		rotations := []uint64{0, 1, 2, 3, 4, 5}
		for _, n := range rotations {
			testContext.kgen.GenRotationKey(RotationLeft, testContext.sk, n, rotKey)
		}

		ciphertexts := testContext.evaluator.RotateHoisted(ciphertext1, rotations, rotKey)

		for _, n := range rotations {

			for i := range values1 {
				values2[i] = values1[(i+int(n))%len(values1)]
			}

			verifyTestVectors(testContext, testContext.decryptor, values2, ciphertexts[n], t, 0)
		}
	})
}

func testMarshaller(testContext *testParams, t *testing.T) {

	ringQP := testContext.ringQP

	t.Run("Marshaller/Ciphertext/", func(t *testing.T) {
		t.Run(testString(testContext, "EndToEnd/"), func(t *testing.T) {

			ciphertextWant := NewCiphertextRandom(testContext.prng, testContext.params, 2, testContext.params.MaxLevel(), testContext.params.Scale())

			marshalledCiphertext, err := ciphertextWant.MarshalBinary()
			require.NoError(t, err)

			ciphertextTest := new(Ciphertext)
			require.NoError(t, ciphertextTest.UnmarshalBinary(marshalledCiphertext))

			require.Equal(t, ciphertextWant.Degree(), ciphertextTest.Degree())
			require.Equal(t, ciphertextWant.Level(), ciphertextTest.Level())
			require.Equal(t, ciphertextWant.Scale(), ciphertextTest.Scale())

			for i := range ciphertextWant.value {
				require.True(t, testContext.ringQ.EqualLvl(ciphertextWant.Level(), ciphertextWant.Value()[i], ciphertextTest.Value()[i]))
			}
		})

		t.Run(testString(testContext, "Minimal/"), func(t *testing.T) {

			ciphertext := NewCiphertextRandom(testContext.prng, testContext.params, 0, testContext.params.MaxLevel(), testContext.params.Scale())

			marshalledCiphertext, err := ciphertext.MarshalBinary()
			require.NoError(t, err)

			ciphertextTest := new(Ciphertext)
			require.Error(t, ciphertextTest.UnmarshalBinary(nil))
			require.NoError(t, ciphertextTest.UnmarshalBinary(marshalledCiphertext))

			require.Equal(t, ciphertext.Degree(), uint64(0))
			require.Equal(t, ciphertext.Level(), testContext.params.MaxLevel())
			require.Equal(t, ciphertext.Scale(), testContext.params.Scale())
			require.Equal(t, len(ciphertext.Value()), 1)
		})
	})

	t.Run(testString(testContext, "Marshaller/Sk/"), func(t *testing.T) {

		marshalledSk, err := testContext.sk.MarshalBinary()
		require.NoError(t, err)

		sk := new(SecretKey)
		err = sk.UnmarshalBinary(marshalledSk)
		require.NoError(t, err)

		require.True(t, ringQP.Equal(sk.sk, testContext.sk.sk))

	})

	t.Run(testString(testContext, "Marshaller/Pk/"), func(t *testing.T) {

		marshalledPk, err := testContext.pk.MarshalBinary()
		require.NoError(t, err)

		pk := new(PublicKey)
		err = pk.UnmarshalBinary(marshalledPk)
		require.NoError(t, err)

		for k := range testContext.pk.pk {
			require.Truef(t, ringQP.Equal(pk.pk[k], testContext.pk.pk[k]), "Marshal PublicKey element [%d]", k)
		}
	})

	t.Run(testString(testContext, "Marshaller/EvaluationKey/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		evalKey := testContext.kgen.GenRelinKey(testContext.sk)
		data, err := evalKey.MarshalBinary()
		require.NoError(t, err)

		resEvalKey := new(EvaluationKey)
		err = resEvalKey.UnmarshalBinary(data)
		require.NoError(t, err)

		evakeyWant := evalKey.evakey.evakey
		evakeyTest := resEvalKey.evakey.evakey

		for j := range evakeyWant {
			for k := range evakeyWant[j] {
				require.Truef(t, ringQP.Equal(evakeyWant[j][k], evakeyTest[j][k]), "Marshal EvaluationKey element [%d][%d]", j, k)
			}
		}
	})

	t.Run(testString(testContext, "Marshaller/SwitchingKey/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		skOut := testContext.kgen.GenSecretKey()

		switchingKey := testContext.kgen.GenSwitchingKey(testContext.sk, skOut)
		data, err := switchingKey.MarshalBinary()
		require.NoError(t, err)

		resSwitchingKey := new(SwitchingKey)
		err = resSwitchingKey.UnmarshalBinary(data)
		require.NoError(t, err)

		evakeyWant := switchingKey.evakey
		evakeyTest := resSwitchingKey.evakey

		for j := range evakeyWant {
			for k := range evakeyWant[j] {
				require.True(t, ringQP.Equal(evakeyWant[j][k], evakeyTest[j][k]))
			}
		}
	})

	t.Run(testString(testContext, "Marshaller/RotationKey/"), func(t *testing.T) {

		if testContext.params.PiCount() == 0 {
			t.Skip("#Pi is empty")
		}

		rotationKey := NewRotationKeys()

		testContext.kgen.GenRotationKey(Conjugate, testContext.sk, 0, rotationKey)
		testContext.kgen.GenRotationKey(RotationLeft, testContext.sk, 1, rotationKey)
		testContext.kgen.GenRotationKey(RotationLeft, testContext.sk, 2, rotationKey)
		testContext.kgen.GenRotationKey(RotationRight, testContext.sk, 3, rotationKey)
		testContext.kgen.GenRotationKey(RotationRight, testContext.sk, 5, rotationKey)

		data, err := rotationKey.MarshalBinary()
		require.NoError(t, err)

		resRotationKey := new(RotationKeys)
		err = resRotationKey.UnmarshalBinary(data)
		require.NoError(t, err)

		for i := uint64(1); i < testContext.ringQ.N>>1; i++ {

			if rotationKey.evakeyRotColLeft[i] != nil {

				evakeyWant := rotationKey.evakeyRotColLeft[i].evakey
				evakeyTest := resRotationKey.evakeyRotColLeft[i].evakey

				evakeyNTTIndexWant := rotationKey.permuteNTTLeftIndex[i]
				evakeyNTTIndexTest := resRotationKey.permuteNTTLeftIndex[i]

				require.True(t, utils.EqualSliceUint64(evakeyNTTIndexWant, evakeyNTTIndexTest))

				for j := range evakeyWant {
					for k := range evakeyWant[j] {
						require.Truef(t, ringQP.Equal(evakeyWant[j][k], evakeyTest[j][k]), "Marshal RotationKey RotateLeft %d element [%d][%d]", i, j, k)
					}
				}
			}

			if rotationKey.evakeyRotColRight[i] != nil {

				evakeyWant := rotationKey.evakeyRotColRight[i].evakey
				evakeyTest := resRotationKey.evakeyRotColRight[i].evakey

				evakeyNTTIndexWant := rotationKey.permuteNTTRightIndex[i]
				evakeyNTTIndexTest := resRotationKey.permuteNTTRightIndex[i]

				require.True(t, utils.EqualSliceUint64(evakeyNTTIndexWant, evakeyNTTIndexTest))

				for j := range evakeyWant {
					for k := range evakeyWant[j] {
						require.Truef(t, ringQP.Equal(evakeyWant[j][k], evakeyTest[j][k]), "Marshal RotationKey RotateRight %d element [%d][%d]", i, j, k)
					}
				}
			}
		}

		if rotationKey.evakeyConjugate != nil {

			evakeyWant := rotationKey.evakeyConjugate.evakey
			evakeyTest := resRotationKey.evakeyConjugate.evakey

			evakeyNTTIndexWant := rotationKey.permuteNTTConjugateIndex
			evakeyNTTIndexTest := resRotationKey.permuteNTTConjugateIndex

			require.True(t, utils.EqualSliceUint64(evakeyNTTIndexWant, evakeyNTTIndexTest))

			for j := range evakeyWant {
				for k := range evakeyWant[j] {
					require.Truef(t, ringQP.Equal(evakeyWant[j][k], evakeyTest[j][k]), "Marshal RotationKey RotateRow element [%d][%d]", j, k)
				}
			}
		}
	})
}
