package ckks

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"math/cmplx"
	"runtime"
	"testing"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters + secure bootstrapping). Overrides -short and requires -timeout=0.")
var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")
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
		testContext.params.PCount(),
		testContext.params.Beta())
}

type testParams struct {
	params      Parameters
	ringQ       *ring.Ring
	ringP       *ring.Ring
	ringQP      *ring.Ring
	prng        utils.PRNG
	encoder     Encoder
	kgen        KeyGenerator
	sk          *rlwe.SecretKey
	pk          *rlwe.PublicKey
	rlk         *rlwe.RelinearizationKey
	encryptorPk *Encryptor
	encryptorSk *Encryptor
	decryptor   *Decryptor
	evaluator   Evaluator
}

func TestCKKS(t *testing.T) {

	defaultParams := DefaultParams[:4] // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15
	if testing.Short() {
		defaultParams = DefaultParams[:2] // the short test suite runs for ring degree N=2^12, 2^13
	}
	if *flagLongTest {
		defaultParams = append(DefaultParams, DefaultPostQuantumParams...) // the long test suite runs for all default parameters
	}
	if *flagParamString != "" {
		var jsonParams ParametersLiteral
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, defaultParam := range defaultParams {
		params, err := NewParametersFromLiteral(defaultParam)
		if err != nil {
			panic(err)
		}
		var testContext *testParams
		if testContext, err = genTestParams(params, 0); err != nil {
			panic(err)
		}

		for _, testSet := range []func(testContext *testParams, t *testing.T){
			testParameters,
			testEncoder,
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
			testAutomorphisms,
			testInnerSum,
			testReplicate,
			testLinearTransform,
			testMarshaller,
		} {
			testSet(testContext, t)
			runtime.GC()
		}
	}
}

func genTestParams(defaultParam Parameters, hw int) (testContext *testParams, err error) {

	testContext = new(testParams)

	testContext.params = defaultParam

	testContext.kgen = NewKeyGenerator(testContext.params)

	if hw == 0 {
		testContext.sk, testContext.pk = testContext.kgen.GenKeyPair()
	} else {
		testContext.sk, testContext.pk = testContext.kgen.GenKeyPairSparse(hw)
	}

	testContext.ringQ = defaultParam.RingQ()
	testContext.ringQP = defaultParam.RingQP()
	if testContext.params.PCount() != 0 {
		testContext.ringP = defaultParam.RingP()
		testContext.rlk = testContext.kgen.GenRelinearizationKey(testContext.sk)
	}

	if testContext.prng, err = utils.NewPRNG(); err != nil {
		return nil, err
	}

	testContext.encoder = NewEncoder(testContext.params)

	testContext.encryptorPk = NewEncryptor(testContext.params, testContext.pk)
	testContext.encryptorSk = NewEncryptor(testContext.params, testContext.sk)
	testContext.decryptor = NewDecryptor(testContext.params, testContext.sk)

	testContext.evaluator = NewEvaluator(testContext.params, rlwe.EvaluationKey{Rlk: testContext.rlk})

	return testContext, nil

}

func newTestVectors(testContext *testParams, encryptor *Encryptor, a, b complex128, t *testing.T) (values []complex128, plaintext *Plaintext, ciphertext *Ciphertext) {

	logSlots := testContext.params.LogSlots()

	values = make([]complex128, 1<<logSlots)

	for i := 0; i < 1<<logSlots; i++ {
		values[i] = complex(utils.RandFloat64(real(a), real(b)), utils.RandFloat64(imag(a), imag(b)))
	}

	values[0] = complex(0.607538, 0)

	plaintext = testContext.encoder.EncodeNTTAtLvlNew(testContext.params.MaxLevel(), values, logSlots)

	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return values, plaintext, ciphertext
}

func verifyTestVectors(testContext *testParams, decryptor *Decryptor, valuesWant []complex128, element interface{}, logSlots int, bound float64, t *testing.T) {

	precStats := GetPrecisionStats(testContext.params, testContext.encoder, decryptor, valuesWant, element, logSlots, bound)

	if *printPrecisionStats {
		t.Log(precStats.String())
	}

	require.GreaterOrEqual(t, real(precStats.MeanPrecision), minPrec)
	require.GreaterOrEqual(t, imag(precStats.MeanPrecision), minPrec)
}

func testParameters(testContext *testParams, t *testing.T) {
	t.Run(testString(testContext, "Parameters/CopyNew/"), func(t *testing.T) {
		params1, params2 := testContext.params.CopyNew(), testContext.params.CopyNew()
		assert.True(t, params1.Equals(testContext.params) && params2.Equals(testContext.params))
		params1.logSlots = 7
		params1.scale = 3.14
		assert.False(t, params1.Equals(testContext.params))
		assert.True(t, params2.Equals(testContext.params))
	})
}

func testEncoder(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Encoder/Encode/"), func(t *testing.T) {

		values, plaintext, _ := newTestVectors(testContext, nil, complex(-1, -1), complex(1, 1), t)

		verifyTestVectors(testContext, testContext.decryptor, values, plaintext, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "Encoder/EncodeCoeffs/"), func(t *testing.T) {

		slots := testContext.params.N()

		valuesWant := make([]float64, slots)

		for i := 0; i < slots; i++ {
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

func testEvaluatorAdd(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Evaluator/Add/CtCt/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		testContext.evaluator.Add(ciphertext1, ciphertext2, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "Evaluator/AddNew/CtCt/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		ciphertext3 := testContext.evaluator.AddNew(ciphertext1, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext3, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "Evaluator/Add/CtPlain/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, plaintext2, _ := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		testContext.evaluator.Add(ciphertext1, plaintext2, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, testContext.params.LogSlots(), 0, t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		testContext.evaluator.Add(plaintext2, ciphertext1, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "Evaluator/AddNew/CtPlain/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, plaintext2, _ := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		ciphertext3 := testContext.evaluator.AddNew(ciphertext1, plaintext2)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext3, testContext.params.LogSlots(), 0, t)

		ciphertext3 = testContext.evaluator.AddNew(plaintext2, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext3, testContext.params.LogSlots(), 0, t)
	})

}

func testEvaluatorSub(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Evaluator/Sub/CtCt/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] -= values2[i]
		}

		testContext.evaluator.Sub(ciphertext1, ciphertext2, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "Evaluator/SubNew/CtCt/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] -= values2[i]
		}

		ciphertext3 := testContext.evaluator.SubNew(ciphertext1, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext3, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "Evaluator/Sub/CtPlain/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, plaintext2, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		valuesTest := make([]complex128, len(values1))
		for i := range values1 {
			valuesTest[i] = values1[i] - values2[i]
		}

		testContext.evaluator.Sub(ciphertext1, plaintext2, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, valuesTest, ciphertext2, testContext.params.LogSlots(), 0, t)

		for i := range values1 {
			valuesTest[i] = values2[i] - values1[i]
		}

		testContext.evaluator.Sub(plaintext2, ciphertext1, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, valuesTest, ciphertext2, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "Evaluator/SubNew/CtPlain/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, plaintext2, _ := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		valuesTest := make([]complex128, len(values1))
		for i := range values1 {
			valuesTest[i] = values1[i] - values2[i]
		}

		ciphertext3 := testContext.evaluator.SubNew(ciphertext1, plaintext2)

		verifyTestVectors(testContext, testContext.decryptor, valuesTest, ciphertext3, testContext.params.LogSlots(), 0, t)

		for i := range values1 {
			valuesTest[i] = values2[i] - values1[i]
		}

		ciphertext3 = testContext.evaluator.SubNew(plaintext2, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, valuesTest, ciphertext3, testContext.params.LogSlots(), 0, t)
	})

}

func testEvaluatorRescale(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Evaluator/Rescale/Single/"), func(t *testing.T) {

		if testContext.params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		constant := testContext.ringQ.Modulus[ciphertext.Level()]

		testContext.evaluator.MultByConst(ciphertext, constant, ciphertext)

		ciphertext.Scale *= float64(constant)

		testContext.evaluator.Rescale(ciphertext, testContext.params.Scale(), ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "Evaluator/Rescale/Many/"), func(t *testing.T) {

		if testContext.params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		nbRescales := testContext.params.MaxLevel()
		if nbRescales > 5 {
			nbRescales = 5
		}

		for i := 0; i < nbRescales; i++ {
			constant := testContext.ringQ.Modulus[ciphertext.Level()-i]
			testContext.evaluator.MultByConst(ciphertext, constant, ciphertext)
			ciphertext.Scale *= float64(constant)
		}

		testContext.evaluator.Rescale(ciphertext, testContext.params.Scale(), ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, testContext.params.LogSlots(), 0, t)
	})
}

func testEvaluatorAddConst(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Evaluator/AddConst/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		constant := complex(3.1415, -1.4142)

		for i := range values {
			values[i] += constant
		}

		testContext.evaluator.AddConst(ciphertext, constant, ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, testContext.params.LogSlots(), 0, t)
	})

}

func testEvaluatorMultByConst(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Evaluator/MultByConst/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		constant := 1.0 / complex(3.1415, -1.4142)

		for i := range values {
			values[i] *= constant
		}

		testContext.evaluator.MultByConst(ciphertext, constant, ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, testContext.params.LogSlots(), 0, t)
	})

}

func testEvaluatorMultByConstAndAdd(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Evaluator/MultByConstAndAdd/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		constant := 1.0 / complex(3.1415, -1.4142)

		for i := range values1 {
			values2[i] += (constant * values1[i])
		}

		testContext.evaluator.MultByConstAndAdd(ciphertext1, constant, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext2, testContext.params.LogSlots(), 0, t)
	})

}

func testEvaluatorMul(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Evaluator/Mul/ct0*pt->ct0/"), func(t *testing.T) {

		values1, plaintext1, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, plaintext1, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "Evaluator/Mul/pt*ct0->ct0/"), func(t *testing.T) {

		values1, plaintext1, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, plaintext1, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "Evaluator/Mul/ct0*pt->ct1/"), func(t *testing.T) {

		values1, plaintext1, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		ciphertext2 := testContext.evaluator.MulRelinNew(ciphertext1, plaintext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext2, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "Evaluator/Mul/ct0*ct1->ct0/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values2[i] *= values1[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, ciphertext2, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext1, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "Evaluator/Mul/ct0*ct1->ct1/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values2[i] *= values1[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, ciphertext2, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext2, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "Evaluator/Mul/ct0*ct1->ct2/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values2[i] *= values1[i]
		}

		ciphertext3 := testContext.evaluator.MulRelinNew(ciphertext1, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext3, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "Evaluator/Mul/ct0*ct0->ct0/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, ciphertext1, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "Evaluator/Mul/ct0*ct0->ct1/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		ciphertext2 := testContext.evaluator.MulRelinNew(ciphertext1, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext2, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "Evaluator/Mul/Relinearize(ct0*ct1->ct0)/"), func(t *testing.T) {

		if testContext.params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values1[i] *= values2[i]
		}

		testContext.evaluator.Mul(ciphertext1, ciphertext2, ciphertext1)
		require.Equal(t, ciphertext1.Degree(), 2)
		testContext.evaluator.Relinearize(ciphertext1, ciphertext1)
		require.Equal(t, ciphertext1.Degree(), 1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "Evaluator/Mul/Relinearize(ct0*ct1->ct1)/"), func(t *testing.T) {

		if testContext.params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values1 {
			values2[i] *= values1[i]
		}

		testContext.evaluator.Mul(ciphertext1, ciphertext2, ciphertext2)
		require.Equal(t, ciphertext2.Degree(), 2)
		testContext.evaluator.Relinearize(ciphertext2, ciphertext2)
		require.Equal(t, ciphertext2.Degree(), 1)

		verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext2, testContext.params.LogSlots(), 0, t)
	})

}

func testFunctions(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Evaluator/PowerOf2/"), func(t *testing.T) {

		if testContext.params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		if testContext.params.MaxLevel() < 3 {
			t.Skip("skipping test for params max level < 3")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		n := 2

		valuesWant := make([]complex128, len(values))
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values[i]
		}

		for i := 0; i < n; i++ {
			for j := 0; j < len(valuesWant); j++ {
				valuesWant[j] *= valuesWant[j]
			}
		}

		testContext.evaluator.PowerOf2(ciphertext, n, ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, valuesWant, ciphertext, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "Evaluator/Power/"), func(t *testing.T) {

		if testContext.params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		if testContext.params.MaxLevel() < 4 {
			t.Skip("skipping test for params max level < 4")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		n := 3

		for i := range values {
			values[i] = cmplx.Pow(values[i], complex(float64(n), 0))
		}

		testContext.evaluator.Power(ciphertext, n, ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "Evaluator/Inverse/"), func(t *testing.T) {

		if testContext.params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		if testContext.params.MaxLevel() < 7 {
			t.Skip("skipping test for params max level < 7")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(0.1, 0), complex(1, 0), t)

		n := 7

		for i := range values {
			values[i] = 1.0 / values[i]
		}

		ciphertext = testContext.evaluator.InverseNew(ciphertext, n)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, testContext.params.LogSlots(), 0, t)
	})
}

func testEvaluatePoly(testContext *testParams, t *testing.T) {

	var err error

	t.Run(testString(testContext, "EvaluatePoly/Exp/"), func(t *testing.T) {

		if testContext.params.PCount() == 0 {
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

		if ciphertext, err = testContext.evaluator.EvaluatePoly(ciphertext, poly, ciphertext.Scale); err != nil {
			t.Error(err)
		}

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, testContext.params.LogSlots(), 0, t)
	})
}

func testChebyshevInterpolator(testContext *testParams, t *testing.T) {

	var err error

	t.Run(testString(testContext, "ChebyshevInterpolator/Sin/"), func(t *testing.T) {

		if testContext.params.PCount() == 0 {
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

		if ciphertext, err = eval.EvaluateCheby(ciphertext, cheby, ciphertext.Scale); err != nil {
			t.Error(err)
		}

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, testContext.params.LogSlots(), 0, t)
	})
}

func testDecryptPublic(testContext *testParams, t *testing.T) {

	var err error

	t.Run(testString(testContext, "DecryptPublic/Sin/"), func(t *testing.T) {

		if testContext.params.PCount() == 0 {
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

		if ciphertext, err = eval.EvaluateCheby(ciphertext, cheby, ciphertext.Scale); err != nil {
			t.Error(err)
		}

		plaintext := testContext.decryptor.DecryptNew(ciphertext)

		valuesHave := testContext.encoder.Decode(plaintext, testContext.params.LogSlots())

		verifyTestVectors(testContext, nil, values, valuesHave, testContext.params.LogSlots(), 0, t)

		sigma := testContext.encoder.GetErrSTDCoeffDomain(values, valuesHave, plaintext.Scale)

		valuesHave = testContext.encoder.DecodePublic(plaintext, testContext.params.LogSlots(), sigma)

		verifyTestVectors(testContext, nil, values, valuesHave, testContext.params.LogSlots(), 0, t)
	})
}

func testSwitchKeys(testContext *testParams, t *testing.T) {

	var sk2 *rlwe.SecretKey
	var decryptorSk2 *Decryptor
	var switchingKey *rlwe.SwitchingKey

	if testContext.params.PCount() != 0 {
		sk2 = testContext.kgen.GenSecretKey()
		decryptorSk2 = NewDecryptor(testContext.params, sk2)
		switchingKey = testContext.kgen.GenSwitchingKey(testContext.sk, sk2)
	}

	t.Run(testString(testContext, "SwitchKeys/"), func(t *testing.T) {

		if testContext.params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		testContext.evaluator.SwitchKeys(ciphertext, switchingKey, ciphertext)

		verifyTestVectors(testContext, decryptorSk2, values, ciphertext, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "SwitchKeysNew/"), func(t *testing.T) {

		if testContext.params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		ciphertext = testContext.evaluator.SwitchKeysNew(ciphertext, switchingKey)

		verifyTestVectors(testContext, decryptorSk2, values, ciphertext, testContext.params.LogSlots(), 0, t)
	})

}

func testAutomorphisms(testContext *testParams, t *testing.T) {

	if testContext.params.PCount() == 0 {
		t.Skip("#Pi is empty")
	}
	rots := []int{0, 1, -1, 4, -4, 63, -63}
	rotKey := testContext.kgen.GenRotationKeysForRotations(rots, true, testContext.sk)
	evaluator := testContext.evaluator.WithKey(rlwe.EvaluationKey{Rlk: testContext.rlk, Rtks: rotKey})

	t.Run(testString(testContext, "Conjugate/"), func(t *testing.T) {

		if testContext.params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values {
			values[i] = complex(real(values[i]), -imag(values[i]))
		}

		evaluator.Conjugate(ciphertext, ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "ConjugateNew/"), func(t *testing.T) {

		if testContext.params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for i := range values {
			values[i] = complex(real(values[i]), -imag(values[i]))
		}

		ciphertext = evaluator.ConjugateNew(ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "Rotate/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		ciphertext2 := NewCiphertext(testContext.params, ciphertext1.Degree(), ciphertext1.Level(), ciphertext1.Scale)

		for _, n := range rots {
			evaluator.Rotate(ciphertext1, n, ciphertext2)
			verifyTestVectors(testContext, testContext.decryptor, utils.RotateComplex128Slice(values1, n), ciphertext2, testContext.params.LogSlots(), 0, t)
		}
	})

	t.Run(testString(testContext, "RotateNew/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		for _, n := range rots {
			verifyTestVectors(testContext, testContext.decryptor, utils.RotateComplex128Slice(values1, n), evaluator.RotateNew(ciphertext1, n), testContext.params.LogSlots(), 0, t)
		}

	})

	t.Run(testString(testContext, "RotateHoisted/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		ciphertexts := evaluator.RotateHoisted(ciphertext1, rots)

		for _, n := range rots {
			verifyTestVectors(testContext, testContext.decryptor, utils.RotateComplex128Slice(values1, n), ciphertexts[n], testContext.params.LogSlots(), 0, t)
		}
	})
}

func testInnerSum(testContext *testParams, t *testing.T) {

	if testContext.params.PCount() == 0 {
		t.Skip("#Pi is empty")
	}

	t.Run(testString(testContext, "InnerSum/"), func(t *testing.T) {
		batch := 2
		n := 35

		rotKey := testContext.kgen.GenRotationKeysForRotations(testContext.params.RotationsForInnerSum(batch, n), false, testContext.sk)
		eval := testContext.evaluator.WithKey(rlwe.EvaluationKey{Rlk: testContext.rlk, Rtks: rotKey})

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		eval.InnerSum(ciphertext1, batch, n, ciphertext1)

		tmp0 := make([]complex128, len(values1))
		copy(tmp0, values1)

		for i := 1; i < n; i++ {

			tmp1 := utils.RotateComplex128Slice(tmp0, i*batch)

			for j := range values1 {
				values1[j] += tmp1[j]
			}
		}

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "InnerSumLog/"), func(t *testing.T) {

		batch := 3
		n := 15

		rotKey := testContext.kgen.GenRotationKeysForRotations(testContext.params.RotationsForInnerSumLog(batch, n), false, testContext.sk)
		eval := testContext.evaluator.WithKey(rlwe.EvaluationKey{Rlk: testContext.rlk, Rtks: rotKey})

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		eval.InnerSumLog(ciphertext1, batch, n, ciphertext1)

		tmp0 := make([]complex128, len(values1))
		copy(tmp0, values1)

		for i := 1; i < n; i++ {

			tmp1 := utils.RotateComplex128Slice(tmp0, i*batch)

			for j := range values1 {
				values1[j] += tmp1[j]
			}
		}

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, testContext.params.LogSlots(), 0, t)

	})
}

func testReplicate(testContext *testParams, t *testing.T) {

	if testContext.params.PCount() == 0 {
		t.Skip("#Pi is empty")
	}

	t.Run(testString(testContext, "Replicate/"), func(t *testing.T) {
		batch := 2
		n := 35

		rotKey := testContext.kgen.GenRotationKeysForRotations(testContext.params.RotationsForReplicate(batch, n), false, testContext.sk)
		eval := testContext.evaluator.WithKey(rlwe.EvaluationKey{Rlk: testContext.rlk, Rtks: rotKey})

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		eval.Replicate(ciphertext1, batch, n, ciphertext1)

		tmp0 := make([]complex128, len(values1))
		copy(tmp0, values1)

		for i := 1; i < n; i++ {

			tmp1 := utils.RotateComplex128Slice(tmp0, i*-batch)

			for j := range values1 {
				values1[j] += tmp1[j]
			}
		}

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "ReplicateLog/"), func(t *testing.T) {

		batch := 3
		n := 15

		rotKey := testContext.kgen.GenRotationKeysForRotations(testContext.params.RotationsForReplicateLog(batch, n), false, testContext.sk)
		eval := testContext.evaluator.WithKey(rlwe.EvaluationKey{Rlk: testContext.rlk, Rtks: rotKey})

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		eval.ReplicateLog(ciphertext1, batch, n, ciphertext1)

		tmp0 := make([]complex128, len(values1))
		copy(tmp0, values1)

		for i := 1; i < n; i++ {

			tmp1 := utils.RotateComplex128Slice(tmp0, i*-batch)

			for j := range values1 {
				values1[j] += tmp1[j]
			}
		}

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, testContext.params.LogSlots(), 0, t)

	})
}

func testLinearTransform(testContext *testParams, t *testing.T) {

	if testContext.params.PCount() == 0 {
		t.Skip("#Pi is empty")
	}

	t.Run(testString(testContext, "LinearTransform/BSGS/"), func(t *testing.T) {

		params := testContext.params

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		diagMatrix := make(map[int][]complex128)

		diagMatrix[-15] = make([]complex128, params.Slots())
		diagMatrix[-4] = make([]complex128, params.Slots())
		diagMatrix[-1] = make([]complex128, params.Slots())
		diagMatrix[0] = make([]complex128, params.Slots())
		diagMatrix[1] = make([]complex128, params.Slots())
		diagMatrix[4] = make([]complex128, params.Slots())
		diagMatrix[15] = make([]complex128, params.Slots())

		for i := 0; i < params.Slots(); i++ {
			diagMatrix[-15][i] = complex(1, 0)
			diagMatrix[-4][i] = complex(1, 0)
			diagMatrix[-1][i] = complex(1, 0)
			diagMatrix[0][i] = complex(1, 0)
			diagMatrix[1][i] = complex(1, 0)
			diagMatrix[4][i] = complex(1, 0)
			diagMatrix[15][i] = complex(1, 0)
		}

		ptDiagMatrix := testContext.encoder.EncodeDiagMatrixBSGSAtLvl(params.MaxLevel(), diagMatrix, params.Scale(), 1.0, params.LogSlots())

		rots := testContext.params.RotationsForDiagMatrixMult(ptDiagMatrix)

		rotKey := testContext.kgen.GenRotationKeysForRotations(rots, false, testContext.sk)

		eval := testContext.evaluator.WithKey(rlwe.EvaluationKey{Rlk: testContext.rlk, Rtks: rotKey})

		res := eval.LinearTransform(ciphertext1, ptDiagMatrix)[0]

		tmp := make([]complex128, params.Slots())
		copy(tmp, values1)

		for i := 0; i < params.Slots(); i++ {
			values1[i] += tmp[(i-15+params.Slots())%params.Slots()]
			values1[i] += tmp[(i-4+params.Slots())%params.Slots()]
			values1[i] += tmp[(i-1+params.Slots())%params.Slots()]
			values1[i] += tmp[(i+1)%params.Slots()]
			values1[i] += tmp[(i+4)%params.Slots()]
			values1[i] += tmp[(i+15)%params.Slots()]
		}

		verifyTestVectors(testContext, testContext.decryptor, values1, res, testContext.params.LogSlots(), 0, t)
	})

	t.Run(testString(testContext, "LinearTransform/Naive/"), func(t *testing.T) {

		params := testContext.params

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, complex(-1, -1), complex(1, 1), t)

		diagMatrix := make(map[int][]complex128)

		diagMatrix[-1] = make([]complex128, params.Slots())
		diagMatrix[0] = make([]complex128, params.Slots())

		for i := 0; i < params.Slots(); i++ {
			diagMatrix[-1][i] = complex(1, 0)
			diagMatrix[0][i] = complex(1, 0)
		}

		ptDiagMatrix := testContext.encoder.EncodeDiagMatrixAtLvl(params.MaxLevel(), diagMatrix, params.Scale(), params.LogSlots())

		rots := testContext.params.RotationsForDiagMatrixMult(ptDiagMatrix)

		rotKey := testContext.kgen.GenRotationKeysForRotations(rots, false, testContext.sk)

		eval := testContext.evaluator.WithKey(rlwe.EvaluationKey{Rlk: testContext.rlk, Rtks: rotKey})

		res := eval.LinearTransform(ciphertext1, ptDiagMatrix)[0]

		tmp := make([]complex128, params.Slots())
		copy(tmp, values1)

		for i := 0; i < params.Slots(); i++ {
			values1[i] += tmp[(i-1+params.Slots())%params.Slots()]
		}

		verifyTestVectors(testContext, testContext.decryptor, values1, res, testContext.params.LogSlots(), 0, t)
	})
}

func testMarshaller(testctx *testParams, t *testing.T) {

	ringQ := testctx.ringQ
	ringP := testctx.ringP

	t.Run("Marshaller/Parameters/Binary", func(t *testing.T) {
		bytes, err := testctx.params.MarshalBinary()
		assert.Nil(t, err)
		var p Parameters
		err = p.UnmarshalBinary(bytes)
		assert.Nil(t, err)
		assert.Equal(t, testctx.params, p)
		assert.Equal(t, testctx.params.RingQ(), p.RingQ())
	})

	t.Run("Marshaller/Parameters/JSON", func(t *testing.T) {
		// checks that parameters can be marshalled without error
		data, err := json.Marshal(testctx.params)
		assert.Nil(t, err)
		assert.NotNil(t, data)

		// checks that bfv.Parameters can be unmarshalled without error
		var paramsRec Parameters
		err = json.Unmarshal(data, &paramsRec)
		assert.Nil(t, err)
		assert.True(t, testctx.params.Equals(paramsRec))

		// checks that rlwe.Parameters can be unmarshalled without error
		var rlweParams rlwe.Parameters
		err = json.Unmarshal(data, &rlweParams)
		assert.Nil(t, err)
		assert.True(t, testctx.params.Parameters.Equals(rlweParams))

		// checks that bfv.Paramters can be unmarshalled with log-moduli definition without error
		dataWithLogModuli := []byte(fmt.Sprintf(`{"LogN":%d,"LogQ":[50,50],"LogP":[60],"Sigma":3.2,"T":65537}`, testctx.params.LogN()))
		var paramsWithLogModuli Parameters
		err = json.Unmarshal(dataWithLogModuli, &paramsWithLogModuli)
		assert.Nil(t, err)
		assert.Equal(t, 2, paramsWithLogModuli.QCount())
		assert.Equal(t, 1, paramsWithLogModuli.PCount())

		// checks that bfv.Paramters can be unmarshalled with log-moduli definition with empty P without error
		dataWithLogModuliNoP := []byte(fmt.Sprintf(`{"LogN":%d,"LogQ":[50,50],"LogP":[],"Sigma":3.2,"T":65537}`, testctx.params.LogN()))
		var paramsWithLogModuliNoP Parameters
		err = json.Unmarshal(dataWithLogModuliNoP, &paramsWithLogModuliNoP)
		assert.Nil(t, err)
		assert.Equal(t, 2, paramsWithLogModuliNoP.QCount())
		assert.Equal(t, 0, paramsWithLogModuliNoP.PCount())
	})

	t.Run("Marshaller/Ciphertext/", func(t *testing.T) {
		t.Run(testString(testctx, "EndToEnd/"), func(t *testing.T) {

			ciphertextWant := NewCiphertextRandom(testctx.prng, testctx.params, 2, testctx.params.MaxLevel(), testctx.params.Scale())

			marshalledCiphertext, err := ciphertextWant.MarshalBinary()
			require.NoError(t, err)

			ciphertextTest := new(Ciphertext)
			require.NoError(t, ciphertextTest.UnmarshalBinary(marshalledCiphertext))

			require.Equal(t, ciphertextWant.Degree(), ciphertextTest.Degree())
			require.Equal(t, ciphertextWant.Level(), ciphertextTest.Level())
			require.Equal(t, ciphertextWant.Scale, ciphertextTest.Scale)

			for i := range ciphertextWant.Value {
				require.True(t, testctx.ringQ.EqualLvl(ciphertextWant.Level(), ciphertextWant.Value[i], ciphertextTest.Value[i]))
			}
		})

		t.Run(testString(testctx, "Minimal/"), func(t *testing.T) {

			ciphertext := NewCiphertextRandom(testctx.prng, testctx.params, 0, testctx.params.MaxLevel(), testctx.params.Scale())

			marshalledCiphertext, err := ciphertext.MarshalBinary()
			require.NoError(t, err)

			ciphertextTest := new(Ciphertext)
			require.Error(t, ciphertextTest.UnmarshalBinary(nil))
			require.NoError(t, ciphertextTest.UnmarshalBinary(marshalledCiphertext))

			require.Equal(t, ciphertext.Degree(), 0)
			require.Equal(t, ciphertext.Level(), testctx.params.MaxLevel())
			require.Equal(t, ciphertext.Scale, testctx.params.Scale())
			require.Equal(t, len(ciphertext.Value), 1)
		})
	})

	t.Run(testString(testctx, "Marshaller/Sk/"), func(t *testing.T) {

		marshalledSk, err := testctx.sk.MarshalBinary()
		require.NoError(t, err)

		sk := new(rlwe.SecretKey)
		err = sk.UnmarshalBinary(marshalledSk)
		require.NoError(t, err)

		require.True(t, ringQ.Equal(sk.Value[0], testctx.sk.Value[0]))
		require.True(t, ringP.Equal(sk.Value[1], testctx.sk.Value[1]))

	})

	t.Run(testString(testctx, "Marshaller/Pk/"), func(t *testing.T) {

		marshalledPk, err := testctx.pk.MarshalBinary()
		require.NoError(t, err)

		pk := new(rlwe.PublicKey)
		err = pk.UnmarshalBinary(marshalledPk)
		require.NoError(t, err)

		for k := range testctx.pk.Value {
			require.Truef(t, ringQ.Equal(pk.Value[k][0], testctx.pk.Value[k][0]), "Marshal PublicKey element [%d][0]", k)
			require.Truef(t, ringP.Equal(pk.Value[k][1], testctx.pk.Value[k][1]), "Marshal PublicKey element [%d][1]", k)
		}
	})

	t.Run(testString(testctx, "Marshaller/EvaluationKey/"), func(t *testing.T) {

		if testctx.params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		evalKey := testctx.kgen.GenRelinearizationKey(testctx.sk)
		data, err := evalKey.MarshalBinary()
		require.NoError(t, err)

		resEvalKey := new(rlwe.RelinearizationKey)
		err = resEvalKey.UnmarshalBinary(data)
		require.NoError(t, err)

		evakeyWant := evalKey.Keys[0].Value
		evakeyTest := resEvalKey.Keys[0].Value

		for j := range evakeyWant {
			for k := range evakeyWant[j] {
				require.Truef(t, ringQ.Equal(evakeyWant[j][k][0], evakeyTest[j][k][0]), "Marshal EvaluationKey element [%d][%d][0]", j, k)
				require.Truef(t, ringP.Equal(evakeyWant[j][k][1], evakeyTest[j][k][1]), "Marshal EvaluationKey element [%d][%d][1]", j, k)
			}
		}
	})

	t.Run(testString(testctx, "Marshaller/SwitchingKey/"), func(t *testing.T) {

		if testctx.params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		skOut := testctx.kgen.GenSecretKey()

		switchingKey := testctx.kgen.GenSwitchingKey(testctx.sk, skOut)
		data, err := switchingKey.MarshalBinary()
		require.NoError(t, err)

		resSwitchingKey := new(rlwe.SwitchingKey)
		err = resSwitchingKey.UnmarshalBinary(data)
		require.NoError(t, err)

		evakeyWant := switchingKey.Value
		evakeyTest := resSwitchingKey.Value

		for j := range evakeyWant {
			for k := range evakeyWant[j] {
				require.True(t, ringQ.Equal(evakeyWant[j][k][0], evakeyTest[j][k][0]))
				require.True(t, ringP.Equal(evakeyWant[j][k][1], evakeyTest[j][k][1]))
			}
		}
	})

	t.Run(testString(testctx, "Marshaller/RotationKey/"), func(t *testing.T) {

		if testctx.params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		rots := []int{1, -1, 63, -63}
		galEls := []uint64{testctx.params.GaloisElementForRowRotation()}
		for _, n := range rots {
			galEls = append(galEls, testctx.params.GaloisElementForColumnRotationBy(n))
		}

		rotationKey := testctx.kgen.GenRotationKeys(galEls, testctx.sk)

		data, err := rotationKey.MarshalBinary()
		require.NoError(t, err)

		resRotationKey := new(rlwe.RotationKeySet)
		err = resRotationKey.UnmarshalBinary(data)
		require.NoError(t, err)

		for _, galEl := range galEls {

			evakeyWant := rotationKey.Keys[galEl].Value
			evakeyTest := resRotationKey.Keys[galEl].Value

			for j := range evakeyWant {
				for k := range evakeyWant[j] {
					require.Truef(t, ringQ.Equal(evakeyWant[j][k][0], evakeyTest[j][k][0]), "Marshal RotationKey RotateLeft %d element [%d][%d][0]", galEl, j, k)
					require.Truef(t, ringP.Equal(evakeyWant[j][k][1], evakeyTest[j][k][1]), "Marshal RotationKey RotateLeft %d element [%d][%d][1]", galEl, j, k)
				}
			}
		}
	})
}
