package bfv

import (
	"encoding/json"
	"flag"
	"fmt"
	"runtime"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters). Overrides -short and requires -timeout=0.")
var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")

func testString(opname string, p Parameters) string {
	return fmt.Sprintf("%sLogN=%d/logQ=%d/alpha=%d/beta=%d", opname, p.LogN(), p.LogQP(), p.PCount(), p.Beta())
}

type testContext struct {
	params      Parameters
	ringQ       *ring.Ring
	ringQP      *ring.Ring
	ringT       *ring.Ring
	prng        utils.PRNG
	uSampler    *ring.UniformSampler
	encoder     Encoder
	kgen        rlwe.KeyGenerator
	sk          *rlwe.SecretKey
	pk          *rlwe.PublicKey
	rlk         *rlwe.RelinearizationKey
	encryptorPk Encryptor
	encryptorSk Encryptor
	decryptor   Decryptor
	evaluator   Evaluator
}

func TestBFV(t *testing.T) {

	defaultParams := DefaultParams // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15
	if testing.Short() {
		defaultParams = DefaultParams[:2] // the short test suite runs for ring degree N=2^12, 2^13
	}
	if *flagLongTest {
		defaultParams = append(defaultParams, DefaultPostQuantumParams...) // the long test suite runs for all default parameters
	}
	if *flagParamString != "" {
		var jsonParams ParametersLiteral
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range defaultParams {

		params, err := NewParametersFromLiteral(p)
		if err != nil {
			panic(err)
		}
		var testctx *testContext
		if testctx, err = genTestParams(params); err != nil {
			panic(err)
		}

		for _, testSet := range []func(testctx *testContext, t *testing.T){
			testParameters,
			testEncoder,
			testEvaluator,
			testEvaluatorKeySwitch,
			testEvaluatorRotate,
			testMarshaller,
		} {
			testSet(testctx, t)
			runtime.GC()
		}
	}
}

func genTestParams(params Parameters) (testctx *testContext, err error) {

	testctx = new(testContext)
	testctx.params = params

	if testctx.prng, err = utils.NewPRNG(); err != nil {
		return nil, err
	}

	testctx.ringQ = params.RingQ()
	testctx.ringQP = params.RingQP()
	testctx.ringT = params.RingT()

	testctx.uSampler = ring.NewUniformSampler(testctx.prng, testctx.ringT)
	testctx.kgen = NewKeyGenerator(testctx.params)
	testctx.sk, testctx.pk = testctx.kgen.GenKeyPair()
	if params.PCount() != 0 {
		testctx.rlk = testctx.kgen.GenRelinearizationKey(testctx.sk, 1)
	}

	testctx.encoder = NewEncoder(testctx.params)
	testctx.encryptorPk = NewEncryptor(testctx.params, testctx.pk)
	testctx.encryptorSk = NewEncryptor(testctx.params, testctx.sk)
	testctx.decryptor = NewDecryptor(testctx.params, testctx.sk)
	testctx.evaluator = NewEvaluator(testctx.params, rlwe.EvaluationKey{Rlk: testctx.rlk})
	return

}

func testParameters(testctx *testContext, t *testing.T) {

	t.Run("Parameters/InverseGaloisElement/", func(t *testing.T) {
		for i := 1; i < int(testctx.params.N()/2); i++ {
			galEl := testctx.params.GaloisElementForColumnRotationBy(i)
			mod := uint64(2 * testctx.params.N())
			inv := testctx.params.InverseGaloisElement(galEl)
			res := (inv * galEl) % mod
			assert.Equal(t, uint64(1), res)
		}
	})

	t.Run(testString("Parameters/CopyNew/", testctx.params), func(t *testing.T) {
		params1, params2 := testctx.params.CopyNew(), testctx.params.CopyNew()
		assert.True(t, params1.Equals(testctx.params) && params2.Equals(testctx.params))
		params1.ringT, _ = ring.NewRing(testctx.params.N(), []uint64{7})
		assert.False(t, params1.Equals(testctx.params))
		assert.True(t, params2.Equals(testctx.params))
	})
}

func newTestVectorsRingQ(testctx *testContext, encryptor Encryptor, t *testing.T) (coeffs *ring.Poly, plaintext *Plaintext, ciphertext *Ciphertext) {

	coeffs = testctx.uSampler.ReadNew()

	plaintext = NewPlaintext(testctx.params)

	testctx.encoder.EncodeUint(coeffs.Coeffs[0], plaintext)

	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return coeffs, plaintext, ciphertext
}

func newTestVectorsRingT(testctx *testContext, t *testing.T) (coeffs *ring.Poly, plaintext *PlaintextRingT) {

	coeffs = testctx.uSampler.ReadNew()

	plaintext = NewPlaintextRingT(testctx.params)

	testctx.encoder.EncodeUintRingT(coeffs.Coeffs[0], plaintext)

	return coeffs, plaintext
}

func newTestVectorsMul(testctx *testContext, t *testing.T) (coeffs *ring.Poly, plaintext *PlaintextMul) {

	coeffs = testctx.uSampler.ReadNew()

	plaintext = NewPlaintextMul(testctx.params)

	testctx.encoder.EncodeUintMul(coeffs.Coeffs[0], plaintext)

	return coeffs, plaintext
}

func verifyTestVectors(testctx *testContext, decryptor Decryptor, coeffs *ring.Poly, element Operand, t *testing.T) {

	var coeffsTest []uint64

	switch el := element.(type) {
	case *Plaintext, *PlaintextMul, *PlaintextRingT:
		coeffsTest = testctx.encoder.DecodeUintNew(el)
	case *Ciphertext:
		coeffsTest = testctx.encoder.DecodeUintNew(decryptor.DecryptNew(el))
	default:
		t.Error("invalid test object to verify")
	}

	require.True(t, utils.EqualSliceUint64(coeffs.Coeffs[0], coeffsTest))
}

func testEncoder(testctx *testContext, t *testing.T) {
	t.Run(testString("Encoder/Encode&Decode/RingT/Uint/", testctx.params), func(t *testing.T) {
		values, plaintext := newTestVectorsRingT(testctx, t)
		verifyTestVectors(testctx, nil, values, plaintext, t)
	})

	t.Run(testString("Encoder/Encode&Decode/RingT/Int/", testctx.params), func(t *testing.T) {

		T := testctx.params.T()
		THalf := T >> 1
		coeffs := testctx.uSampler.ReadNew()
		coeffsInt := make([]int64, len(coeffs.Coeffs[0]))
		for i, c := range coeffs.Coeffs[0] {
			c %= T
			if c >= THalf {
				coeffsInt[i] = -int64(T - c)
			} else {
				coeffsInt[i] = int64(c)
			}
		}
		plaintext := NewPlaintextRingT(testctx.params)
		testctx.encoder.EncodeIntRingT(coeffsInt, plaintext)
		coeffsTest := testctx.encoder.DecodeIntNew(plaintext)

		require.True(t, utils.EqualSliceInt64(coeffsInt, coeffsTest))
	})

	t.Run(testString("Encoder/Encode&Decode/RingQ/Uint/", testctx.params), func(t *testing.T) {
		values, plaintext, _ := newTestVectorsRingQ(testctx, nil, t)
		verifyTestVectors(testctx, nil, values, plaintext, t)
	})

	t.Run(testString("Encoder/Encode&Decode/RingQ/Int/", testctx.params), func(t *testing.T) {

		T := testctx.params.T()
		THalf := T >> 1
		coeffs := testctx.uSampler.ReadNew()
		coeffsInt := make([]int64, len(coeffs.Coeffs[0]))
		for i, c := range coeffs.Coeffs[0] {
			c %= T
			if c >= THalf {
				coeffsInt[i] = -int64(T - c)
			} else {
				coeffsInt[i] = int64(c)
			}
		}
		plaintext := NewPlaintext(testctx.params)
		testctx.encoder.EncodeInt(coeffsInt, plaintext)
		coeffsTest := testctx.encoder.DecodeIntNew(plaintext)

		require.True(t, utils.EqualSliceInt64(coeffsInt, coeffsTest))
	})

	t.Run(testString("Encoder/Encode&Decode/PlaintextMul/", testctx.params), func(t *testing.T) {
		values, plaintext := newTestVectorsMul(testctx, t)
		verifyTestVectors(testctx, nil, values, plaintext, t)
	})
}

func testEvaluator(testctx *testContext, t *testing.T) {

	t.Run(testString("Evaluator/Add/op1=Ciphertext/op2=Ciphertext/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		testctx.evaluator.Add(ciphertext1, ciphertext2, ciphertext1)
		testctx.ringT.Add(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/AddNoMod/op1=Ciphertext/op2=Ciphertext/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		testctx.evaluator.AddNoMod(ciphertext1, ciphertext2, ciphertext1)
		testctx.evaluator.Reduce(ciphertext1, ciphertext1)
		testctx.ringT.Add(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/AddNew/op1=Ciphertext/op2=Ciphertext/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		ciphertext1 = testctx.evaluator.AddNew(ciphertext1, ciphertext2)
		testctx.ringT.Add(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/AddNoModNew/op1=Ciphertext/op2=Ciphertext/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		ciphertext1 = testctx.evaluator.AddNoModNew(ciphertext1, ciphertext2)
		ciphertext1 = testctx.evaluator.ReduceNew(ciphertext1)
		testctx.ringT.Add(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/Add/op1=Ciphertext/op2=PlaintextRingT/", testctx.params), func(t *testing.T) {

		values1, plaintextRingT := newTestVectorsRingT(testctx, t)
		values2, _, ciphertext := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		ciphertextOut := NewCiphertext(testctx.params, 1)

		testctx.evaluator.Add(ciphertext, plaintextRingT, ciphertextOut)
		testctx.ringT.Add(values1, values2, values2)

		verifyTestVectors(testctx, testctx.decryptor, values2, ciphertextOut, t)

		testctx.evaluator.Add(plaintextRingT, ciphertext, ciphertextOut)

		verifyTestVectors(testctx, testctx.decryptor, values2, ciphertextOut, t)
	})

	t.Run(testString("Evaluator/Add/op1=Ciphertext/op2=Plaintext/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, plaintext2, ciphertext2 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		testctx.evaluator.Add(ciphertext1, plaintext2, ciphertext2)
		testctx.ringT.Add(values1, values2, values2)

		verifyTestVectors(testctx, testctx.decryptor, values2, ciphertext2, t)

		testctx.evaluator.Add(plaintext2, ciphertext1, ciphertext2)

		verifyTestVectors(testctx, testctx.decryptor, values2, ciphertext2, t)
	})

	t.Run(testString("Evaluator/Sub/op1=Ciphertext/op2=Ciphertext/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		testctx.evaluator.Sub(ciphertext1, ciphertext2, ciphertext1)
		testctx.ringT.Sub(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/SubNoMod/op1=Ciphertext/op2=Ciphertext/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		testctx.evaluator.SubNoMod(ciphertext1, ciphertext2, ciphertext1)
		testctx.evaluator.Reduce(ciphertext1, ciphertext1)
		testctx.ringT.Sub(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/SubNew/op1=Ciphertext/op2=Ciphertext/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		ciphertext1 = testctx.evaluator.SubNew(ciphertext1, ciphertext2)
		testctx.ringT.Sub(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/Sub/op1=Ciphertext/op2=PlaintextRingT/", testctx.params), func(t *testing.T) {

		values1, plaintextRingT := newTestVectorsRingT(testctx, t)
		values2, _, ciphertext := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		ciphertextOut := NewCiphertext(testctx.params, 1)
		plaintextWant := NewPlaintextRingT(testctx.params)

		testctx.evaluator.Sub(ciphertext, plaintextRingT, ciphertextOut)
		testctx.ringT.Sub(values2, values1, plaintextWant.Value)
		verifyTestVectors(testctx, testctx.decryptor, plaintextWant.Value, ciphertextOut, t)

		testctx.evaluator.Sub(plaintextRingT, ciphertext, ciphertextOut)
		testctx.ringT.Sub(values1, values2, plaintextWant.Value)
		verifyTestVectors(testctx, testctx.decryptor, plaintextWant.Value, ciphertextOut, t)
	})

	t.Run(testString("Evaluator/Sub/op1=Ciphertext/op2=Plaintext/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, plaintext2, ciphertext2 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		valuesWant := testctx.ringT.NewPoly()

		testctx.evaluator.Sub(ciphertext1, plaintext2, ciphertext2)
		testctx.ringT.Sub(values1, values2, valuesWant)
		verifyTestVectors(testctx, testctx.decryptor, valuesWant, ciphertext2, t)

		testctx.evaluator.Sub(plaintext2, ciphertext1, ciphertext2)
		testctx.ringT.Sub(values2, values1, valuesWant)
		verifyTestVectors(testctx, testctx.decryptor, valuesWant, ciphertext2, t)
	})

	t.Run(testString("Evaluator/SubNoModNew/op1=Ciphertext/op2=Ciphertext/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		ciphertext1 = testctx.evaluator.SubNoModNew(ciphertext1, ciphertext2)
		ciphertext1 = testctx.evaluator.ReduceNew(ciphertext1)
		testctx.ringT.Sub(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/SubNoMod/op1=Ciphertext/op2=Plaintext/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, plaintext2, ciphertext2 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		valuesWant := testctx.ringT.NewPoly()

		testctx.evaluator.SubNoMod(ciphertext1, plaintext2, ciphertext2)
		testctx.evaluator.Reduce(ciphertext2, ciphertext2)
		testctx.ringT.Sub(values1, values2, valuesWant)
		verifyTestVectors(testctx, testctx.decryptor, valuesWant, ciphertext2, t)

		testctx.evaluator.SubNoMod(plaintext2, ciphertext1, ciphertext2)
		testctx.evaluator.Reduce(ciphertext2, ciphertext2)
		testctx.ringT.Sub(values2, values1, valuesWant)
		verifyTestVectors(testctx, testctx.decryptor, valuesWant, ciphertext2, t)
	})

	t.Run(testString("Evaluator/Neg/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		testctx.evaluator.Neg(ciphertext1, ciphertext1)
		testctx.ringT.Neg(values1, values1)
		testctx.ringT.Reduce(values1, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/NegNew/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		ciphertext1 = testctx.evaluator.NegNew(ciphertext1)
		testctx.ringT.Neg(values1, values1)
		testctx.ringT.Reduce(values1, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/MulScalar/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		testctx.evaluator.MulScalar(ciphertext1, 37, ciphertext1)
		testctx.ringT.MulScalar(values1, 37, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/MulScalarNew/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		ciphertext1 = testctx.evaluator.MulScalarNew(ciphertext1, 37)
		testctx.ringT.MulScalar(values1, 37, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=Ciphertext/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		receiver := NewCiphertext(testctx.params, ciphertext1.Degree()+ciphertext2.Degree())
		testctx.evaluator.Mul(ciphertext1, ciphertext2, receiver)
		testctx.ringT.MulCoeffs(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, receiver, t)
	})

	t.Run(testString("Evaluator/MulNew/op1=Ciphertext/op2=Ciphertext/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		receiver := testctx.evaluator.MulNew(ciphertext1, ciphertext2)
		testctx.ringT.MulCoeffs(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, receiver, t)
	})

	t.Run(testString("Evaluator/MulSquare/op1=Ciphertext/op2=Ciphertext/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		receiver := NewCiphertext(testctx.params, ciphertext1.Degree()+ciphertext1.Degree())
		testctx.evaluator.Mul(ciphertext1, ciphertext1, receiver)
		testctx.ringT.MulCoeffs(values1, values1, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, receiver, t)

		if testctx.params.LogN() < 13 {
			t.Skip()
		}
		receiver2 := NewCiphertext(testctx.params, receiver.Degree()+receiver.Degree())
		testctx.evaluator.Mul(receiver, receiver, receiver2)
		testctx.ringT.MulCoeffs(values1, values1, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, receiver2, t)
	})

	t.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=Plaintext/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, plaintext2, _ := newTestVectorsRingQ(testctx, nil, t)

		testctx.evaluator.Mul(ciphertext1, plaintext2, ciphertext1)
		testctx.ringT.MulCoeffs(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=PlaintextRingT/", testctx.params), func(t *testing.T) {

		values1, plaintextRingT := newTestVectorsRingT(testctx, t)
		values2, _, ciphertext := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		ciphertextOut := NewCiphertext(testctx.params, 1)

		testctx.evaluator.Mul(ciphertext, plaintextRingT, ciphertextOut)
		testctx.ringT.MulCoeffs(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertextOut, t)
	})

	t.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=PlaintextMul/", testctx.params), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, plaintext2 := newTestVectorsMul(testctx, t)

		testctx.evaluator.Mul(ciphertext1, plaintext2, ciphertext1)
		testctx.ringT.MulCoeffs(values1, values2, values1)

		verifyTestVectors(testctx, testctx.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString("Evaluator/Mul/Relinearize/", testctx.params), func(t *testing.T) {

		if testctx.params.PCount() == 0 {
			t.Skip("#Pi is empty")
		}

		values1, _, ciphertext1 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		values2, _, ciphertext2 := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		receiver := NewCiphertext(testctx.params, ciphertext1.Degree()+ciphertext2.Degree())
		testctx.evaluator.Mul(ciphertext1, ciphertext2, receiver)
		testctx.ringT.MulCoeffs(values1, values2, values1)

		receiver2 := testctx.evaluator.RelinearizeNew(receiver)
		verifyTestVectors(testctx, testctx.decryptor, values1, receiver2, t)

		testctx.evaluator.Relinearize(receiver, receiver)
		verifyTestVectors(testctx, testctx.decryptor, values1, receiver, t)
	})
}

func testEvaluatorKeySwitch(testctx *testContext, t *testing.T) {

	if testctx.params.PCount() == 0 {
		t.Skip("#Pi is empty")
	}

	sk2 := testctx.kgen.GenSecretKey()
	decryptorSk2 := NewDecryptor(testctx.params, sk2)
	switchKey := testctx.kgen.GenSwitchingKey(testctx.sk, sk2)

	t.Run(testString("Evaluator/KeySwitch/InPlace/", testctx.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		testctx.evaluator.SwitchKeys(ciphertext, switchKey, ciphertext)
		verifyTestVectors(testctx, decryptorSk2, values, ciphertext, t)
	})

	t.Run(testString("Evaluator/KeySwitch/New/", testctx.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		ciphertext = testctx.evaluator.SwitchKeysNew(ciphertext, switchKey)
		verifyTestVectors(testctx, decryptorSk2, values, ciphertext, t)
	})
}

func testEvaluatorRotate(testctx *testContext, t *testing.T) {

	if testctx.params.PCount() == 0 {
		t.Skip("#Pi is empty")
	}

	rots := []int{1, -1, 4, -4, 63, -63}
	rotkey := testctx.kgen.GenRotationKeysForRotations(rots, true, testctx.sk)
	evaluator := testctx.evaluator.WithKey(rlwe.EvaluationKey{Rlk: testctx.rlk, Rtks: rotkey})

	t.Run(testString("Evaluator/RotateRows/", testctx.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		evaluator.RotateRows(ciphertext, ciphertext)
		values.Coeffs[0] = append(values.Coeffs[0][testctx.params.N()>>1:], values.Coeffs[0][:testctx.params.N()>>1]...)
		verifyTestVectors(testctx, testctx.decryptor, values, ciphertext, t)
	})

	t.Run(testString("Evaluator/RotateRowsNew/", testctx.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)
		ciphertext = evaluator.RotateRowsNew(ciphertext)
		values.Coeffs[0] = append(values.Coeffs[0][testctx.params.N()>>1:], values.Coeffs[0][:testctx.params.N()>>1]...)
		verifyTestVectors(testctx, testctx.decryptor, values, ciphertext, t)
	})

	t.Run(testString("Evaluator/RotateColumns", testctx.params), func(t *testing.T) {

		values, _, ciphertext := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		receiver := NewCiphertext(testctx.params, 1)
		for _, n := range rots {

			evaluator.RotateColumns(ciphertext, n, receiver)
			valuesWant := utils.RotateUint64Slots(values.Coeffs[0], n)

			verifyTestVectors(testctx, testctx.decryptor, &ring.Poly{Coeffs: [][]uint64{valuesWant}}, receiver, t)
		}
	})

	t.Run(testString("Evaluator/RotateColumnsNew", testctx.params), func(t *testing.T) {

		values, _, ciphertext := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		for _, n := range rots {

			receiver := evaluator.RotateColumnsNew(ciphertext, n)
			valuesWant := utils.RotateUint64Slots(values.Coeffs[0], n)

			verifyTestVectors(testctx, testctx.decryptor, &ring.Poly{Coeffs: [][]uint64{valuesWant}}, receiver, t)
		}
	})

	rotkey = testctx.kgen.GenRotationKeysForInnerSum(testctx.sk)
	evaluator = evaluator.WithKey(rlwe.EvaluationKey{Rlk: testctx.rlk, Rtks: rotkey})

	t.Run(testString("Evaluator/Rotate/InnerSum/", testctx.params), func(t *testing.T) {
		values, _, ciphertext := newTestVectorsRingQ(testctx, testctx.encryptorPk, t)

		evaluator.InnerSum(ciphertext, ciphertext)

		var sum uint64
		for _, c := range values.Coeffs[0] {
			sum += c
		}

		sum %= testctx.params.T()

		for i := range values.Coeffs[0] {
			values.Coeffs[0][i] = sum
		}
		verifyTestVectors(testctx, testctx.decryptor, values, ciphertext, t)
	})
}

func testMarshaller(testctx *testContext, t *testing.T) {

	t.Run("Marshaller/Parameters/Binary", func(t *testing.T) {
		bytes, err := testctx.params.MarshalBinary()
		assert.Nil(t, err)
		var p Parameters
		err = p.UnmarshalBinary(bytes)
		assert.Nil(t, err)
		assert.Equal(t, testctx.params, p)
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

	t.Run(testString("Marshaller/Ciphertext/", testctx.params), func(t *testing.T) {

		ciphertextWant := NewCiphertextRandom(testctx.prng, testctx.params, 2)

		marshalledCiphertext, err := ciphertextWant.MarshalBinary()
		require.NoError(t, err)

		ciphertextTest := new(Ciphertext)
		err = ciphertextTest.UnmarshalBinary(marshalledCiphertext)
		require.NoError(t, err)

		for i := range ciphertextWant.Value {
			require.True(t, testctx.ringQ.Equal(ciphertextWant.Value[i], ciphertextTest.Value[i]))
		}
	})
}
