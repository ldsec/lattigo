package rckks

import (
	"flag"
	"fmt"
	"math"
	"math/rand"
	"runtime"
	"testing"
	"time"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")
var minPrec float64 = 10.0

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

func TestRCKKS(t *testing.T) {

	rand.Seed(time.Now().UnixNano())

	var err error
	var testContext = new(testParams)

	var defaultParams []*Parameters

	if testing.Short() {
		defaultParams = DefaultParams[PN12QP109 : PN12QP109+2]
	} else {
		defaultParams = DefaultParams
	}

	for _, defaultParam := range defaultParams {

		if testContext, err = genTestParams(defaultParam, 0); err != nil {
			panic(err)
		}

		for _, testSet := range []func(testContext *testParams, t *testing.T){
			testNTTRCKKS,
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
			testChebyshevInterpolator,
			testSwitchKeys,
			testRotate,
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
	}

	if testContext.prng, err = utils.NewPRNG(); err != nil {
		return nil, err
	}

	testContext.encoder = NewEncoder(testContext.params)

	testContext.rlk = testContext.kgen.GenRelinKey(testContext.sk)

	testContext.encryptorPk = NewEncryptorFromPk(testContext.params, testContext.pk)
	testContext.encryptorSk = NewEncryptorFromSk(testContext.params, testContext.sk)
	testContext.decryptor = NewDecryptor(testContext.params, testContext.sk)

	testContext.evaluator = NewEvaluator(testContext.params)

	return testContext, nil

}

func newTestVectors(testContext *testParams, encryptor Encryptor, a, b float64, t *testing.T) (values []float64, plaintext *Plaintext, ciphertext *Ciphertext) {

	slots := testContext.params.Slots()

	values = make([]float64, slots)

	for i := uint64(0); i < slots; i++ {
		values[i] = randomFloat(a, b)
	}

	values[0] = 0.607538

	plaintext = NewPlaintext(testContext.params, testContext.params.MaxLevel(), testContext.params.Scale())

	testContext.encoder.EncodeNTT(plaintext, values, slots)

	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return values, plaintext, ciphertext
}

func verifyTestVectors(testContext *testParams, decryptor Decryptor, valuesWant []float64, element interface{}, t *testing.T) {

	precStats := GetPrecisionStats(testContext.params, testContext.encoder, decryptor, valuesWant, element)

	if *printPrecisionStats {
		t.Log(precStats.String())
	}

	require.GreaterOrEqual(t, precStats.MeanPrecision, minPrec)
}

func testNTTRCKKS(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "NTTRCKKS/"), func(t *testing.T) {

		ringQ := testContext.ringQ
		ringQ2N, _ := ring.NewRing(ringQ.N<<1, ringQ.Modulus)
		ringQ4NthRoot, _ := ring.NewRingWithNthRoot(ringQ.N, ringQ.N<<2, ringQ.Modulus)

		sampler := ring.NewUniformSampler(testContext.prng, ringQ)
		p1 := sampler.ReadNew()
		p2 := p1.CopyNew()

		for i, qi := range ringQ.Modulus {
			p2.Coeffs[i] = append(p2.Coeffs[i], make([]uint64, ringQ.N)...)
			p2.Coeffs[i][ringQ.N] = 0
			for j := uint64(1); j < ringQ.N; j++ {
				p2.Coeffs[i][ringQ.N*2-j] = qi - p2.Coeffs[i][j]
			}
		}

		ringQ2N.NTT(p2, p2)
		ringQ2N.MForm(p2, p2)
		ringQ2N.MulCoeffsMontgomery(p2, p2, p2)
		ringQ2N.InvMForm(p2, p2)
		ringQ2N.InvNTT(p2, p2)

		NTTRCKKS(ringQ4NthRoot, p1, p1)
		ringQ4NthRoot.MForm(p1, p1)
		ringQ4NthRoot.MulCoeffsMontgomery(p1, p1, p1)
		ringQ4NthRoot.InvMForm(p1, p1)
		InvNTTRCKKS(ringQ4NthRoot, p1, p1)

		for j := range ringQ.Modulus {
			for i := uint64(0); i < ringQ.N; i++ {
				require.Equal(t, p1.Coeffs[j][i], p2.Coeffs[j][i])
			}
		}
	})
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

	t.Run(testString(testContext, "Encoder/EncodeBatch/"), func(t *testing.T) {

		values, plaintext, _ := newTestVectors(testContext, nil, -1, 1, t)

		verifyTestVectors(testContext, testContext.decryptor, values, plaintext, t)
	})

	t.Run(testString(testContext, "Encoder/EncodeCoeffs/"), func(t *testing.T) {

		slots := testContext.params.N()

		valuesWant := make([]float64, slots)

		for i := uint64(0); i < slots; i++ {
			valuesWant[i] = randomFloat(-1, 1)
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

		require.GreaterOrEqual(t, math.Log2(1/meanprec), minPrec)
	})

}

func testEncryptor(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Encryptor/EncryptFromPk/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorPk, -1, 1, t)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t)
	})

	t.Run(testString(testContext, "Encryptor/EncryptFromPkFast/"), func(t *testing.T) {

		slots := testContext.params.Slots()

		values := make([]float64, slots)

		for i := uint64(0); i < slots; i++ {
			values[i] = 0 // randomFloat(-1, 1)
		}

		values[0] = 0 // 0.607538

		plaintext := NewPlaintext(testContext.params, testContext.params.MaxLevel(), testContext.params.Scale())

		testContext.encoder.Encode(plaintext, values, slots)

		verifyTestVectors(testContext, testContext.decryptor, values, testContext.encryptorPk.EncryptFastNew(plaintext), t)
	})

	t.Run(testString(testContext, "Encryptor/EncryptFromSk/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t)
	})

}

func testEvaluatorAdd(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Evaluator/Add/CtCtInPlace/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		testContext.evaluator.Add(ciphertext1, ciphertext2, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString(testContext, "Evaluator/Add/CtCtNew/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		ciphertext3 := testContext.evaluator.AddNew(ciphertext1, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext3, t)
	})

	t.Run(testString(testContext, "Evaluator/Add/CtPlainInPlace/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)
		values2, plaintext2, _ := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		testContext.evaluator.Add(ciphertext1, plaintext2, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		testContext.evaluator.Add(plaintext2, ciphertext1, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString(testContext, "Evaluator/Add/CtPlainInPlaceNew/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)
		values2, plaintext2, _ := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		for i := range values1 {
			values1[i] += values2[i]
		}

		ciphertext3 := testContext.evaluator.AddNew(ciphertext1, plaintext2)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext3, t)

		ciphertext3 = testContext.evaluator.AddNew(plaintext2, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext3, t)
	})

}

func testEvaluatorSub(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Evaluator/Sub/CtCtInPlace/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		for i := range values1 {
			values1[i] -= values2[i]
		}

		testContext.evaluator.Sub(ciphertext1, ciphertext2, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString(testContext, "Evaluator/Sub/CtCtNew/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		for i := range values1 {
			values1[i] -= values2[i]
		}

		ciphertext3 := testContext.evaluator.SubNew(ciphertext1, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext3, t)
	})

	t.Run(testString(testContext, "Evaluator/Sub/CtPlainInPlace/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)
		values2, plaintext2, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		valuesTest := make([]float64, len(values1))
		for i := range values1 {
			valuesTest[i] = values1[i] - values2[i]
		}

		testContext.evaluator.Sub(ciphertext1, plaintext2, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, valuesTest, ciphertext2, t)

		for i := range values1 {
			valuesTest[i] = values2[i] - values1[i]
		}

		testContext.evaluator.Sub(plaintext2, ciphertext1, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, valuesTest, ciphertext2, t)
	})

	t.Run(testString(testContext, "Evaluator/Sub/CtPlainNew/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)
		values2, plaintext2, _ := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		valuesTest := make([]float64, len(values1))
		for i := range values1 {
			valuesTest[i] = values1[i] - values2[i]
		}

		ciphertext3 := testContext.evaluator.SubNew(ciphertext1, plaintext2)

		verifyTestVectors(testContext, testContext.decryptor, valuesTest, ciphertext3, t)

		for i := range values1 {
			valuesTest[i] = values2[i] - values1[i]
		}

		ciphertext3 = testContext.evaluator.SubNew(plaintext2, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, valuesTest, ciphertext3, t)
	})

}

func testEvaluatorRescale(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Evaluator/Rescale/Single/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		constant := testContext.ringQ.Modulus[ciphertext.Level()]

		testContext.evaluator.MultByConst(ciphertext, constant, ciphertext)

		ciphertext.MulScale(float64(constant))

		testContext.evaluator.Rescale(ciphertext, testContext.params.Scale(), ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t)
	})

	t.Run(testString(testContext, "Evaluator/Rescale/Many/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

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

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t)
	})
}

func testEvaluatorAddConst(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Evaluator/AddConst/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		constant := 3.1415

		for i := range values {
			values[i] += constant
		}

		testContext.evaluator.AddConst(ciphertext, constant, ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t)
	})

}

func testEvaluatorMultByConst(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Evaluator/MultByConst/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		constant := 1.0 / 3.1415

		for i := range values {
			values[i] *= constant
		}

		testContext.evaluator.MultByConst(ciphertext, constant, ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t)
	})

}

func testEvaluatorMultByConstAndAdd(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Evaluator/MultByConstAndAdd/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		constant := 1.0 / 3.1415

		for i := range values1 {
			values2[i] += (constant * values1[i])
		}

		testContext.evaluator.MultByConstAndAdd(ciphertext1, constant, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext2, t)
	})

}

func testEvaluatorMul(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Evaluator/Mul/ct0*pt->ct0/"), func(t *testing.T) {

		values1, plaintext1, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, plaintext1, nil, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString(testContext, "Evaluator/Mul/pt*ct0->ct0/"), func(t *testing.T) {

		values1, plaintext1, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, plaintext1, nil, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString(testContext, "Evaluator/Mul/ct0*pt->ct1/"), func(t *testing.T) {

		values1, plaintext1, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		ciphertext2 := testContext.evaluator.MulRelinNew(ciphertext1, plaintext1, nil)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext2, t)
	})

	t.Run(testString(testContext, "Evaluator/Mul/ct0*ct1->ct0/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		for i := range values1 {
			values2[i] *= values1[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, ciphertext2, nil, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext1, t)
	})

	t.Run(testString(testContext, "Evaluator/Mul/ct0*ct1->ct1/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		for i := range values1 {
			values2[i] *= values1[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, ciphertext2, nil, ciphertext2)

		verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext2, t)
	})

	t.Run(testString(testContext, "Evaluator/Mul/ct0*ct1->ct2/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		for i := range values1 {
			values2[i] *= values1[i]
		}

		ciphertext3 := testContext.evaluator.MulRelinNew(ciphertext1, ciphertext2, nil)

		verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext3, t)
	})

	t.Run(testString(testContext, "Evaluator/Mul/ct0*ct0->ct0/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, ciphertext1, nil, ciphertext1)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString(testContext, "Evaluator/Mul/ct0*ct0->ct1/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		for i := range values1 {
			values1[i] *= values1[i]
		}

		ciphertext2 := testContext.evaluator.MulRelinNew(ciphertext1, ciphertext1, nil)

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext2, t)
	})

	t.Run(testString(testContext, "Evaluator/Mul/Relinearize(ct0*ct1->ct0)/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		for i := range values1 {
			values1[i] *= values2[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, ciphertext2, nil, ciphertext1)

		testContext.evaluator.Relinearize(ciphertext1, testContext.rlk, ciphertext1)

		require.Equal(t, ciphertext1.Degree(), uint64(1))

		verifyTestVectors(testContext, testContext.decryptor, values1, ciphertext1, t)
	})

	t.Run(testString(testContext, "Evaluator/Mul/Relinearize(ct0*ct1->ct1)/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)
		values2, _, ciphertext2 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		for i := range values1 {
			values2[i] *= values1[i]
		}

		testContext.evaluator.MulRelin(ciphertext1, ciphertext2, nil, ciphertext2)

		testContext.evaluator.Relinearize(ciphertext2, testContext.rlk, ciphertext2)

		testContext.evaluator.Rescale(ciphertext2, testContext.params.scale, ciphertext2)

		require.Equal(t, ciphertext1.Degree(), uint64(1))

		verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext2, t)
	})

}

func testFunctions(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "Functions/PowerOf2/"), func(t *testing.T) {

		if testContext.params.MaxLevel() < 3 {
			t.Skip()
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		n := uint64(2)

		valuesWant := make([]float64, len(values))
		for i := 0; i < len(valuesWant); i++ {
			valuesWant[i] = values[i]
		}

		for i := uint64(0); i < n; i++ {
			for j := 0; j < len(valuesWant); j++ {
				valuesWant[j] *= valuesWant[j]
			}
		}

		testContext.evaluator.PowerOf2(ciphertext, 2, testContext.rlk, ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, valuesWant, ciphertext, t)
	})

	t.Run(testString(testContext, "Functions/Power/"), func(t *testing.T) {

		if testContext.params.MaxLevel() < 4 {
			t.Skip()
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		n := uint64(3)

		for i := range values {
			values[i] = math.Pow(values[i], float64(n))
		}

		testContext.evaluator.Power(ciphertext, n, testContext.rlk, ciphertext)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t)
	})

	t.Run(testString(testContext, "Functions/Inverse/"), func(t *testing.T) {

		if testContext.params.MaxLevel() < 7 {
			t.Skip()
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, 0.1, 1, t)

		n := uint64(7)

		for i := range values {
			values[i] = 1.0 / values[i]
		}

		ciphertext = testContext.evaluator.InverseNew(ciphertext, n, testContext.rlk)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t)
	})
}

func testEvaluatePoly(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "EvaluatePoly/Exp/"), func(t *testing.T) {

		if testContext.params.MaxLevel() < 3 {
			t.Skip()
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		coeffs := []float64{
			1.0,
			1.0,
			1.0 / 2,
			1.0 / 6,
			1.0 / 24,
			1.0 / 120,
			1.0 / 720,
			1.0 / 5040,
		}

		poly := NewPoly(coeffs)

		for i := range values {
			values[i] = math.Exp(values[i])
		}

		ciphertext = testContext.evaluator.EvaluatePoly(ciphertext, poly, testContext.rlk)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t)
	})
}

func testChebyshevInterpolator(testContext *testParams, t *testing.T) {

	t.Run(testString(testContext, "ChebyshevInterpolator/Sin/"), func(t *testing.T) {

		if testContext.params.MaxLevel() < 5 {
			t.Skip()
		}

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		cheby := Approximate(math.Sin, -1.5, 1.5, 15)

		for i := range values {
			values[i] = math.Sin(values[i])
		}

		ciphertext = testContext.evaluator.EvaluateCheby(ciphertext, cheby, testContext.rlk)

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, t)
	})
}

func testSwitchKeys(testContext *testParams, t *testing.T) {

	sk2 := testContext.kgen.GenSecretKey()
	decryptorSk2 := NewDecryptor(testContext.params, sk2)
	switchingKey := testContext.kgen.GenSwitchingKey(testContext.sk, sk2)

	t.Run(testString(testContext, "SwitchKeys/InPlace/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		testContext.evaluator.SwitchKeys(ciphertext, switchingKey, ciphertext)

		verifyTestVectors(testContext, decryptorSk2, values, ciphertext, t)
	})

	t.Run(testString(testContext, "SwitchKeys/New/"), func(t *testing.T) {

		values, _, ciphertext := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		ciphertext = testContext.evaluator.SwitchKeysNew(ciphertext, switchingKey)

		verifyTestVectors(testContext, decryptorSk2, values, ciphertext, t)
	})

}

func testRotate(testContext *testParams, t *testing.T) {

	rotKey := testContext.kgen.GenRotationKeysPow2(testContext.sk)

	t.Run(testString(testContext, "Rotate/InPlace/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		values2 := make([]float64, len(values1))
		ciphertext2 := NewCiphertext(testContext.params, ciphertext1.Degree(), ciphertext1.Level(), ciphertext1.Scale())

		for n := 1; n < len(values1); n <<= 1 {

			// Applies the slots rotation to the values
			for i := range values1 {
				values2[i] = values1[(i+n)%len(values1)]
			}

			testContext.evaluator.Rotate(ciphertext1, uint64(n), rotKey, ciphertext2)

			verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext2, t)
		}

	})

	t.Run(testString(testContext, "Rotate/New/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		values2 := make([]float64, len(values1))
		ciphertext2 := NewCiphertext(testContext.params, ciphertext1.Degree(), ciphertext1.Level(), ciphertext1.Scale())

		for n := 1; n < len(values1); n <<= 1 {

			// Applies the slots rotation to the values
			for i := range values1 {
				values2[i] = values1[(i+n)%len(values1)]
			}

			ciphertext2 = testContext.evaluator.RotateNew(ciphertext1, uint64(n), rotKey)

			verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext2, t)
		}

	})

	t.Run(testString(testContext, "Rotate/Random/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		values2 := make([]float64, len(values1))
		ciphertext2 := NewCiphertext(testContext.params, ciphertext1.Degree(), ciphertext1.Level(), ciphertext1.Scale())

		for n := 1; n < 4; n++ {

			rand := rand.Uint64() % uint64(len(values1))

			// Applies the slots rotation to the values
			for i := range values1 {
				values2[i] = values1[(i+int(rand))%len(values1)]
			}

			testContext.evaluator.Rotate(ciphertext1, rand, rotKey, ciphertext2)

			verifyTestVectors(testContext, testContext.decryptor, values2, ciphertext2, t)
		}

	})

	t.Run(testString(testContext, "Rotate/Hoisted/"), func(t *testing.T) {

		values1, _, ciphertext1 := newTestVectors(testContext, testContext.encryptorSk, -1, 1, t)

		values2 := make([]float64, len(values1))
		rotations := []uint64{0, 1, 2, 3, 4, 5}
		for _, n := range rotations {
			testContext.kgen.GenRotationKey(RotationLeft, testContext.sk, n, rotKey)
		}

		ciphertexts := testContext.evaluator.RotateHoisted(ciphertext1, rotations, rotKey)

		for _, n := range rotations {

			for i := range values1 {
				values2[i] = values1[(i+int(n))%len(values1)]
			}

			verifyTestVectors(testContext, testContext.decryptor, values2, ciphertexts[n], t)
		}

	})
}

func testMarshaller(testContext *testParams, t *testing.T) {

	ringQP := testContext.ringQP

	t.Run("Marshaller/Ciphertext/", func(t *testing.T) {
		t.Run(testString(testContext, "EndToEnd/"), func(t *testing.T) {
			t.Parallel()

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
			t.Parallel()

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

		rotationKey := NewRotationKeys()

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
	})
}
