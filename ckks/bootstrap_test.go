package ckks

import (
	"github.com/ldsec/lattigo/v2/ckks/bettersine"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
	"math"
	"math/cmplx"
	"runtime"
	"testing"
)

func TestBootstrap(t *testing.T) {

	if !*testBootstrapping {
		t.Skip("skipping bootstrapping test")
	}

	if runtime.GOARCH == "wasm" {
		t.Skip("skipping bootstrapping tests for GOARCH=wasm")
	}

	var testContext = new(testParams)

	paramSet := 0

	bootstrapParams := DefaultBootstrapParams[paramSet : paramSet+1]

	for paramSet := range bootstrapParams {

		btpParams := bootstrapParams[paramSet]

		// Insecure params for fast testing only
		if !*flagLongTest {
			btpParams.LogN = 14
			btpParams.LogSlots = 13
		}

		// Tests homomorphic modular reduction encoding and bootstrapping on sparse slots
		params, err := btpParams.Params()
		if err != nil {
			panic(err)
		}

		if testContext, err = genTestParams(params, btpParams.H); err != nil { // TODO: setting the param.scale field is not something the user can do
			panic(err)
		}

		for _, testSet := range []func(testContext *testParams, btpParams *BootstrappingParameters, t *testing.T){
			testEvalSine,
		} {
			testSet(testContext, btpParams, t)
			runtime.GC()
		}

		for _, testSet := range []func(testContext *testParams, btpParams *BootstrappingParameters, t *testing.T){
			testCoeffsToSlots,
			testSlotsToCoeffs,
			testbootstrap,
		} {
			testSet(testContext, btpParams, t)
			runtime.GC()
		}

		if !*flagLongTest {
			btpParams.LogSlots = 12
		}

		// Tests homomorphic encoding and bootstrapping on full slots
		params, err = btpParams.Params()
		if err != nil {
			panic(err)
		}

		if testContext, err = genTestParams(params, btpParams.H); err != nil { // TODO: setting the param.scale field is not something the user can do
			panic(err)
		}

		for _, testSet := range []func(testContext *testParams, btpParams *BootstrappingParameters, t *testing.T){
			testCoeffsToSlots,
			testSlotsToCoeffs,
			testbootstrap,
		} {
			testSet(testContext, btpParams, t)
			runtime.GC()
		}
	}
}

func testEvalSine(testContext *testParams, btpParams *BootstrappingParameters, t *testing.T) {

	t.Run(testString(testContext, "Sin/"), func(t *testing.T) {

		var err error

		eval := testContext.evaluator

		DefaultScale := testContext.params.Scale()

		SineScale := btpParams.SineEvalModuli.ScalingFactor

		testContext.params.scale = SineScale
		eval.(*evaluator).scale = SineScale

		deg := 127
		K := float64(15)

		values, _, ciphertext := newTestVectorsSineBootstrapp(testContext, btpParams, testContext.encryptorSk, -K+1, K-1, t)
		eval.DropLevel(ciphertext, btpParams.CtSDepth(true)-1)

		cheby := Approximate(sin2pi2pi, -complex(K, 0), complex(K, 0), deg)

		for i := range values {
			values[i] = sin2pi2pi(values[i])
		}

		eval.MultByConst(ciphertext, 2/(cheby.b-cheby.a), ciphertext)
		eval.AddConst(ciphertext, (-cheby.a-cheby.b)/(cheby.b-cheby.a), ciphertext)
		eval.Rescale(ciphertext, eval.(*evaluator).scale, ciphertext)

		if ciphertext, err = eval.EvaluateCheby(ciphertext, cheby, ciphertext.Scale); err != nil {
			t.Error(err)
		}

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, testContext.params.LogSlots(), 0, t)

		testContext.params.scale = DefaultScale
		eval.(*evaluator).scale = DefaultScale
	})

	t.Run(testString(testContext, "Cos1/"), func(t *testing.T) {

		var err error

		eval := testContext.evaluator

		DefaultScale := testContext.params.Scale()

		SineScale := btpParams.SineEvalModuli.ScalingFactor

		testContext.params.scale = SineScale
		eval.(*evaluator).scale = SineScale

		K := 25
		deg := 63
		dev := btpParams.MessageRatio
		scNum := 2

		scFac := complex(float64(int(1<<scNum)), 0)

		values, _, ciphertext := newTestVectorsSineBootstrapp(testContext, btpParams, testContext.encryptorSk, float64(-K+1), float64(K-1), t)
		eval.DropLevel(ciphertext, btpParams.CtSDepth(true)-1)

		cheby := new(ChebyshevInterpolation)
		cheby.coeffs = bettersine.Approximate(K, deg, dev, scNum)
		cheby.maxDeg = cheby.Degree()
		cheby.a = complex(float64(-K), 0) / scFac
		cheby.b = complex(float64(K), 0) / scFac
		cheby.lead = true

		var sqrt2pi float64
		if btpParams.ArcSineDeg > 0 {
			sqrt2pi = math.Pow(1, 1.0/real(scFac))
		} else {
			sqrt2pi = math.Pow(0.15915494309189535, 1.0/real(scFac))
		}

		for i := range cheby.coeffs {
			cheby.coeffs[i] *= complex(sqrt2pi, 0)
		}

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, testContext.params.LogSlots(), 0, t)

		for i := range values {

			values[i] = cmplx.Cos(6.283185307179586 * (1 / scFac) * (values[i] - 0.25))

			for j := 0; j < scNum; j++ {
				values[i] = 2*values[i]*values[i] - 1
			}

			if btpParams.ArcSineDeg == 0 {
				values[i] /= 6.283185307179586
			}
		}

		eval.AddConst(ciphertext, -0.25, ciphertext)

		eval.MultByConst(ciphertext, 2/((cheby.b-cheby.a)*scFac), ciphertext)
		eval.AddConst(ciphertext, (-cheby.a-cheby.b)/(cheby.b-cheby.a), ciphertext)
		eval.Rescale(ciphertext, eval.(*evaluator).scale, ciphertext)

		if ciphertext, err = eval.EvaluateCheby(ciphertext, cheby, ciphertext.Scale); err != nil {
			t.Error(err)
		}

		for i := 0; i < scNum; i++ {
			sqrt2pi *= sqrt2pi
			eval.MulRelin(ciphertext, ciphertext, ciphertext)
			eval.Add(ciphertext, ciphertext, ciphertext)
			eval.AddConst(ciphertext, -sqrt2pi, ciphertext)
			eval.Rescale(ciphertext, eval.(*evaluator).scale, ciphertext)
		}

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, testContext.params.LogSlots(), 0, t)

		testContext.params.scale = DefaultScale
		eval.(*evaluator).scale = DefaultScale

	})

	t.Run(testString(testContext, "Cos2/"), func(t *testing.T) {

		if len(btpParams.SineEvalModuli.Qi) < 12 {
			t.Skip()
		}

		var err error

		eval := testContext.evaluator

		DefaultScale := testContext.params.Scale()

		SineScale := btpParams.SineEvalModuli.ScalingFactor

		testContext.params.scale = SineScale
		eval.(*evaluator).scale = SineScale

		K := 325
		deg := 255
		scNum := 4

		scFac := complex(float64(int(1<<scNum)), 0)

		values, _, ciphertext := newTestVectorsSineBootstrapp(testContext, btpParams, testContext.encryptorSk, float64(-K+1), float64(K-1), t)
		eval.DropLevel(ciphertext, btpParams.CtSDepth(true)-1)

		cheby := Approximate(cos2pi, -complex(float64(K), 0)/scFac, complex(float64(K), 0)/scFac, deg)

		sqrt2pi := math.Pow(0.15915494309189535, 1.0/real(scFac))

		for i := range cheby.coeffs {
			cheby.coeffs[i] *= complex(sqrt2pi, 0)
		}

		for i := range values {

			values[i] = cmplx.Cos(6.283185307179586 * (1 / scFac) * (values[i] - 0.25))

			for j := 0; j < scNum; j++ {
				values[i] = 2*values[i]*values[i] - 1
			}

			values[i] /= 6.283185307179586
		}

		testContext.evaluator.AddConst(ciphertext, -0.25, ciphertext)

		eval.MultByConst(ciphertext, 2/((cheby.b-cheby.a)*scFac), ciphertext)
		eval.AddConst(ciphertext, (-cheby.a-cheby.b)/(cheby.b-cheby.a), ciphertext)
		eval.Rescale(ciphertext, eval.(*evaluator).scale, ciphertext)

		if ciphertext, err = eval.EvaluateCheby(ciphertext, cheby, ciphertext.Scale); err != nil {
			t.Error(err)
		}

		for i := 0; i < scNum; i++ {
			sqrt2pi *= sqrt2pi
			eval.MulRelin(ciphertext, ciphertext, ciphertext)
			eval.Add(ciphertext, ciphertext, ciphertext)
			eval.AddConst(ciphertext, -sqrt2pi, ciphertext)
			eval.Rescale(ciphertext, eval.(*evaluator).scale, ciphertext)
		}

		verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, testContext.params.LogSlots(), 0, t)

		testContext.params.scale = DefaultScale
		eval.(*evaluator).scale = DefaultScale

	})
}

func testCoeffsToSlots(testContext *testParams, btpParams *BootstrappingParameters, t *testing.T) {
	t.Run(testString(testContext, "CoeffsToSlots/"), func(t *testing.T) {

		// This test tests the homomorphic encoding
		// It first generates a vector of complex values of size params.Slots()
		//
		// vReal + i*vImag
		//
		// Then encode coefficient-wise and encrypts the vectors :
		//
		// Enc(bitReverse(vReal)||bitReverse(vImg))
		//
		// And applies the homomorphic Encoding (will merge both vectors if there was two)
		//
		// Enc(iFFT(vReal+ i*vImag))
		//
		// And returns the result in one ciphextext if the ciphertext can store it else in two ciphertexts
		//
		// Enc(Ecd(vReal) || Ecd(vImag)) or Enc(Ecd(vReal)) and Enc(Ecd(vImag))
		//
		// Then checks that Dcd(Dec(Enc(Ecd(vReal)))) = vReal and Dcd(Dec(Enc(Ecd(vImag)))) = vImag

		params := testContext.params

		n := math.Pow(1.0/float64(2*params.Slots()), 1.0/float64(btpParams.CtSDepth(true)))

		// Generates the encoding matrices
		CoeffsToSlotMatrices := btpParams.GenCoeffsToSlotsMatrix(complex(n, 0), testContext.encoder)

		// Gets the rotations indexes for CoeffsToSlots
		rotations := btpParams.RotationsForCoeffsToSlots(params.LogSlots())

		// Generates the rotation keys
		rotKey := testContext.kgen.GenRotationKeysForRotations(rotations, true, testContext.sk)

		// Generates the vector of random complex values
		values := make([]complex128, params.Slots())
		for i := range values {
			values[i] = complex(utils.RandFloat64(-1, 1), utils.RandFloat64(-1, 1))
		}

		// Splits between real and imaginary
		valuesReal := make([]complex128, params.Slots())
		for i := range valuesReal {
			valuesReal[i] = complex(real(values[i]), 0)
		}

		valuesImag := make([]complex128, params.Slots())
		for i := range valuesImag {
			valuesImag[i] = complex(imag(values[i]), 0)
		}

		// Applies bit-reverse on the original complex vector
		sliceBitReverseInPlaceComplex128(values, params.Slots())

		// Maps to a float vector
		// Add gaps if sparse packing
		valuesFloat := make([]float64, params.N())
		gap := params.N() / (2 * params.Slots())
		for i, jdx, idx := 0, params.N()>>1, 0; i < params.Slots(); i, jdx, idx = i+1, jdx+gap, idx+gap {
			valuesFloat[idx] = real(values[i])
			valuesFloat[jdx] = imag(values[i])
		}

		// Encodes coefficient-wise and encrypts the test vector
		plaintext := NewPlaintext(params, params.MaxLevel(), params.Scale())
		testContext.encoder.EncodeCoeffs(valuesFloat, plaintext)
		ciphertext := testContext.encryptorPk.EncryptNew(plaintext)

		// Creates an evaluator with the rotation keys
		eval := testContext.evaluator.WithKey(rlwe.EvaluationKey{Rlk: testContext.rlk, Rtks: rotKey})

		// Applies the homomorphic DFT
		ct0, ct1 := CoeffsToSlots(ciphertext, CoeffsToSlotMatrices, eval)

		// Checks against the original coefficients
		var coeffsReal, coeffsImag []complex128
		if params.LogSlots() < params.LogN()-1 {
			coeffsRealImag := testContext.encoder.DecodePublic(testContext.decryptor.DecryptNew(ct0), params.LogSlots()+1, 0)
			coeffsReal = coeffsRealImag[:params.Slots()]
			coeffsImag = coeffsRealImag[params.Slots():]
		} else {
			coeffsReal = testContext.encoder.DecodePublic(testContext.decryptor.DecryptNew(ct0), params.LogSlots(), 0)
			coeffsImag = testContext.encoder.DecodePublic(testContext.decryptor.DecryptNew(ct1), params.LogSlots(), 0)
		}

		verifyTestVectors(testContext, nil, valuesReal, coeffsReal, params.LogSlots(), 0, t)
		verifyTestVectors(testContext, nil, valuesImag, coeffsImag, params.LogSlots(), 0, t)
	})
}

func testSlotsToCoeffs(testContext *testParams, btpParams *BootstrappingParameters, t *testing.T) {
	t.Run(testString(testContext, "SlotsToCoeffs/"), func(t *testing.T) {

		// This test tests the homomorphic decoding
		// It first generates a complex vector of size 2*slots
		// if 2*slots == N, then two vectors are generated, one for the real part, one for the imaginary part :
		//
		// vReal and vReal (both floating point vectors because the encoding always result in a real vector)
		//
		// Then encode and encrypts the vectors :
		//
		// Enc(Ecd(vReal)) and Enc(Ecd(vImag))
		//
		// And applies the homomorphic decoding (will merge both vectors if there was two)
		//
		// Enc(FFT(Ecd(vReal) + i*Ecd(vImag)))
		//
		// The result should be the decoding of the initial vectors bit-reversed
		//
		// Enc(FFT(Ecd(vReal) + i*Ecd(vImag))) = Enc(BitReverse(Dcd(Ecd(vReal + i*vImag))))
		//
		// The first N/2 slots of the plaintext will be the real part while the last N/2 the imaginary part
		// In case of 2*slots < N, then there is a gap of N/(2*slots) between each values

		params := testContext.params

		// Generates the encoding matrices
		SlotsToCoeffsMatrix := btpParams.GenSlotsToCoeffsMatrix(1.0, testContext.encoder)

		// Gets the rotations indexes for SlotsToCoeffs
		rotations := btpParams.RotationsForSlotsToCoeffs(params.LogSlots())

		// Generates the rotation keys
		rotKey := testContext.kgen.GenRotationKeysForRotations(rotations, true, testContext.sk)

		// Creates an evaluator with the rotation keys
		eval := testContext.evaluator.WithKey(rlwe.EvaluationKey{Rlk: testContext.rlk, Rtks: rotKey})

		// Generates the n first slots of the test vector (real part to encode)
		valuesReal := make([]complex128, params.Slots())
		for i := range valuesReal {
			valuesReal[i] = complex(float64(i+1)/float64(params.Slots()), 0)
		}

		// Generates the n first slots of the test vector (imaginary part to encode)
		valuesImag := make([]complex128, params.Slots())
		for i := range valuesImag {
			valuesImag[i] = complex(1+float64(i+1)/float64(params.Slots()), 0)
		}

		// If sparse, there there is the space to store both vectors in one
		if params.LogSlots() < params.LogN()-1 {
			for i := range valuesReal {
				valuesReal[i] += complex(0, real(valuesImag[i]))
			}
		}

		// Encodes and encrypts the test vectors
		logSlots := params.LogSlots()
		if params.LogSlots() < params.LogN()-1 {
			logSlots++
		}
		encoder := testContext.encoder.(*encoderComplex128)
		plaintext := NewPlaintext(params, params.MaxLevel(), params.Scale())
		encoder.Encode(plaintext, valuesReal, logSlots)
		ct0 := testContext.encryptorPk.EncryptNew(plaintext)
		var ct1 *Ciphertext
		if params.LogSlots() == params.LogN()-1 {
			testContext.encoder.Encode(plaintext, valuesImag, logSlots)
			ct1 = testContext.encryptorPk.EncryptNew(plaintext)
		}

		// Applies the homomorphic DFT
		res := SlotsToCoeffs(ct0, ct1, SlotsToCoeffsMatrix, eval)

		// Decrypt and decode in the coefficient domain
		coeffsFloat := testContext.encoder.DecodeCoeffsPublic(testContext.decryptor.DecryptNew(res), 0)

		// Extracts the coefficients and construct the complex vector
		// This is simply coefficient ordering
		valuesTest := make([]complex128, params.Slots())
		gap := params.N() / (2 * params.Slots())
		for i, idx := 0, 0; i < params.Slots(); i, idx = i+1, idx+gap {
			valuesTest[i] = complex(coeffsFloat[idx], coeffsFloat[idx+(params.N()>>1)])
		}

		// The result is always returned as a single complex vector, so if full-packing (2 initial vectors)
		// then repacks both vectors together
		if params.LogSlots() == params.LogN()-1 {
			for i := range valuesReal {
				valuesReal[i] += complex(0, real(valuesImag[i]))
			}
		}

		// Result is bit-reversed, so applies the bit-reverse permutation on the reference vector
		sliceBitReverseInPlaceComplex128(valuesReal, params.Slots())

		verifyTestVectors(testContext, testContext.decryptor, valuesReal, valuesTest, params.LogSlots(), 0, t)
	})
}

func testbootstrap(testContext *testParams, btpParams *BootstrappingParameters, t *testing.T) {

	t.Run(testString(testContext, "Bootstrapping/FullCircuit/"), func(t *testing.T) {

		params := testContext.params

		rotations := btpParams.RotationsForBootstrapping(testContext.params.LogSlots())
		rotkeys := testContext.kgen.GenRotationKeysForRotations(rotations, true, testContext.sk)
		btpKey := BootstrappingKey{testContext.rlk, rotkeys}

		btp, err := NewBootstrapper(testContext.params, btpParams, btpKey)
		if err != nil {
			panic(err)
		}

		values := make([]complex128, 1<<params.LogSlots())
		for i := range values {
			values[i] = utils.RandComplex128(-1, 1)
		}

		values[0] = complex(0.9238795325112867, 0.3826834323650898)
		values[1] = complex(0.9238795325112867, 0.3826834323650898)
		if 1<<params.LogSlots() > 2 {
			values[2] = complex(0.9238795325112867, 0.3826834323650898)
			values[3] = complex(0.9238795325112867, 0.3826834323650898)
		}

		plaintext := NewPlaintext(params, params.MaxLevel(), params.Scale())
		testContext.encoder.Encode(plaintext, values, params.LogSlots())

		ciphertext := testContext.encryptorPk.EncryptNew(plaintext)

		eval := testContext.evaluator
		for ciphertext.Level() != 0 {
			eval.DropLevel(ciphertext, 1)
		}

		for i := 0; i < 1; i++ {

			ciphertext = btp.Bootstrapp(ciphertext)
			//testContext.evaluator.SetScale(ciphertext, testContext.params.scale)
			verifyTestVectors(testContext, testContext.decryptor, values, ciphertext, testContext.params.LogSlots(), 0, t)
		}

	})
}

func newTestVectorsSineBootstrapp(testContext *testParams, btpParams *BootstrappingParameters, encryptor Encryptor, a, b float64, t *testing.T) (values []complex128, plaintext *Plaintext, ciphertext *Ciphertext) {

	logSlots := testContext.params.LogSlots()

	values = make([]complex128, 1<<logSlots)

	ratio := btpParams.MessageRatio

	for i := uint64(0); i < 1<<logSlots; i++ {
		values[i] = complex(math.Round(utils.RandFloat64(a, b))+utils.RandFloat64(-1, 1)/ratio, 0)
	}

	plaintext = NewPlaintext(testContext.params, testContext.params.MaxLevel(), testContext.params.Scale())

	testContext.encoder.EncodeNTT(plaintext, values, logSlots)

	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return values, plaintext, ciphertext
}
