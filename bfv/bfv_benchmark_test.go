package bfv

import (
	"encoding/json"
	"testing"

	"github.com/ldsec/lattigo/v2/rlwe"
)

func BenchmarkBFV(b *testing.B) {

	defaultParams := DefaultParams
	if testing.Short() {
		defaultParams = DefaultParams[:2]
	} else {

	}

	if *flagParamString != "" {
		var jsonParams ParametersLiteral
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range defaultParams {

		params, err := NewParametersFromLiteral(p)
		var testctx *testContext
		if testctx, err = genTestParams(params); err != nil {
			panic(err)
		}

		benchEncoder(testctx, b)
		benchKeyGen(testctx, b)
		benchEncrypt(testctx, b)
		benchDecrypt(testctx, b)
		benchEvaluator(testctx, b)
	}
}

func benchEncoder(testctx *testContext, b *testing.B) {

	encoder := testctx.encoder
	coeffs := testctx.uSampler.ReadNew()
	coeffsOut := make([]uint64, testctx.params.N())

	plaintext := NewPlaintext(testctx.params)
	plaintextRingT := NewPlaintextRingT(testctx.params)
	plaintextMul := NewPlaintextMul(testctx.params)

	b.Run(testString("Encoder/EncodeUint/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encoder.EncodeUint(coeffs.Coeffs[0], plaintext)
		}
	})

	b.Run(testString("Encoder/DecodeUint/pt=Plaintext/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			testctx.encoder.DecodeUint(plaintext, coeffsOut)
		}
	})

	b.Run(testString("Encoder/EncodeUintRingT/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encoder.EncodeUintRingT(coeffs.Coeffs[0], plaintextRingT)
		}
	})

	b.Run(testString("Encoder/DecodeUint/pt=PlaintextRingT/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			testctx.encoder.DecodeUint(plaintextRingT, coeffsOut)
		}
	})

	b.Run(testString("Encoder/EncodeUintMul/", testctx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			encoder.EncodeUintMul(coeffs.Coeffs[0], plaintextMul)
		}
	})
}

func benchKeyGen(testctx *testContext, b *testing.B) {

	kgen := testctx.kgen
	sk := testctx.sk

	b.Run(testString("KeyGen/KeyPairGen/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			kgen.GenKeyPair()
		}
	})

	b.Run(testString("KeyGen/SwitchKeyGen/", testctx.params), func(b *testing.B) {

		if testctx.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		for i := 0; i < b.N; i++ {
			kgen.GenRelinearizationKey(sk, 1)
		}
	})
}

func benchEncrypt(testctx *testContext, b *testing.B) {

	encryptorPk := testctx.encryptorPk
	encryptorPkFast := NewFastEncryptor(testctx.params, testctx.pk)
	encryptorSk := testctx.encryptorSk

	plaintext := NewPlaintext(testctx.params)
	ciphertext := NewCiphertextRandom(testctx.prng, testctx.params, 1)

	b.Run(testString("Encrypt/key=Pk/", testctx.params), func(b *testing.B) {

		if testctx.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		for i := 0; i < b.N; i++ {
			encryptorPk.Encrypt(plaintext, ciphertext)
		}
	})

	b.Run(testString("EncryptFast/key=Pk/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorPkFast.Encrypt(plaintext, ciphertext)
		}
	})

	b.Run(testString("Encrypt/key=Sk/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorSk.Encrypt(plaintext, ciphertext)
		}
	})
}

func benchDecrypt(testctx *testContext, b *testing.B) {

	decryptor := testctx.decryptor
	ciphertext := NewCiphertextRandom(testctx.prng, testctx.params, 1)

	b.Run(testString("Decrypt/", testctx.params), func(b *testing.B) {
		plaintext := NewPlaintext(testctx.params)
		for i := 0; i < b.N; i++ {
			decryptor.Decrypt(ciphertext, plaintext)
		}
	})
}

func benchEvaluator(testctx *testContext, b *testing.B) {

	encoder := testctx.encoder

	plaintext := NewPlaintext(testctx.params)
	plaintextRingT := NewPlaintextRingT(testctx.params)
	plaintextMul := NewPlaintextMul(testctx.params)

	coeffs := testctx.uSampler.ReadNew()
	encoder.EncodeUintRingT(coeffs.Coeffs[0], plaintextRingT)
	encoder.EncodeUint(coeffs.Coeffs[0], plaintext)
	encoder.EncodeUintMul(coeffs.Coeffs[0], plaintextMul)

	ciphertext1 := NewCiphertextRandom(testctx.prng, testctx.params, 1)
	ciphertext2 := NewCiphertextRandom(testctx.prng, testctx.params, 1)
	receiver := NewCiphertextRandom(testctx.prng, testctx.params, 2)

	var rotkey *rlwe.RotationKeySet
	if testctx.params.PCount() != 0 {
		rotkey = testctx.kgen.GenRotationKeysForRotations([]int{1}, true, testctx.sk)
	}
	evaluator := testctx.evaluator.WithKey(rlwe.EvaluationKey{Rlk: testctx.rlk, Rtks: rotkey})

	b.Run(testString("Evaluator/Add/Ct/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Add(ciphertext1, ciphertext2, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/Add/op1=Ciphertext/op2=PlaintextRingT/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Add(ciphertext1, plaintextRingT, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/Add/op1=Ciphertext/op2=Plaintext/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Add(ciphertext1, plaintext, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/MulScalar/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.MulScalar(ciphertext1, 5, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=Ciphertext/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Mul(ciphertext1, ciphertext2, receiver)
		}
	})

	b.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=Plaintext/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Mul(ciphertext1, plaintext, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=PlaintextRingT/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Mul(ciphertext1, plaintextRingT, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/Mul/op1=Ciphertext/op2=PlaintextMul/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Mul(ciphertext1, plaintextMul, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/Square/", testctx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Mul(ciphertext1, ciphertext1, receiver)
		}
	})

	b.Run(testString("Evaluator/Relin/", testctx.params), func(b *testing.B) {

		if testctx.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		for i := 0; i < b.N; i++ {
			evaluator.Relinearize(receiver, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/RotateRows/", testctx.params), func(b *testing.B) {

		if testctx.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		for i := 0; i < b.N; i++ {
			evaluator.RotateRows(ciphertext1, ciphertext1)
		}
	})

	b.Run(testString("Evaluator/RotateCols/", testctx.params), func(b *testing.B) {

		if testctx.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		for i := 0; i < b.N; i++ {
			evaluator.RotateColumns(ciphertext1, 1, ciphertext1)
		}
	})
}
