package ckks

import (
	"encoding/json"
	"testing"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

func BenchmarkCKKSScheme(b *testing.B) {

	defaultParams := DefaultParams
	if testing.Short() {
		defaultParams = DefaultParams[:2]
	}
	if *flagParamString != "" {
		var jsonParams ParametersLiteral
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, defaultParams := range defaultParams {
		params, err := NewParametersFromLiteral(defaultParams)
		if err != nil {
			panic(err)
		}
		var testContext *testParams
		if testContext, err = genTestParams(params, 0); err != nil {
			panic(err)
		}

		benchEncoder(testContext, b)
		benchKeyGen(testContext, b)
		benchEncrypt(testContext, b)
		benchDecrypt(testContext, b)
		benchEvaluator(testContext, b)
		benchInnerSum(testContext, b)
	}
}

func benchEncoder(testContext *testParams, b *testing.B) {

	encoder := testContext.encoder
	logSlots := testContext.params.LogSlots()

	b.Run(testString(testContext, "Encoder/Encode/"), func(b *testing.B) {

		values := make([]complex128, 1<<logSlots)
		for i := 0; i < 1<<logSlots; i++ {
			values[i] = utils.RandComplex128(-1, 1)
		}

		plaintext := NewPlaintext(testContext.params, testContext.params.MaxLevel(), testContext.params.Scale())

		for i := 0; i < b.N; i++ {
			encoder.Encode(plaintext, values, logSlots)
		}
	})

	b.Run(testString(testContext, "Encoder/Decode/"), func(b *testing.B) {

		values := make([]complex128, 1<<logSlots)
		for i := 0; i < 1<<logSlots; i++ {
			values[i] = utils.RandComplex128(-1, 1)
		}

		plaintext := NewPlaintext(testContext.params, testContext.params.MaxLevel(), testContext.params.Scale())
		encoder.Encode(plaintext, values, logSlots)

		for i := 0; i < b.N; i++ {
			encoder.Decode(plaintext, logSlots)
		}
	})
}

func benchKeyGen(testContext *testParams, b *testing.B) {

	kgen := testContext.kgen
	sk := testContext.sk

	b.Run(testString(testContext, "KeyGen/KeyPairGen/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			kgen.GenKeyPair()
		}
	})

	b.Run(testString(testContext, "KeyGen/SwitchKeyGen/"), func(b *testing.B) {

		if testContext.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		for i := 0; i < b.N; i++ {
			kgen.GenRelinearizationKey(sk, 2)
		}
	})
}

func benchEncrypt(testContext *testParams, b *testing.B) {

	encryptorPk := testContext.encryptorPk
	encryptorPkFast := NewFastEncryptor(testContext.params, testContext.pk)
	encryptorSk := testContext.encryptorSk

	plaintext := NewPlaintext(testContext.params, testContext.params.MaxLevel(), testContext.params.Scale())
	ciphertext := NewCiphertext(testContext.params, 1, testContext.params.MaxLevel(), testContext.params.Scale())

	b.Run(testString(testContext, "Encrypt/key=Pk/"), func(b *testing.B) {

		if testContext.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		for i := 0; i < b.N; i++ {
			encryptorPk.Encrypt(plaintext, ciphertext)
		}
	})

	b.Run(testString(testContext, "EncryptFast/key=Pk/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorPkFast.Encrypt(plaintext, ciphertext)
		}
	})

	b.Run(testString(testContext, "Encrypt/key=Sk/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorSk.Encrypt(plaintext, ciphertext)
		}
	})
}

func benchDecrypt(testContext *testParams, b *testing.B) {

	decryptor := testContext.decryptor

	plaintext := NewPlaintext(testContext.params, testContext.params.MaxLevel(), testContext.params.Scale())
	ciphertext := NewCiphertextRandom(testContext.prng, testContext.params, 1, testContext.params.MaxLevel(), testContext.params.Scale())

	b.Run(testString(testContext, "Decrypt/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			decryptor.Decrypt(ciphertext, plaintext)
		}
	})
}

func benchEvaluator(testContext *testParams, b *testing.B) {

	plaintext := NewPlaintext(testContext.params, testContext.params.MaxLevel(), testContext.params.Scale())
	ciphertext1 := NewCiphertextRandom(testContext.prng, testContext.params, 1, testContext.params.MaxLevel(), testContext.params.Scale())
	ciphertext2 := NewCiphertextRandom(testContext.prng, testContext.params, 1, testContext.params.MaxLevel(), testContext.params.Scale())
	receiver := NewCiphertextRandom(testContext.prng, testContext.params, 2, testContext.params.MaxLevel(), testContext.params.Scale())

	var rlk *rlwe.RelinearizationKey
	var rotkey *rlwe.RotationKeySet
	if testContext.params.PCount() != 0 {
		rlk = testContext.kgen.GenRelinearizationKey(testContext.sk, 2)
		rotkey = testContext.kgen.GenRotationKeysForRotations([]int{1}, true, testContext.sk)
	}

	eval := testContext.evaluator.WithKey(rlwe.EvaluationKey{Rlk: rlk, Rtks: rotkey})

	b.Run(testString(testContext, "Evaluator/Add/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Add(ciphertext1, ciphertext2, ciphertext1)
		}
	})

	b.Run(testString(testContext, "Evaluator/AddScalar/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.AddConst(ciphertext1, ciphertext2, ciphertext1)
		}
	})

	b.Run(testString(testContext, "Evaluator/MulScalar/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.MultByConst(ciphertext1, complex(3.1415, -1.4142), ciphertext1)
		}
	})

	b.Run(testString(testContext, "Evaluator/MulPlain/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Mul(ciphertext1, plaintext, ciphertext1)
		}
	})

	b.Run(testString(testContext, "Evaluator/Mul/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Mul(ciphertext1, ciphertext2, receiver)
		}
	})

	b.Run(testString(testContext, "Evaluator/Square/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			eval.Mul(ciphertext1, ciphertext1, receiver)
		}
	})

	b.Run(testString(testContext, "Evaluator/Rescale/"), func(b *testing.B) {

		if testContext.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		ciphertext1.Scale = testContext.params.Scale() * testContext.params.Scale()

		for i := 0; i < b.N; i++ {
			if err := eval.Rescale(ciphertext1, testContext.params.Scale(), ciphertext2); err != nil {
				panic(err)
			}
		}
	})

	b.Run(testString(testContext, "Evaluator/PermuteNTT/"), func(b *testing.B) {

		if testContext.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		galEL := testContext.params.GaloisElementForColumnRotationBy(1)
		for i := 0; i < b.N; i++ {
			ring.PermuteNTTWithIndexLvl(ciphertext1.Level(), ciphertext1.Value[0], eval.(*evaluator).permuteNTTIndex[galEL], ciphertext1.Value[0])
			ring.PermuteNTTWithIndexLvl(ciphertext1.Level(), ciphertext1.Value[1], eval.(*evaluator).permuteNTTIndex[galEL], ciphertext1.Value[1])
		}
	})

	b.Run(testString(testContext, "Evaluator/Conjugate/"), func(b *testing.B) {

		if testContext.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		for i := 0; i < b.N; i++ {
			eval.Conjugate(ciphertext1, ciphertext1)
		}
	})

	b.Run(testString(testContext, "Evaluator/Relin/"), func(b *testing.B) {

		if testContext.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		for i := 0; i < b.N; i++ {
			eval.Relinearize(receiver, ciphertext1)
		}
	})

	b.Run(testString(testContext, "Evaluator/Rotate/"), func(b *testing.B) {

		if testContext.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		for i := 0; i < b.N; i++ {
			eval.Rotate(ciphertext1, 1, ciphertext1)
		}
	})
}

func benchInnerSum(testContext *testParams, b *testing.B) {

	ciphertext1 := NewCiphertextRandom(testContext.prng, testContext.params, 1, testContext.params.MaxLevel(), testContext.params.Scale())

	batch := 1
	n := 4

	b.Run(testString(testContext, "InnerSum/"), func(b *testing.B) {

		if testContext.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		rotKey := testContext.kgen.GenRotationKeysForRotations(testContext.params.RotationsForInnerSum(batch, n), false, testContext.sk)
		eval := testContext.evaluator.WithKey(rlwe.EvaluationKey{Rlk: testContext.rlk, Rtks: rotKey})

		b.ResetTimer()

		for i := 0; i < b.N; i++ {
			eval.InnerSum(ciphertext1, batch, n, ciphertext1)
		}
	})

	b.Run(testString(testContext, "InnerSumLog/"), func(b *testing.B) {

		if testContext.params.PCount() == 0 {
			b.Skip("#Pi is empty")
		}

		rotKey := testContext.kgen.GenRotationKeysForRotations(testContext.params.RotationsForInnerSumLog(batch, n), false, testContext.sk)
		eval := testContext.evaluator.WithKey(rlwe.EvaluationKey{Rlk: testContext.rlk, Rtks: rotKey})

		b.ResetTimer()

		for i := 0; i < b.N; i++ {
			eval.InnerSumLog(ciphertext1, batch, n, ciphertext1)
		}
	})

}
