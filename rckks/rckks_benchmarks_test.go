package rckks

import (
	"testing"

	"github.com/ldsec/lattigo/v2/ring"
)

func BenchmarkRCKKSScheme(b *testing.B) {

	var err error
	var testContext = new(testParams)
	var defaultParams []*Parameters

	if testing.Short() {
		defaultParams = DefaultParams[PN12QP109+3 : PN12QP109+4]
	} else {
		defaultParams = DefaultParams
	}

	for _, defaultParams := range defaultParams {

		if testContext, err = genTestParams(defaultParams, 0); err != nil {
			panic(err)
		}
		benchNTT(testContext, b)
		benchEncoder(testContext, b)
		benchKeyGen(testContext, b)
		benchEncrypt(testContext, b)
		benchDecrypt(testContext, b)
		benchEvaluator(testContext, b)
		benchHoistedRotations(testContext, b)
	}
}

func benchNTT(testContext *testParams, b *testing.B) {

	ringQ := testContext.ringQ
	ringQ4NthRoot, _ := ring.NewRingWithNthRoot(ringQ.N, ringQ.N<<2, ringQ.Modulus)

	sampler := ring.NewUniformSampler(testContext.prng, ringQ)
	p1 := sampler.ReadNew()

	b.Run(testString(testContext, "NTTRCKKS/NTT"), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			NTTRCKKS(ringQ4NthRoot, p1, p1)
		}

	})

	b.Run(testString(testContext, "NTTRCKKS/InvNTT"), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			InvNTTRCKKS(ringQ4NthRoot, p1, p1)
		}
	})
}

func benchEncoder(testContext *testParams, b *testing.B) {

	encoder := testContext.encoder
	slots := testContext.params.Slots()

	b.Run(testString(testContext, "Encoder/Encode/"), func(b *testing.B) {

		values := make([]float64, slots)
		for i := uint64(0); i < slots; i++ {
			values[i] = randomFloat(-1, 1)
		}

		plaintext := NewPlaintext(testContext.params, testContext.params.MaxLevel(), testContext.params.Scale())

		for i := 0; i < b.N; i++ {
			encoder.Encode(plaintext, values, slots)
		}
	})

	b.Run(testString(testContext, "Encoder/Decode/"), func(b *testing.B) {

		values := make([]float64, slots)
		for i := uint64(0); i < slots; i++ {
			values[i] = randomFloat(-1, 1)
		}

		plaintext := NewPlaintext(testContext.params, testContext.params.MaxLevel(), testContext.params.Scale())
		encoder.Encode(plaintext, values, slots)

		for i := 0; i < b.N; i++ {
			encoder.Decode(plaintext, slots)
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
		for i := 0; i < b.N; i++ {
			kgen.GenRelinKey(sk)
		}
	})
}

func benchEncrypt(testContext *testParams, b *testing.B) {

	encryptorPk := testContext.encryptorPk
	encryptorSk := testContext.encryptorSk

	plaintext := NewPlaintext(testContext.params, testContext.params.MaxLevel(), testContext.params.Scale())
	ciphertext := NewCiphertext(testContext.params, 1, testContext.params.MaxLevel(), testContext.params.Scale())

	b.Run(testString(testContext, "Encrypt/Pk/Slow"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorPk.Encrypt(plaintext, ciphertext)
		}
	})

	b.Run(testString(testContext, "Encrypt/Pk/Fast"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			encryptorPk.EncryptFast(plaintext, ciphertext)
		}
	})

	b.Run(testString(testContext, "Encrypt/Sk/"), func(b *testing.B) {
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

	evaluator := testContext.evaluator

	plaintext := NewPlaintext(testContext.params, testContext.params.MaxLevel(), testContext.params.scale)
	ciphertext1 := NewCiphertextRandom(testContext.prng, testContext.params, 1, testContext.params.MaxLevel(), testContext.params.Scale())
	ciphertext2 := NewCiphertextRandom(testContext.prng, testContext.params, 1, testContext.params.MaxLevel(), testContext.params.Scale())
	receiver := NewCiphertextRandom(testContext.prng, testContext.params, 2, testContext.params.MaxLevel(), testContext.params.Scale())

	rlk := testContext.kgen.GenRelinKey(testContext.sk)
	rotkey := NewRotationKeys()
	testContext.kgen.GenRotationKey(RotationLeft, testContext.sk, 1, rotkey)

	b.Run(testString(testContext, "Evaluator/Add/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Add(ciphertext1, ciphertext2, ciphertext1)
		}
	})

	b.Run(testString(testContext, "Evaluator/AddScalar/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.AddConst(ciphertext1, ciphertext2, ciphertext1)
		}
	})

	b.Run(testString(testContext, "Evaluator/MulScalar/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.MultByConst(ciphertext1, complex(3.1415, -1.4142), ciphertext1)
		}
	})

	b.Run(testString(testContext, "Evaluator/MulPlain/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.MulRelin(ciphertext1, plaintext, nil, receiver)
		}
	})

	b.Run(testString(testContext, "Evaluator/Mul/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.MulRelin(ciphertext1, ciphertext2, nil, receiver)
		}
	})

	b.Run(testString(testContext, "Evaluator/Square/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.MulRelin(ciphertext1, ciphertext1, nil, receiver)
		}
	})

	b.Run(testString(testContext, "Evaluator/Rescale/"), func(b *testing.B) {

		ringQ := testContext.ringQ

		for i := 0; i < b.N; i++ {
			ringQ.DivRoundByLastModulusNTT(ciphertext1.Value()[0])
			ringQ.DivRoundByLastModulusNTT(ciphertext1.Value()[1])

			b.StopTimer()
			ciphertext1 = NewCiphertextRandom(testContext.prng, testContext.params, 1, testContext.params.MaxLevel(), testContext.params.Scale())
			b.StartTimer()
		}
	})

	b.Run(testString(testContext, "Evaluator/PermuteNTT/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ring.PermuteNTTWithIndexLvl(ciphertext1.Level(), ciphertext1.value[0], rotkey.permuteNTTLeftIndex[1], ciphertext1.value[0])
			ring.PermuteNTTWithIndexLvl(ciphertext1.Level(), ciphertext1.value[1], rotkey.permuteNTTLeftIndex[1], ciphertext1.value[1])
		}
	})

	b.Run(testString(testContext, "Evaluator/Relin/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Relinearize(receiver, rlk, ciphertext1)
		}
	})

	b.Run(testString(testContext, "Evaluator/Rotate/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.Rotate(ciphertext1, 1, rotkey, ciphertext1)
		}
	})
}

func benchHoistedRotations(testContext *testParams, b *testing.B) {

	evaluator := testContext.evaluator.(*evaluator)

	rotkey := NewRotationKeys()
	testContext.kgen.GenRotationKey(RotationLeft, testContext.sk, 5, rotkey)

	ciphertext := NewCiphertextRandom(testContext.prng, testContext.params, 1, testContext.params.MaxLevel(), testContext.params.Scale())

	ringQ := testContext.ringQ
	ringP := testContext.ringP

	c2NTT := ciphertext.value[1]
	c2InvNTT := ringQ.NewPoly()
	ringQ.InvNTTLvl(ciphertext.Level(), c2NTT, c2InvNTT)

	c2QiQDecomp := make([]*ring.Poly, testContext.params.Beta())
	c2QiPDecomp := make([]*ring.Poly, testContext.params.Beta())

	for i := uint64(0); i < testContext.params.Beta(); i++ {
		c2QiQDecomp[i] = ringQ.NewPoly()
		c2QiPDecomp[i] = ringP.NewPoly()
	}

	b.Run(testString(testContext, "HoistedRotations/DecomposeNTT/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			for j := uint64(0); j < testContext.params.Beta(); j++ {
				evaluator.decomposeAndSplitNTT(ciphertext.Level(), j, c2NTT, c2InvNTT, c2QiQDecomp[j], c2QiPDecomp[j])
			}
		}
	})

	b.Run(testString(testContext, "HoistedRotations/RotateHoisted/"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			evaluator.permuteNTTHoisted(ciphertext, c2QiQDecomp, c2QiPDecomp, 5, rotkey, ciphertext)
		}
	})
}
