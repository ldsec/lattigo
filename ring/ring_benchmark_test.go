package ring

import (
	"fmt"
	"math/big"
	"math/bits"
	"math/rand"
	"testing"
)

func BenchmarkRing(b *testing.B) {

	var err error

	var defaultParams []*Parameters

	if testing.Short() {
		defaultParams = DefaultParams[:3]
	} else {
		defaultParams = DefaultParams
	}

	var testContext = new(testParams)

	for _, defaultParam := range defaultParams {

		if testContext, err = genTestParams(defaultParam); err != nil {
			panic(err)
		}

		benchGenRing(testContext, b)
		benchMarshalling(testContext, b)
		benchSampling(testContext, b)
		benchMontgomery(testContext, b)
		benchNTT(testContext, b)
		benchMurakami(testContext, b)
		benchMulCoeffs(testContext, b)
		benchAddCoeffs(testContext, b)
		benchSubCoeffs(testContext, b)
		benchNegCoeffs(testContext, b)
		benchMulScalar(testContext, b)
		benchExtendBasis(testContext, b)
		benchDivByLastModulus(testContext, b)
		benchDivByRNSBasis(testContext, b)
		benchMRed(testContext, b)
		benchBRed(testContext, b)
		benchBRedAdd(testContext, b)
	}
}

func benchGenRing(testContext *testParams, b *testing.B) {

	b.Run(testString("GenRing/", testContext.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			NewRing(testContext.ringQ.N, testContext.ringQ.Modulus)
		}
	})
}

func benchMarshalling(testContext *testParams, b *testing.B) {

	p := testContext.uniformSamplerQ.ReadNew()

	b.Run(testString("Marshalling/MarshalPoly/", testContext.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.MarshalBinary()
		}
	})

	data, _ := p.MarshalBinary()

	b.Run(testString("Marshalling/UnmarshalPoly/", testContext.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.UnmarshalBinary(data)
		}
	})
}

func benchSampling(testContext *testParams, b *testing.B) {

	pol := testContext.ringQ.NewPoly()

	b.Run(testString("Sampling/Gaussian/", testContext.ringQ), func(b *testing.B) {

		gaussianSampler := NewGaussianSampler(testContext.prng, testContext.ringQ, DefaultSigma, DefaultBound)

		for i := 0; i < b.N; i++ {
			gaussianSampler.ReadLvl(uint64(len(testContext.ringQ.Modulus)-1), pol)
		}
	})

	b.Run(testString("Sampling/Ternary/0.3/", testContext.ringQ), func(b *testing.B) {

		ternarySampler := NewTernarySampler(testContext.prng, testContext.ringQ, 1.0/3, true)

		for i := 0; i < b.N; i++ {
			ternarySampler.Read(pol)
		}
	})

	b.Run(testString("Sampling/Ternary/0.5/", testContext.ringQ), func(b *testing.B) {

		ternarySampler := NewTernarySampler(testContext.prng, testContext.ringQ, 0.5, true)

		for i := 0; i < b.N; i++ {
			ternarySampler.Read(pol)
		}
	})

	b.Run(testString("Sampling/Ternary/sparse128/", testContext.ringQ), func(b *testing.B) {

		ternarySampler := NewTernarySamplerSparse(testContext.prng, testContext.ringQ, 128, true)

		for i := 0; i < b.N; i++ {
			ternarySampler.Read(pol)
		}
	})

	b.Run(testString("Sampling/Uniform/", testContext.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			testContext.uniformSamplerQ.Read(pol)
		}
	})
}

func benchMontgomery(testContext *testParams, b *testing.B) {

	p := testContext.uniformSamplerQ.ReadNew()

	b.Run(testString("Montgomery/MForm/", testContext.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			testContext.ringQ.MForm(p, p)
		}
	})

	b.Run(testString("Montgomery/InvMForm/", testContext.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			testContext.ringQ.InvMForm(p, p)
		}
	})
}

func benchNTT(testContext *testParams, b *testing.B) {

	p := testContext.uniformSamplerQ.ReadNew()

	b.Run(testString("NTT/NTT/Montgomery/", testContext.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			testContext.ringQ.NTT(p, p)
		}
	})

	b.Run(testString("NTT/InvNTT/Montgomery/", testContext.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			testContext.ringQ.InvNTT(p, p)
		}
	})

	b.Run(testString("NTT/NTT/Barrett/", testContext.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			testContext.ringQ.NTTBarrett(p, p)
		}
	})

	b.Run(testString("NTT/InvNTT/Barrett/", testContext.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			testContext.ringQ.InvNTTBarrett(p, p)
		}
	})
}

func benchMurakami(testContext *testParams, b *testing.B) {

	ringQ := testContext.ringQ
	if err := ringQ.GenMurakamiParams(); err != nil {
		panic(err)
	}

	p := testContext.uniformSamplerQ.ReadNew()

	b.Run(testString("Murakami/Forward/", testContext.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ringQ.MapXX2NToXNAndMurakami(uint64(len(ringQ.Modulus)-1), p)
		}
	})

	b.Run(testString("Murakami/Backward/", testContext.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			ringQ.MapXNToXX2NAndMurakami(uint64(len(ringQ.Modulus)-1), p)
		}
	})
}

func benchMulCoeffs(testContext *testParams, b *testing.B) {

	p0 := testContext.uniformSamplerQ.ReadNew()
	p1 := testContext.uniformSamplerQ.ReadNew()

	b.Run(testString("MulCoeffs/Montgomery/", testContext.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			testContext.ringQ.MulCoeffsMontgomery(p0, p1, p0)
		}
	})

	b.Run(testString("MulCoeffs/MontgomeryConstant/", testContext.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			testContext.ringQ.MulCoeffsMontgomeryConstant(p0, p1, p0)
		}
	})

	b.Run(testString("MulCoeffs/Barrett/", testContext.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			testContext.ringQ.MulCoeffs(p0, p1, p0)
		}
	})

	b.Run(testString("MulCoeffs/BarrettConstant/", testContext.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			testContext.ringQ.MulCoeffsConstant(p0, p1, p0)
		}
	})
}

func benchAddCoeffs(testContext *testParams, b *testing.B) {

	p0 := testContext.uniformSamplerQ.ReadNew()
	p1 := testContext.uniformSamplerQ.ReadNew()

	b.Run(testString("AddCoeffs/Add/", testContext.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			testContext.ringQ.Add(p0, p1, p0)
		}
	})

	b.Run(testString("AddCoeffs/AddConstant/", testContext.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			testContext.ringQ.AddNoMod(p0, p1, p0)
		}
	})
}

func benchSubCoeffs(testContext *testParams, b *testing.B) {

	p0 := testContext.uniformSamplerQ.ReadNew()
	p1 := testContext.uniformSamplerQ.ReadNew()

	b.Run(testString("SubCoeffs/Sub/", testContext.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			testContext.ringQ.Sub(p0, p1, p0)
		}
	})

	b.Run(testString("SubCoeffs/SubConstant/", testContext.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			testContext.ringQ.SubNoMod(p0, p1, p0)
		}
	})
}

func benchNegCoeffs(testContext *testParams, b *testing.B) {

	p0 := testContext.uniformSamplerQ.ReadNew()

	b.Run(testString("NegCoeffs", testContext.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			testContext.ringQ.Neg(p0, p0)
		}
	})
}

func benchMulScalar(testContext *testParams, b *testing.B) {

	p := testContext.uniformSamplerQ.ReadNew()

	rand1 := RandUniform(testContext.prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)
	rand2 := RandUniform(testContext.prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)

	scalarBigint := NewUint(rand1)
	scalarBigint.Mul(scalarBigint, NewUint(rand2))

	b.Run(testString("MulScalar/uint64/", testContext.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			testContext.ringQ.MulScalar(p, rand1, p)
		}
	})

	b.Run(testString("MulScalar/big.Int/", testContext.ringQ), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			testContext.ringQ.MulScalarBigint(p, scalarBigint, p)
		}
	})
}

func benchExtendBasis(testContext *testParams, b *testing.B) {

	rescaleParams := make([]uint64, len(testContext.ringP.Modulus))

	for i, pi := range testContext.ringP.Modulus {
		rescaleParams[i] = RandUniform(testContext.prng, pi, (1<<uint64(bits.Len64(pi)-1) - 1))
	}

	basisExtender := NewFastBasisExtender(testContext.ringQ, testContext.ringP)

	p0 := testContext.uniformSamplerQ.ReadNew()
	p1 := testContext.uniformSamplerP.ReadNew()

	level := uint64(len(testContext.ringQ.Modulus) - 1)

	b.Run(fmt.Sprintf("ExtendBasis/ModUp/N=%d/limbsQ=%d/limbsP=%d", testContext.ringQ.N, len(testContext.ringQ.Modulus), len(testContext.ringP.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			basisExtender.ModUpSplitQP(level, p0, p1)
		}
	})

	b.Run(fmt.Sprintf("ExtendBasis/ModDown/N=%d/limbsQ=%d/limbsP=%d", testContext.ringQ.N, len(testContext.ringQ.Modulus), len(testContext.ringP.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			basisExtender.ModDownSplitPQ(level, p0, p1, p0)
		}
	})

	b.Run(fmt.Sprintf("ExtendBasis/ModDownNTT/N=%d/limbsQ=%d/limbsP=%d", testContext.ringQ.N, len(testContext.ringQ.Modulus), len(testContext.ringP.Modulus)), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			basisExtender.ModDownSplitNTTPQ(level, p0, p1, p0)
		}
	})
}

func benchDivByLastModulus(testContext *testParams, b *testing.B) {

	var p0 *Poly

	b.Run(testString("DivByLastModulus/Floor/", testContext.ringQ), func(b *testing.B) {

		for i := 0; i < b.N; i++ {

			b.StopTimer()
			p0 = testContext.uniformSamplerQ.ReadNew()
			b.StartTimer()

			testContext.ringQ.DivFloorByLastModulus(p0)
		}
	})

	b.Run(testString("DivByLastModulus/FloorNTT/", testContext.ringQ), func(b *testing.B) {

		for i := 0; i < b.N; i++ {

			b.StopTimer()
			p0 = testContext.uniformSamplerQ.ReadNew()
			b.StartTimer()

			testContext.ringQ.DivFloorByLastModulusNTT(p0)
		}
	})

	b.Run(testString("DivByLastModulus/Round/", testContext.ringQ), func(b *testing.B) {

		for i := 0; i < b.N; i++ {

			b.StopTimer()
			p0 = testContext.uniformSamplerQ.ReadNew()
			b.StartTimer()

			testContext.ringQ.DivRoundByLastModulus(p0)
		}
	})

	b.Run(testString("DivByLastModulus/RoundNTT/", testContext.ringQ), func(b *testing.B) {

		for i := 0; i < b.N; i++ {

			b.StopTimer()
			p0 = testContext.uniformSamplerQ.ReadNew()
			b.StartTimer()

			testContext.ringQ.DivRoundByLastModulusNTT(p0)
		}
	})
}

func benchDivByRNSBasis(testContext *testParams, b *testing.B) {

	b.Run(testString("DivByRNSBasis/Simple/DivByQOverTRounded/reconstructAndScale/", testContext.ringQ), func(b *testing.B) {

		rescaler := NewSimpleScaler(T, testContext.ringQ)

		coeffs := make([]*big.Int, testContext.ringQ.N)
		for i := uint64(0); i < testContext.ringQ.N; i++ {
			coeffs[i] = RandInt(testContext.ringQ.ModulusBigint)
		}

		tmp0 := testContext.ringQ.NewPoly()
		tmp1 := testContext.ringQ.NewPoly()

		testContext.ringQ.SetCoefficientsBigint(coeffs, tmp0)

		for i := 0; i < b.N; i++ {
			rescaler.reconstructAndScale(tmp0, tmp1)
		}
	})

	b.Run(testString("DivByRNSBasis/Simple/DivByQOverTRounded/reconstructThenScale/", testContext.ringQ), func(b *testing.B) {

		rescaler := NewSimpleScaler(T, testContext.ringQ)

		coeffs := make([]*big.Int, testContext.ringQ.N)
		for i := uint64(0); i < testContext.ringQ.N; i++ {
			coeffs[i] = RandInt(testContext.ringQ.ModulusBigint)
		}

		tmp0 := testContext.ringQ.NewPoly()
		tmp1 := testContext.ringQ.NewPoly()

		testContext.ringQ.SetCoefficientsBigint(coeffs, tmp0)

		for i := 0; i < b.N; i++ {
			rescaler.reconstructThenScale(tmp0, tmp1)
		}
	})

	b.Run(testString("DivByRNSBasis/RNS/DivByQOverTRounded/", testContext.ringQ), func(b *testing.B) {

		coeffs := make([]*big.Int, testContext.ringQ.N)
		for i := uint64(0); i < testContext.ringQ.N; i++ {
			coeffs[i] = RandInt(testContext.ringQ.ModulusBigint)
		}

		scaler := NewRNSScaler(T, testContext.ringQ)
		polyQ := testContext.ringQ.NewPoly()
		polyT := NewPoly(testContext.ringQ.N, 1)
		testContext.ringQ.SetCoefficientsBigint(coeffs, polyQ)

		for i := 0; i < b.N; i++ {
			scaler.DivByQOverTRounded(polyQ, polyT)
		}
	})
}

func benchBRed(testContext *testParams, b *testing.B) {

	q := uint64(1033576114481528833)
	u := BRedParams(q)

	x := rand.Uint64() % q
	y := rand.Uint64() % q

	b.ResetTimer()

	b.Run(fmt.Sprintf("BRed"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			x = BRed(x, y, q, u)
		}
	})
}

func benchMRed(testContext *testParams, b *testing.B) {

	q := uint64(1033576114481528833)

	x := rand.Uint64() % q
	y := rand.Uint64() % q

	u := BRedParams(q)

	y = MForm(y, q, u)

	m := MRedParams(q)

	b.ResetTimer()

	b.Run(fmt.Sprintf("MRed"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			x = MRed(x, y, q, m)
		}
	})
}

func benchBRedAdd(testContext *testParams, b *testing.B) {

	q := uint64(1033576114481528833)
	u := BRedParams(q)

	x := rand.Uint64()

	b.ResetTimer()

	b.Run(fmt.Sprintf("BRedAdd"), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			BRedAdd(x, q, u)
		}
	})
}
