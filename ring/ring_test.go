package ring

import (
	"flag"
	"fmt"
	"math/big"
	"testing"

	"github.com/ldsec/lattigo/v2/utils"

	"github.com/stretchr/testify/require"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters). Overrides -short and requires -timeout=0.")

var T = uint64(0x3ee0001)
var DefaultSigma = float64(3.2)
var DefaultBound = int(6 * DefaultSigma)

func testString(opname string, ringQ *Ring) string {
	return fmt.Sprintf("%sN=%d/limbs=%d", opname, ringQ.N, len(ringQ.Modulus))
}

type testParams struct {
	ringQ           *Ring
	ringP           *Ring
	prng            utils.PRNG
	uniformSamplerQ *UniformSampler
	uniformSamplerP *UniformSampler
}

func genTestParams(defaultParams Parameters) (testContext *testParams, err error) {

	testContext = new(testParams)

	if testContext.ringQ, err = NewRing(1<<defaultParams.logN, defaultParams.qi); err != nil {
		return nil, err
	}
	if testContext.ringP, err = NewRing(1<<defaultParams.logN, defaultParams.pi); err != nil {
		return nil, err
	}
	if testContext.prng, err = utils.NewPRNG(); err != nil {
		return nil, err
	}
	testContext.uniformSamplerQ = NewUniformSampler(testContext.prng, testContext.ringQ)
	testContext.uniformSamplerP = NewUniformSampler(testContext.prng, testContext.ringP)
	return
}

func TestRing(t *testing.T) {

	var err error

	var defaultParams = DefaultParams[0:4] // the default test
	if testing.Short() {
		defaultParams = DefaultParams[0:2] // the short test suite
	}
	if *flagLongTest {
		defaultParams = DefaultParams // the long test suite
	}

	testNewRing(t)
	for _, defaultParam := range defaultParams {

		var testContext *testParams
		if testContext, err = genTestParams(defaultParam); err != nil {
			t.Error(err)
		}
		testPRNG(testContext, t)
		testGenerateNTTPrimes(testContext, t)
		testImportExportPolyString(testContext, t)
		testDivFloorByLastModulusMany(testContext, t)
		testDivRoundByLastModulusMany(testContext, t)
		testMarshalBinary(testContext, t)
		testUniformSampler(testContext, t)
		testGaussianSampler(testContext, t)
		testTernarySampler(testContext, t)
		testGaloisShift(testContext, t)
		testModularReduction(testContext, t)
		testMForm(testContext, t)
		testMulScalarBigint(testContext, t)
		testExtendBasis(testContext, t)
		testScaling(testContext, t)
		testMultByMonomial(testContext, t)
	}
}

func testNewRing(t *testing.T) {
	t.Run("NewRing/", func(t *testing.T) {
		r, err := NewRing(0, nil)
		require.Nil(t, r)
		require.Error(t, err)

		r, err = NewRing(0, []uint64{})
		require.Nil(t, r)
		require.Error(t, err)

		r, err = NewRing(4, []uint64{})
		require.Nil(t, r)
		require.Error(t, err)

		r, err = NewRing(8, []uint64{})
		require.Nil(t, r)
		require.Error(t, err)

		r, err = NewRing(16, []uint64{7}) // Passing non NTT-enabling coeff modulus
		require.NotNil(t, r)              // Should still return a Ring instance
		require.Error(t, err)             // Should also return an error due to non NTT

		r, err = NewRing(16, []uint64{4}) // Passing non prime moduli
		require.NotNil(t, r)              // Should still return a Ring instance
		require.Error(t, err)             // Should also return an error due to non NTT

		r, err = NewRing(16, []uint64{97, 7}) // Passing a NTT-enabling and a non NTT-enabling coeff modulus
		require.NotNil(t, r)                  // Should still return a Ring instance
		require.Error(t, err)                 // Should also return an error due to non NTT

		r, err = NewRing(16, []uint64{97, 97}) // Passing non CRT-enabling coeff modulus
		require.Nil(t, r)                      // Should not return a Ring instance
		require.Error(t, err)

		r, err = NewRing(16, []uint64{97}) // Passing NTT-enabling coeff modulus
		require.NotNil(t, r)
		require.NoError(t, err)

	})
}

func testPRNG(testContext *testParams, t *testing.T) {

	sum := make([]byte, testContext.ringQ.N)
	t.Run(testString("PRNG/", testContext.ringQ), func(t *testing.T) {
		prng1, err := utils.NewKeyedPRNG(nil)
		if err != nil {
			panic(err)
		}
		prng2, err := utils.NewKeyedPRNG(nil)
		if err != nil {
			panic(err)
		}

		prng1.SetClock(sum, 256)
		prng2.SetClock(sum, 256)

		crsGenerator1 := NewUniformSampler(prng1, testContext.ringQ)
		crsGenerator2 := NewUniformSampler(prng2, testContext.ringQ)

		p0 := crsGenerator1.ReadNew()
		p1 := crsGenerator2.ReadNew()

		require.True(t, testContext.ringQ.Equal(p0, p1))
	})

}

func testGenerateNTTPrimes(testContext *testParams, t *testing.T) {

	t.Run(testString("GenerateNTTPrimes/", testContext.ringQ), func(t *testing.T) {

		primes := GenerateNTTPrimes(55, testContext.ringQ.N<<1, len(testContext.ringQ.Modulus))

		for _, q := range primes {
			require.Equal(t, q&uint64((testContext.ringQ.N<<1)-1), uint64(1))
			require.True(t, IsPrime(q), q)
		}
	})
}

func testImportExportPolyString(testContext *testParams, t *testing.T) {

	t.Run(testString("ImportExportPolyString/", testContext.ringQ), func(t *testing.T) {

		p0 := testContext.uniformSamplerQ.ReadNew()
		p1 := testContext.ringQ.NewPoly()

		testContext.ringQ.SetCoefficientsString(testContext.ringQ.PolyToString(p0), p1)

		require.True(t, testContext.ringQ.Equal(p0, p1))
	})
}

func testDivFloorByLastModulusMany(testContext *testParams, t *testing.T) {

	t.Run(testString("DivFloorByLastModulusMany/", testContext.ringQ), func(t *testing.T) {

		coeffs := make([]*big.Int, testContext.ringQ.N)
		for i := 0; i < testContext.ringQ.N; i++ {
			coeffs[i] = RandInt(testContext.ringQ.ModulusBigint)
			coeffs[i].Quo(coeffs[i], NewUint(10))
		}

		nbRescals := len(testContext.ringQ.Modulus) - 1

		coeffsWant := make([]*big.Int, testContext.ringQ.N)
		for i := range coeffs {
			coeffsWant[i] = new(big.Int).Set(coeffs[i])
			for j := 0; j < nbRescals; j++ {
				coeffsWant[i].Quo(coeffsWant[i], NewUint(testContext.ringQ.Modulus[len(testContext.ringQ.Modulus)-1-j]))
			}
		}

		polTest0 := testContext.ringQ.NewPoly()
		polTest1 := testContext.ringQ.NewPoly()
		polWant := testContext.ringQ.NewPoly()

		testContext.ringQ.SetCoefficientsBigint(coeffs, polTest0)
		testContext.ringQ.SetCoefficientsBigint(coeffsWant, polWant)

		testContext.ringQ.DivFloorByLastModulusMany(polTest0, polTest1, nbRescals)
		for i := 0; i < testContext.ringQ.N; i++ {
			for j := 0; j < len(testContext.ringQ.Modulus)-nbRescals; j++ {
				require.Equalf(t, polWant.Coeffs[j][i], polTest1.Coeffs[j][i], "coeff %v Qi%v = %s", i, j, coeffs[i].String())
			}
		}
	})
}

func testDivRoundByLastModulusMany(testContext *testParams, t *testing.T) {

	t.Run(testString("DivRoundByLastModulusMany/", testContext.ringQ), func(t *testing.T) {

		coeffs := make([]*big.Int, testContext.ringQ.N)
		for i := 0; i < testContext.ringQ.N; i++ {
			coeffs[i] = RandInt(testContext.ringQ.ModulusBigint)
			coeffs[i].Quo(coeffs[i], NewUint(10))
		}

		nbRescals := len(testContext.ringQ.Modulus) - 1

		coeffsWant := make([]*big.Int, testContext.ringQ.N)
		for i := range coeffs {
			coeffsWant[i] = new(big.Int).Set(coeffs[i])
			for j := 0; j < nbRescals; j++ {
				DivRound(coeffsWant[i], NewUint(testContext.ringQ.Modulus[len(testContext.ringQ.Modulus)-1-j]), coeffsWant[i])
			}
		}

		polTest0 := testContext.ringQ.NewPoly()
		polTest1 := testContext.ringQ.NewPoly()
		polWant := testContext.ringQ.NewPoly()

		testContext.ringQ.SetCoefficientsBigint(coeffs, polTest0)
		testContext.ringQ.SetCoefficientsBigint(coeffsWant, polWant)

		testContext.ringQ.DivRoundByLastModulusMany(polTest0, polTest1, nbRescals)
		for i := 0; i < testContext.ringQ.N; i++ {
			for j := 0; j < len(testContext.ringQ.Modulus)-nbRescals; j++ {
				require.Equalf(t, polWant.Coeffs[j][i], polTest1.Coeffs[j][i], "coeff %v Qi%v = %s", i, j, coeffs[i].String())
			}
		}
	})
}

func testMarshalBinary(testContext *testParams, t *testing.T) {

	t.Run(testString("MarshalBinary/Ring/", testContext.ringQ), func(t *testing.T) {

		data, _ := testContext.ringQ.MarshalBinary()

		ringQTest := new(Ring)
		ringQTest.UnmarshalBinary(data)

		require.Equal(t, ringQTest.N, testContext.ringQ.N)
		require.Equal(t, ringQTest.Modulus, testContext.ringQ.Modulus)
	})

	t.Run(testString("MarshalBinary/Poly/", testContext.ringQ), func(t *testing.T) {

		p := testContext.uniformSamplerQ.ReadNew()
		pTest := testContext.ringQ.NewPoly()

		data, _ := p.MarshalBinary()

		_ = pTest.UnmarshalBinary(data)

		for i := range testContext.ringQ.Modulus {
			require.Equal(t, p.Coeffs[i][:testContext.ringQ.N], pTest.Coeffs[i][:testContext.ringQ.N])
		}
	})
}

func testUniformSampler(testContext *testParams, t *testing.T) {

	t.Run(testString("UniformSampler/Read/", testContext.ringQ), func(t *testing.T) {
		pol := testContext.ringQ.NewPoly()
		testContext.uniformSamplerQ.Read(pol)
		for i := 0; i < testContext.ringQ.N; i++ {
			for j, qi := range testContext.ringQ.Modulus {
				require.False(t, pol.Coeffs[j][i] > qi)
			}
		}
	})

	t.Run(testString("UniformSampler/ReadNew/", testContext.ringQ), func(t *testing.T) {
		pol := testContext.uniformSamplerQ.ReadNew()
		for i := 0; i < testContext.ringQ.N; i++ {
			for j, qi := range testContext.ringQ.Modulus {
				require.False(t, pol.Coeffs[j][i] > qi)
			}
		}
	})
}

func testGaussianSampler(testContext *testParams, t *testing.T) {

	t.Run(testString("GaussianSampler/", testContext.ringQ), func(t *testing.T) {
		gaussianSampler := NewGaussianSampler(testContext.prng, testContext.ringQ, DefaultSigma, DefaultBound)
		pol := gaussianSampler.ReadNew()

		for i := 0; i < testContext.ringQ.N; i++ {
			for j, qi := range testContext.ringQ.Modulus {
				require.False(t, uint64(DefaultBound) < pol.Coeffs[j][i] && pol.Coeffs[j][i] < (qi-uint64(DefaultBound)))
			}
		}
	})
}

func testTernarySampler(testContext *testParams, t *testing.T) {

	for _, p := range []float64{.5, 1. / 3., 128. / 65536.} {
		t.Run(testString(fmt.Sprintf("TernarySampler/p=%1.2f/", p), testContext.ringQ), func(t *testing.T) {

			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}
			ternarySampler := NewTernarySampler(prng, testContext.ringQ, p, false)

			pol := ternarySampler.ReadNew()
			for i, mod := range testContext.ringQ.Modulus {
				minOne := mod - 1
				for _, c := range pol.Coeffs[i] {
					require.True(t, c == 0 || c == minOne || c == 1)
				}
			}
		})
	}

	for _, p := range []int{0, 64, 96, 128, 256} {
		t.Run(testString(fmt.Sprintf("TernarySampler/hw=%d/", p), testContext.ringQ), func(t *testing.T) {

			prng, err := utils.NewPRNG()
			if err != nil {
				panic(err)
			}

			ternarySampler := NewTernarySamplerSparse(prng, testContext.ringQ, p, false)

			pol := ternarySampler.ReadNew()

			for i := range testContext.ringQ.Modulus {
				hw := 0
				for _, c := range pol.Coeffs[i] {
					if c != 0 {
						hw++
					}
				}
				require.True(t, hw == p)
			}
		})
	}
}

func testModularReduction(testContext *testParams, t *testing.T) {

	t.Run(testString("ModularReduction/BRed/", testContext.ringQ), func(t *testing.T) {

		var x, y uint64
		var bigQ, result *big.Int

		for j, q := range testContext.ringQ.Modulus {

			bigQ = NewUint(q)

			bredParams := testContext.ringQ.BredParams[j]

			x = 1
			y = 1

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, BRed(x, y, q, bredParams), result.Uint64(), "x = %v, y=%v", x, y)

			x = 1
			y = q - 1

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, BRed(x, y, q, bredParams), result.Uint64(), "x = %v, y=%v", x, y)

			x = 1
			y = 0xFFFFFFFFFFFFFFFF

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, BRed(x, y, q, bredParams), result.Uint64(), "x = %v, y=%v", x, y)

			x = q - 1
			y = q - 1

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, BRed(x, y, q, bredParams), result.Uint64(), "x = %v, y=%v", x, y)

			x = q - 1
			y = 0xFFFFFFFFFFFFFFFF

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, BRed(x, y, q, bredParams), result.Uint64(), "x = %v, y=%v", x, y)

			x = 0xFFFFFFFFFFFFFFFF
			y = 0xFFFFFFFFFFFFFFFF

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, BRed(x, y, q, bredParams), result.Uint64(), "x = %v, y=%v", x, y)
		}
	})

	t.Run(testString("ModularReduction/MRed/", testContext.ringQ), func(t *testing.T) {

		var x, y uint64
		var bigQ, result *big.Int

		for j, q := range testContext.ringQ.Modulus {

			bigQ = NewUint(q)

			bredParams := testContext.ringQ.BredParams[j]
			mredparams := testContext.ringQ.MredParams[j]

			x = 1
			y = 1

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, MRed(x, MForm(y, q, bredParams), q, mredparams), result.Uint64(), "x = %v, y=%v", x, y)

			x = 1
			y = q - 1

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, MRed(x, MForm(y, q, bredParams), q, mredparams), result.Uint64(), "x = %v, y=%v", x, y)

			x = 1
			y = 0xFFFFFFFFFFFFFFFF

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, MRed(x, MForm(y, q, bredParams), q, mredparams), result.Uint64(), "x = %v, y=%v", x, y)

			x = q - 1
			y = q - 1

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, MRed(x, MForm(y, q, bredParams), q, mredparams), result.Uint64(), "x = %v, y=%v", x, y)

			x = q - 1
			y = 0xFFFFFFFFFFFFFFFF

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, MRed(x, MForm(y, q, bredParams), q, mredparams), result.Uint64(), "x = %v, y=%v", x, y)

			x = 0xFFFFFFFFFFFFFFFF
			y = 0xFFFFFFFFFFFFFFFF

			result = NewUint(x)
			result.Mul(result, NewUint(y))
			result.Mod(result, bigQ)

			require.Equalf(t, MRed(x, MForm(y, q, bredParams), q, mredparams), result.Uint64(), "x = %v, y=%v", x, y)
		}
	})
}

func testGaloisShift(testContext *testParams, t *testing.T) {

	t.Run(testString("GaloisShift/", testContext.ringQ), func(t *testing.T) {

		pWant := testContext.uniformSamplerQ.ReadNew()
		pTest := pWant.CopyNew()

		testContext.ringQ.BitReverse(pTest, pTest)
		testContext.ringQ.InvNTT(pTest, pTest)
		testContext.ringQ.Rotate(pTest, 1, pTest)
		testContext.ringQ.NTT(pTest, pTest)
		testContext.ringQ.BitReverse(pTest, pTest)
		testContext.ringQ.Reduce(pTest, pTest)

		testContext.ringQ.Shift(pWant, 1, pWant)

		for i := range testContext.ringQ.Modulus {
			require.Equal(t, pTest.Coeffs[i][:testContext.ringQ.N], pWant.Coeffs[i][:testContext.ringQ.N])
		}
	})
}

func testMForm(testContext *testParams, t *testing.T) {

	t.Run(testString("MForm/", testContext.ringQ), func(t *testing.T) {

		polWant := testContext.uniformSamplerQ.ReadNew()
		polTest := testContext.ringQ.NewPoly()

		testContext.ringQ.MForm(polWant, polTest)
		testContext.ringQ.InvMForm(polTest, polTest)

		require.True(t, testContext.ringQ.Equal(polWant, polTest))
	})
}

func testMulScalarBigint(testContext *testParams, t *testing.T) {

	t.Run(testString("MulScalarBigint/", testContext.ringQ), func(t *testing.T) {

		polWant := testContext.uniformSamplerQ.ReadNew()
		polTest := polWant.CopyNew()

		rand1 := RandUniform(testContext.prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)
		rand2 := RandUniform(testContext.prng, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF)

		scalarBigint := NewUint(rand1)
		scalarBigint.Mul(scalarBigint, NewUint(rand2))

		testContext.ringQ.MulScalar(polWant, rand1, polWant)
		testContext.ringQ.MulScalar(polWant, rand2, polWant)
		testContext.ringQ.MulScalarBigint(polTest, scalarBigint, polTest)

		require.True(t, testContext.ringQ.Equal(polWant, polTest))
	})
}

func testExtendBasis(testContext *testParams, t *testing.T) {

	t.Run(testString("ExtendBasis/", testContext.ringQ), func(t *testing.T) {

		basisextender := NewFastBasisExtender(testContext.ringQ, testContext.ringP)

		coeffs := make([]*big.Int, testContext.ringQ.N)
		for i := 0; i < testContext.ringQ.N; i++ {
			coeffs[i] = RandInt(testContext.ringQ.ModulusBigint)
		}

		Pol := testContext.ringQ.NewPoly()
		PolTest := testContext.ringP.NewPoly()
		PolWant := testContext.ringP.NewPoly()

		testContext.ringQ.SetCoefficientsBigint(coeffs, Pol)
		testContext.ringP.SetCoefficientsBigint(coeffs, PolWant)

		basisextender.ModUpSplitQP(len(testContext.ringQ.Modulus)-1, Pol, PolTest)

		testContext.ringP.Reduce(PolTest, PolTest)

		for i := range testContext.ringP.Modulus {
			require.Equal(t, PolTest.Coeffs[i][:testContext.ringQ.N], PolWant.Coeffs[i][:testContext.ringQ.N])
		}
	})
}

func testScaling(testContext *testParams, t *testing.T) {

	t.Run(testString("Scaling/Simple/", testContext.ringQ), func(t *testing.T) {

		rescaler := NewSimpleScaler(T, testContext.ringQ)

		coeffs := make([]*big.Int, testContext.ringQ.N)
		for i := 0; i < testContext.ringQ.N; i++ {
			coeffs[i] = RandInt(testContext.ringQ.ModulusBigint)
		}

		coeffsWant := make([]*big.Int, testContext.ringQ.N)
		for i := range coeffs {
			coeffsWant[i] = new(big.Int).Set(coeffs[i])
			coeffsWant[i].Mul(coeffsWant[i], NewUint(T))
			DivRound(coeffsWant[i], testContext.ringQ.ModulusBigint, coeffsWant[i])
			coeffsWant[i].Mod(coeffsWant[i], NewUint(T))
		}

		PolTest := testContext.ringQ.NewPoly()

		testContext.ringQ.SetCoefficientsBigint(coeffs, PolTest)

		rescaler.DivByQOverTRounded(PolTest, PolTest)

		for i := 0; i < testContext.ringQ.N; i++ {
			require.Equal(t, PolTest.Coeffs[0][i], coeffsWant[i].Uint64())
		}
	})

	t.Run(testString("Scaling/RNS", testContext.ringQ), func(t *testing.T) {

		scaler := NewRNSScaler(T, testContext.ringQ)

		coeffs := make([]*big.Int, testContext.ringQ.N)
		for i := 0; i < testContext.ringQ.N; i++ {
			coeffs[i] = RandInt(testContext.ringQ.ModulusBigint)
		}

		coeffsWant := make([]*big.Int, testContext.ringQ.N)
		for i := range coeffs {
			coeffsWant[i] = new(big.Int).Set(coeffs[i])
			coeffsWant[i].Mul(coeffsWant[i], NewUint(T))
			DivRound(coeffsWant[i], testContext.ringQ.ModulusBigint, coeffsWant[i])
			coeffsWant[i].Mod(coeffsWant[i], NewUint(T))
		}

		polyQ := testContext.ringQ.NewPoly()
		polyT := NewPoly(testContext.ringQ.N, 1)
		testContext.ringQ.SetCoefficientsBigint(coeffs, polyQ)

		scaler.DivByQOverTRounded(polyQ, polyT)

		for i := 0; i < testContext.ringQ.N; i++ {
			require.Equal(t, polyT.Coeffs[0][i], coeffsWant[i].Uint64())
		}
	})
}

func testMultByMonomial(testContext *testParams, t *testing.T) {

	t.Run(testString("MultByMonomial/", testContext.ringQ), func(t *testing.T) {

		p1 := testContext.uniformSamplerQ.ReadNew()

		p3Test := testContext.ringQ.NewPoly()
		p3Want := testContext.ringQ.NewPoly()

		testContext.ringQ.MultByMonomial(p1, 1, p3Test)
		testContext.ringQ.MultByMonomial(p3Test, 8, p3Test)

		testContext.ringQ.MultByMonomial(p1, 9, p3Want)

		require.Equal(t, p3Want.Coeffs[0][:testContext.ringQ.N], p3Test.Coeffs[0][:testContext.ringQ.N])
	})
}
