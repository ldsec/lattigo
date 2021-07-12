package rlwe

import (
	"encoding/json"
	"flag"
	"fmt"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/stretchr/testify/require"
	"math"
	"math/big"
	"math/bits"
	"runtime"
	"testing"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters + secure bootstrapping). Overrides -short and requires -timeout=0.")
var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")

var (
	// PN12QP109 is a set of default parameters with logN=12 and logQP=109
	TestPN12QP109 = ParametersLiteral{
		LogN:  12,
		Q:     []uint64{0x7ffffec001, 0x40002001}, // 39 + 39 bits
		P:     []uint64{0x8000016001},             // 30 bits
		Sigma: DefaultSigma,
	}
	// PN13QP218 is a set of default parameters with logN=13 and logQP=218
	TestPN13QP218 = ParametersLiteral{
		LogN:  13,
		Q:     []uint64{0x3fffffffef8001, 0x4000000011c001, 0x40000000120001}, // 54 + 54 + 54 bits
		P:     []uint64{0x7ffffffffb4001},                                     // 55 bits
		Sigma: DefaultSigma,
	}

	// PN14QP438 is a set of default parameters with logN=14 and logQP=438
	TestPN14QP438 = ParametersLiteral{
		LogN: 14,
		Q: []uint64{0x100000000060001, 0x80000000068001, 0x80000000080001,
			0x3fffffffef8001, 0x40000000120001, 0x3fffffffeb8001}, // 56 + 55 + 55 + 54 + 54 + 54 bits
		P:     []uint64{0x80000000130001, 0x7fffffffe90001}, // 55 + 55 bits
		Sigma: DefaultSigma,
	}

	// PN15QP880 is a set of default parameters with logN=15 and logQP=880
	TestPN15QP880 = ParametersLiteral{
		LogN: 15,
		Q: []uint64{0x7ffffffffe70001, 0x7ffffffffe10001, 0x7ffffffffcc0001, // 59 + 59 + 59 bits
			0x400000000270001, 0x400000000350001, 0x400000000360001, // 58 + 58 + 58 bits
			0x3ffffffffc10001, 0x3ffffffffbe0001, 0x3ffffffffbd0001, // 58 + 58 + 58 bits
			0x4000000004d0001, 0x400000000570001, 0x400000000660001}, // 58 + 58 + 58 bits
		P:     []uint64{0xffffffffffc0001, 0x10000000001d0001, 0x10000000006e0001}, // 60 + 60 + 60 bits
		Sigma: DefaultSigma,
	}
)

// TestParams is a set of test parameters for the correctness of the rlwe pacakge.
var TestParams = []ParametersLiteral{TestPN12QP109, TestPN13QP218, TestPN14QP438, TestPN15QP880}

func testString(params Parameters, opname string) string {
	return fmt.Sprintf("%slogN=%d/logQ=%d/logP=%d/#Qi=%d/#Pi=%d",
		opname,
		params.LogN(),
		params.LogQ(),
		params.LogP(),
		params.QCount(),
		params.PCount())
}

func TestRLWE(t *testing.T) {
	defaultParams := TestParams // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15
	if testing.Short() {
		defaultParams = TestParams[:1] // the short test suite runs for ring degree N=2^12, 2^13
	}
	if *flagLongTest {
		defaultParams = append(TestParams) //, DefaultPostQuantumParams...) // the long test suite runs for all default parameters
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

		kgen := NewKeyGenerator(params)

		for _, testSet := range []func(kgen KeyGenerator, t *testing.T){
			testGenKeyPair,
			testSwitchKeyGen,
			testEncryptor,
			testDecryptor,
			testKeySwitcher,
		} {
			testSet(kgen, t)
			runtime.GC()
		}
	}
}

// Returns the ceil(log2) of the sum of the absolute value of all the coefficients
func log2OfInnerSum(level int, ringQ, ringP *ring.Ring, polQ, polP *ring.Poly) (logSum int) {
	sumRNSQ := make([]uint64, level+1)
	var sum uint64
	for i := 0; i < level+1; i++ {

		qi := ringQ.Modulus[i]
		qiHalf := qi >> 1
		coeffs := polQ.Coeffs[i]
		sum = 0

		for j := 0; j < ringQ.N; j++ {

			v := coeffs[j]

			if v >= qiHalf {
				sum = ring.CRed(sum+qi-v, qi)
			} else {
				sum = ring.CRed(sum+v, qi)
			}
		}

		sumRNSQ[i] = sum
	}

	var sumRNSP []uint64
	if polP != nil {

		sumRNSP = make([]uint64, len(ringP.Modulus))

		for i := 0; i < len(ringP.Modulus); i++ {
			pi := ringP.Modulus[i]
			piHalf := pi >> 1
			coeffs := polP.Coeffs[i]
			sum = 0

			for j := 0; j < ringQ.N; j++ {

				v := coeffs[j]

				if v >= piHalf {
					sum = ring.CRed(sum+pi-v, pi)
				} else {
					sum = ring.CRed(sum+v, pi)
				}
			}

			sumRNSP[i] = sum
		}
	}

	var smallNorm = true
	for i := 1; i < level+1; i++ {
		smallNorm = smallNorm && (sumRNSQ[0] == sumRNSQ[i])
	}

	if polP != nil {
		for i := 0; i < len(ringP.Modulus); i++ {
			smallNorm = smallNorm && (sumRNSQ[0] == sumRNSP[i])
		}
	}

	if !smallNorm {
		var qi, pi uint64
		var crtReconstruction *big.Int

		QPBig := new(big.Int).Set(ringQ.ModulusBigint)

		if polP != nil {
			QPBig.Mul(QPBig, ringP.ModulusBigint)
		}

		sumBigInt := ring.NewUint(0)
		QiB := new(big.Int)
		PiB := new(big.Int)
		tmp := new(big.Int)
		modulusBigint := ring.NewUint(1)

		for i := 0; i < level+1; i++ {

			qi = ringQ.Modulus[i]
			QiB.SetUint64(qi)

			modulusBigint.Mul(modulusBigint, QiB)

			crtReconstruction = new(big.Int)
			crtReconstruction.Quo(QPBig, QiB)
			tmp.ModInverse(crtReconstruction, QiB)
			tmp.Mod(tmp, QiB)
			crtReconstruction.Mul(crtReconstruction, tmp)

			sumBigInt.Add(sumBigInt, tmp.Mul(ring.NewUint(sumRNSQ[i]), crtReconstruction))
		}

		if polP != nil {
			for i := 0; i < len(ringP.Modulus); i++ {

				pi = ringP.Modulus[i]
				PiB.SetUint64(pi)

				modulusBigint.Mul(modulusBigint, PiB)

				crtReconstruction = new(big.Int)
				crtReconstruction.Quo(QPBig, PiB)
				tmp.ModInverse(crtReconstruction, PiB)
				tmp.Mod(tmp, PiB)
				crtReconstruction.Mul(crtReconstruction, tmp)

				sumBigInt.Add(sumBigInt, tmp.Mul(ring.NewUint(sumRNSP[i]), crtReconstruction))
			}
		}

		sumBigInt.Mod(sumBigInt, modulusBigint)

		logSum = sumBigInt.BitLen()
	} else {
		logSum = bits.Len64(sumRNSQ[0])
	}

	return

}

func testGenKeyPair(kgen KeyGenerator, t *testing.T) {

	params := kgen.(*keyGenerator).params

	// Checks that sum([-as + e, a] + [as])) <= N * 6 * sigma
	t.Run(testString(params, "PKGen/"), func(t *testing.T) {

		ringQ := params.RingQ()
		ringP := params.RingP()

		sk, pk := kgen.GenKeyPair()

		// [-as + e] + [as]
		ringQ.MulCoeffsMontgomeryAndAdd(sk.Value[0], pk.Value[1][0], pk.Value[0][0])
		ringP.MulCoeffsMontgomeryAndAdd(sk.Value[1], pk.Value[1][1], pk.Value[0][1])
		ringQ.InvNTT(pk.Value[0][0], pk.Value[0][0])
		ringP.InvNTT(pk.Value[0][1], pk.Value[0][1])

		log2Bound := bits.Len64(uint64(math.Floor(DefaultSigma*6)) * uint64(params.N()))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(pk.Value[0][0].Level(), ringQ, ringP, pk.Value[0][0], pk.Value[0][1]))
	})
}

func testSwitchKeyGen(kgen KeyGenerator, t *testing.T) {

	params := kgen.(*keyGenerator).params

	// Checks that switching keys are en encryption under the output key
	// of the RNS decomposition of the input key by
	// 1) Decrypting the RNS decomposed input key
	// 2) Reconstructing the key
	// 3) Checking that the difference with the input key has a small norm
	t.Run(testString(params, "SWKGen/"), func(t *testing.T) {

		ringQ := params.RingQ()
		ringP := params.RingP()
		skIn := kgen.GenSecretKey()
		skOut := kgen.GenSecretKey()

		// Generates Decomp([-asIn + w*P*sOut + e, a])
		swk := NewSwitchingKey(params)
		kgen.(*keyGenerator).newSwitchingKey(skIn.Value, skOut.Value, swk)

		// Decrypts
		// [-asIn + w*P*sOut + e, a] + [asIn]
		for j := range swk.Value {
			ringQ.MulCoeffsMontgomeryAndAdd(swk.Value[j][1][0], skOut.Value[0], swk.Value[j][0][0])
			ringP.MulCoeffsMontgomeryAndAdd(swk.Value[j][1][1], skOut.Value[1], swk.Value[j][0][1])
		}

		polyQ := swk.Value[0][0][0]
		polyP := swk.Value[0][0][1]

		// Sums all basis together (equivalent to multiplying with CRT decomposition of 1)
		// sum([1]_w * [w*P*sOut + e]) = P*sOut + sum(e)
		for j := range swk.Value {
			if j > 0 {
				ringQ.Add(polyQ, swk.Value[j][0][0], polyQ)
				ringP.Add(polyP, swk.Value[j][0][1], polyP)
			}
		}

		// sOut * P
		ringQ.MulScalarBigint(skIn.Value[0], ringP.ModulusBigint, skIn.Value[0])
		ringP.MulScalarBigint(skIn.Value[1], ringP.ModulusBigint, skIn.Value[1])

		// P*s^i + sum(e) - P*s^i = sum(e)
		ringQ.Sub(polyQ, skIn.Value[0], polyQ)
		ringP.Sub(polyP, skIn.Value[1], polyP)

		// Checks that the error is below the bound
		ringQ.InvNTT(polyQ, polyQ)
		ringP.InvNTT(polyP, polyP)
		ringQ.InvMForm(polyQ, polyQ)
		ringP.InvMForm(polyP, polyP)

		log2Bound := bits.Len64(uint64(math.Floor(DefaultSigma*6)) * uint64(params.N()*len(swk.Value)))
		require.GreaterOrEqual(t, log2Bound, log2OfInnerSum(len(ringQ.Modulus)-1, ringQ, ringP, polyQ, polyP))
	})
}

func testEncryptor(kgen KeyGenerator, t *testing.T) {

	params := kgen.(*keyGenerator).params

	sk, pk := kgen.GenKeyPair()

	ringQ := params.RingQ()

	t.Run(testString(params, "Encrypt/Pk/Fast/MaxLevel/"), func(t *testing.T) {
		plaintext := NewPlaintext(params, params.MaxLevel())
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptorFromPk(params, pk)
		ciphertext := encryptor.EncryptFastNTTNew(plaintext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], sk.Value[0], ciphertext.Value[0])
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value[0], ciphertext.Value[0])
		require.GreaterOrEqual(t, 12+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, nil, ciphertext.Value[0], nil))
	})

	t.Run(testString(params, "Encrypt/Pk/Fast/MinLevel/"), func(t *testing.T) {
		plaintext := NewPlaintext(params, 0)
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptorFromPk(params, pk)
		ciphertext := encryptor.EncryptFastNTTNew(plaintext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], sk.Value[0], ciphertext.Value[0])
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value[0], ciphertext.Value[0])
		require.GreaterOrEqual(t, 12+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, nil, ciphertext.Value[0], nil))
	})

	t.Run(testString(params, "Encrypt/Pk/Slow/MaxLevel/"), func(t *testing.T) {
		if params.PCount() == 0 {
			t.Skip()
		}
		plaintext := NewPlaintext(params, params.MaxLevel())
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptorFromPk(params, pk)
		ciphertext := encryptor.EncryptNTTNew(plaintext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], sk.Value[0], ciphertext.Value[0])
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value[0], ciphertext.Value[0])
		require.GreaterOrEqual(t, 9+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, nil, ciphertext.Value[0], nil))
	})

	t.Run(testString(params, "Encrypt/Pk/Slow/MinLevel/"), func(t *testing.T) {
		if params.PCount() == 0 {
			t.Skip()
		}
		plaintext := NewPlaintext(params, 0)
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptorFromPk(params, pk)
		ciphertext := encryptor.EncryptNTTNew(plaintext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], sk.Value[0], ciphertext.Value[0])
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value[0], ciphertext.Value[0])
		require.GreaterOrEqual(t, 9+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, nil, ciphertext.Value[0], nil))
	})

	t.Run(testString(params, "Encrypt/Sk/MaxLevel/"), func(t *testing.T) {
		plaintext := NewPlaintext(params, params.MaxLevel())
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptorFromSk(params, sk)
		ciphertext := encryptor.EncryptNTTNew(plaintext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], sk.Value[0], ciphertext.Value[0])
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value[0], ciphertext.Value[0])
		require.GreaterOrEqual(t, 5+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, nil, ciphertext.Value[0], nil))
	})

	t.Run(testString(params, "Encrypt/Sk/MinLevel/"), func(t *testing.T) {
		plaintext := NewPlaintext(params, 0)
		plaintext.Value.IsNTT = true
		encryptor := NewEncryptorFromSk(params, sk)
		ciphertext := encryptor.EncryptNTTNew(plaintext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], sk.Value[0], ciphertext.Value[0])
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value[0], ciphertext.Value[0])
		require.GreaterOrEqual(t, 5+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, nil, ciphertext.Value[0], nil))
	})
}

func testDecryptor(kgen KeyGenerator, t *testing.T) {
	params := kgen.(*keyGenerator).params
	sk := kgen.GenSecretKey()
	ringQ := params.RingQ()
	encryptor := NewEncryptorFromSk(params, sk)
	decryptor := NewDecryptor(params, sk)

	t.Run(testString(params, "Decrypt/MaxLevel/"), func(t *testing.T) {
		plaintext := NewPlaintext(params, params.MaxLevel())
		plaintext.Value.IsNTT = true
		ciphertext := encryptor.EncryptNTTNew(plaintext)
		plaintext = decryptor.DecryptNTTNew(ciphertext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.InvNTTLvl(plaintext.Level(), plaintext.Value, plaintext.Value)
		require.GreaterOrEqual(t, 5+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, nil, plaintext.Value, nil))
	})

	t.Run(testString(params, "Encrypt/MinLevel/"), func(t *testing.T) {
		plaintext := NewPlaintext(params, 0)
		plaintext.Value.IsNTT = true
		ciphertext := encryptor.EncryptNTTNew(plaintext)
		plaintext = decryptor.DecryptNTTNew(ciphertext)
		require.Equal(t, plaintext.Level(), ciphertext.Level())
		ringQ.InvNTTLvl(plaintext.Level(), plaintext.Value, plaintext.Value)
		require.GreaterOrEqual(t, 5+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, nil, plaintext.Value, nil))
	})
}

func testKeySwitcher(kgen KeyGenerator, t *testing.T) {

	params := kgen.(*keyGenerator).params
	sk := kgen.GenSecretKey()
	skOut := kgen.GenSecretKey()
	ks := NewKeySwitcher(params)

	ringQP := params.RingQP()
	ringQ := params.RingQ()

	plaintext := NewPlaintext(params, params.MaxLevel())
	plaintext.Value.IsNTT = true
	encryptor := NewEncryptorFromSk(params, sk)
	ciphertext := encryptor.EncryptNTTNew(plaintext)

	// Tests that a random polynomial decomposed is equal to its
	// reconstruction mod each RNS
	t.Run(testString(params, "DecomposeNTT/"), func(t *testing.T) {

		c2 := ciphertext.Value[1]

		ks.DecomposeNTT(ciphertext.Level(), c2, ks.PoolDecompQ, ks.PoolDecompP)

		coeffsBigintHave := make([]*big.Int, ringQ.N)
		coeffsBigintRef := make([]*big.Int, ringQ.N)
		coeffsBigintWant := make([]*big.Int, ringQ.N)

		for i := range coeffsBigintRef {
			coeffsBigintHave[i] = new(big.Int)
			coeffsBigintRef[i] = new(big.Int)
			coeffsBigintWant[i] = new(big.Int)
		}

		ringQ.PolyToBigintCenteredLvl(len(ringQ.Modulus)-1, c2, coeffsBigintRef)

		for i := 0; i < len(ks.PoolDecompQ); i++ {

			// Compute q_alpha_i in bigInt
			modulus := ring.NewInt(1)

			for j := 0; j < params.PCount(); j++ {
				idx := i*params.PCount() + j
				if idx > params.QCount()-1 {
					break
				}
				modulus.Mul(modulus, ring.NewUint(ringQ.Modulus[idx]))
			}

			// Reconstruct the decomposed polynomial
			polyQP := new(ring.Poly)
			polyQP.Coeffs = append(ks.PoolDecompQ[i].Coeffs, ks.PoolDecompP[i].Coeffs...)
			ringQP.PolyToBigintCenteredLvl(len(ringQP.Modulus)-1, polyQP, coeffsBigintHave)

			// Checks that Reconstruct(NTT(c2 mod Q)) mod q_alpha_i == Reconstruct(NTT(Decomp(c2 mod Q, q_alpha-i) mod QP))
			for i := range coeffsBigintWant {
				coeffsBigintHave[i].Mod(coeffsBigintHave[i], modulus)
				coeffsBigintWant[i].Mod(coeffsBigintRef[i], modulus)
				require.Equal(t, coeffsBigintHave[i].Cmp(coeffsBigintWant[i]), 0)
			}
		}
	})

	// Test that Dec(KS(Enc(ct, sk), skOut), skOut) has a small norm
	t.Run(testString(params, "KeySwitch/"), func(t *testing.T) {
		swk := kgen.GenSwitchingKey(sk, skOut)
		ks.SwitchKeysInPlace(ciphertext.Value[1].Level(), ciphertext.Value[1], swk, ks.PoolQ[1], ks.PoolQ[2])
		ringQ.Add(ciphertext.Value[0], ks.PoolQ[1], ciphertext.Value[0])
		ring.CopyValues(ks.PoolQ[2], ciphertext.Value[1])
		ringQ.MulCoeffsMontgomeryAndAddLvl(ciphertext.Level(), ciphertext.Value[1], skOut.Value[0], ciphertext.Value[0])
		ringQ.InvNTTLvl(ciphertext.Level(), ciphertext.Value[0], ciphertext.Value[0])
		require.GreaterOrEqual(t, 10+params.LogN(), log2OfInnerSum(ciphertext.Level(), ringQ, nil, ciphertext.Value[0], nil))
	})
}
