package dckks

import (
	"encoding/json"
	"flag"
	"fmt"
	"testing"

	"github.com/stretchr/testify/require"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters). Overrides -short and requires -timeout=0.")
var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")
var printPrecisionStats = flag.Bool("print-precision", false, "print precision stats")
var minPrec float64 = 15.0
var parties int = 2

func testString(opname string, parties int, params ckks.Parameters) string {
	return fmt.Sprintf("%sparties=%d/logN=%d/logSlots=%d/logQ=%d/levels=%d/alpha=%d/beta=%d",
		opname,
		parties,
		params.LogN(),
		params.LogSlots(),
		params.LogQP(),
		params.MaxLevel()+1,
		params.PCount(),
		params.Beta())
}

type testContext struct {
	params ckks.Parameters

	ringQ  *ring.Ring
	ringQP *ring.Ring

	prng utils.PRNG

	encoder   ckks.Encoder
	evaluator ckks.Evaluator

	encryptorPk0 ckks.Encryptor
	decryptorSk0 ckks.Decryptor
	decryptorSk1 ckks.Decryptor

	pk0 *rlwe.PublicKey
	pk1 *rlwe.PublicKey

	sk0 *rlwe.SecretKey
	sk1 *rlwe.SecretKey

	sk0Shards []*rlwe.SecretKey
	sk1Shards []*rlwe.SecretKey
}

func TestDCKKS(t *testing.T) {

	var defaultParams = ckks.DefaultParams[:4] // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15
	if testing.Short() {
		defaultParams = ckks.DefaultParams[:2] // the short test runs for ring degree N=2^12, 2^13
	}
	if *flagLongTest {
		defaultParams = ckks.DefaultParams // the long test suite runs for all default parameters
	}
	if *flagParamString != "" {
		var jsonParams ckks.ParametersLiteral
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []ckks.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range defaultParams {

		params, err := ckks.NewParametersFromLiteral(p)
		if err != nil {
			panic(err)
		}

		var testCtx *testContext
		if testCtx, err = genTestParams(params); err != nil {
			panic(err)
		}

		testPublicKeyGen(testCtx, t)
		testRelinKeyGen(testCtx, t)
		testKeyswitching(testCtx, t)
		testPublicKeySwitching(testCtx, t)
		testRotKeyGenConjugate(testCtx, t)
		testRotKeyGenCols(testCtx, t)
		testE2SProtocol(testCtx, t)
		testRefresh(testCtx, t)
		testRefreshAndTransform(testCtx, t)
		testMarshalling(testCtx, t)
	}
}

func genTestParams(defaultParams ckks.Parameters) (testCtx *testContext, err error) {

	testCtx = new(testContext)

	testCtx.params = defaultParams

	testCtx.ringQ = defaultParams.RingQ()

	testCtx.ringQP = defaultParams.RingQP()

	if testCtx.prng, err = utils.NewPRNG(); err != nil {
		return nil, err
	}

	testCtx.encoder = ckks.NewEncoder(testCtx.params)
	testCtx.evaluator = ckks.NewEvaluator(testCtx.params, rlwe.EvaluationKey{})

	kgen := ckks.NewKeyGenerator(testCtx.params)

	// SecretKeys
	testCtx.sk0Shards = make([]*rlwe.SecretKey, parties)
	testCtx.sk1Shards = make([]*rlwe.SecretKey, parties)
	tmp0 := testCtx.ringQP.NewPoly()
	tmp1 := testCtx.ringQP.NewPoly()

	for j := 0; j < parties; j++ {
		testCtx.sk0Shards[j] = kgen.GenSecretKey()
		testCtx.sk1Shards[j] = kgen.GenSecretKey()
		testCtx.ringQP.Add(tmp0, testCtx.sk0Shards[j].Value, tmp0)
		testCtx.ringQP.Add(tmp1, testCtx.sk1Shards[j].Value, tmp1)
	}

	testCtx.sk0 = ckks.NewSecretKey(testCtx.params)
	testCtx.sk1 = ckks.NewSecretKey(testCtx.params)
	testCtx.sk0.Value.Copy(tmp0)
	testCtx.sk1.Value.Copy(tmp1)

	// Publickeys
	testCtx.pk0 = kgen.GenPublicKey(testCtx.sk0)
	testCtx.pk1 = kgen.GenPublicKey(testCtx.sk1)

	testCtx.encryptorPk0 = ckks.NewEncryptor(testCtx.params, testCtx.pk0)
	testCtx.decryptorSk0 = ckks.NewDecryptor(testCtx.params, testCtx.sk0)
	testCtx.decryptorSk1 = ckks.NewDecryptor(testCtx.params, testCtx.sk1)

	return
}

func testPublicKeyGen(testCtx *testContext, t *testing.T) {

	decryptorSk0 := testCtx.decryptorSk0
	sk0Shards := testCtx.sk0Shards
	params := testCtx.params

	t.Run(testString("PublicKeyGen/", parties, params), func(t *testing.T) {
		crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.ringQP)
		crp := crpGenerator.ReadNew()

		type Party struct {
			*CKGProtocol
			s  *rlwe.SecretKey
			s1 *drlwe.CKGShare
		}

		ckgParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			p.CKGProtocol = NewCKGProtocol(params)
			p.s = sk0Shards[i]
			p.s1 = p.AllocateShares()
			ckgParties[i] = p
		}
		P0 := ckgParties[0]

		var _ drlwe.CollectivePublicKeyGenerator = P0.CKGProtocol

		// Each party creates a new CKGProtocol instance
		for i, p := range ckgParties {
			p.GenShare(p.s, crp, p.s1)
			if i > 0 {
				P0.AggregateShares(p.s1, P0.s1, P0.s1)
			}
		}

		pk := ckks.NewPublicKey(params)
		P0.GenPublicKey(P0.s1, crp, pk)

		// Verifies that decrypt((encryptp(collectiveSk, m), collectivePk) = m
		encryptorTest := ckks.NewEncryptor(params, pk)

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorTest, -1, 1, t)

		verifyTestVectors(testCtx, decryptorSk0, coeffs, ciphertext, t)
	})

}

func testRelinKeyGen(testCtx *testContext, t *testing.T) {

	encryptorPk0 := testCtx.encryptorPk0
	decryptorSk0 := testCtx.decryptorSk0
	sk0Shards := testCtx.sk0Shards
	params := testCtx.params

	t.Run(testString("RelinKeyGen/", parties, params), func(t *testing.T) {

		type Party struct {
			*RKGProtocol
			ephSk  *rlwe.SecretKey
			sk     *rlwe.SecretKey
			share1 *drlwe.RKGShare
			share2 *drlwe.RKGShare
		}

		rkgParties := make([]*Party, parties)

		for i := range rkgParties {
			p := new(Party)
			p.RKGProtocol = NewRKGProtocol(params)
			p.sk = sk0Shards[i]
			p.ephSk, p.share1, p.share2 = p.AllocateShares()
			rkgParties[i] = p
		}

		P0 := rkgParties[0]

		// Checks that ckks.RKGProtocol complies to the drlwe.RelinearizationKeyGenerator interface
		var _ drlwe.RelinearizationKeyGenerator = P0.RKGProtocol

		crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.ringQP)
		crp := make([]*ring.Poly, params.Beta())

		for i := 0; i < params.Beta(); i++ {
			crp[i] = crpGenerator.ReadNew()
		}

		// ROUND 1
		for i, p := range rkgParties {
			p.GenShareRoundOne(p.sk, crp, p.ephSk, p.share1)
			if i > 0 {
				P0.AggregateShares(p.share1, P0.share1, P0.share1)
			}
		}

		//ROUND 2
		for i, p := range rkgParties {
			p.GenShareRoundTwo(p.ephSk, p.sk, P0.share1, crp, p.share2)
			if i > 0 {
				P0.AggregateShares(p.share2, P0.share2, P0.share2)
			}
		}

		rlk := ckks.NewRelinearizationKey(params)
		P0.GenRelinearizationKey(P0.share1, P0.share2, rlk)

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, -1, 1, t)

		for i := range coeffs {
			coeffs[i] *= coeffs[i]
		}

		evaluator := testCtx.evaluator.WithKey(rlwe.EvaluationKey{Rlk: rlk, Rtks: nil})
		evaluator.MulRelin(ciphertext, ciphertext, ciphertext)

		evaluator.Rescale(ciphertext, params.Scale(), ciphertext)

		require.Equal(t, ciphertext.Degree(), 1)

		verifyTestVectors(testCtx, decryptorSk0, coeffs, ciphertext, t)

	})

}

func testKeyswitching(testCtx *testContext, t *testing.T) {

	encryptorPk0 := testCtx.encryptorPk0
	decryptorSk1 := testCtx.decryptorSk1
	sk0Shards := testCtx.sk0Shards
	sk1Shards := testCtx.sk1Shards
	params := testCtx.params

	t.Run(testString("Keyswitching/", parties, params), func(t *testing.T) {

		coeffs, _, ciphertextFullLevels := newTestVectors(testCtx, encryptorPk0, -1, 1, t)

		for _, dropped := range []int{0, ciphertextFullLevels.Level()} { // runs the test for full and level zero
			ciphertext := testCtx.evaluator.DropLevelNew(ciphertextFullLevels, dropped)

			t.Run(fmt.Sprintf("atLevel=%d", ciphertext.Level()), func(t *testing.T) {

				type Party struct {
					cks   *CKSProtocol
					s0    *rlwe.SecretKey
					s1    *rlwe.SecretKey
					share *drlwe.CKSShare
				}

				cksParties := make([]*Party, parties)
				for i := 0; i < parties; i++ {
					p := new(Party)
					p.cks = NewCKSProtocol(params, 3.2)
					p.s0 = sk0Shards[i]
					p.s1 = sk1Shards[i]
					p.share = p.cks.AllocateShare(ciphertext.Level())
					cksParties[i] = p
				}
				P0 := cksParties[0]

				// Checks that the protocol complies to the drlwe.KeySwitchingProtocol interface
				var _ drlwe.KeySwitchingProtocol = P0.cks

				// Each party creates its CKSProtocol instance with tmp = si-si'
				for i, p := range cksParties {
					p.cks.GenShare(p.s0, p.s1, ciphertext.Ciphertext, p.share)
					if i > 0 {
						P0.cks.AggregateShares(p.share, P0.share, P0.share)
					}
				}

				ksCiphertext := ckks.NewCiphertext(params, 1, ciphertext.Level(), ciphertext.Scale/2)

				P0.cks.KeySwitchCKKS(P0.share, ciphertext, ksCiphertext)

				verifyTestVectors(testCtx, decryptorSk1, coeffs, ksCiphertext, t)

				P0.cks.KeySwitchCKKS(P0.share, ciphertext, ciphertext)

				verifyTestVectors(testCtx, decryptorSk1, coeffs, ksCiphertext, t)

			})
		}
	})
}

func testPublicKeySwitching(testCtx *testContext, t *testing.T) {

	encryptorPk0 := testCtx.encryptorPk0
	decryptorSk1 := testCtx.decryptorSk1
	sk0Shards := testCtx.sk0Shards
	pk1 := testCtx.pk1
	params := testCtx.params

	t.Run(testString("PublicKeySwitching/", parties, params), func(t *testing.T) {

		coeffs, _, ciphertextFullLevels := newTestVectors(testCtx, encryptorPk0, -1, 1, t)

		for _, dropped := range []int{0, ciphertextFullLevels.Level()} { // runs the test for full and level zero
			ciphertext := testCtx.evaluator.DropLevelNew(ciphertextFullLevels, dropped)

			t.Run(fmt.Sprintf("atLevel=%d", ciphertext.Level()), func(t *testing.T) {

				type Party struct {
					*PCKSProtocol
					s     *rlwe.SecretKey
					share *drlwe.PCKSShare
				}

				pcksParties := make([]*Party, parties)
				for i := 0; i < parties; i++ {
					p := new(Party)
					p.PCKSProtocol = NewPCKSProtocol(params, 3.2)
					p.s = sk0Shards[i]
					p.share = p.AllocateShare(ciphertext.Level())
					pcksParties[i] = p
				}
				P0 := pcksParties[0]

				// Checks that the protocol complies to the drlwe.KeySwitchingProtocol interface
				var _ drlwe.PublicKeySwitchingProtocol = P0.PCKSProtocol

				ciphertextSwitched := ckks.NewCiphertext(params, 1, ciphertext.Level(), ciphertext.Scale)

				for i, p := range pcksParties {
					p.GenShare(p.s, pk1, ciphertext.Ciphertext, p.share)
					if i > 0 {
						P0.AggregateShares(p.share, P0.share, P0.share)
					}
				}

				P0.KeySwitchCKKS(P0.share, ciphertext, ciphertextSwitched)

				verifyTestVectors(testCtx, decryptorSk1, coeffs, ciphertextSwitched, t)
			})
		}

	})
}

func testRotKeyGenConjugate(testCtx *testContext, t *testing.T) {

	ringQP := testCtx.ringQP
	encryptorPk0 := testCtx.encryptorPk0
	decryptorSk0 := testCtx.decryptorSk0
	sk0Shards := testCtx.sk0Shards
	params := testCtx.params

	t.Run(testString("RotKeyGenConjugate/", parties, params), func(t *testing.T) {

		type Party struct {
			*RTGProtocol
			s     *rlwe.SecretKey
			share *drlwe.RTGShare
		}

		pcksParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			p.RTGProtocol = NewRotKGProtocol(params)
			p.s = sk0Shards[i]
			p.share = p.AllocateShares()
			pcksParties[i] = p
		}
		P0 := pcksParties[0]

		// checks that ckks.RTGProtocol complies to the drlwe.RotationKeyGenerator interface
		var _ drlwe.RotationKeyGenerator = P0.RTGProtocol

		crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.ringQP)
		crp := make([]*ring.Poly, params.Beta())

		for i := 0; i < params.Beta(); i++ {
			crp[i] = crpGenerator.ReadNew()
		}

		galEl := params.GaloisElementForRowRotation()
		rotKeySet := ckks.NewRotationKeySet(params, []uint64{galEl})

		for i, p := range pcksParties {
			p.GenShare(p.s, galEl, crp, p.share)
			if i > 0 {
				P0.Aggregate(p.share, P0.share, P0.share)
			}
		}

		P0.GenRotationKey(P0.share, crp, rotKeySet.Keys[galEl])

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, -1, 1, t)

		evaluator := testCtx.evaluator.WithKey(rlwe.EvaluationKey{Rlk: nil, Rtks: rotKeySet})
		evaluator.Conjugate(ciphertext, ciphertext)

		coeffsWant := make([]complex128, ringQP.N>>1)

		for i := 0; i < ringQP.N>>1; i++ {
			coeffsWant[i] = complex(real(coeffs[i]), -imag(coeffs[i]))
		}

		verifyTestVectors(testCtx, decryptorSk0, coeffsWant, ciphertext, t)

	})
}

func testRotKeyGenCols(testCtx *testContext, t *testing.T) {

	ringQP := testCtx.ringQP
	encryptorPk0 := testCtx.encryptorPk0
	decryptorSk0 := testCtx.decryptorSk0
	sk0Shards := testCtx.sk0Shards
	params := testCtx.params

	t.Run(testString("RotKeyGenCols/", parties, params), func(t *testing.T) {

		type Party struct {
			*RTGProtocol
			s     *rlwe.SecretKey
			share *drlwe.RTGShare
		}

		pcksParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			p.RTGProtocol = NewRotKGProtocol(params)
			p.s = sk0Shards[i]
			p.share = p.AllocateShares()
			pcksParties[i] = p
		}

		P0 := pcksParties[0]

		crpGenerator := ring.NewUniformSampler(testCtx.prng, ringQP)
		crp := make([]*ring.Poly, params.Beta())

		for i := 0; i < params.Beta(); i++ {
			crp[i] = crpGenerator.ReadNew()
		}

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, -1, 1, t)

		receiver := ckks.NewCiphertext(params, ciphertext.Degree(), ciphertext.Level(), ciphertext.Scale)

		galEls := params.GaloisElementsForRowInnerSum()
		rotKeySet := ckks.NewRotationKeySet(params, galEls)

		for _, galEl := range galEls {
			for i, p := range pcksParties {
				p.GenShare(p.s, galEl, crp, p.share)
				if i > 0 {
					P0.Aggregate(p.share, P0.share, P0.share)
				}
			}
			P0.GenRotationKey(P0.share, crp, rotKeySet.Keys[galEl])
		}

		evaluator := testCtx.evaluator.WithKey(rlwe.EvaluationKey{Rlk: nil, Rtks: rotKeySet})

		for k := 1; k < ringQP.N>>1; k <<= 1 {
			evaluator.Rotate(ciphertext, int(k), receiver)

			coeffsWant := utils.RotateComplex128Slice(coeffs, int(k))

			verifyTestVectors(testCtx, decryptorSk0, coeffsWant, receiver, t)
		}
	})
}

func testE2SProtocol(testCtx *testContext, t *testing.T) {

	params := testCtx.params

	t.Run(testString("E2SProtocol/", parties, params), func(t *testing.T) {

		var minLevel, logBound int
		var ok bool
		if minLevel, logBound, ok = GetMinimumLevelForBootstrapping(128, params.Scale(), parties, params.Q()); ok != true {
			t.Skip("Not enough levels to ensure correcness and 128 security")
		}

		type Party struct {
			e2s            *E2SProtocol
			s2e            *S2EProtocol
			sk             *rlwe.SecretKey
			publicShareE2S *drlwe.CKSShare
			publicShareS2E *drlwe.CKSShare
			secretShare    *rlwe.AdditiveShareBigint
		}

		coeffs, _, ciphertext := newTestVectors(testCtx, testCtx.encryptorPk0, -1, 1, t)

		testCtx.evaluator.DropLevel(ciphertext, ciphertext.Level()-minLevel)

		params := testCtx.params
		P := make([]Party, parties)
		for i := range P {
			P[i].e2s = NewE2SProtocol(params, 3.2)
			P[i].s2e = NewS2EProtocol(params, 3.2)
			P[i].sk = testCtx.sk0Shards[i]
			P[i].publicShareE2S = P[i].e2s.AllocateShare(ciphertext.Level())
			P[i].publicShareS2E = P[i].s2e.AllocateShare(params.Parameters.MaxLevel())
			P[i].secretShare = rlwe.NewAdditiveShareBigint(params.Parameters)
		}

		if testCtx.params.MaxLevel() < minLevel {
			t.Skip("Not enough levels to ensure correcness and 128 security")
		}

		for i, p := range P {
			// Enc(-M_i)
			p.e2s.GenShare(p.sk, logBound, params.LogSlots(), ciphertext, p.secretShare, p.publicShareE2S)
			if i > 0 {
				// Enc(sum(-M_i))
				p.e2s.AggregateShares(P[0].publicShareE2S, p.publicShareE2S, P[0].publicShareE2S)
			}
		}

		// sum(-M_i) + x
		P[0].e2s.GetShare(P[0].secretShare, P[0].publicShareE2S, ciphertext, P[0].secretShare)

		// sum(-M_i) + x + sum(M_i) = x
		rec := rlwe.NewAdditiveShareBigint(params.Parameters)
		for _, p := range P {
			a := rec.Value
			b := p.secretShare.Value

			for i := range a {
				a[i].Add(a[i], b[i])
			}
		}

		pt := ckks.NewPlaintext(params, ciphertext.Level(), ciphertext.Scale)
		pt.Value.IsNTT = false
		testCtx.ringQ.SetCoefficientsBigintLvl(pt.Level(), rec.Value, pt.Value)

		verifyTestVectors(testCtx, nil, coeffs, pt, t)

		crs := ring.NewUniformSampler(testCtx.prng, testCtx.ringQ)
		c1 := crs.ReadLvlNew(params.Parameters.MaxLevel())

		for i, p := range P {

			p.s2e.GenShare(p.sk, c1, p.secretShare, p.publicShareS2E)

			if i > 0 {
				p.s2e.AggregateShares(P[0].publicShareS2E, p.publicShareS2E, P[0].publicShareS2E)
			}
		}

		ctRec := ckks.NewCiphertext(params, 1, params.Parameters.MaxLevel(), ciphertext.Scale)
		P[0].s2e.GetEncryption(P[0].publicShareS2E, c1, ctRec)

		verifyTestVectors(testCtx, testCtx.decryptorSk0, coeffs, ctRec, t)

	})
}

func testRefresh(testCtx *testContext, t *testing.T) {

	encryptorPk0 := testCtx.encryptorPk0
	sk0Shards := testCtx.sk0Shards
	decryptorSk0 := testCtx.decryptorSk0
	params := testCtx.params

	t.Run(testString("Refresh/", parties, params), func(t *testing.T) {

		var minLevel, logBound int
		var ok bool
		if minLevel, logBound, ok = GetMinimumLevelForBootstrapping(128, params.Scale(), parties, params.Q()); ok != true {
			t.Skip("Not enough levels to ensure correcness and 128 security")
		}

		type Party struct {
			*RefreshProtocol
			s     *rlwe.SecretKey
			share *RefreshShare
		}

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, -1, 1, t)

		// Brings ciphertext to level 2
		testCtx.evaluator.DropLevel(ciphertext, ciphertext.Level()-minLevel)

		levelMin := ciphertext.Level()
		levelMax := params.MaxLevel()

		RefreshParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			p.RefreshProtocol = NewRefreshProtocol(params, logBound, 3.2)
			p.s = sk0Shards[i]
			p.share = p.AllocateShare(levelMin, levelMax)
			RefreshParties[i] = p
		}

		P0 := RefreshParties[0]

		crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.ringQ)
		crp := crpGenerator.ReadLvlNew(levelMax)

		for i, p := range RefreshParties {
			p.GenShares(p.s, logBound, params.LogSlots(), ciphertext, crp, p.share)
			if i > 0 {
				P0.Aggregate(p.share, P0.share, P0.share)
			}
		}

		P0.Finalize(ciphertext, params.LogSlots(), crp, P0.share, ciphertext)

		verifyTestVectors(testCtx, decryptorSk0, coeffs, ciphertext, t)
	})
}

func testRefreshAndTransform(testCtx *testContext, t *testing.T) {

	encryptorPk0 := testCtx.encryptorPk0
	sk0Shards := testCtx.sk0Shards
	params := testCtx.params
	decryptorSk0 := testCtx.decryptorSk0

	t.Run(testString("RefreshAndTransform/", parties, params), func(t *testing.T) {

		var minLevel, logBound int
		var ok bool
		if minLevel, logBound, ok = GetMinimumLevelForBootstrapping(128, params.Scale(), parties, params.Q()); ok != true {
			t.Skip("Not enough levels to ensure correcness and 128 security")
		}

		type Party struct {
			*MaskedTransformProtocol
			s     *rlwe.SecretKey
			share *MaskedTransformShare
		}

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, -1, 1, t)

		// Drops the ciphertext to the minimum level that ensures correctness and 128-bit security
		testCtx.evaluator.DropLevel(ciphertext, ciphertext.Level()-minLevel)

		levelMin := ciphertext.Level()
		levelMax := params.MaxLevel()

		RefreshParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			p.MaskedTransformProtocol = NewMaskedTransformProtocol(params, logBound, 3.2)
			p.s = sk0Shards[i]
			p.share = p.AllocateShare(levelMin, levelMax)
			RefreshParties[i] = p
		}

		P0 := RefreshParties[0]

		crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.ringQ)
		crp := crpGenerator.ReadLvlNew(levelMax)

		permute := func(ptIn, ptOut []*ring.Complex) {
			for i := range ptIn {
				ptOut[i][0].Mul(ptIn[i][0], ring.NewFloat(0.9238795325112867, logBound))
				ptOut[i][1].Mul(ptIn[i][1], ring.NewFloat(0.7071067811865476, logBound))
			}
		}

		for i, p := range RefreshParties {
			p.GenShares(p.s, logBound, params.LogSlots(), ciphertext, crp, permute, p.share)
			if i > 0 {
				P0.Aggregate(p.share, P0.share, P0.share)
			}
		}

		P0.Transform(ciphertext, testCtx.params.LogSlots(), permute, crp, P0.share, ciphertext)

		for i := range coeffs {
			coeffs[i] = complex(real(coeffs[i])*0.9238795325112867, imag(coeffs[i])*0.7071067811865476)
		}

		verifyTestVectors(testCtx, decryptorSk0, coeffs, ciphertext, t)
	})
}

func testMarshalling(testCtx *testContext, t *testing.T) {
	crsGen := ring.NewUniformSampler(testCtx.prng, testCtx.ringQP)
	params := testCtx.params

	t.Run(testString("Marshalling/Refresh/", parties, params), func(t *testing.T) {

		var minLevel, logBound int
		var ok bool
		if minLevel, logBound, ok = GetMinimumLevelForBootstrapping(128, params.Scale(), parties, params.Q()); ok != true {
			t.Skip("Not enough levels to ensure correcness and 128 security")
		}

		ciphertext := ckks.NewCiphertextRandom(testCtx.prng, testCtx.params, 1, minLevel, testCtx.params.Scale())

		crsLevel := crsGen.ReadLvlNew(minLevel)

		//testing refresh shares
		refreshproto := NewRefreshProtocol(testCtx.params, logBound, 3.2)
		refreshshare := refreshproto.AllocateShare(ciphertext.Level(), ciphertext.Level())
		refreshproto.GenShares(testCtx.sk0, logBound, params.LogSlots(), ciphertext, crsLevel, refreshshare)

		data, err := refreshshare.MarshalBinary()

		if err != nil {
			t.Fatal("Could not marshal RefreshShare", err)
		}

		resRefreshShare := new(MaskedTransformShare)
		err = resRefreshShare.UnmarshalBinary(data)

		if err != nil {
			t.Fatal("Could not unmarshal RefreshShare", err)
		}

		for i, r := range refreshshare.e2sShare.Value.Coeffs {
			if !utils.EqualSliceUint64(resRefreshShare.e2sShare.Value.Coeffs[i], r) {
				t.Fatal("Resulting of marshalling not the same as original : RefreshShare")
			}

		}
		for i, r := range refreshshare.s2eShare.Value.Coeffs {
			if !utils.EqualSliceUint64(resRefreshShare.s2eShare.Value.Coeffs[i], r) {
				t.Fatal("Resulting of marshalling not the same as original : RefreshShare")
			}

		}
	})
}

func newTestVectors(testContext *testContext, encryptor ckks.Encryptor, a, b complex128, t *testing.T) (values []complex128, plaintext *ckks.Plaintext, ciphertext *ckks.Ciphertext) {

	params := testContext.params

	logSlots := params.LogSlots()

	values = make([]complex128, 1<<logSlots)

	for i := 0; i < 1<<logSlots; i++ {
		values[i] = complex(utils.RandFloat64(real(a), real(b)), utils.RandFloat64(imag(a), imag(b)))
	}

	plaintext = testContext.encoder.EncodeNTTAtLvlNew(params.MaxLevel(), values, logSlots)

	if encryptor != nil {
		ciphertext = encryptor.EncryptNew(plaintext)
	}

	return values, plaintext, ciphertext
}

func verifyTestVectors(testCtx *testContext, decryptor ckks.Decryptor, valuesWant []complex128, element interface{}, t *testing.T) {

	precStats := ckks.GetPrecisionStats(testCtx.params, testCtx.encoder, decryptor, valuesWant, element, testCtx.params.LogSlots(), 0)

	if *printPrecisionStats {
		t.Log(precStats.String())
	}

	require.GreaterOrEqual(t, real(precStats.MeanPrecision), minPrec)
	require.GreaterOrEqual(t, imag(precStats.MeanPrecision), minPrec)
}
