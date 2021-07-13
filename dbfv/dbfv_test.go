package dbfv

import (
	"encoding/json"
	"flag"
	"fmt"
	"math/big"
	"testing"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

var flagLongTest = flag.Bool("long", false, "run the long test suite (all parameters). Overrides -short and requires -timeout=0.")
var flagParamString = flag.String("params", "", "specify the test cryptographic parameters as a JSON string. Overrides -short and -long.")
var parties int = 3

func testString(opname string, parties int, params bfv.Parameters) string {
	return fmt.Sprintf("%sparties=%d/LogN=%d/logQ=%d", opname, parties, params.LogN(), params.LogQP())
}

type testContext struct {
	params bfv.Parameters

	// Polynomial degree
	n int

	// floor(Q/T) mod each Qi in Montgomery form
	deltaMont []uint64

	// Polynomial contexts
	ringT  *ring.Ring
	ringQ  *ring.Ring
	ringQP *ring.Ring

	prng utils.PRNG

	encoder bfv.Encoder

	sk0Shards []*rlwe.SecretKey
	sk0       *rlwe.SecretKey

	sk1       *rlwe.SecretKey
	sk1Shards []*rlwe.SecretKey

	pk0 *rlwe.PublicKey
	pk1 *rlwe.PublicKey

	encryptorPk0 bfv.Encryptor
	decryptorSk0 bfv.Decryptor
	decryptorSk1 bfv.Decryptor
	evaluator    bfv.Evaluator
}

func Test_DBFV(t *testing.T) {

	defaultParams := bfv.DefaultParams // the default test runs for ring degree N=2^12, 2^13, 2^14, 2^15
	if testing.Short() {
		defaultParams = bfv.DefaultParams[:2] // the short test runs for ring degree N=2^12, 2^13
	}
	if *flagLongTest {
		defaultParams = bfv.DefaultParams // the long test suite runs for all default parameters
	}
	if *flagParamString != "" {
		var jsonParams bfv.ParametersLiteral
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []bfv.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range defaultParams {

		params, err := bfv.NewParametersFromLiteral(p)
		var testCtx *testContext
		if testCtx, err = gentestContext(params); err != nil {
			panic(err)
		}

		testPublicKeyGen(testCtx, t)
		testRelinKeyGen(testCtx, t)
		testKeyswitching(testCtx, t)
		testPublicKeySwitching(testCtx, t)
		testRotKeyGenRotRows(testCtx, t)
		testRotKeyGenRotCols(testCtx, t)
		testEncToShares(testCtx, t)
		testRefresh(testCtx, t)
		testRefreshAndPermutation(testCtx, t)
		testMarshalling(testCtx, t)
	}
}

func gentestContext(params bfv.Parameters) (testCtx *testContext, err error) {

	testCtx = new(testContext)

	testCtx.params = params

	testCtx.n = params.N()

	testCtx.ringT = params.RingT()
	testCtx.ringQ = params.RingQ()
	testCtx.ringQP = params.RingQP()

	testCtx.deltaMont = bfv.GenLiftParams(testCtx.ringQ, params.T())

	if testCtx.prng, err = utils.NewPRNG(); err != nil {
		return nil, err
	}

	testCtx.encoder = bfv.NewEncoder(testCtx.params)
	testCtx.evaluator = bfv.NewEvaluator(testCtx.params, rlwe.EvaluationKey{})

	kgen := bfv.NewKeyGenerator(testCtx.params)

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

	testCtx.sk0 = bfv.NewSecretKey(testCtx.params)
	testCtx.sk1 = bfv.NewSecretKey(testCtx.params)

	testCtx.sk0.Value.Copy(tmp0)
	testCtx.sk1.Value.Copy(tmp1)

	// Publickeys
	testCtx.pk0 = kgen.GenPublicKey(testCtx.sk0)
	testCtx.pk1 = kgen.GenPublicKey(testCtx.sk1)

	testCtx.encryptorPk0 = bfv.NewEncryptor(testCtx.params, testCtx.pk0)
	testCtx.decryptorSk0 = bfv.NewDecryptor(testCtx.params, testCtx.sk0)
	testCtx.decryptorSk1 = bfv.NewDecryptor(testCtx.params, testCtx.sk1)

	return
}

func testPublicKeyGen(testCtx *testContext, t *testing.T) {

	sk0Shards := testCtx.sk0Shards
	decryptorSk0 := testCtx.decryptorSk0

	t.Run(testString("PublicKeyGen/", parties, testCtx.params), func(t *testing.T) {

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
			p.CKGProtocol = NewCKGProtocol(testCtx.params)
			p.s = sk0Shards[i]
			p.s1 = p.AllocateShares()
			ckgParties[i] = p
		}
		P0 := ckgParties[0]

		// Checks that bfv.CKGProtocol complies to the drlwe.CollectivePublicKeyGenerator interface
		var _ drlwe.CollectivePublicKeyGenerator = P0.CKGProtocol

		// Each party creates a new CKGProtocol instance
		for i, p := range ckgParties {
			p.GenShare(p.s, crp, p.s1)
			if i > 0 {
				P0.AggregateShares(p.s1, P0.s1, P0.s1)
			}
		}

		pk := bfv.NewPublicKey(testCtx.params)
		P0.GenPublicKey(P0.s1, crp, pk)

		// Verifies that decrypt((encryptp(collectiveSk, m), collectivePk) = m
		encryptorTest := bfv.NewEncryptor(testCtx.params, pk)

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorTest, t)

		verifyTestVectors(testCtx, decryptorSk0, coeffs, ciphertext, t)
	})
}

func testRelinKeyGen(testCtx *testContext, t *testing.T) {

	sk0Shards := testCtx.sk0Shards
	encryptorPk0 := testCtx.encryptorPk0
	decryptorSk0 := testCtx.decryptorSk0

	t.Run(testString("RelinKeyGen/", parties, testCtx.params), func(t *testing.T) {

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
			p.RKGProtocol = NewRKGProtocol(testCtx.params)
			p.sk = sk0Shards[i]
			p.ephSk, p.share1, p.share2 = p.AllocateShares()
			rkgParties[i] = p
		}

		P0 := rkgParties[0]

		// checks that bfv.RKGProtocol complies to the drlwe.RelinearizationKeyGenerator interface
		var _ drlwe.RelinearizationKeyGenerator = P0.RKGProtocol

		crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.ringQP)
		crp := make([]*ring.Poly, testCtx.params.Beta())

		for i := 0; i < testCtx.params.Beta(); i++ {
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

		evk := bfv.NewRelinearizationKey(testCtx.params, 1)
		P0.GenRelinearizationKey(P0.share1, P0.share2, evk)

		evaluator := testCtx.evaluator.WithKey(rlwe.EvaluationKey{Rlk: evk, Rtks: nil})

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)
		for i := range coeffs {
			coeffs[i] *= coeffs[i]
			coeffs[i] %= testCtx.ringT.Modulus[0]
		}

		ciphertextMul := bfv.NewCiphertext(testCtx.params, ciphertext.Degree()*2)
		evaluator.Mul(ciphertext, ciphertext, ciphertextMul)

		res := bfv.NewCiphertext(testCtx.params, 1)
		evaluator.Relinearize(ciphertextMul, res)

		verifyTestVectors(testCtx, decryptorSk0, coeffs, res, t)
	})

}

func testKeyswitching(testCtx *testContext, t *testing.T) {

	sk0Shards := testCtx.sk0Shards
	sk1Shards := testCtx.sk1Shards
	encryptorPk0 := testCtx.encryptorPk0
	decryptorSk1 := testCtx.decryptorSk1

	t.Run(testString("Keyswitching/", parties, testCtx.params), func(t *testing.T) {

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

		type Party struct {
			cks   *CKSProtocol
			s0    *rlwe.SecretKey
			s1    *rlwe.SecretKey
			share *drlwe.CKSShare
		}

		cksParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			p.cks = NewCKSProtocol(testCtx.params, 6.36)
			p.s0 = sk0Shards[i]
			p.s1 = sk1Shards[i]
			p.share = p.cks.AllocateShare(ciphertext.Level())
			cksParties[i] = p
		}
		P0 := cksParties[0]

		// checks that the protocol complies to the drlwe.PublicKeySwitchingProtocol interface
		var _ drlwe.KeySwitchingProtocol = P0.cks

		// Each party creates its CKSProtocol instance with tmp = si-si'
		for i, p := range cksParties {
			p.cks.GenShare(p.s0, p.s1, ciphertext.Ciphertext, p.share)
			if i > 0 {
				P0.cks.AggregateShares(p.share, P0.share, P0.share)
			}
		}

		ksCiphertext := bfv.NewCiphertext(testCtx.params, 1)
		P0.cks.KeySwitch(P0.share, ciphertext.Ciphertext, ksCiphertext.Ciphertext)

		verifyTestVectors(testCtx, decryptorSk1, coeffs, ksCiphertext, t)

		P0.cks.KeySwitch(P0.share, ciphertext.Ciphertext, ciphertext.Ciphertext)

		verifyTestVectors(testCtx, decryptorSk1, coeffs, ciphertext, t)

	})
}

func testPublicKeySwitching(testCtx *testContext, t *testing.T) {

	sk0Shards := testCtx.sk0Shards
	pk1 := testCtx.pk1
	encryptorPk0 := testCtx.encryptorPk0
	decryptorSk1 := testCtx.decryptorSk1

	t.Run(testString("PublicKeySwitching/", parties, testCtx.params), func(t *testing.T) {

		type Party struct {
			*PCKSProtocol
			s     *rlwe.SecretKey
			share *drlwe.PCKSShare
		}

		pcksParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			p.PCKSProtocol = NewPCKSProtocol(testCtx.params, 6.36)
			p.s = sk0Shards[i]
			p.share = p.AllocateShare(testCtx.params.MaxLevel())
			pcksParties[i] = p
		}
		P0 := pcksParties[0]

		// checks that the protocol complies to the drlwe.PublicKeySwitchingProtocol interface
		var _ drlwe.PublicKeySwitchingProtocol = P0.PCKSProtocol

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

		ciphertextSwitched := bfv.NewCiphertext(testCtx.params, 1)

		for i, p := range pcksParties {
			p.GenShare(p.s, pk1, ciphertext.Ciphertext, p.share)
			if i > 0 {
				P0.AggregateShares(p.share, P0.share, P0.share)
			}
		}

		P0.KeySwitch(P0.share, ciphertext.Ciphertext, ciphertextSwitched.Ciphertext)

		verifyTestVectors(testCtx, decryptorSk1, coeffs, ciphertextSwitched, t)
	})
}

func testRotKeyGenRotRows(testCtx *testContext, t *testing.T) {

	encryptorPk0 := testCtx.encryptorPk0
	decryptorSk0 := testCtx.decryptorSk0
	sk0Shards := testCtx.sk0Shards

	t.Run(testString("RotKeyGenRotRows/", parties, testCtx.params), func(t *testing.T) {

		type Party struct {
			*RTGProtocol
			s     *rlwe.SecretKey
			share *drlwe.RTGShare
		}

		pcksParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			p.RTGProtocol = NewRotKGProtocol(testCtx.params)
			p.s = sk0Shards[i]
			p.share = p.AllocateShares()
			pcksParties[i] = p
		}
		P0 := pcksParties[0]

		// Checks that bfv.RTGProtocol complies to the drlwe.RotationKeyGenerator interface
		var _ drlwe.RotationKeyGenerator = P0.RTGProtocol

		crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.ringQP)
		crp := make([]*ring.Poly, testCtx.params.Beta())

		for i := 0; i < testCtx.params.Beta(); i++ {
			crp[i] = crpGenerator.ReadNew()
		}

		galEl := testCtx.params.GaloisElementForRowRotation()
		rotKeySet := bfv.NewRotationKeySet(testCtx.params, []uint64{galEl})

		for i, p := range pcksParties {
			p.GenShare(p.s, galEl, crp, p.share)
			if i > 0 {
				P0.Aggregate(p.share, P0.share, P0.share)
			}
		}

		P0.GenRotationKey(P0.share, crp, rotKeySet.Keys[galEl])

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

		evaluator := testCtx.evaluator.WithKey(rlwe.EvaluationKey{Rlk: nil, Rtks: rotKeySet})
		result := evaluator.RotateRowsNew(ciphertext)
		coeffsWant := append(coeffs[testCtx.params.N()>>1:], coeffs[:testCtx.params.N()>>1]...)

		verifyTestVectors(testCtx, decryptorSk0, coeffsWant, result, t)

	})
}

func testRotKeyGenRotCols(testCtx *testContext, t *testing.T) {

	encryptorPk0 := testCtx.encryptorPk0
	decryptorSk0 := testCtx.decryptorSk0
	sk0Shards := testCtx.sk0Shards

	t.Run(testString("RotKeyGenRotCols/", parties, testCtx.params), func(t *testing.T) {

		type Party struct {
			*RTGProtocol
			s     *rlwe.SecretKey
			share *drlwe.RTGShare
		}

		pcksParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			p.RTGProtocol = NewRotKGProtocol(testCtx.params)
			p.s = sk0Shards[i]
			p.share = p.AllocateShares()
			pcksParties[i] = p
		}

		P0 := pcksParties[0]

		// Checks that bfv.RTGProtocol complies to the drlwe.RotationKeyGenerator interface
		var _ drlwe.RotationKeyGenerator = P0.RTGProtocol

		crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.ringQP)
		crp := make([]*ring.Poly, testCtx.params.Beta())

		for i := 0; i < testCtx.params.Beta(); i++ {
			crp[i] = crpGenerator.ReadNew()
		}

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

		galEls := testCtx.params.GaloisElementsForRowInnerSum()
		rotKeySet := bfv.NewRotationKeySet(testCtx.params, galEls)

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
		for k := 1; k < testCtx.params.N()>>1; k <<= 1 {
			result := evaluator.RotateColumnsNew(ciphertext, int(k))
			coeffsWant := utils.RotateUint64Slots(coeffs, int(k))
			verifyTestVectors(testCtx, decryptorSk0, coeffsWant, result, t)
		}
	})
}

func testEncToShares(testCtx *testContext, t *testing.T) {

	coeffs, _, ciphertext := newTestVectors(testCtx, testCtx.encryptorPk0, t)

	type Party struct {
		e2s         *E2SProtocol
		s2e         *S2EProtocol
		sk          *rlwe.SecretKey
		publicShare *drlwe.CKSShare
		secretShare *rlwe.AdditiveShare
	}

	params := testCtx.params
	P := make([]Party, parties)

	for i := range P {
		P[i].e2s = NewE2SProtocol(params, 3.2)
		P[i].s2e = NewS2EProtocol(params, 3.2)
		P[i].sk = testCtx.sk0Shards[i]
		P[i].publicShare = P[i].e2s.AllocateShare(ciphertext.Level())
		P[i].secretShare = rlwe.NewAdditiveShare(params.Parameters)
	}

	// The E2S protocol is run in all tests, as a setup to the S2E test.
	for i, p := range P {
		p.e2s.GenShare(p.sk, ciphertext, p.secretShare, p.publicShare)
		if i > 0 {
			p.e2s.AggregateShares(P[0].publicShare, p.publicShare, P[0].publicShare)
		}
	}

	P[0].e2s.GetShare(P[0].secretShare, P[0].publicShare, ciphertext, P[0].secretShare)

	t.Run(testString("E2SProtocol/", parties, testCtx.params), func(t *testing.T) {

		rec := rlwe.NewAdditiveShare(params.Parameters)
		for _, p := range P {
			//fmt.Println("P[", i, "] share:", p.secretShare.Value.Coeffs[0][:see])
			testCtx.ringT.Add(&rec.Value, &p.secretShare.Value, &rec.Value)
		}

		ptRt := bfv.NewPlaintextRingT(testCtx.params)
		ptRt.Value.Copy(&rec.Value)

		assert.True(t, utils.EqualSliceUint64(coeffs, testCtx.encoder.DecodeUintNew(ptRt)))

	})

	crs := ring.NewUniformSampler(testCtx.prng, testCtx.ringQ)
	c1 := crs.ReadNew()

	t.Run(testString("S2EProtocol/", parties, testCtx.params), func(t *testing.T) {
		for i, p := range P {
			p.s2e.GenShare(p.sk, c1, p.secretShare, p.publicShare)
			if i > 0 {
				p.s2e.AggregateShares(P[0].publicShare, p.publicShare, P[0].publicShare)
			}
		}

		ctRec := bfv.NewCiphertext(testCtx.params, 1)
		P[0].s2e.GetEncryption(P[0].publicShare, c1, ctRec)

		verifyTestVectors(testCtx, testCtx.decryptorSk0, coeffs, ctRec, t)
	})
}

func testRefresh(testCtx *testContext, t *testing.T) {

	encryptorPk0 := testCtx.encryptorPk0
	sk0Shards := testCtx.sk0Shards
	encoder := testCtx.encoder
	decryptorSk0 := testCtx.decryptorSk0

	kgen := bfv.NewKeyGenerator(testCtx.params)

	rlk := kgen.GenRelinearizationKey(testCtx.sk0, 2)

	t.Run(testString("Refresh/", parties, testCtx.params), func(t *testing.T) {

		type Party struct {
			*RefreshProtocol
			s       *rlwe.SecretKey
			share   *RefreshShare
			ptShare *bfv.Plaintext
		}

		RefreshParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			p.RefreshProtocol = NewRefreshProtocol(testCtx.params, 3.2)
			p.s = sk0Shards[i]
			p.share = p.AllocateShare()
			p.ptShare = bfv.NewPlaintext(testCtx.params)
			RefreshParties[i] = p
		}

		P0 := RefreshParties[0]

		crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.ringQ)
		crp := crpGenerator.ReadNew()

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

		maxDepth := 0

		ciphertextTmp := ciphertext.CopyNew()
		coeffsTmp := make([]uint64, len(coeffs))

		copy(coeffsTmp, coeffs)

		evaluator := testCtx.evaluator.WithKey(rlwe.EvaluationKey{Rlk: rlk, Rtks: nil})
		// Finds the maximum multiplicative depth
		for {

			evaluator.Relinearize(testCtx.evaluator.MulNew(ciphertextTmp, ciphertextTmp), ciphertextTmp)

			for j := range coeffsTmp {
				coeffsTmp[j] = ring.BRed(coeffsTmp[j], coeffsTmp[j], testCtx.ringT.Modulus[0], testCtx.ringT.BredParams[0])
			}

			if utils.EqualSliceUint64(coeffsTmp, encoder.DecodeUintNew(decryptorSk0.DecryptNew(ciphertextTmp))) {
				maxDepth++
			} else {
				break
			}
		}

		// Simulated added error of size Q/(T^2) and add it to the fresh ciphertext
		coeffsBigint := make([]*big.Int, testCtx.params.N())
		testCtx.ringQ.PolyToBigint(ciphertext.Value[0], coeffsBigint)

		errorRange := new(big.Int).Set(testCtx.ringQ.ModulusBigint)
		errorRange.Quo(errorRange, testCtx.ringT.ModulusBigint)
		errorRange.Quo(errorRange, testCtx.ringT.ModulusBigint)

		for i := 0; i < testCtx.params.N(); i++ {
			coeffsBigint[i].Add(coeffsBigint[i], ring.RandInt(errorRange))
		}

		testCtx.ringQ.SetCoefficientsBigint(coeffsBigint, ciphertext.Value[0])

		for i, p := range RefreshParties {
			p.GenShares(p.s, ciphertext, crp, p.share)
			if i > 0 {
				P0.Aggregate(p.share, P0.share, P0.share)
			}
		}

		ctRes := bfv.NewCiphertext(testCtx.params, 1)
		P0.Finalize(ciphertext, crp, P0.share, ctRes)

		// Square the refreshed ciphertext up to the maximum depth-1
		for i := 0; i < maxDepth-1; i++ {

			evaluator.Relinearize(testCtx.evaluator.MulNew(ctRes, ctRes), ctRes)

			for j := range coeffs {
				coeffs[j] = ring.BRed(coeffs[j], coeffs[j], testCtx.ringT.Modulus[0], testCtx.ringT.BredParams[0])
			}
		}

		//Decrypts and compare
		require.True(t, utils.EqualSliceUint64(coeffs, encoder.DecodeUintNew(decryptorSk0.DecryptNew(ctRes))))
	})
}

func testRefreshAndPermutation(testCtx *testContext, t *testing.T) {

	encryptorPk0 := testCtx.encryptorPk0
	sk0Shards := testCtx.sk0Shards
	encoder := testCtx.encoder
	decryptorSk0 := testCtx.decryptorSk0

	t.Run(testString("RefreshAndPermutation/", parties, testCtx.params), func(t *testing.T) {

		type Party struct {
			*MaskedTransformProtocol
			s       *rlwe.SecretKey
			share   *MaskedTransformShare
			ptShare *bfv.Plaintext
		}

		RefreshParties := make([]*Party, parties)
		for i := 0; i < parties; i++ {
			p := new(Party)
			p.MaskedTransformProtocol = NewMaskedTransformProtocol(testCtx.params, 3.2)
			p.s = sk0Shards[i]
			p.share = p.AllocateShare()
			p.ptShare = bfv.NewPlaintext(testCtx.params)
			RefreshParties[i] = p
		}

		P0 := RefreshParties[0]

		crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.ringQ)
		crp := crpGenerator.ReadNew()

		coeffs, _, ciphertext := newTestVectors(testCtx, encryptorPk0, t)

		permutation := make([]uint64, len(coeffs))
		N := uint64(testCtx.params.N())
		for i := range permutation {
			permutation[i] = ring.RandUniform(testCtx.prng, N, N-1)
		}

		permute := func(ptIn, ptOut bfv.PlaintextRingT) {
			coeffs := testCtx.encoder.DecodeUintNew(&ptIn)
			coeffsPerm := make([]uint64, len(coeffs))
			for i := range coeffs {
				coeffsPerm[i] = coeffs[permutation[i]]
			}
			testCtx.encoder.EncodeUintRingT(coeffsPerm, &ptOut)
		}

		for i, p := range RefreshParties {
			p.GenShares(p.s, ciphertext, crp, permute, p.share)
			if i > 0 {
				P0.Aggregate(p.share, P0.share, P0.share)
			}
		}

		P0.Transform(ciphertext, permute, crp, P0.share, ciphertext)

		coeffsPermute := make([]uint64, len(coeffs))
		for i := range coeffsPermute {
			coeffsPermute[i] = coeffs[permutation[i]]
		}

		coeffsHave := encoder.DecodeUintNew(decryptorSk0.DecryptNew(ciphertext))

		//Decrypts and compare
		require.True(t, utils.EqualSliceUint64(coeffsPermute, coeffsHave))
	})
}

func newTestVectors(testCtx *testContext, encryptor bfv.Encryptor, t *testing.T) (coeffs []uint64, plaintext *bfv.Plaintext, ciphertext *bfv.Ciphertext) {

	uniformSampler := ring.NewUniformSampler(testCtx.prng, testCtx.ringT)

	coeffsPol := uniformSampler.ReadNew()
	plaintext = bfv.NewPlaintext(testCtx.params)
	testCtx.encoder.EncodeUint(coeffsPol.Coeffs[0], plaintext)
	ciphertext = encryptor.EncryptNew(plaintext)
	return coeffsPol.Coeffs[0], plaintext, ciphertext
}

func verifyTestVectors(testCtx *testContext, decryptor bfv.Decryptor, coeffs []uint64, ciphertext *bfv.Ciphertext, t *testing.T) {
	require.True(t, utils.EqualSliceUint64(coeffs, testCtx.encoder.DecodeUintNew(decryptor.DecryptNew(ciphertext))))
}

func testMarshalling(testCtx *testContext, t *testing.T) {

	//verify if the un.marshalling works properly

	crsGen := ring.NewUniformSampler(testCtx.prng, testCtx.ringQP)
	crs := crsGen.ReadNew()
	ringQ := testCtx.ringQ

	ciphertext := bfv.NewCiphertextRandom(testCtx.prng, testCtx.params, 1)

	t.Run(fmt.Sprintf("Marshalling/Refresh/N=%d/limbQ=%d/limbsP=%d", ringQ.N, testCtx.params.QCount(), testCtx.params.PCount()), func(t *testing.T) {

		//testing refresh shares
		refreshproto := NewRefreshProtocol(testCtx.params, 3.2)
		refreshshare := refreshproto.AllocateShare()
		refreshproto.GenShares(testCtx.sk0, ciphertext, crs, refreshshare)

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
