package dbfv

import (
	"encoding/json"
	"testing"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

func Benchmark_DBFV(b *testing.B) {

	defaultParams := bfv.DefaultParams
	if testing.Short() {
		defaultParams = bfv.DefaultParams[:2]
	}
	if *flagParamString != "" {
		var jsonParams bfv.ParametersLiteral
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []bfv.ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, p := range defaultParams {
		params, err := bfv.NewParametersFromLiteral(p)
		if err != nil {
			panic(err)
		}
		var testCtx *testContext
		if testCtx, err = gentestContext(params); err != nil {
			panic(err)
		}

		benchPublicKeyGen(testCtx, b)
		benchRelinKeyGen(testCtx, b)
		benchKeyswitching(testCtx, b)
		benchPublicKeySwitching(testCtx, b)
		benchRotKeyGen(testCtx, b)
		benchRefresh(testCtx, b)
	}
}

func benchPublicKeyGen(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards

	crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.ringQP)

	crp := crpGenerator.ReadNew()

	type Party struct {
		*CKGProtocol
		s  *rlwe.SecretKey
		s1 *drlwe.CKGShare
	}

	p := new(Party)
	p.CKGProtocol = NewCKGProtocol(testCtx.params)
	p.s = sk0Shards[0]
	p.s1 = p.AllocateShares()

	b.Run(testString("PublicKeyGen/Round1/Gen", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, crp, p.s1)
		}
	})

	b.Run(testString("PublicKeyGen/Round1/Agg", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.s1, p.s1, p.s1)
		}
	})

	pk := bfv.NewPublicKey(testCtx.params)
	b.Run(testString("PublicKeyGen/Finalize", parties, testCtx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.GenPublicKey(p.s1, crp, pk)
		}
	})
}

func benchRelinKeyGen(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards

	type Party struct {
		*RKGProtocol
		ephSk  *rlwe.SecretKey
		sk     *rlwe.SecretKey
		share1 *drlwe.RKGShare
		share2 *drlwe.RKGShare

		rlk *rlwe.RelinearizationKey
	}

	p := new(Party)
	p.RKGProtocol = NewRKGProtocol(testCtx.params)
	p.sk = sk0Shards[0]
	p.ephSk, p.share1, p.share2 = p.RKGProtocol.AllocateShares()
	p.rlk = bfv.NewRelinearizationKey(testCtx.params, 2)

	crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.ringQP)

	crp := make([]*ring.Poly, testCtx.params.Beta())

	for i := 0; i < testCtx.params.Beta(); i++ {
		crp[i] = crpGenerator.ReadNew()
	}

	b.Run(testString("RelinKeyGen/Round1/Gen", parties, testCtx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.GenShareRoundOne(p.sk, crp, p.ephSk, p.share1)
		}
	})

	b.Run(testString("RelinKeyGen/Round1/Agg", parties, testCtx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share1, p.share1, p.share1)
		}
	})

	b.Run(testString("RelinKeyGen/Round2/Gen", parties, testCtx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.GenShareRoundTwo(p.ephSk, p.sk, p.share1, crp, p.share2)
		}
	})

	b.Run(testString("RelinKeyGen/Round2/Agg", parties, testCtx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share2, p.share2, p.share2)
		}
	})

	b.Run(testString("RelinKeyGen/Finalize", parties, testCtx.params), func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			p.GenRelinearizationKey(p.share1, p.share2, p.rlk)
		}
	})
}

func benchKeyswitching(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards
	sk1Shards := testCtx.sk1Shards

	type Party struct {
		*CKSProtocol
		s0    *rlwe.SecretKey
		s1    *rlwe.SecretKey
		share *drlwe.CKSShare
	}

	ciphertext := bfv.NewCiphertextRandom(testCtx.prng, testCtx.params, 1)

	p := new(Party)
	p.CKSProtocol = NewCKSProtocol(testCtx.params, 6.36)
	p.s0 = sk0Shards[0]
	p.s1 = sk1Shards[0]
	p.share = p.AllocateShare(ciphertext.Level())

	b.Run(testString("Keyswitching/Round1/Gen", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s0, p.s1, ciphertext, p.share)
		}
	})

	b.Run(testString("Keyswitching/Round1/Agg", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share, p.share, p.share)
		}
	})

	b.Run(testString("Keyswitching/Finalize", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.KeySwitch(p.share, ciphertext, ciphertext)
		}
	})
}

func benchPublicKeySwitching(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards
	pk1 := testCtx.pk1

	ciphertext := bfv.NewCiphertextRandom(testCtx.prng, testCtx.params, 1)

	type Party struct {
		*PCKSProtocol
		s     *rlwe.SecretKey
		share *drlwe.PCKSShare
	}

	p := new(Party)
	p.PCKSProtocol = NewPCKSProtocol(testCtx.params, 6.36)
	p.s = sk0Shards[0]
	p.share = p.AllocateShares(ciphertext.Level())

	b.Run(testString("PublicKeySwitching/Round1/Gen", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, pk1, ciphertext, p.share)

		}
	})

	b.Run(testString("PublicKeySwitching/Round1/Agg", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share, p.share, p.share)
		}
	})

	b.Run(testString("PublicKeySwitching/Finalize", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.KeySwitch(p.share, ciphertext, ciphertext)
		}
	})
}

func benchRotKeyGen(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards

	type Party struct {
		*RTGProtocol
		s     *rlwe.SecretKey
		share *drlwe.RTGShare
	}

	p := new(Party)
	p.RTGProtocol = NewRotKGProtocol(testCtx.params)
	p.s = sk0Shards[0]
	p.share = p.AllocateShares()

	crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.ringQP)
	crp := make([]*ring.Poly, testCtx.params.Beta())

	for i := 0; i < testCtx.params.Beta(); i++ {
		crp[i] = crpGenerator.ReadNew()
	}

	b.Run(testString("RotKeyGen/Round1/Gen", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, testCtx.params.GaloisElementForRowRotation(), crp, p.share)
		}
	})

	b.Run(testString("RotKeyGen/Round1/Agg", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.Aggregate(p.share, p.share, p.share)
		}
	})

	rotKey := bfv.NewSwitchingKey(testCtx.params)
	b.Run(testString("RotKeyGen/Finalize", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenRotationKey(p.share, crp, rotKey)
		}
	})
}

func benchRefresh(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards

	type Party struct {
		*RefreshProtocol
		s     *rlwe.SecretKey
		share RefreshShare
	}

	p := new(Party)
	p.RefreshProtocol = NewRefreshProtocol(testCtx.params, 3.2)
	p.s = sk0Shards[0]
	p.share = p.AllocateShares()

	crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.ringQ)
	crp := crpGenerator.ReadNew()

	ciphertext := bfv.NewCiphertextRandom(testCtx.prng, testCtx.params, 1)

	b.Run(testString("Refresh/Round1/Gen", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShares(p.s, ciphertext, crp, p.share)
		}
	})

	b.Run(testString("Refresh/Round1/Agg", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.Aggregate(p.share, p.share, p.share)
		}
	})

	b.Run(testString("Refresh/Finalize", parties, testCtx.params), func(b *testing.B) {
		ctOut := bfv.NewCiphertext(testCtx.params, 1)
		for i := 0; i < b.N; i++ {
			p.Finalize(ciphertext, crp, p.share, ctOut)
		}
	})
}
