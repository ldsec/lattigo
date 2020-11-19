package drckks

import (
	"testing"

	"github.com/ldsec/lattigo/v2/rckks"
	"github.com/ldsec/lattigo/v2/ring"
)

func BenchmarkDrckks(b *testing.B) {

	var err error
	var testCtx = new(testContext)
	var defaultParams []*rckks.Parameters

	if testing.Short() {
		defaultParams = rckks.DefaultParams[rckks.PN12QP109 : rckks.PN12QP109+3]
	} else {
		defaultParams = rckks.DefaultParams
	}

	for _, p := range defaultParams {

		if testCtx, err = genTestParams(p); err != nil {
			panic(err)
		}

		benchPublicKeyGen(testCtx, b)
		benchRelinKeyGen(testCtx, b)
		benchKeySwitching(testCtx, b)
		benchPublicKeySwitching(testCtx, b)
		benchRotKeyGen(testCtx, b)
		benchRefresh(testCtx, b)
	}
}

func benchPublicKeyGen(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards

	crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.drckksContext.ringQP)
	crp := crpGenerator.ReadNew()

	type Party struct {
		*CKGProtocol
		s  *ring.Poly
		s1 CKGShare
	}

	p := new(Party)
	p.CKGProtocol = NewCKGProtocol(testCtx.params)
	p.s = sk0Shards[0].Get()
	p.s1 = p.AllocateShares()

	b.Run(testString("PublicKeyGen/Gen/", parties, testCtx.params), func(b *testing.B) {

		// Each party creates a new CKGProtocol instance
		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, crp, p.s1)
		}
	})

	b.Run(testString("PublicKeyGen/Agg/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.s1, p.s1, p.s1)
		}
	})

}

func benchRelinKeyGen(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards

	type Party struct {
		*RKGProtocol
		u      *ring.Poly
		s      *ring.Poly
		share1 RKGShare
		share2 RKGShare
	}

	p := new(Party)
	p.RKGProtocol = NewEkgProtocol(testCtx.params)
	p.u = p.RKGProtocol.NewEphemeralKey()
	p.s = sk0Shards[0].Get()
	p.share1, p.share2 = p.RKGProtocol.AllocateShares()

	crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.drckksContext.ringQP)
	crp := make([]*ring.Poly, testCtx.params.Beta())

	for i := uint64(0); i < testCtx.params.Beta(); i++ {
		crp[i] = crpGenerator.ReadNew()
	}

	b.Run(testString("RelinKeyGen/Round1Gen/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShareRoundOne(p.u, p.s, crp, p.share1)
		}
	})

	b.Run(testString("RelinKeyGen/Round1Agg/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShareRoundOne(p.share1, p.share1, p.share1)
		}
	})

	b.Run(testString("RelinKeyGen/Round2Gen/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShareRoundTwo(p.share1, p.u, p.s, crp, p.share2)
		}
	})

	b.Run(testString("RelinKeyGen/Round2Agg/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShareRoundTwo(p.share2, p.share2, p.share2)
		}
	})

}

func benchKeySwitching(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards
	sk1Shards := testCtx.sk1Shards

	type Party struct {
		*CKSProtocol
		s0    *ring.Poly
		s1    *ring.Poly
		share CKSShare
	}

	p := new(Party)
	p.CKSProtocol = NewCKSProtocol(testCtx.params, 6.36)
	p.s0 = sk0Shards[0].Get()
	p.s1 = sk1Shards[0].Get()
	p.share = p.AllocateShare()

	ciphertext := rckks.NewCiphertextRandom(testCtx.prng, testCtx.params, 1, testCtx.params.MaxLevel(), testCtx.params.Scale())

	b.Run(testString("KeySwitching/Gen/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s0, p.s1, ciphertext, p.share)
		}
	})

	b.Run(testString("KeySwitching/Agg/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share, p.share, p.share)
		}
	})

	b.Run(testString("KeySwitching/KS/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.KeySwitch(p.share, ciphertext, ciphertext)
		}
	})
}

func benchPublicKeySwitching(testCtx *testContext, b *testing.B) {

	sk0Shards := testCtx.sk0Shards
	pk1 := testCtx.pk1

	ciphertext := rckks.NewCiphertextRandom(testCtx.prng, testCtx.params, 1, testCtx.params.MaxLevel(), testCtx.params.Scale())

	type Party struct {
		*PCKSProtocol
		s     *ring.Poly
		share PCKSShare
	}

	p := new(Party)
	p.PCKSProtocol = NewPCKSProtocol(testCtx.params, 6.36)
	p.s = sk0Shards[0].Get()
	p.share = p.AllocateShares(ciphertext.Level())

	b.Run(testString("PublicKeySwitching/Gen/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(p.s, pk1, ciphertext, p.share)
		}
	})

	b.Run(testString("PublicKeySwitching/Agg/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.AggregateShares(p.share, p.share, p.share)
		}
	})

	b.Run(testString("PublicKeySwitching/KS/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.KeySwitch(p.share, ciphertext, ciphertext)
		}
	})
}

func benchRotKeyGen(testCtx *testContext, b *testing.B) {

	ringQP := testCtx.drckksContext.ringQP
	sk0Shards := testCtx.sk0Shards

	type Party struct {
		*RTGProtocol
		s     *ring.Poly
		share RTGShare
	}

	p := new(Party)
	p.RTGProtocol = NewRotKGProtocol(testCtx.params)
	p.s = sk0Shards[0].Get()
	p.share = p.AllocateShare()

	crpGenerator := ring.NewUniformSampler(testCtx.prng, ringQP)
	crp := make([]*ring.Poly, testCtx.params.Beta())

	for i := uint64(0); i < testCtx.params.Beta(); i++ {
		crp[i] = crpGenerator.ReadNew()
	}

	mask := uint64((ringQP.N >> 1) - 1)

	b.Run(testString("RotKeyGen/Round1/Gen/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShare(rckks.RotationRight, uint64(i)&mask, sk0Shards[0].Get(), crp, &p.share)
		}
	})

	b.Run(testString("RotKeyGen/Round1/Agg/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.Aggregate(p.share, p.share, p.share)
		}
	})

	rotKey := rckks.NewRotationKeys()
	b.Run(testString("RotKeyGen/Finalize/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.Finalize(testCtx.params, p.share, crp, rotKey)
		}
	})
}

func benchRefresh(testCtx *testContext, b *testing.B) {

	if testCtx.params.MaxLevel() < 3 {
		b.Skip()
	}

	sk0Shards := testCtx.sk0Shards
	ringQ := testCtx.drckksContext.ringQ

	levelStart := uint64(3)

	type Party struct {
		*RefreshProtocol
		s      *ring.Poly
		share1 RefreshShareDecrypt
		share2 RefreshShareRecrypt
	}

	p := new(Party)
	p.RefreshProtocol = NewRefreshProtocol(testCtx.params)
	p.s = sk0Shards[0].Get()
	p.share1, p.share2 = p.AllocateShares(levelStart)

	crpGenerator := ring.NewUniformSampler(testCtx.prng, ringQ)
	crp := crpGenerator.ReadNew()

	ciphertext := rckks.NewCiphertextRandom(testCtx.prng, testCtx.params, 1, levelStart, testCtx.params.Scale())

	b.Run(testString("Refresh/Gen/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShares(p.s, levelStart, parties, ciphertext, crp, p.share1, p.share2)
		}
	})

	b.Run(testString("Refresh/Agg/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.Aggregate(p.share1, p.share1, p.share1)
		}
	})

	b.Run(testString("Refresh/Decrypt/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.Decrypt(ciphertext, p.share1)
		}
	})

	b.Run(testString("Refresh/Recode/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.Recode(ciphertext)
		}
	})

	b.Run(testString("Refresh/Recrypt/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.Recrypt(ciphertext, crp, p.share2)
		}
	})
}

func benchRefreshAndPermute(testCtx *testContext, b *testing.B) {

	if testCtx.params.MaxLevel() < 3 {
		b.Skip()
	}

	sk0Shards := testCtx.sk0Shards

	levelStart := uint64(2)

	type Party struct {
		*PermuteProtocol
		s      *ring.Poly
		share1 RefreshShareDecrypt
		share2 RefreshShareRecrypt
	}

	p := new(Party)
	p.PermuteProtocol = NewPermuteProtocol(testCtx.params)
	p.s = sk0Shards[0].Get()
	p.share1, p.share2 = p.AllocateShares(levelStart)

	crpGenerator := ring.NewUniformSampler(testCtx.prng, testCtx.drckksContext.ringQP)
	crp := crpGenerator.ReadNew()

	ciphertext := rckks.NewCiphertextRandom(testCtx.prng, testCtx.params, 1, levelStart, testCtx.params.Scale())

	crpGenerator.Readlvl(levelStart, ciphertext.Value()[0])
	crpGenerator.Readlvl(levelStart, ciphertext.Value()[1])

	permutation := make([]uint64, testCtx.params.Slots())

	for i := range permutation {
		permutation[i] = ring.RandUniform(testCtx.prng, testCtx.params.Slots(), testCtx.params.Slots()-1)
	}

	b.Run(testString("RefreshAndPermute/Gen/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.GenShares(p.s, levelStart, parties, ciphertext, crp, testCtx.params.Slots(), permutation, p.share1, p.share2)
		}
	})

	b.Run(testString("RefreshAndPermute/Agg/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.Aggregate(p.share1, p.share1, p.share1)
		}
	})

	b.Run(testString("RefreshAndPermute/Decrypt/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.Decrypt(ciphertext, p.share1)
		}
	})

	b.Run(testString("RefreshAndPermute/Permute/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.Permute(ciphertext, permutation, testCtx.params.Slots())
		}
	})

	b.Run(testString("RefreshAndPermute/Recrypt/", parties, testCtx.params), func(b *testing.B) {

		for i := 0; i < b.N; i++ {
			p.Recrypt(ciphertext, crp, p.share2)
		}
	})
}
