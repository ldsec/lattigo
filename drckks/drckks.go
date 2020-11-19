package drckks

import (
	"github.com/ldsec/lattigo/v2/rckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

type drckksContext struct {
	params *rckks.Parameters

	n uint64

	ringQ  *ring.Ring
	ringP  *ring.Ring
	ringQP *ring.Ring

	alpha uint64
	beta  uint64
}

func newDrckksContext(params *rckks.Parameters) (context *drckksContext) {

	context = new(drckksContext)

	context.params = params.Copy()

	context.n = params.N()

	context.alpha = params.Alpha()
	context.beta = params.Beta()

	var err error
	if context.ringQ, err = ring.NewRingWithNthRoot(params.N(), params.N()<<2, params.Qi()); err != nil {
		panic(err)
	}

	if context.ringP, err = ring.NewRingWithNthRoot(params.N(), params.N()<<2, params.Pi()); err != nil {
		panic(err)
	}

	if context.ringQP, err = ring.NewRingWithNthRoot(params.N(), params.N()<<2, append(params.Qi(), params.Pi()...)); err != nil {
		panic(err)
	}

	return
}

// NewCRPGenerator creates a new deterministic random polynomial generator.
func NewCRPGenerator(params *rckks.Parameters, key []byte) *ring.UniformSampler {
	ctx := newDrckksContext(params)
	prng, err := utils.NewKeyedPRNG(key)
	if err != nil {
		panic(err)
	}
	return ring.NewUniformSampler(prng, ctx.ringQP)
}
