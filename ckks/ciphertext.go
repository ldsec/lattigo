package ckks

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// Ciphertext is *ring.Poly array representing a polynomial of degree > 0 with coefficients in R_Q.
type Ciphertext struct {
	*rlwe.Ciphertext
	Scale float64
}

// NewCiphertext creates a new Ciphertext parameterized by degree, level and scale.
func NewCiphertext(params Parameters, degree, level int, scale float64) (ciphertext *Ciphertext) {
	ciphertext = &Ciphertext{Ciphertext: rlwe.NewCiphertext(params.Parameters, degree, level)}
	for _, pol := range ciphertext.Value {
		pol.IsNTT = true
	}
	ciphertext.Scale = scale
	return ciphertext
}

// NewCiphertextRandom generates a new uniformly distributed Ciphertext of degree, level and scale.
func NewCiphertextRandom(prng utils.PRNG, params Parameters, degree, level int, scale float64) (ciphertext *Ciphertext) {

	ringQ, err := ring.NewRing(params.N(), params.Q()[:level+1])
	if err != nil {
		panic(err)
	}

	sampler := ring.NewUniformSampler(prng, ringQ)
	ciphertext = NewCiphertext(params, degree, level, scale)
	for i := 0; i < degree+1; i++ {
		sampler.Read(ciphertext.Value[i])
	}

	return ciphertext
}

// ScalingFactor returns the scaling factor of the ciphertext
func (ct *Ciphertext) ScalingFactor() float64 {
	return ct.Scale
}

// SetScalingFactor sets the scaling factor of the ciphertext
func (ct *Ciphertext) SetScalingFactor(scale float64) {
	ct.Scale = scale
}

// Copy copies the given ciphertext ctp into the receiver ciphertext.
func (ct *Ciphertext) Copy(ctp *Ciphertext) {
	ct.Ciphertext.Copy(ctp.Ciphertext)
	ct.Scale = ctp.Scale
}

// CopyNew makes a deep copy of the receiver ciphertext and returns it.
func (ct *Ciphertext) CopyNew() (ctc *Ciphertext) {
	ctc = &Ciphertext{Ciphertext: ct.Ciphertext.CopyNew(), Scale: ct.Scale}
	return
}
