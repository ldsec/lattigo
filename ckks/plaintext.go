package ckks

import (
	"github.com/ldsec/lattigo/v2/rlwe"
)

// Plaintext is is a Element with only one Poly.
type Plaintext struct {
	*rlwe.Plaintext
	Scale float64
}

// NewPlaintext creates a new Plaintext of level level and scale scale.
func NewPlaintext(params Parameters, level int, scale float64) *Plaintext {
	pt := &Plaintext{Plaintext: rlwe.NewPlaintext(params.Parameters, level), Scale: scale}
	pt.Value.IsNTT = true
	return pt
}

// ScalingFactor returns the scaling factor of the plaintext
func (p *Plaintext) ScalingFactor() float64 {
	return p.Scale
}

// SetScalingFactor sets the scaling factor of the target plaintext
func (p *Plaintext) SetScalingFactor(scale float64) {
	p.Scale = scale
}
