package rckks

import (
	"github.com/ldsec/lattigo/v2/ring"
)

// Plaintext is is a Element with only one Poly.
type Plaintext struct {
	*Element
	value *ring.Poly
}

// NewPlaintext creates a new Plaintext of level level and scale scale.
func NewPlaintext(params *Parameters, level uint64, scale float64) *Plaintext {

	plaintext := &Plaintext{&Element{}, nil}

	plaintext.Element.value = []*ring.Poly{ring.NewPoly(params.N(), level+1)}

	plaintext.value = plaintext.Element.value[0]

	plaintext.scale = scale
	plaintext.isNTT = true

	return plaintext
}
