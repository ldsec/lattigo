package rlwe

import (
	"math/big"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// Plaintext is a common base type for RLWE plaintexts.
type Plaintext struct {
	Value *ring.Poly
}

// Ciphertext is a generic type for RLWE ciphertext.
type Ciphertext struct {
	Value []*ring.Poly
}

// AdditiveShare is a type for storing additively shared values in Z_Q[X] (RNS domain)
type AdditiveShare struct {
	Value ring.Poly
}

// AdditiveShareBigint is a type for storing additively shared values
// in Z (positional domain)
type AdditiveShareBigint struct {
	Value []*big.Int
}

// NewAdditiveShare instantiate a new additive share struct for the ring defined
// by the given parameters at maximum level.
func NewAdditiveShare(params Parameters) *AdditiveShare {
	return &AdditiveShare{Value: *ring.NewPoly(params.N(), 1)}
}

// NewAdditiveShareAtLevel instantiate a new additive share struct for the ring defined
// by the given parameters at level `level`.
func NewAdditiveShareAtLevel(params Parameters, level int) *AdditiveShare {
	return &AdditiveShare{Value: *ring.NewPoly(params.N(), level+1)}
}

// NewAdditiveShareBigint instantiate a new additive share struct composed of big.Int elements
func NewAdditiveShareBigint(params Parameters) *AdditiveShareBigint {
	v := make([]*big.Int, params.N())
	for i := range v {
		v[i] = new(big.Int)
	}
	return &AdditiveShareBigint{Value: v}
}

// NewPlaintext creates a new Plaintext at level `level` from the parameters.
func NewPlaintext(params Parameters, level int) *Plaintext {
	return &Plaintext{Value: ring.NewPoly(params.N(), level+1)}
}

// Degree returns the degree of the target element.
func (pt Plaintext) Degree() int {
	return 0
}

// Level returns the level of the target element.
func (pt Plaintext) Level() int {
	return len(pt.Value.Coeffs) - 1
}

// El returns the plaintext as a new `Element` for which the value points
// to the receiver `Value` field.
func (pt Plaintext) El() *Ciphertext {
	return &Ciphertext{Value: []*ring.Poly{pt.Value}}
}

// Copy copies the `other` plaintext value into the reciever plaintext.
func (pt *Plaintext) Copy(other *Plaintext) {
	if other != nil && other.Value != nil {
		pt.Value.Copy(other.Value)
	}
}

// NewCiphertext returns a new Element with zero values.
func NewCiphertext(params Parameters, degree, level int) *Ciphertext {
	el := new(Ciphertext)
	el.Value = make([]*ring.Poly, degree+1)
	for i := 0; i < degree+1; i++ {
		el.Value[i] = ring.NewPoly(params.N(), level+1)
	}
	return el
}

// NewCiphertextNTT returns a new Element with zero values and the NTT flags set.
func NewCiphertextNTT(params Parameters, degree, level int) *Ciphertext {
	el := new(Ciphertext)
	el.Value = make([]*ring.Poly, degree+1)
	for i := 0; i < degree+1; i++ {
		el.Value[i] = ring.NewPoly(params.N(), level+1)
		el.Value[i].IsNTT = true
	}
	return el
}

// SetValue sets the input slice of polynomials as the value of the target element.
func (el *Ciphertext) SetValue(value []*ring.Poly) {
	el.Value = value
}

// Degree returns the degree of the target element.
func (el *Ciphertext) Degree() int {
	return len(el.Value) - 1
}

// Level returns the level of the target element.
func (el *Ciphertext) Level() int {
	return len(el.Value[0].Coeffs) - 1
}

// Resize resizes the degree of the target element.
// Sets the NTT flag of the added poly equal to the NTT flag
// to the poly at degree zero.
func (el *Ciphertext) Resize(params Parameters, degree int) {
	if el.Degree() > degree {
		el.Value = el.Value[:degree+1]
	} else if el.Degree() < degree {
		for el.Degree() < degree {
			el.Value = append(el.Value, []*ring.Poly{new(ring.Poly)}...)
			el.Value[el.Degree()].Coeffs = make([][]uint64, el.Level()+1)
			for i := 0; i < el.Level()+1; i++ {
				el.Value[el.Degree()].Coeffs[i] = make([]uint64, params.N())
				el.Value[el.Degree()].IsNTT = el.Value[0].IsNTT
			}
		}
	}
}

// CopyNew creates a new element as a copy of the target element.
func (el *Ciphertext) CopyNew() *Ciphertext {

	ctxCopy := new(Ciphertext)

	ctxCopy.Value = make([]*ring.Poly, el.Degree()+1)
	for i := range el.Value {
		ctxCopy.Value[i] = el.Value[i].CopyNew()
	}

	return ctxCopy
}

// Copy copies the input element and its parameters on the target element.
func (el *Ciphertext) Copy(ctxCopy *Ciphertext) {

	if el != ctxCopy {
		for i := range ctxCopy.Value {
			el.Value[i].Copy(ctxCopy.Value[i])
		}
	}
}

// El returns a pointer to this Element
func (el *Ciphertext) El() *Ciphertext {
	return el
}

// RLWEElement returns a pointer to this Element
func (el *Ciphertext) RLWEElement() *Ciphertext {
	return el
}

// GetSmallestLargest returns the provided element that has the smallest degree as a first
// returned value and the largest degree as second return value. If the degree match, the
// order is the same as for the input.
func GetSmallestLargest(el0, el1 *Ciphertext) (smallest, largest *Ciphertext, sameDegree bool) {
	switch {
	case el0.Degree() > el1.Degree():
		return el1, el0, false
	case el0.Degree() < el1.Degree():
		return el0, el1, false
	}
	return el0, el1, true
}

// PopulateElementRandom creates a new rlwe.Element with random coefficients
func PopulateElementRandom(prng utils.PRNG, params Parameters, el *Ciphertext) {

	ringQ, err := ring.NewRing(params.N(), params.Q())
	if err != nil {
		panic(err)
	}
	sampler := ring.NewUniformSampler(prng, ringQ)
	for i := range el.Value {
		sampler.Read(el.Value[i])
	}
}
