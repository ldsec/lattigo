// Package bfv implements a RNS-accelerated Fan-Vercauteren version of Brakerski's scale invariant homomorphic encryption scheme. It provides modular arithmetic over the integers.
package bfv

import (
	"fmt"
	"math/big"
	"unsafe"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// GaloisGen is an integer of order N=2^d modulo M=2N and that spans Z_M with the integer -1.
// The j-th ring automorphism takes the root zeta to zeta^(5j).
const GaloisGen uint64 = 5

// Encoder is an interface for plaintext encoding and decoding operations. It provides methods to embed []uint64 and []int64 types into
// the various plaintext types and the inverse operations. It also provides methodes to convert between the different plaintext types.
// The different plaintext types represent different embedings of the message in the polynomial space. This relation is illustrated in
// The figure below:
//
// []uint64 --- Encoder.EncodeUintRingT(.) -┬-> PlaintextRingT -┬-> Encoder.ScaleUp(.) -----> Plaintext
// []uint64 --- Encoder.EncodeIntRingT(.) --┘                   └-> Encoder.RingTToMul(.) ---> PlaintextMul
//
//
// The different plaintext types have different efficiency-related caracteristics that we summarize in the Table below. For more information
// about the different plaintext types, see plaintext.go.
//
// Relative efficiency of operation
//  -----------------------------------------------------------------------
// |                      |  PlaintextRingT  |  Plaintext  | PlaintextMul  |
//  -----------------------------------------------------------------------
// | Encoding/Decoding    |    Faster      |    Slower   |    Slower       |
// | Memory size          |    Smaller     |    Larger   |    Larger       |
// | Ct-Pt Add / Sub      |    Slower      |    Faster   |    N/A          |
// | Ct-Pt Mul            |    Faster      |    Slower   |    Much Faster  |
//  -----------------------------------------------------------------------
//
type Encoder interface {
	EncodeUint(coeffs []uint64, pt *Plaintext)
	EncodeUintRingT(coeffs []uint64, pt *PlaintextRingT)
	EncodeUintMul(coeffs []uint64, pt *PlaintextMul)
	EncodeInt(coeffs []int64, pt *Plaintext)
	EncodeIntRingT(coeffs []int64, pt *PlaintextRingT)
	EncodeIntMul(coeffs []int64, pt *PlaintextMul)

	ScaleUp(*PlaintextRingT, *Plaintext)
	ScaleDown(pt *Plaintext, ptRt *PlaintextRingT)
	RingTToMul(ptRt *PlaintextRingT, ptmul *PlaintextMul)
	MulToRingT(pt *PlaintextMul, ptRt *PlaintextRingT)

	DecodeRingT(pt interface{}, ptRt *PlaintextRingT)
	DecodeUint(pt interface{}, coeffs []uint64)
	DecodeInt(pt interface{}, coeffs []int64)
	DecodeUintNew(pt interface{}) (coeffs []uint64)
	DecodeIntNew(pt interface{}) (coeffs []int64)
}

// Encoder is a structure that stores the parameters to encode values on a plaintext in a SIMD (Single-Instruction Multiple-Data) fashion.
type encoder struct {
	params *Parameters

	ringQ *ring.Ring
	ringT *ring.Ring

	indexMatrix []uint64
	scaler      ring.Scaler
	deltaMont   []ring.FastBRedOperand

	tmpPoly *ring.Poly
	tmpPtRt *PlaintextRingT
}

// NewEncoder creates a new encoder from the provided parameters.
func NewEncoder(params *Parameters) Encoder {

	var ringQ, ringT *ring.Ring
	var err error

	if ringQ, err = ring.NewRing(params.N(), params.qi); err != nil {
		panic(err)
	}

	if ringT, err = ring.NewRing(params.N(), []uint64{params.t}); err != nil {
		panic(err)
	}

	var m, pos, index1, index2 uint64

	slots := params.N()

	indexMatrix := make([]uint64, slots)

	logN := params.LogN()

	rowSize := params.N() >> 1
	m = (params.N() << 1)
	pos = 1

	for i := uint64(0); i < rowSize; i++ {

		index1 = (pos - 1) >> 1
		index2 = (m - pos - 1) >> 1

		indexMatrix[i] = utils.BitReverse64(index1, logN)
		indexMatrix[i|rowSize] = utils.BitReverse64(index2, logN)

		pos *= GaloisGen
		pos &= (m - 1)
	}

	return &encoder{
		params:      params.Copy(),
		ringQ:       ringQ,
		ringT:       ringT,
		indexMatrix: indexMatrix,
		deltaMont:   GenLiftParams(ringQ, params.t),
		scaler:      ring.NewRNSScaler(params.t, ringQ),
		tmpPoly:     ringT.NewPoly(),
		tmpPtRt:     NewPlaintextRingT(params),
	}
}

// GenLiftParams generates the lifting parameters.
func GenLiftParams(ringQ *ring.Ring, t uint64) (deltaMont []ring.FastBRedOperand) {

	delta := new(big.Int).Quo(ringQ.ModulusBigint, ring.NewUint(t))

	deltaMont = make([]ring.FastBRedOperand, len(ringQ.Modulus))

	tmp := new(big.Int)
	for i, Qi := range ringQ.Modulus {
		deltaMont[i] = ring.NewFastBRedOperand(tmp.Mod(delta, ring.NewUint(Qi)).Uint64(), Qi)
	}

	return
}

// EncodeUintRingT encodes a slice of uint64 into a Plaintext in R_t
func (encoder *encoder) EncodeUintRingT(coeffs []uint64, p *PlaintextRingT) {
	if len(coeffs) > len(encoder.indexMatrix) {
		panic("invalid input to encode: number of coefficients must be smaller or equal to the ring degree")
	}

	if len(p.value.Coeffs[0]) != len(encoder.indexMatrix) {
		panic("invalid plaintext to receive encoding: number of coefficients does not match the ring degree")
	}

	for i := 0; i < len(coeffs); i++ {
		p.value.Coeffs[0][encoder.indexMatrix[i]] = coeffs[i]
	}

	for i := len(coeffs); i < len(encoder.indexMatrix); i++ {
		p.value.Coeffs[0][encoder.indexMatrix[i]] = 0
	}

	encoder.ringT.InvNTT(p.value, p.value)
}

// EncodeUint encodes an uint64 slice of size at most N on a plaintext.
func (encoder *encoder) EncodeUint(coeffs []uint64, p *Plaintext) {
	ptRt := &PlaintextRingT{p.Element, p.Element.value[0]}

	// Encodes the values in RingT
	encoder.EncodeUintRingT(coeffs, ptRt)

	// Scales by Q/t
	encoder.ScaleUp(ptRt, p)
}

func (encoder *encoder) EncodeUintMul(coeffs []uint64, p *PlaintextMul) {

	ptRt := &PlaintextRingT{p.Element, p.Element.value[0]}

	// Encodes the values in RingT
	encoder.EncodeUintRingT(coeffs, ptRt)

	// Puts in NTT+Montgomerry domains of ringQ
	encoder.RingTToMul(ptRt, p)
}

// EncodeInt encodes an int64 slice of size at most N on a plaintext. It also encodes the sign of the given integer (as its inverse modulo the plaintext modulus).
// The sign will correctly decode as long as the absolute value of the coefficient does not exceed half of the plaintext modulus.
func (encoder *encoder) EncodeIntRingT(coeffs []int64, p *PlaintextRingT) {

	if len(coeffs) > len(encoder.indexMatrix) {
		panic("invalid input to encode: number of coefficients must be smaller or equal to the ring degree")
	}

	if len(p.value.Coeffs[0]) != len(encoder.indexMatrix) {
		panic("invalid plaintext to receive encoding: number of coefficients does not match the ring degree")
	}

	for i := 0; i < len(coeffs); i++ {

		if coeffs[i] < 0 {
			p.value.Coeffs[0][encoder.indexMatrix[i]] = uint64(int64(encoder.params.t) + coeffs[i])
		} else {
			p.value.Coeffs[0][encoder.indexMatrix[i]] = uint64(coeffs[i])
		}
	}

	for i := len(coeffs); i < len(encoder.indexMatrix); i++ {
		p.value.Coeffs[0][encoder.indexMatrix[i]] = 0
	}

	encoder.ringT.InvNTTLazy(p.value, p.value)
}

func (encoder *encoder) EncodeInt(coeffs []int64, p *Plaintext) {
	ptRt := &PlaintextRingT{p.Element, p.value}

	// Encodes the values in RingT
	encoder.EncodeIntRingT(coeffs, ptRt)

	// Scales by Q/t
	encoder.ScaleUp(ptRt, p)
}

func (encoder *encoder) EncodeIntMul(coeffs []int64, p *PlaintextMul) {
	ptRt := &PlaintextRingT{p.Element, p.value}

	// Encodes the values in RingT
	encoder.EncodeIntRingT(coeffs, ptRt)

	// Puts in NTT+Montgomerry domains of ringQ
	encoder.RingTToMul(ptRt, p)
}

// ScaleUp transforms a PlaintextRingT (R_t) into a Plaintext (R_q) by scaling up the coefficient by Q/t.
func (encoder *encoder) ScaleUp(ptRt *PlaintextRingT, pt *Plaintext) {
	scaleUp(encoder.ringQ, encoder.deltaMont, ptRt.value, pt.value)
}

func scaleUp(ringQ *ring.Ring, deltaMont []ring.FastBRedOperand, pIn, pOut *ring.Poly) {

	for i := len(ringQ.Modulus) - 1; i >= 0; i-- {
		out := pOut.Coeffs[i]
		in := pIn.Coeffs[0]
		d := deltaMont[i]
		qi := ringQ.Modulus[i]

		for j := uint64(0); j < ringQ.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&in[j]))
			z := (*[8]uint64)(unsafe.Pointer(&out[j]))

			z[0] = ring.FastBRed(x[0], d, qi)
			z[1] = ring.FastBRed(x[1], d, qi)
			z[2] = ring.FastBRed(x[2], d, qi)
			z[3] = ring.FastBRed(x[3], d, qi)
			z[4] = ring.FastBRed(x[4], d, qi)
			z[5] = ring.FastBRed(x[5], d, qi)
			z[6] = ring.FastBRed(x[6], d, qi)
			z[7] = ring.FastBRed(x[7], d, qi)
		}
	}
}

// ScaleDown transforms a Plaintext (R_q) into a PlaintextRingT (R_t) by scaling down the coefficient by t/Q and rounding.
func (encoder *encoder) ScaleDown(pt *Plaintext, ptRt *PlaintextRingT) {
	encoder.scaler.DivByQOverTRounded(pt.value, ptRt.value)
}

// RingTToMul transforms a PlaintextRingT into a PlaintextMul by operating the NTT transform
// of R_q and putting the coefficients in Montgommery form.
func (encoder *encoder) RingTToMul(ptRt *PlaintextRingT, ptMul *PlaintextMul) {
	if ptRt.value != ptMul.value {
		copy(ptMul.value.Coeffs[0], ptRt.value.Coeffs[0])
	}
	for i := 1; i < len(encoder.ringQ.Modulus); i++ {
		copy(ptMul.value.Coeffs[i], ptRt.value.Coeffs[0])
	}

	encoder.ringQ.NTTLazy(ptMul.value, ptMul.value)
	encoder.ringQ.MForm(ptMul.value, ptMul.value)
}

// MulToRingT transforms a PlaintextMul into PlaintextRingT by operating the inverse NTT transform of R_q and
// putting the coefficients out of the Montgommery form.
func (encoder *encoder) MulToRingT(pt *PlaintextMul, ptRt *PlaintextRingT) {
	encoder.ringQ.InvNTTLvl(0, pt.value, ptRt.value)
	encoder.ringQ.InvMFormLvl(0, ptRt.value, ptRt.value)
}

// DecodeRingT decodes any plaintext type into a PlaintextRingT. It panics if p is not PlaintextRingT, Plaintext or PlaintextMul.
func (encoder *encoder) DecodeRingT(p interface{}, ptRt *PlaintextRingT) {
	switch pt := p.(type) {
	case *Plaintext:
		encoder.ScaleDown(pt, ptRt)
	case *PlaintextMul:
		encoder.MulToRingT(pt, ptRt)
	case *PlaintextRingT:
		ptRt.Copy(pt.Element)
	default:
		panic(fmt.Errorf("unsuported plaintext type (%T)", pt))
	}
}

// DecodeUint decodes a any plaintext type and write the coefficients in coeffs. It panics if p is not PlaintextRingT, Plaintext or PlaintextMul.
func (encoder *encoder) DecodeUint(p interface{}, coeffs []uint64) {

	var ptRt *PlaintextRingT
	var isInRingT bool
	if ptRt, isInRingT = p.(*PlaintextRingT); !isInRingT {
		encoder.DecodeRingT(p, encoder.tmpPtRt)
		ptRt = encoder.tmpPtRt
	}

	encoder.ringT.NTT(ptRt.value, encoder.tmpPoly)

	for i := uint64(0); i < encoder.ringQ.N; i++ {
		coeffs[i] = encoder.tmpPoly.Coeffs[0][encoder.indexMatrix[i]]
	}
}

// DecodeUintNew decodes any plaintext type and returns the coefficients in a new []uint64.
// It panics if p is not PlaintextRingT, Plaintext or PlaintextMul.
func (encoder *encoder) DecodeUintNew(p interface{}) (coeffs []uint64) {
	coeffs = make([]uint64, encoder.ringQ.N)
	encoder.DecodeUint(p, coeffs)
	return
}

// DecodeInt decodes a any plaintext type and write the coefficients in coeffs. It also decodes the sign
// modulus (by centering the values around the plaintext). It panics if p is not PlaintextRingT, Plaintext or PlaintextMul.
func (encoder *encoder) DecodeInt(p interface{}, coeffs []int64) {

	encoder.DecodeRingT(p, encoder.tmpPtRt)

	encoder.ringT.NTT(encoder.tmpPtRt.value, encoder.tmpPoly)

	modulus := int64(encoder.params.t)
	modulusHalf := modulus >> 1
	var value int64
	for i := uint64(0); i < encoder.ringQ.N; i++ {

		value = int64(encoder.tmpPoly.Coeffs[0][encoder.indexMatrix[i]])
		coeffs[i] = value
		if value >= modulusHalf {
			coeffs[i] -= modulus
		}
	}
}

// DecodeIntNew decodes any plaintext type and returns the coefficients in a new []int64. It also decodes the sign
// modulus (by centering the values around the plaintext). It panics if p is not PlaintextRingT, Plaintext or PlaintextMul.
func (encoder *encoder) DecodeIntNew(p interface{}) (coeffs []int64) {
	coeffs = make([]int64, encoder.ringQ.N)
	encoder.DecodeInt(p, coeffs)
	return
}
