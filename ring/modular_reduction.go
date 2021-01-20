package ring

import (
	"math/big"
	"math/bits"
)

// MForm switches a to the Montgomery domain by computing
// a*2^64 mod q.
func MForm(a, q uint64, u []uint64) (r uint64) {
	mhi, _ := bits.Mul64(a, u[1])
	r = -(a*u[0] + mhi) * q
	if r >= q {
		r -= q
	}
	return
}

// MFormConstant switches a to the Montgomery domain by computing
// a*2^64 mod q in constant time.
// The result is between 0 and 2*q-1.
func MFormConstant(a, q uint64, u []uint64) (r uint64) {
	mhi, _ := bits.Mul64(a, u[1])
	r = -(a*u[0] + mhi) * q
	return
}

// InvMForm switches a from the Montgomery domain back to the
// standard domain by computing a*(1/2^64) mod q.
func InvMForm(a, q, qInv uint64) (r uint64) {
	r, _ = bits.Mul64(a*qInv, q)
	r = q - r
	if r >= q {
		r -= q
	}
	return
}

// InvMFormConstant switches a from the Montgomery domain back to the
// standard domain by computing a*(1/2^64) mod q in constant time.
// The result is between 0 and 2*q-1.
func InvMFormConstant(a, q, qInv uint64) (r uint64) {
	r, _ = bits.Mul64(a*qInv, q)
	r = q - r
	return
}

// MRedParams computes the parameter qInv = (q^-1) mod 2^64,
// required for MRed.
func MRedParams(q uint64) (qInv uint64) {
	qInv = 1
	for i := 0; i < 63; i++ {
		qInv *= q
		q *= q
	}
	return
}

// MRed computes x * y * (1/2^64) mod q.
func MRed(x, y, q, qInv uint64) (r uint64) {
	mhi, mlo := bits.Mul64(x, y)
	hhi, _ := bits.Mul64(mlo*qInv, q)
	r = mhi - hhi + q
	if r >= q {
		r -= q
	}
	return
}

// MRedConstant computes x * y * (1/2^64) mod q in constant time.
// The result is between 0 and 2*q-1.
func MRedConstant(x, y, q, qInv uint64) (r uint64) {
	ahi, alo := bits.Mul64(x, y)
	H, _ := bits.Mul64(alo*qInv, q)
	r = ahi - H + q
	return
}

// FastBRedOperand is a struct storing y and floor((y * 2^64)/q).
type FastBRedOperand struct {
	Operand  uint64
	Quotient uint64
}

// NewFastBRedOperand returns a new FastBRedOperand storing y and floor((y * 2^64)/q).
func NewFastBRedOperand(x, q uint64) (t FastBRedOperand) {
	return FastBRedOperand{x, FastBRedQuotient(x, q)}
}

// FastBRedQuotient returns floor((y * 2^64)/q).
func FastBRedQuotient(y, q uint64) (quotient uint64) {
	quotient, _ = bits.Div64(y, 0, q)
	return quotient
}

// FastBRed returns (x*y) % q, y must be smaller than q.
// y is an operand holding y and floor((y * 2^64)/q).
func FastBRed(x uint64, y FastBRedOperand, q uint64) (r uint64) {
	hw, _ := bits.Mul64(x, y.Quotient)
	r = y.Operand*x - hw*q
	if r >= q {
		r -= q
	}
	return
}

// FastBRedConstant returns (x*y) % q, y must be smaller than q.
// y is an operand holding y and floor((y * 2^64)/q).
// The result is between 0 and 2*q-1.
func FastBRedConstant(x uint64, y FastBRedOperand, q uint64) uint64 {
	hw, _ := bits.Mul64(x, y.Quotient)
	return y.Operand*x - hw*q
}

// BRedParams computes the parameters for the BRed algorithm.
// Returns ((2^128)/q)/(2^64) and (2^128)/q mod 2^64.
func BRedParams(q uint64) (params []uint64) {
	bigR := new(big.Int).Lsh(NewUint(1), 128)
	bigR.Quo(bigR, NewUint(q))

	mhi := new(big.Int).Rsh(bigR, 64).Uint64()
	mlo := bigR.Uint64()

	return []uint64{mhi, mlo}
}

// BRedAdd computes a mod q.
func BRedAdd(a, q uint64, u []uint64) (r uint64) {
	mhi, _ := bits.Mul64(a, u[0])
	r = a - mhi*q
	if r >= q {
		r -= q
	}
	return
}

// BRedAddConstant computes a mod q in constant time.
// The result is between 0 and 2*q-1.
func BRedAddConstant(x, q uint64, u []uint64) uint64 {
	s0, _ := bits.Mul64(x, u[0])
	return x - s0*q
}

// BRed computes x*y mod q.
func BRed(x, y, q uint64, u []uint64) (r uint64) {

	var mhi, mlo, lhi, hhi, hlo, s0, carry uint64

	mhi, mlo = bits.Mul64(x, y)

	// computes r = mhi * uhi + (mlo * uhi + mhi * ulo)<<64 + (mlo * ulo)) >> 128

	r = mhi * u[0] // r = mhi * uhi

	hhi, hlo = bits.Mul64(mlo, u[0]) // mlo * uhi

	r += hhi

	lhi, _ = bits.Mul64(mlo, u[1]) // mlo * ulo

	s0, carry = bits.Add64(hlo, lhi, 0)

	r += carry

	hhi, hlo = bits.Mul64(mhi, u[1]) // mhi * ulo

	r += hhi

	_, carry = bits.Add64(hlo, s0, 0)

	r += carry

	r = mlo - r*q

	if r >= q {
		r -= q
	}

	return
}

// BRedConstant computes x*y mod q in constant time.
// The result is between 0 and 2*q-1.
func BRedConstant(x, y, q uint64, u []uint64) (r uint64) {

	var mhi, mlo, lhi, hhi, hlo, s0, carry uint64

	mhi, mlo = bits.Mul64(x, y)

	// computes r = mhi * uhi + (mlo * uhi + mhi * ulo)<<64 + (mlo * ulo)) >> 128

	r = mhi * u[0] // r = mhi * uhi

	hhi, hlo = bits.Mul64(mlo, u[0]) // mlo * uhi

	r += hhi

	lhi, _ = bits.Mul64(mlo, u[1]) // mlo * ulo

	s0, carry = bits.Add64(hlo, lhi, 0)

	r += carry

	hhi, hlo = bits.Mul64(mhi, u[1]) // mhi * ulo

	r += hhi

	_, carry = bits.Add64(hlo, s0, 0)

	r += carry

	r = mlo - r*q

	return
}

// CRed reduce returns a mod q where a is between 0 and 2*q-1.
func CRed(a, q uint64) uint64 {
	if a >= q {
		return a - q
	}
	return a
}
