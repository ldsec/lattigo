//Package rckks implements a RNS-accelerated version of the Approximate Homomorphic Encryption over the Conjugate-invariant Ring
//(Real-HEAAN) scheme. It provides approximate arithmetic over the real numbers.
package rckks

import (
	"math"
	"math/big"

	"github.com/ldsec/lattigo/v2/ring"
)

// GaloisGen is an integer of order N/2 modulo M and that spans Z_M with the integer -1.
// The j-th ring automorphism takes the root zeta to zeta^(5j).
const GaloisGen uint64 = 5

var pi = "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491298336733624406566430860213949463952247371907021798609437027705392171762931767523846748184676694051320005681271452635608277857713427577896091736371787214684409012249534301465495853710507922796892589235420199561121290219608640344181598136297747713099605187072113499999983729780499510597317328160963185950244594553469083026425223082533446850352619311881710100031378387528865875332083814206171776691473035982534904287554687311595628638823537875937519577818577805321712268066130019278766111959092164201989"

// Encoder is an interface implenting the encoding algorithms.
type Encoder interface {
	Encode(plaintext *Plaintext, values []float64, slots uint64)
	EncodeNew(values []float64, slots uint64) (plaintext *Plaintext)
	EncodeNTT(plaintext *Plaintext, values []float64, slots uint64)
	EncodeNTTNew(values []float64, slots uint64) (plaintext *Plaintext)
	Decode(plaintext *Plaintext, slots uint64) (res []float64)
	EncodeCoeffs(values []float64, plaintext *Plaintext)
	DecodeCoeffs(plaintext *Plaintext) (res []float64)
}

// EncoderBigComplex is an interface implenting the encoding algorithms with arbitrary precision.
type EncoderBigComplex interface {
	Encode(plaintext *Plaintext, values []*ring.Complex, slots uint64)
	EncodeNew(values []*ring.Complex, slots uint64) (plaintext *Plaintext)
	EncodeNTT(plaintext *Plaintext, values []*ring.Complex, slots uint64)
	EncodeNTTNew(values []*ring.Complex, slots uint64) (plaintext *Plaintext)
	Decode(plaintext *Plaintext, slots uint64) (res []*ring.Complex)
	FFT(values []*ring.Complex, N uint64)
	InvFFT(values []*ring.Complex, N uint64)

	//EncodeCoeffs(values []*big.Float, plaintext *Plaintext)
	//DecodeCoeffs(plaintext *Plaintext) (res []*big.Float)
}

// encoder is a struct storing the necessary parameters to encode a slice of complex number on a Plaintext.
type encoder struct {
	params       *Parameters
	ringQ        *ring.Ring
	bigintChain  []*big.Int
	bigintCoeffs []*big.Int
	qHalf        *big.Int
	polypool     *ring.Poly
	m            uint64
	rotGroup     []uint64
}

type encoderComplex128 struct {
	encoder
	values      []complex128
	valuesfloat []float64
	roots       []complex128
}

func newEncoder(params *Parameters) encoder {

	m := 4 * params.N()

	var q *ring.Ring
	var err error
	if q, err = ring.NewRingWithNthRoot(params.N(), params.N()<<2, params.qi); err != nil {
		panic(err)
	}

	rotGroup := make([]uint64, m>>1)
	fivePows := uint64(1)
	for i := uint64(0); i < m>>2; i++ {
		rotGroup[i] = fivePows
		fivePows *= GaloisGen
		fivePows &= (m - 1)
	}

	return encoder{
		params:       params.Copy(),
		ringQ:        q,
		bigintChain:  genBigIntChain(params.qi),
		bigintCoeffs: make([]*big.Int, m>>1),
		qHalf:        ring.NewUint(0),
		polypool:     q.NewPoly(),
		m:            m,
		rotGroup:     rotGroup,
	}
}

// NewEncoder creates a new Encoder that is used to encode a slice of complex values of size at most N/2 (the number of slots) on a Plaintext.
func NewEncoder(params *Parameters) Encoder {

	encoder := newEncoder(params)

	var angle float64
	roots := make([]complex128, encoder.m+1)
	for i := uint64(0); i < encoder.m; i++ {
		angle = 2 * 3.141592653589793 * float64(i) / float64(encoder.m)

		roots[i] = complex(math.Cos(angle), math.Sin(angle))
	}
	roots[encoder.m] = roots[0]

	return &encoderComplex128{
		encoder:     encoder,
		roots:       roots,
		values:      make([]complex128, encoder.m>>2),
		valuesfloat: make([]float64, encoder.m>>2),
	}
}

func (encoder *encoderComplex128) EncodeNew(values []float64, slots uint64) (plaintext *Plaintext) {
	plaintext = NewPlaintext(encoder.params, encoder.params.MaxLevel(), encoder.params.scale)
	encoder.Encode(plaintext, values, slots)
	return
}

func (encoder *encoderComplex128) embed(values []float64, slots uint64) {

	if uint64(len(values)) > encoder.params.N() || uint64(len(values)) > slots {
		panic("cannot Encode: too many values for the given number of slots")
	}

	if slots == 0 && slots&(slots-1) != 0 {
		panic("cannot Encode: slots must be a power of two between 1 and N/2")
	}

	for i := range values {
		encoder.values[i] = complex(values[i], 0)
	}

	encoder.invfft(encoder.values, slots)

	gap := encoder.ringQ.N / slots

	for i, jdx, idx := uint64(0), encoder.ringQ.N, uint64(0); i < slots; i, jdx, idx = i+1, jdx+gap, idx+gap {
		encoder.valuesfloat[idx] = real(encoder.values[i])
	}
}

func (encoder *encoderComplex128) scaleUp(pol *ring.Poly, scale float64, moduli []uint64) {
	scaleUpVecExact(encoder.valuesfloat, scale, moduli, pol.Coeffs)
}

func (encoder *encoderComplex128) wipeInternalMemory() {
	for i := range encoder.values {
		encoder.values[i] = 0
	}

	for i := range encoder.valuesfloat {
		encoder.valuesfloat[i] = 0
	}
}

// Encode takes a slice of complex128 values of size at most N/2 (the number of slots) and encodes it in the receiver Plaintext.
func (encoder *encoderComplex128) Encode(plaintext *Plaintext, values []float64, slots uint64) {
	encoder.embed(values, slots)
	encoder.scaleUp(plaintext.value, plaintext.scale, encoder.ringQ.Modulus[:plaintext.Level()+1])
	encoder.wipeInternalMemory()
	plaintext.isNTT = false
}

func (encoder *encoderComplex128) EncodeNTTNew(values []float64, slots uint64) (plaintext *Plaintext) {
	plaintext = NewPlaintext(encoder.params, encoder.params.MaxLevel(), encoder.params.scale)
	encoder.EncodeNTT(plaintext, values, slots)
	return
}

func (encoder *encoderComplex128) EncodeNTT(plaintext *Plaintext, values []float64, slots uint64) {
	encoder.Encode(plaintext, values, slots)
	NTTRCKKSLvl(encoder.ringQ, plaintext.Level(), plaintext.value, plaintext.value)
	plaintext.isNTT = true
}

// EncodeCoefficients takes as input a polynomial a0 + a1x + a2x^2 + ... + an-1x^n-1 with float coefficient
// and returns a scaled integer plaintext polynomial in NTT.
func (encoder *encoderComplex128) EncodeCoeffs(values []float64, plaintext *Plaintext) {

	if uint64(len(values)) > encoder.params.N() {
		panic("cannot EncodeCoeffs : too many values (maximum is N)")
	}

	scaleUpVecExact(values, plaintext.scale, encoder.ringQ.Modulus[:plaintext.Level()+1], plaintext.value.Coeffs)

	plaintext.isNTT = false
}

// EncodeCoefficients takes as input a polynomial a0 + a1x + a2x^2 + ... + an-1x^n-1 with float coefficient
// and returns a scaled integer plaintext polynomial in NTT.
func (encoder *encoderComplex128) EncodeCoeffsNTT(values []float64, plaintext *Plaintext) {
	encoder.EncodeCoeffs(values, plaintext)
	NTTRCKKSLvl(encoder.ringQ, plaintext.Level(), plaintext.value, plaintext.value)
	plaintext.isNTT = true
}

// DecodeCoeffs takes as input a plaintext and returns the scaled down coefficient of the plaintext in flaot64
func (encoder *encoderComplex128) DecodeCoeffs(plaintext *Plaintext) (res []float64) {

	if plaintext.isNTT {
		InvNTTRCKKSLvl(encoder.ringQ, plaintext.Level(), plaintext.value, encoder.polypool)
	} else {
		encoder.ringQ.CopyLvl(plaintext.Level(), plaintext.value, encoder.polypool)
	}

	res = make([]float64, encoder.params.N())

	// We have more than one moduli and need the CRT reconstruction
	if plaintext.Level() > 0 {

		encoder.ringQ.PolyToBigint(encoder.polypool, encoder.bigintCoeffs)

		Q := encoder.bigintChain[plaintext.Level()]

		encoder.qHalf.Set(Q)
		encoder.qHalf.Rsh(encoder.qHalf, 1)

		var sign int

		for i := range res {

			// Centers the value around the current modulus
			encoder.bigintCoeffs[i].Mod(encoder.bigintCoeffs[i], Q)

			sign = encoder.bigintCoeffs[i].Cmp(encoder.qHalf)
			if sign == 1 || sign == 0 {
				encoder.bigintCoeffs[i].Sub(encoder.bigintCoeffs[i], Q)
			}

			res[i] = scaleDown(encoder.bigintCoeffs[i], plaintext.scale)
		}
		// We can directly get the coefficients
	} else {

		Q := encoder.ringQ.Modulus[0]
		coeffs := encoder.polypool.Coeffs[0]

		for i := range res {

			if coeffs[i] >= Q>>1 {
				res[i] = -float64(Q - coeffs[i])
			} else {
				res[i] = float64(coeffs[i])
			}

			res[i] /= plaintext.scale
		}
	}

	return
}

// Decode decodes the Plaintext values to a slice of complex128 values of size at most N/2.
func (encoder *encoderComplex128) Decode(plaintext *Plaintext, slots uint64) (res []float64) {

	if plaintext.isNTT {
		InvNTTRCKKSLvl(encoder.ringQ, plaintext.Level(), plaintext.value, encoder.polypool)
	} else {
		encoder.ringQ.CopyLvl(plaintext.Level(), plaintext.value, encoder.polypool)
	}

	maxSlots := encoder.ringQ.N
	gap := maxSlots / slots

	// We have more than one moduli and need the CRT reconstruction
	if plaintext.Level() > 0 {

		encoder.ringQ.PolyToBigint(encoder.polypool, encoder.bigintCoeffs)

		Q := encoder.bigintChain[plaintext.Level()]

		encoder.qHalf.Set(Q)
		encoder.qHalf.Rsh(encoder.qHalf, 1)

		var sign int

		for i, idx := uint64(0), uint64(0); i < slots; i, idx = i+1, idx+gap {

			// Centers the value around the current modulus
			encoder.bigintCoeffs[idx].Mod(encoder.bigintCoeffs[idx], Q)
			sign = encoder.bigintCoeffs[idx].Cmp(encoder.qHalf)
			if sign == 1 || sign == 0 {
				encoder.bigintCoeffs[idx].Sub(encoder.bigintCoeffs[idx], Q)
			}

			encoder.values[i] = complex(scaleDown(encoder.bigintCoeffs[idx], plaintext.scale), 0)
		}
		// We can directly get the coefficients
	} else {

		Q := encoder.ringQ.Modulus[0]
		coeffs := encoder.polypool.Coeffs[0]

		var real float64
		for i, idx := uint64(0), uint64(0); i < slots; i, idx = i+1, idx+gap {

			if coeffs[idx] >= Q>>1 {
				real = -float64(Q - coeffs[idx])
			} else {
				real = float64(coeffs[idx])
			}

			encoder.values[i] = complex(real, 0) / complex(plaintext.scale, 0)
		}
	}

	for i := uint64(1); i < slots; i++ {
		encoder.values[i] -= complex(0, real(encoder.values[slots-i]))
	}

	encoder.fft(encoder.values, slots)

	res = make([]float64, slots)

	for i := range res {
		res[i] = real(encoder.values[i])
	}

	for i := uint64(0); i < encoder.ringQ.N>>1; i++ {
		encoder.values[i] = 0
	}

	return
}

func (encoder *encoderComplex128) invfft(values []complex128, N uint64) {

	var lenh, lenq, gap, idx uint64
	var u, v complex128

	for len := N; len >= 1; len >>= 1 {
		for i := uint64(0); i < N; i += len {
			lenh = len >> 1
			lenq = len << 2
			gap = encoder.m / lenq
			for j := uint64(0); j < lenh; j++ {
				idx = (lenq - (encoder.rotGroup[j] % lenq)) * gap
				u = values[i+j] + values[i+j+lenh]
				v = values[i+j] - values[i+j+lenh]
				v *= encoder.roots[idx]
				values[i+j] = u
				values[i+j+lenh] = v

			}
		}
	}

	for i := uint64(0); i < N; i++ {
		values[i] /= complex(float64(N), 0)
	}

	sliceBitReverseInPlaceComplex128(values, N)
}

func (encoder *encoderComplex128) fft(values []complex128, N uint64) {

	var lenh, lenq, gap, idx uint64
	var u, v complex128

	sliceBitReverseInPlaceComplex128(values, N)

	for len := uint64(2); len <= N; len <<= 1 {
		for i := uint64(0); i < N; i += len {
			lenh = len >> 1
			lenq = len << 2
			gap = encoder.m / lenq
			for j := uint64(0); j < lenh; j++ {
				idx = (encoder.rotGroup[j] % lenq) * gap
				u = values[i+j]
				v = values[i+j+lenh]
				v *= encoder.roots[idx]
				values[i+j] = u + v
				values[i+j+lenh] = u - v
			}
		}
	}
}

type encoderBigComplex struct {
	encoder
	zero         *big.Float
	cMul         *ring.ComplexMultiplier
	logPrecision uint64
	values       []*ring.Complex
	valuesfloat  []*big.Float
	roots        []*ring.Complex
}

// NewEncoderBigComplex creates a new encoder using arbitrary precision complex arithmetic
func NewEncoderBigComplex(params *Parameters, logPrecision uint64) EncoderBigComplex {
	encoder := newEncoder(params)

	var PI = new(big.Float)
	PI.SetPrec(uint(logPrecision))
	PI.SetString(pi)

	var PIHalf = new(big.Float)
	PIHalf.SetPrec(uint(logPrecision))
	PIHalf.SetString(pi)
	PIHalf.Quo(PIHalf, ring.NewFloat(2, logPrecision))

	var angle *big.Float
	roots := make([]*ring.Complex, encoder.m+1)
	for i := uint64(0); i < encoder.m; i++ {
		angle = ring.NewFloat(2, logPrecision)
		angle.Mul(angle, PI)
		angle.Mul(angle, ring.NewFloat(float64(i), logPrecision))
		angle.Quo(angle, ring.NewFloat(float64(encoder.m), logPrecision))

		real := ring.Cos(angle)
		angle.Sub(PIHalf, angle)
		imag := ring.Cos(angle)

		roots[i] = ring.NewComplex(real, imag)
	}

	roots[encoder.m] = roots[0].Copy()

	values := make([]*ring.Complex, encoder.m>>2)
	valuesfloat := make([]*big.Float, encoder.m>>2)

	for i := uint64(0); i < encoder.m>>2; i++ {

		values[i] = ring.NewComplex(ring.NewFloat(0, logPrecision), ring.NewFloat(0, logPrecision))
		valuesfloat[i] = ring.NewFloat(0, logPrecision)
	}

	return &encoderBigComplex{
		encoder:      encoder,
		zero:         ring.NewFloat(0, logPrecision),
		cMul:         ring.NewComplexMultiplier(),
		logPrecision: logPrecision,
		roots:        roots,
		values:       values,
		valuesfloat:  valuesfloat,
	}
}

func (encoder *encoderBigComplex) EncodeNew(values []*ring.Complex, slots uint64) (plaintext *Plaintext) {
	plaintext = NewPlaintext(encoder.params, encoder.params.MaxLevel(), encoder.params.scale)
	encoder.Encode(plaintext, values, slots)
	return
}

func (encoder *encoderBigComplex) Encode(plaintext *Plaintext, values []*ring.Complex, slots uint64) {

	if uint64(len(values)) > encoder.ringQ.N>>1 || uint64(len(values)) > slots {
		panic("cannot Encode: too many values for the given number of slots")
	}

	if slots == 0 && slots&(slots-1) == 0 {
		panic("cannot Encode: slots must be a power of two between 1 and N/2")
	}

	if uint64(len(values)) != slots {
		panic("cannot Encode: number of values must be equal to slots")
	}

	for i := uint64(0); i < slots; i++ {
		encoder.values[i].Set(values[i])
	}

	encoder.InvFFT(encoder.values, slots)

	gap := encoder.ringQ.N / slots

	for i, jdx, idx := uint64(0), encoder.ringQ.N, uint64(0); i < slots; i, jdx, idx = i+1, jdx+gap, idx+gap {
		encoder.valuesfloat[idx].Set(encoder.values[i].Real())
	}

	scaleUpVecExactBigFloat(encoder.valuesfloat, plaintext.scale, encoder.ringQ.Modulus[:plaintext.Level()+1], plaintext.value.Coeffs)

	coeffsBigInt := make([]*big.Int, encoder.params.N())

	encoder.ringQ.PolyToBigint(plaintext.value, coeffsBigInt)

	for i := range encoder.values {
		encoder.values[i].Real().Set(encoder.zero)
		encoder.values[i].Imag().Set(encoder.zero)
	}

	for i := range encoder.valuesfloat {
		encoder.valuesfloat[i].Set(encoder.zero)
	}
}

func (encoder *encoderBigComplex) EncodeNTTNew(values []*ring.Complex, slots uint64) (plaintext *Plaintext) {
	plaintext = NewPlaintext(encoder.params, encoder.params.MaxLevel(), encoder.params.scale)
	encoder.EncodeNTT(plaintext, values, slots)
	return
}

func (encoder *encoderBigComplex) EncodeNTT(plaintext *Plaintext, values []*ring.Complex, slots uint64) {

	encoder.Encode(plaintext, values, slots)

	encoder.ringQ.NTTLvl(plaintext.Level(), plaintext.value, plaintext.value)

	plaintext.isNTT = true
}

func (encoder *encoderBigComplex) Decode(plaintext *Plaintext, slots uint64) (res []*ring.Complex) {

	encoder.ringQ.InvNTTLvl(plaintext.Level(), plaintext.value, encoder.polypool)
	encoder.ringQ.PolyToBigint(encoder.polypool, encoder.bigintCoeffs)

	Q := encoder.bigintChain[plaintext.Level()]

	maxSlots := encoder.ringQ.N

	scaleFlo := ring.NewFloat(plaintext.Scale(), encoder.logPrecision)

	encoder.qHalf.Set(Q)
	encoder.qHalf.Rsh(encoder.qHalf, 1)

	gap := maxSlots / slots

	var sign int

	for i, idx := uint64(0), uint64(0); i < slots; i, idx = i+1, idx+gap {

		// Centers the value around the current modulus
		encoder.bigintCoeffs[idx].Mod(encoder.bigintCoeffs[idx], Q)
		sign = encoder.bigintCoeffs[idx].Cmp(encoder.qHalf)
		if sign == 1 || sign == 0 {
			encoder.bigintCoeffs[idx].Sub(encoder.bigintCoeffs[idx], Q)
		}

		// Centers the value around the current modulus
		encoder.bigintCoeffs[idx+maxSlots].Mod(encoder.bigintCoeffs[idx+maxSlots], Q)
		sign = encoder.bigintCoeffs[idx+maxSlots].Cmp(encoder.qHalf)
		if sign == 1 || sign == 0 {
			encoder.bigintCoeffs[idx+maxSlots].Sub(encoder.bigintCoeffs[idx+maxSlots], Q)
		}

		encoder.values[i].Real().SetInt(encoder.bigintCoeffs[idx])
		encoder.values[i].Real().Quo(encoder.values[i].Real(), scaleFlo)
	}

	encoder.FFT(encoder.values, slots)

	res = make([]*ring.Complex, slots)

	for i := range res {
		res[i] = encoder.values[i].Copy()
	}

	for i := range encoder.values {
		encoder.values[i].Real().Set(encoder.zero)
		encoder.values[i].Imag().Set(encoder.zero)
	}

	return
}

func (encoder *encoderBigComplex) InvFFT(values []*ring.Complex, N uint64) {

	var lenh, lenq, gap, idx uint64
	u := ring.NewComplex(nil, nil)
	v := ring.NewComplex(nil, nil)

	for len := N; len >= 1; len >>= 1 {
		for i := uint64(0); i < N; i += len {
			lenh = len >> 1
			lenq = len << 2
			gap = encoder.m / lenq
			for j := uint64(0); j < lenh; j++ {
				idx = (lenq - (encoder.rotGroup[j] % lenq)) * gap
				u.Add(values[i+j], values[i+j+lenh])
				v.Sub(values[i+j], values[i+j+lenh])
				encoder.cMul.Mul(v, encoder.roots[idx], v)
				values[i+j].Set(u)
				values[i+j+lenh].Set(v)
			}
		}
	}

	NBig := ring.NewFloat(float64(N), encoder.logPrecision)
	for i := range values {
		values[i][0].Quo(values[i][0], NBig)
		values[i][1].Quo(values[i][1], NBig)
	}

	sliceBitReverseInPlaceRingComplex(values, N)
}

func (encoder *encoderBigComplex) FFT(values []*ring.Complex, N uint64) {

	var lenh, lenq, gap, idx uint64

	u := ring.NewComplex(nil, nil)
	v := ring.NewComplex(nil, nil)

	sliceBitReverseInPlaceRingComplex(values, N)

	for len := uint64(2); len <= N; len <<= 1 {
		for i := uint64(0); i < N; i += len {
			lenh = len >> 1
			lenq = len << 2
			gap = encoder.m / lenq
			for j := uint64(0); j < lenh; j++ {
				idx = (encoder.rotGroup[j] % lenq) * gap
				u.Set(values[i+j])
				v.Set(values[i+j+lenh])
				encoder.cMul.Mul(v, encoder.roots[idx], v)
				values[i+j].Add(u, v)
				values[i+j+lenh].Sub(u, v)
			}
		}
	}
}
