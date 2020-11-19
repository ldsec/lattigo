//Package ckks implements a RNS-accelerated version of the Homomorphic Encryption for Arithmetic for Approximate Numbers
//(HEAAN, a.k.a. CKKS) scheme. It provides approximate arithmetic over the complex numbers.package ckks
package ckks

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
	Encode(plaintext *Plaintext, values []complex128, logSlots uint64)
	EncodeNew(values []complex128, logSlots uint64) (plaintext *Plaintext)
	EncodeAtLvlNew(level uint64, values []complex128, logSlots uint64) (plaintext *Plaintext)
	EncodeNTT(plaintext *Plaintext, values []complex128, logSlots uint64)
	EncodeNTTAtLvlNew(level uint64, values []complex128, logSlots uint64) (plaintext *Plaintext)
	EncodeDiagMatrixAtLvl(level uint64, vector map[uint64][]complex128, scale, maxM1N2Ratio float64, logSlots uint64) (matrix *PtDiagMatrix)
	Decode(plaintext *Plaintext, logSlots uint64) (res []complex128)
	Embed(values []complex128, logSlots uint64)
	ScaleUp(pol *ring.Poly, scale float64, moduli []uint64)
	WipeInternalMemory()
	EncodeCoeffs(values []float64, plaintext *Plaintext)
	DecodeCoeffs(plaintext *Plaintext) (res []float64)
}

// EncoderBigComplex is an interface implenting the encoding algorithms with arbitrary precision.
type EncoderBigComplex interface {
	Encode(plaintext *Plaintext, values []*ring.Complex, logSlots uint64)
	EncodeNew(values []*ring.Complex, logSlots uint64) (plaintext *Plaintext)
	EncodeAtLvlNew(level uint64, values []*ring.Complex, logSlots uint64) (plaintext *Plaintext)
	EncodeNTT(plaintext *Plaintext, values []*ring.Complex, logSlots uint64)
	EncodeNTTAtLvlNew(level uint64, values []*ring.Complex, logSlots uint64) (plaintext *Plaintext)
	Decode(plaintext *Plaintext, logSlots uint64) (res []*ring.Complex)
	FFT(values []*ring.Complex, N uint64)
	InvFFT(values []*ring.Complex, N uint64)

	//EncodeCoeffs(values []*big.Float, plaintext *Plaintext)
	//DecodeCoeffs(plaintext *Plaintext) (res []*big.Float)
}

// encoder is a struct storing the necessary parameters to encode a slice of complex number on a Plaintext.
type encoder struct {
	params       *Parameters
	ringQ        *ring.Ring
	ringP        *ring.Ring
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

	m := 2 * params.N()

	var q *ring.Ring
	var err error
	if q, err = ring.NewRing(params.N(), params.qi); err != nil {
		panic(err)
	}

	var p *ring.Ring
	if params.PiCount() != 0 {
		if p, err = ring.NewRing(params.N(), params.pi); err != nil {
			panic(err)
		}
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
		ringP:        p,
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
		valuesfloat: make([]float64, encoder.m>>1),
	}
}

// EncodeNew encodes a slice of complex128 of length slots = 2^{logSlots} on new plaintext at the maximum level.
func (encoder *encoderComplex128) EncodeNew(values []complex128, logSlots uint64) (plaintext *Plaintext) {
	return encoder.EncodeAtLvlNew(encoder.params.MaxLevel(), values, logSlots)
}

// EncodeAtLvlNew encodes a slice of complex128 of length slots = 2^{logSlots} on new plaintext at the desired level.
func (encoder *encoderComplex128) EncodeAtLvlNew(level uint64, values []complex128, logSlots uint64) (plaintext *Plaintext) {
	plaintext = NewPlaintext(encoder.params, level, encoder.params.scale)
	encoder.Encode(plaintext, values, logSlots)
	return
}

// Encode encodes a slice of complex128 of length slots = 2^{logSlots} on the input plaintext.
func (encoder *encoderComplex128) Encode(plaintext *Plaintext, values []complex128, logSlots uint64) {
	encoder.Embed(values, logSlots)
	encoder.ScaleUp(plaintext.value, plaintext.scale, encoder.ringQ.Modulus[:plaintext.Level()+1])
	encoder.WipeInternalMemory()
	plaintext.isNTT = false
}

// EncodeNTTNew encodes a slice of complex128 of length slots = 2^{logSlots} on new plaintext at the maximum level.
// Returns a plaintext in the NTT domain.
func (encoder *encoderComplex128) EncodeNTTNew(values []complex128, logSlots uint64) (plaintext *Plaintext) {
	return encoder.EncodeNTTAtLvlNew(encoder.params.MaxLevel(), values, logSlots)
}

// EncodeNTTAtLvlNew encodes a slice of complex128 of length slots = 2^{logSlots} on new plaintext at the desired level.
// Returns a plaintext in the NTT domain.
func (encoder *encoderComplex128) EncodeNTTAtLvlNew(level uint64, values []complex128, logSlots uint64) (plaintext *Plaintext) {
	plaintext = NewPlaintext(encoder.params, encoder.params.MaxLevel(), encoder.params.scale)
	encoder.EncodeNTT(plaintext, values, logSlots)
	return
}

// EncodeNTT encodes a slice of complex128 of length slots = 2^{logSlots} on the input plaintext.
// Returns a plaintext in the NTT domain.
func (encoder *encoderComplex128) EncodeNTT(plaintext *Plaintext, values []complex128, logSlots uint64) {
	encoder.Encode(plaintext, values, logSlots)
	encoder.ringQ.NTTLvl(plaintext.Level(), plaintext.value, plaintext.value)
	plaintext.isNTT = true
}

// PtDiagMatrix is a struct storing a plaintext diagonalized matrix
// ready to be evaluated on a ciphertext using evaluator.MultiplyByDiagMatrice.
type PtDiagMatrix struct {
	LogSlots uint64                   // Log of the number of slots of the plaintext (needed to compute the appropriate rotation keys)
	N1       uint64                   // N1 is the number of inner loops of the baby-step giant-step algo used in the evaluation.
	Level    uint64                   // Level is the level at which the matrix is encoded (can be circuit dependant)
	Scale    float64                  // Scale is the scale at which the matrix is encoded (can be circuit dependant)
	Vec      map[uint64][2]*ring.Poly // Vec is the matrix, in diagonal form, where each entry of vec is an indexed non zero diagonal.
}

// EncodeDiagMatrixAtLvl encodes a diagonalized plaintext matrix into PtDiagMatrix struct.
// It can then be evaluated on a ciphertext using evaluator.MultiplyByDiagMatrice.
// maxM1N2Ratio is the maximum ratio between the inner and outer loop of the baby-step giant-step algorithm used in evaluator.MultiplyByDiagMatrice.
// Optimal maxM1N2Ratio value is between 4 and 16 depending on the sparsity of the matrix.
func (encoder *encoderComplex128) EncodeDiagMatrixAtLvl(level uint64, vector map[uint64][]complex128, scale, maxM1N2Ratio float64, logSlots uint64) (matrix *PtDiagMatrix) {

	matrix = new(PtDiagMatrix)
	matrix.LogSlots = logSlots
	matrix.N1 = findbestbabygiantstepsplit(vector, 1<<logSlots, maxM1N2Ratio)

	var N, N1 uint64

	ringQ := encoder.ringQ
	ringP := encoder.ringP

	// N1*N2 = N
	N = ringQ.N
	N1 = matrix.N1

	if len(vector) > 2 {

		index := make(map[uint64][]uint64)
		for key := range vector {
			idx1 := key / N1
			idx2 := key & (N1 - 1)
			if index[idx1] == nil {
				index[idx1] = []uint64{idx2}
			} else {
				index[idx1] = append(index[idx1], idx2)
			}
		}

		matrix.Vec = make(map[uint64][2]*ring.Poly)

		matrix.Level = level
		matrix.Scale = scale

		for j := range index {

			for _, i := range index[j] {

				encoder.Embed(rotate(vector[N1*j+uint64(i)], (N>>1)-(N1*j)), logSlots)

				plaintextQ := ring.NewPoly(N, level+1)
				encoder.ScaleUp(plaintextQ, scale, ringQ.Modulus[:level+1])
				ringQ.NTTLvl(level, plaintextQ, plaintextQ)
				ringQ.MFormLvl(level, plaintextQ, plaintextQ)

				plaintextP := ring.NewPoly(N, level+1)
				encoder.ScaleUp(plaintextP, scale, ringP.Modulus)
				ringP.NTT(plaintextP, plaintextP)
				ringP.MForm(plaintextP, plaintextP)

				matrix.Vec[N1*j+uint64(i)] = [2]*ring.Poly{plaintextQ, plaintextP}

				encoder.WipeInternalMemory()
			}
		}
	} else {

		matrix.Vec = make(map[uint64][2]*ring.Poly)

		matrix.Level = level
		matrix.Scale = scale

		for i := range vector {

			encoder.Embed(rotate(vector[i], i), logSlots)

			plaintextQ := ring.NewPoly(N, level+1)
			encoder.ScaleUp(plaintextQ, scale, ringQ.Modulus[:level+1])
			ringQ.NTTLvl(level, plaintextQ, plaintextQ)
			ringQ.MFormLvl(level, plaintextQ, plaintextQ)

			plaintextP := ring.NewPoly(N, level+1)
			encoder.ScaleUp(plaintextP, scale, ringP.Modulus)
			ringP.NTT(plaintextP, plaintextP)
			ringP.MForm(plaintextP, plaintextP)

			matrix.Vec[i] = [2]*ring.Poly{plaintextQ, plaintextP}

			encoder.WipeInternalMemory()
		}
	}

	return
}

// Finds the best N1*N2 = N for the baby-step giant-step algorithm for matrix multiplication.
func findbestbabygiantstepsplit(vector map[uint64][]complex128, maxN uint64, maxRatio float64) (minN uint64) {

	for N1 := uint64(1); N1 < maxN; N1 <<= 1 {

		index := make(map[uint64][]uint64)

		for key := range vector {

			idx1 := key / N1
			idx2 := key & (N1 - 1)

			if index[idx1] == nil {
				index[idx1] = []uint64{idx2}
			} else {
				index[idx1] = append(index[idx1], idx2)
			}
		}

		if len(index[0]) > 0 {

			hoisted := len(index[0]) - 1
			normal := len(index) - 1

			// The matrice is very sparse already
			if normal == 0 {
				return N1 / 2
			}

			if hoisted > normal {
				// Finds the next split that has a ratio hoisted/normal greater or equal to maxRatio
				for float64(hoisted)/float64(normal) < maxRatio {

					if normal/2 == 0 {
						break
					}
					N1 *= 2
					hoisted = hoisted*2 + 1
					normal = normal / 2
				}
				return N1
			}
		}
	}

	return 1
}

// Embed encodes a vector and stores internaly the encoded values.
// To be used in conjonction with ScaleUp.
func (encoder *encoderComplex128) Embed(values []complex128, logSlots uint64) {

	slots := uint64(1 << logSlots)

	if uint64(len(values)) > encoder.params.N()/2 || uint64(len(values)) > slots || logSlots > encoder.params.LogN()-1 {
		panic("cannot Encode: too many values/slots for the given ring degree")
	}

	for i := range values {
		encoder.values[i] = values[i]
	}

	encoder.invfft(encoder.values, slots)

	gap := (encoder.ringQ.N >> 1) / slots

	for i, jdx, idx := uint64(0), encoder.ringQ.N>>1, uint64(0); i < slots; i, jdx, idx = i+1, jdx+gap, idx+gap {
		encoder.valuesfloat[idx] = real(encoder.values[i])
		encoder.valuesfloat[jdx] = imag(encoder.values[i])
	}
}

// ScaleUp writes the internaly stored encoded values on a polynomial.
func (encoder *encoderComplex128) ScaleUp(pol *ring.Poly, scale float64, moduli []uint64) {
	scaleUpVecExact(encoder.valuesfloat, scale, moduli, pol.Coeffs)
}

// WipeInternalMemory sets the internaly stored encoded values of the encoder to zero.
func (encoder *encoderComplex128) WipeInternalMemory() {
	for i := range encoder.values {
		encoder.values[i] = 0
	}

	for i := range encoder.valuesfloat {
		encoder.valuesfloat[i] = 0
	}
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
	encoder.ringQ.NTTLvl(plaintext.Level(), plaintext.value, plaintext.value)
	plaintext.isNTT = true
}

// DecodeCoeffs takes as input a plaintext and returns the scaled down coefficient of the plaintext in float64.
func (encoder *encoderComplex128) DecodeCoeffs(plaintext *Plaintext) (res []float64) {

	if plaintext.isNTT {
		encoder.ringQ.InvNTTLvl(plaintext.Level(), plaintext.value, encoder.polypool)
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
func (encoder *encoderComplex128) Decode(plaintext *Plaintext, logSlots uint64) (res []complex128) {

	if logSlots > encoder.params.LogN()-1 {
		panic("cannot Decode: too many slots for the given ring degree")
	}

	slots := uint64(1 << logSlots)

	if plaintext.isNTT {
		encoder.ringQ.InvNTTLvl(plaintext.Level(), plaintext.value, encoder.polypool)
	} else {
		encoder.ringQ.CopyLvl(plaintext.Level(), plaintext.value, encoder.polypool)
	}

	maxSlots := encoder.ringQ.N >> 1
	gap := maxSlots / slots

	// We have more than one moduli and need the CRT reconstruction

	if plaintext.Scale() < float64(encoder.ringQ.Modulus[0]) || plaintext.Level() == 0 {

		Q := encoder.ringQ.Modulus[0]
		coeffs := encoder.polypool.Coeffs[0]

		var real, imag float64
		for i, idx := uint64(0), uint64(0); i < slots; i, idx = i+1, idx+gap {

			if coeffs[idx] >= Q>>1 {
				real = -float64(Q - coeffs[idx])
			} else {
				real = float64(coeffs[idx])
			}

			if coeffs[idx+maxSlots] >= Q>>1 {
				imag = -float64(Q - coeffs[idx+maxSlots])
			} else {
				imag = float64(coeffs[idx+maxSlots])
			}

			encoder.values[i] = complex(real, imag) / complex(plaintext.scale, 0)
		}
	} else {

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

			// Centers the value around the current modulus
			encoder.bigintCoeffs[idx+maxSlots].Mod(encoder.bigintCoeffs[idx+maxSlots], Q)
			sign = encoder.bigintCoeffs[idx+maxSlots].Cmp(encoder.qHalf)
			if sign == 1 || sign == 0 {
				encoder.bigintCoeffs[idx+maxSlots].Sub(encoder.bigintCoeffs[idx+maxSlots], Q)
			}

			encoder.values[i] = complex(scaleDown(encoder.bigintCoeffs[idx], plaintext.scale), scaleDown(encoder.bigintCoeffs[idx+maxSlots], plaintext.scale))
		}
		// We can directly get the coefficients
	}

	encoder.fft(encoder.values, slots)

	res = make([]complex128, slots)

	for i := range res {
		res[i] = encoder.values[i]
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

// NewEncoderBigComplex creates a new encoder using arbitrary precision complex arithmetic.
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
	valuesfloat := make([]*big.Float, encoder.m>>1)

	for i := uint64(0); i < encoder.m>>2; i++ {

		values[i] = ring.NewComplex(ring.NewFloat(0, logPrecision), ring.NewFloat(0, logPrecision))
		valuesfloat[i*2] = ring.NewFloat(0, logPrecision)
		valuesfloat[(i*2)+1] = ring.NewFloat(0, logPrecision)
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

// EncodeNew encodes a slice of ring.Complex of length slots = 2^{logSlots} on a new plaintext at the maximum level.
func (encoder *encoderBigComplex) EncodeNew(values []*ring.Complex, logSlots uint64) (plaintext *Plaintext) {
	return encoder.EncodeAtLvlNew(encoder.params.MaxLevel(), values, logSlots)
}

// EncodeAtLvlNew encodes a slice of ring.Complex of length slots = 2^{logSlots} on a new plaintext at the desired level.
func (encoder *encoderBigComplex) EncodeAtLvlNew(level uint64, values []*ring.Complex, logSlots uint64) (plaintext *Plaintext) {
	plaintext = NewPlaintext(encoder.params, level, encoder.params.scale)
	encoder.Encode(plaintext, values, logSlots)
	return
}

// EncodeNTTNew encodes a slice of ring.Complex of length slots = 2^{logSlots} on a plaintext at the maximum level.
// Returns a plaintext in the NTT domain.
func (encoder *encoderBigComplex) EncodeNTTNew(values []*ring.Complex, logSlots uint64) (plaintext *Plaintext) {
	return encoder.EncodeNTTAtLvlNew(encoder.params.MaxLevel(), values, logSlots)
}

// EncodeNTTAtLvlNew encodes a slice of ring.Complex of length slots = 2^{logSlots} on a plaintext at the desired level.
// Returns a plaintext in the NTT domain.
func (encoder *encoderBigComplex) EncodeNTTAtLvlNew(level uint64, values []*ring.Complex, logSlots uint64) (plaintext *Plaintext) {
	plaintext = NewPlaintext(encoder.params, encoder.params.MaxLevel(), encoder.params.scale)
	encoder.EncodeNTT(plaintext, values, logSlots)
	return
}

// Encode encodes a slice of ring.Complex of length slots = 2^{logSlots} on a plaintext at the input plaintext level.
// Returns a plaintext in the NTT domain.
func (encoder *encoderBigComplex) EncodeNTT(plaintext *Plaintext, values []*ring.Complex, logSlots uint64) {
	encoder.Encode(plaintext, values, logSlots)
	encoder.ringQ.NTTLvl(plaintext.Level(), plaintext.value, plaintext.value)
	plaintext.isNTT = true
}

// Encode encodes a slice of ring.Complex of length slots = 2^{logSlots} on a plaintext at the input plaintext level.
func (encoder *encoderBigComplex) Encode(plaintext *Plaintext, values []*ring.Complex, logSlots uint64) {

	slots := uint64(1 << logSlots)

	if uint64(len(values)) > encoder.params.N()/2 || uint64(len(values)) > slots || logSlots > encoder.params.LogN()-1 {
		panic("cannot Encode: too many values/slots for the given ring degree")
	}

	if uint64(len(values)) != slots {
		panic("cannot Encode: number of values must be equal to slots")
	}

	for i := uint64(0); i < slots; i++ {
		encoder.values[i].Set(values[i])
	}

	encoder.InvFFT(encoder.values, slots)

	gap := (encoder.ringQ.N >> 1) / slots

	for i, jdx, idx := uint64(0), (encoder.ringQ.N >> 1), uint64(0); i < slots; i, jdx, idx = i+1, jdx+gap, idx+gap {
		encoder.valuesfloat[idx].Set(encoder.values[i].Real())
		encoder.valuesfloat[jdx].Set(encoder.values[i].Imag())
	}

	scaleUpVecExactBigFloat(encoder.valuesfloat, plaintext.scale, encoder.ringQ.Modulus[:plaintext.Level()+1], plaintext.value.Coeffs)

	coeffsBigInt := make([]*big.Int, encoder.params.N())

	encoder.ringQ.PolyToBigint(plaintext.value, coeffsBigInt)

	for i := uint64(0); i < (encoder.ringQ.N >> 1); i++ {
		encoder.values[i].Real().Set(encoder.zero)
		encoder.values[i].Imag().Set(encoder.zero)
	}

	for i := uint64(0); i < encoder.ringQ.N; i++ {
		encoder.valuesfloat[i].Set(encoder.zero)
	}
}

// Decode decodes the Plaintext values to a slice of complex128 values of size at most N/2.
func (encoder *encoderBigComplex) Decode(plaintext *Plaintext, logSlots uint64) (res []*ring.Complex) {

	slots := uint64(1 << logSlots)

	if logSlots > encoder.params.LogN()-1 {
		panic("cannot Decode: too many slots for the given ring degree")
	}

	encoder.ringQ.InvNTTLvl(plaintext.Level(), plaintext.value, encoder.polypool)
	encoder.ringQ.PolyToBigint(encoder.polypool, encoder.bigintCoeffs)

	Q := encoder.bigintChain[plaintext.Level()]

	maxSlots := encoder.ringQ.N >> 1

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

		encoder.values[i].Imag().SetInt(encoder.bigintCoeffs[idx+maxSlots])
		encoder.values[i].Imag().Quo(encoder.values[i].Imag(), scaleFlo)
	}

	encoder.FFT(encoder.values, slots)

	res = make([]*ring.Complex, slots)

	for i := range res {
		res[i] = encoder.values[i].Copy()
	}

	for i := uint64(0); i < maxSlots; i++ {
		encoder.values[i].Real().Set(encoder.zero)
		encoder.values[i].Imag().Set(encoder.zero)
	}

	return
}

// InvFFT evaluates the encoding matrix on a slice fo ring.Complex values.
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

// FFT evaluates the decoding matrix on a slice fo ring.Complex values.
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
