package ring

import (
	"math/big"
	"math/bits"
	"unsafe"
)

// Scaler is an interface that rescales polynomial coefficients by a fraction t/Q.
type Scaler interface {
	// DivByQOverTRounded returns p1 scaled by a factor t/Q and mod t on the receiver p2.
	DivByQOverTRounded(p1, p2 *Poly)
}

// RNSScaler implements the Scaler interface by performing a scaling by t/Q in the RNS domain.
// This implementation of the Scaler interface is preferred over the SimpleScaler implementation.
type RNSScaler struct {
	ringQ     *Ring
	polypoolT *Poly

	qHalf     *big.Int // (q-1)/2
	qHalfModT uint64   // (q-1)/2 mod t

	t    uint64
	qInv uint64 //(q mod t)^-1 mod t

	mredParamsT uint64

	paramsQP *modupParams
}

// NewRNSScaler creates a new SimpleScaler from t, the modulus under which the reconstruction is returned, the Ring in which the polynomial to reconstruct is represented.
func NewRNSScaler(t uint64, ringQ *Ring) (rnss *RNSScaler) {

	rnss = new(RNSScaler)

	rnss.ringQ = ringQ

	rnss.mredParamsT = MRedParams(t)

	rnss.polypoolT = NewPoly(ringQ.N, 1)

	rnss.t = t
	rnss.qHalf = new(big.Int)
	rnss.qInv = rnss.qHalf.Mod(ringQ.ModulusBigint, NewUint(t)).Uint64()
	rnss.qInv = ModExp(rnss.qInv, t-2, t)
	rnss.qInv = MForm(rnss.qInv, t, BRedParams(t))

	rnss.qHalf.Set(ringQ.ModulusBigint)
	rnss.qHalf.Rsh(rnss.qHalf, 1)
	rnss.qHalfModT = rnss.qHalf.Mod(rnss.qHalf, NewUint(t)).Uint64()

	rnss.qHalf.Set(ringQ.ModulusBigint)
	rnss.qHalf.Rsh(rnss.qHalf, 1)

	rnss.paramsQP = basisextenderparameters(ringQ.Modulus, []uint64{t})

	return
}

// DivByQOverTRounded returns p1 scaled by a factor t/Q and mod t on the receiver p2.
func (rnss *RNSScaler) DivByQOverTRounded(p1Q, p2T *Poly) {

	ringQ := rnss.ringQ

	T := rnss.t
	p2tmp := p2T.Coeffs[0]
	p3tmp := rnss.polypoolT.Coeffs[0]
	mredParams := rnss.mredParamsT
	qInv := T - rnss.qInv
	qHalfModT := T - rnss.qHalfModT

	// Multiply P_{Q} by t and extend the basis from P_{Q} to t*(P_{Q}||P_{t})
	// Since the coefficients of P_{t} are multiplied by t, they are all zero,
	// hence the basis extension can be omitted
	ringQ.MulScalar(p1Q, T, p1Q)

	// Center t*P_{Q} around (Q-1)/2 to round instead of floor during the division
	ringQ.AddScalarBigint(p1Q, rnss.qHalf, p1Q)

	// Extend the basis of (t*P_{Q} + (Q-1)/2) to (t*P_{t} + (Q-1)/2)
	modUpExact(p1Q.Coeffs, rnss.polypoolT.Coeffs, rnss.paramsQP)

	// Compute [Q^{-1} * (t*P_{t} -   (t*P_{Q} - ((Q-1)/2 mod t)))] mod t which returns round(t/Q * P_{Q}) mod t
	for j := uint64(0); j < ringQ.N; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))
		z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

		z[0] = MRed(qHalfModT+x[0], qInv, T, mredParams)
		z[1] = MRed(qHalfModT+x[1], qInv, T, mredParams)
		z[2] = MRed(qHalfModT+x[2], qInv, T, mredParams)
		z[3] = MRed(qHalfModT+x[3], qInv, T, mredParams)
		z[4] = MRed(qHalfModT+x[4], qInv, T, mredParams)
		z[5] = MRed(qHalfModT+x[5], qInv, T, mredParams)
		z[6] = MRed(qHalfModT+x[6], qInv, T, mredParams)
		z[7] = MRed(qHalfModT+x[7], qInv, T, mredParams)
	}
}

// SimpleScaler implements the Scaler interface by performing an RNS reconstruction and scaling by t/Q.
// This implementation of the Scaler interface is less efficient than the RNSScaler, but uses simple
// multi-precision arithmetic of the math/big package.
type SimpleScaler struct {
	ringQ *Ring

	r, one *big.Int

	t   uint64 // Plaintext modulus
	tBI *big.Int
	wi  []uint64 // Integer parts of   ([Q/Qi]^(-1))_{Qi} * t/Qi
	ti  []*big.Float

	p1BI []*big.Int

	reducealgoMul func(x, y uint64) uint64
	reducealgoAdd func(x uint64) uint64

	reducealgoAddParam uint64
	reducealgoMulParam uint64
}

// SimpleScalerFloatPrecision is the precision in bits for the big.Float in the scaling by t/Q.
const SimpleScalerFloatPrecision = 80

// NewSimpleScaler creates a new SimpleScaler from t, the modulus under which the reconstruction is returned, and
// ringQ, the Ring in which the polynomial to reconstruct is represented.
func NewSimpleScaler(t uint64, ringQ *Ring) (ss *SimpleScaler) {

	ss = new(SimpleScaler)

	ss.r = new(big.Int)
	ss.one = NewInt(1)

	var tmp *big.Float
	QiB := new(big.Int)         // Qi
	QiBF := new(big.Float)      // Qi
	QiStar := new(big.Int)      // Q/Qi
	QiBarre := new(big.Int)     // (Q/Qi)^(-1) mod Qi
	QiBarreBF := new(big.Float) // (Q/Qi)^(-1) mod Qi
	tBF := new(big.Float).SetUint64(t)

	ss.t = t
	ss.tBI = new(big.Int).SetUint64(t)
	ss.ringQ = ringQ

	// Assign the correct reduction algorithm depending on the provided t
	// If t is a power of 2
	if (t&(t-1)) == 0 && t != 0 {

		ss.reducealgoAddParam = t - 1
		ss.reducealgoMulParam = t - 1

		ss.reducealgoMul = func(x, y uint64) uint64 {
			return (x * y) & ss.reducealgoAddParam
		}

		ss.reducealgoAdd = func(x uint64) uint64 {
			return x & ss.reducealgoMulParam
		}

		// Otherwise, we can use Montgomery reduction
	} else {
		ss.reducealgoAddParam = BRedParams(t)[0]
		ss.reducealgoMulParam = MRedParams(t)

		ss.reducealgoMul = func(x, y uint64) uint64 {
			ahi, alo := bits.Mul64(x, y)
			R := alo * ss.reducealgoMulParam
			H, _ := bits.Mul64(R, ss.t)
			r := ahi - H + ss.t

			if r >= ss.t {
				r -= ss.t
			}

			return r
		}

		ss.reducealgoAdd = func(x uint64) uint64 {

			s0, _ := bits.Mul64(x, ss.reducealgoAddParam)

			r := x - s0*ss.t

			if r >= ss.t {
				r -= ss.t
			}

			return r
		}
	}

	// Integer and rational parts of QiBarre * t/Qi
	ss.wi = make([]uint64, len(ringQ.Modulus))
	ss.ti = make([]*big.Float, len(ringQ.Modulus))

	ss.p1BI = make([]*big.Int, ringQ.N)
	for i := range ss.p1BI {
		ss.p1BI[i] = new(big.Int)
		ss.p1BI[i].Mul(ss.ringQ.ModulusBigint, ss.ringQ.ModulusBigint) // Extend to Q^2
	}

	for i, qi := range ringQ.Modulus {

		QiB.SetUint64(qi)
		QiBF.SetUint64(qi)
		QiStar.Quo(ringQ.ModulusBigint, QiB)
		QiBarre.ModInverse(QiStar, QiB)
		QiBarre.Mod(QiBarre, QiB)
		QiBarreBF.SetInt(QiBarre)

		tmp = new(big.Float).Quo(tBF, QiBF)
		tmp.Mul(tmp, QiBarreBF)

		// floor( ([Q/Qi]^(-1))_{Qi} * t/Qi )
		ss.wi[i], _ = tmp.Uint64()

		// If t is not a power of 2, convert to Montgomery form
		if (t&(t-1)) != 0 && t != 0 {
			ss.wi[i] = MForm(ss.wi[i], t, BRedParams(t))
		}

		QiBarre.Mul(QiBarre, NewUint(t))
		QiBarre.Mod(QiBarre, QiB)
		QiBarreBF.SetInt(QiBarre)

		ss.ti[i] = new(big.Float).SetPrec(SimpleScalerFloatPrecision).Quo(QiBarreBF, QiBF)
	}

	return
}

// DivByQOverTRounded returns p1 scaled by a factor t/Q and mod t on the receiver p2.
func (ss *SimpleScaler) DivByQOverTRounded(p1, p2 *Poly) {
	ss.reconstructThenScale(p1, p2)
}

// reconstructThenScale performs the RNS reconstruction and scaling sequentially.
func (ss *SimpleScaler) reconstructThenScale(p1, p2 *Poly) {

	// reconstruction
	ss.ringQ.PolyToBigintNoAlloc(p1, ss.p1BI)

	// scaling
	for i, coeff := range ss.p1BI {
		coeff.Mul(coeff, ss.tBI)
		coeff.QuoRem(coeff, ss.ringQ.ModulusBigint, ss.r)
		if ss.r.Lsh(ss.r, 1).CmpAbs(ss.ringQ.ModulusBigint) != -1 {
			coeff.Add(coeff, ss.one)
		}

		for j := range p2.Coeffs {
			p2.Coeffs[j][i] = coeff.Uint64()
		}
	}
}

// reconstructAndScale performs the RNS reconstruction and scaling operations in one single pass over the coefficients.
// Algorithm from https://eprint.iacr.org/2018/117.pdf.
func (ss *SimpleScaler) reconstructAndScale(p1, p2 *Poly) {

	for i := uint64(0); i < ss.ringQ.N; i++ {

		var a uint64
		var bBF big.Float
		var p1j big.Float
		for j := range ss.ringQ.Modulus {
			// round(xi*wi + xi*ti)%t
			a += ss.reducealgoMul(ss.wi[j], p1.Coeffs[j][i])

			p1j.SetPrec(SimpleScalerFloatPrecision).SetUint64(p1.Coeffs[j][i])
			bBF.SetPrec(SimpleScalerFloatPrecision).Add(&bBF, p1j.Mul(&p1j, ss.ti[j]))
		}

		bBF.Add(&bBF, new(big.Float).SetFloat64(0.5))
		buint64, _ := bBF.Uint64()

		abf := a + buint64

		a = ss.reducealgoAdd(abf)

		for j := 0; j < len(p2.Coeffs); j++ {
			p2.Coeffs[j][i] = a
		}
	}
}

// ============== Scaling-related methods ==============

// DivFloorByLastModulusNTT divides (floored) the polynomial by its last modulus. The input must be in the NTT domain.
func (r *Ring) DivFloorByLastModulusNTT(p0 *Poly) {

	level := len(p0.Coeffs) - 1

	pTmp := make([]uint64, r.N)

	InvNTTLazy(p0.Coeffs[level], p0.Coeffs[level], r.N, r.NttPsiInv[level], r.NttNInv[level], r.Modulus[level])

	for i := 0; i < level; i++ {

		NTTLazy(p0.Coeffs[level], pTmp, r.N, r.NttPsi[i], r.Modulus[i])

		p0tmp := p0.Coeffs[i]

		qi := r.Modulus[i]
		twoqi := qi << 1
		rescaleParams := r.RescaleParams[level-1][i]

		// (x[i] - x[-1]) * InvQ
		for j := uint64(0); j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&pTmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))

			z[0] = FastBRed(twoqi-z[0]+x[0], rescaleParams, qi)
			z[1] = FastBRed(twoqi-z[1]+x[1], rescaleParams, qi)
			z[2] = FastBRed(twoqi-z[2]+x[2], rescaleParams, qi)
			z[3] = FastBRed(twoqi-z[3]+x[3], rescaleParams, qi)
			z[4] = FastBRed(twoqi-z[4]+x[4], rescaleParams, qi)
			z[5] = FastBRed(twoqi-z[5]+x[5], rescaleParams, qi)
			z[6] = FastBRed(twoqi-z[6]+x[6], rescaleParams, qi)
			z[7] = FastBRed(twoqi-z[7]+x[7], rescaleParams, qi)

		}
	}

	p0.Coeffs = p0.Coeffs[:level]
}

// DivFloorByLastModulus divides (floored) the polynomial by its last modulus.
func (r *Ring) DivFloorByLastModulus(p0 *Poly) {

	level := len(p0.Coeffs) - 1

	for i := 0; i < level; i++ {
		p0tmp := p0.Coeffs[level]
		p1tmp := p0.Coeffs[i]
		qi := r.Modulus[i]
		twoqi := qi << 1
		bredParams := r.BredParams[i]
		rescaleParams := r.RescaleParams[level-1][i]

		// (x[i] - x[-1]) * InvQ
		for j := uint64(0); j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))

			z[0] = FastBRed(twoqi-z[0]+BRedAdd(x[0], qi, bredParams), rescaleParams, qi)
			z[1] = FastBRed(twoqi-z[1]+BRedAdd(x[1], qi, bredParams), rescaleParams, qi)
			z[2] = FastBRed(twoqi-z[2]+BRedAdd(x[2], qi, bredParams), rescaleParams, qi)
			z[3] = FastBRed(twoqi-z[3]+BRedAdd(x[3], qi, bredParams), rescaleParams, qi)
			z[4] = FastBRed(twoqi-z[4]+BRedAdd(x[4], qi, bredParams), rescaleParams, qi)
			z[5] = FastBRed(twoqi-z[5]+BRedAdd(x[5], qi, bredParams), rescaleParams, qi)
			z[6] = FastBRed(twoqi-z[6]+BRedAdd(x[6], qi, bredParams), rescaleParams, qi)
			z[7] = FastBRed(twoqi-z[7]+BRedAdd(x[7], qi, bredParams), rescaleParams, qi)
		}
	}

	p0.Coeffs = p0.Coeffs[:level]
}

// DivFloorByLastModulusManyNTT divides (floored) sequentially nbRescales times the polynomial by its last modulus. Input must be in the NTT domain.
func (r *Ring) DivFloorByLastModulusManyNTT(p0 *Poly, nbRescales uint64) {
	r.InvNTTLvl(uint64(len(p0.Coeffs)-1), p0, p0)
	r.DivFloorByLastModulusMany(p0, nbRescales)
	r.NTTLvl(uint64(len(p0.Coeffs)-1), p0, p0)
}

// DivFloorByLastModulusMany divides (floored) sequentially nbRescales times the polynomial by its last modulus.
func (r *Ring) DivFloorByLastModulusMany(p0 *Poly, nbRescales uint64) {
	for k := uint64(0); k < nbRescales; k++ {
		r.DivFloorByLastModulus(p0)
	}
}

// DivRoundByLastModulusNTT divides (rounded) the polynomial by its last modulus. The input must be in the NTT domain.
func (r *Ring) DivRoundByLastModulusNTT(p0 *Poly) {

	var pHalf, pHalfNegQi uint64

	level := len(p0.Coeffs) - 1

	pTmp := make([]uint64, r.N)

	InvNTT(p0.Coeffs[level], p0.Coeffs[level], r.N, r.NttPsiInv[level], r.NttNInv[level], r.Modulus[level])

	// Center by (p-1)/2
	pj := r.Modulus[level]
	pHalf = (pj - 1) >> 1
	p0tmp := p0.Coeffs[level]

	for i := uint64(0); i < r.N; i = i + 8 {

		z := (*[8]uint64)(unsafe.Pointer(&p0tmp[i]))

		z[0] = CRed(z[0]+pHalf, pj)
		z[1] = CRed(z[1]+pHalf, pj)
		z[2] = CRed(z[2]+pHalf, pj)
		z[3] = CRed(z[3]+pHalf, pj)
		z[4] = CRed(z[4]+pHalf, pj)
		z[5] = CRed(z[5]+pHalf, pj)
		z[6] = CRed(z[6]+pHalf, pj)
		z[7] = CRed(z[7]+pHalf, pj)
	}

	for i := 0; i < level; i++ {

		p1tmp := p0.Coeffs[i]
		qi := r.Modulus[i]
		twoqi := qi << 1
		bredParams := r.BredParams[i]
		nttPsi := r.NttPsi[i]
		rescaleParams := r.RescaleParams[level-1][i]

		pHalfNegQi = r.Modulus[i] - BRedAdd(pHalf, qi, bredParams)

		for j := uint64(0); j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&pTmp[j]))

			z[0] = x[0] + pHalfNegQi
			z[1] = x[1] + pHalfNegQi
			z[2] = x[2] + pHalfNegQi
			z[3] = x[3] + pHalfNegQi
			z[4] = x[4] + pHalfNegQi
			z[5] = x[5] + pHalfNegQi
			z[6] = x[6] + pHalfNegQi
			z[7] = x[7] + pHalfNegQi
		}

		NTTLazy(pTmp, pTmp, r.N, nttPsi, qi)

		// (x[i] - x[-1]) * InvQ
		for j := uint64(0); j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&pTmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))

			z[0] = FastBRed(twoqi+x[0]-z[0], rescaleParams, qi)
			z[1] = FastBRed(twoqi+x[1]-z[1], rescaleParams, qi)
			z[2] = FastBRed(twoqi+x[2]-z[2], rescaleParams, qi)
			z[3] = FastBRed(twoqi+x[3]-z[3], rescaleParams, qi)
			z[4] = FastBRed(twoqi+x[4]-z[4], rescaleParams, qi)
			z[5] = FastBRed(twoqi+x[5]-z[5], rescaleParams, qi)
			z[6] = FastBRed(twoqi+x[6]-z[6], rescaleParams, qi)
			z[7] = FastBRed(twoqi+x[7]-z[7], rescaleParams, qi)
		}
	}

	p0.Coeffs = p0.Coeffs[:level]
}

// DivRoundByLastModulus divides (rounded) the polynomial by its last modulus. The input must be in the NTT domain.
func (r *Ring) DivRoundByLastModulus(p0 *Poly) {

	var pHalf, pHalfNegQi uint64

	level := len(p0.Coeffs) - 1

	// Center by (p-1)/2
	pHalf = (r.Modulus[level] - 1) >> 1
	p0tmp := p0.Coeffs[level]
	pj := r.Modulus[level]

	for i := uint64(0); i < r.N; i = i + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p0tmp[i]))

		x[0] = CRed(x[0]+pHalf, pj)
		x[1] = CRed(x[1]+pHalf, pj)
		x[2] = CRed(x[2]+pHalf, pj)
		x[3] = CRed(x[3]+pHalf, pj)
		x[4] = CRed(x[4]+pHalf, pj)
		x[5] = CRed(x[5]+pHalf, pj)
		x[6] = CRed(x[6]+pHalf, pj)
		x[7] = CRed(x[7]+pHalf, pj)
	}

	for i := 0; i < level; i++ {

		p1tmp := p0.Coeffs[i]

		qi := r.Modulus[i]
		twoqi := qi << 1
		bredParams := r.BredParams[i]
		rescaleParams := r.RescaleParams[level-1][i]

		pHalfNegQi = r.Modulus[i] - BRedAdd(pHalf, qi, bredParams)

		// (x[i] - x[-1]) * InvQ
		for j := uint64(0); j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))

			z[0] = FastBRed(x[0]+pHalfNegQi+twoqi-z[0], rescaleParams, qi)
			z[1] = FastBRed(x[1]+pHalfNegQi+twoqi-z[1], rescaleParams, qi)
			z[2] = FastBRed(x[2]+pHalfNegQi+twoqi-z[2], rescaleParams, qi)
			z[3] = FastBRed(x[3]+pHalfNegQi+twoqi-z[3], rescaleParams, qi)
			z[4] = FastBRed(x[4]+pHalfNegQi+twoqi-z[4], rescaleParams, qi)
			z[5] = FastBRed(x[5]+pHalfNegQi+twoqi-z[5], rescaleParams, qi)
			z[6] = FastBRed(x[6]+pHalfNegQi+twoqi-z[6], rescaleParams, qi)
			z[7] = FastBRed(x[7]+pHalfNegQi+twoqi-z[7], rescaleParams, qi)
		}
	}

	p0.Coeffs = p0.Coeffs[:level]
}

// DivRoundByLastModulusManyNTT divides (rounded) sequentially nbRescales times the polynomial by its last modulus. The input must be in the NTT domain.
func (r *Ring) DivRoundByLastModulusManyNTT(p0 *Poly, nbRescales uint64) {
	r.InvNTTLvl(uint64(len(p0.Coeffs)-1), p0, p0)
	r.DivRoundByLastModulusMany(p0, nbRescales)
	r.NTTLvl(uint64(len(p0.Coeffs)-1), p0, p0)
}

// DivRoundByLastModulusMany divides (rounded) sequentially nbRescales times the polynomial by its last modulus.
func (r *Ring) DivRoundByLastModulusMany(p0 *Poly, nbRescales uint64) {
	for k := uint64(0); k < nbRescales; k++ {
		r.DivRoundByLastModulus(p0)
	}
}
