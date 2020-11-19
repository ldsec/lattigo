package rckks

import (
	"math"
	"math/bits"
)

// Poly is a struct storing the coeffients of a polynomial
// that then can be evaluated on the ciphertext
type Poly struct {
	maxDeg uint64
	coeffs []float64
	lead   bool
}

// NewPoly creates a new Poly from the input coefficients
func NewPoly(coeffs []float64) (p *Poly) {

	p = new(Poly)
	p.coeffs = make([]float64, len(coeffs), len(coeffs))
	copy(p.coeffs, coeffs)
	p.maxDeg = uint64(len(coeffs) - 1)
	p.lead = true

	return
}

// Degree returns the degree of the polynomial
func (p *Poly) Degree() uint64 {
	return uint64(len(p.coeffs) - 1)
}

func optimalSplit(logDegree uint64) (logSplit uint64) {
	logSplit = logDegree >> 1
	a := (1 << logSplit) + (1 << (logDegree - logSplit)) + logDegree - logSplit - 3
	b := (1 << (logSplit + 1)) + (1 << (logDegree - logSplit - 1)) + logDegree - logSplit - 4
	if a > b {
		logSplit++
	}

	return
}

func computeSmallPoly(split uint64, coeffs *Poly) (polyList []*Poly) {

	if coeffs.Degree() < (1 << split) {
		return []*Poly{coeffs}
	}

	var nextPower = uint64(1 << split)
	for nextPower < (coeffs.Degree()>>1)+1 {
		nextPower <<= 1
	}

	coeffsq, coeffsr := splitCoeffsCheby(coeffs, nextPower)

	a := computeSmallPoly(split, coeffsq)
	b := computeSmallPoly(split, coeffsr)

	return append(a, b...)
}

// EvaluatePoly evaluates a polynomial in standard basis on the input Ciphertext in ceil(log2(deg+1)) levels.
func (eval *evaluator) EvaluatePoly(ct0 *Ciphertext, pol *Poly, evakey *EvaluationKey) (res *Ciphertext) {

	C := make(map[uint64]*Ciphertext)

	C[1] = ct0.CopyNew().Ciphertext()

	logDegree := uint64(bits.Len64(pol.Degree()))
	logSplit := (logDegree >> 1) //optimalSplit(logDegree) //

	for i := uint64(2); i < (1 << logSplit); i++ {
		computePowerBasis(i, C, eval, evakey)
	}

	for i := logSplit; i < logDegree; i++ {
		computePowerBasis(1<<i, C, eval, evakey)
	}

	res = recurse(eval.scale, logSplit, logDegree, pol, C, eval, evakey)

	C = nil

	return
}

// EvaluateCheby evaluates a polynomial in Chebyshev basis on the input Ciphertext in ceil(log2(deg+1))+1 levels.
func (eval *evaluator) EvaluateCheby(op *Ciphertext, cheby *ChebyshevInterpolation, evakey *EvaluationKey) (opOut *Ciphertext) {

	C := make(map[uint64]*Ciphertext)

	C[1] = op.CopyNew().Ciphertext()

	eval.MultByConst(C[1], 2/(cheby.b-cheby.a), C[1])
	eval.AddConst(C[1], (-cheby.a-cheby.b)/(cheby.b-cheby.a), C[1])
	eval.Rescale(C[1], eval.scale, C[1])

	return eval.evalCheby(cheby, C, evakey)
}

// EvaluateChebyFastSpecial evaluates the input Chebyshev polynomial with the input ciphertext.
// Slower than EvaluateChebyFast but consumes ceil(log(deg)) + 1 levels.
func (eval *evaluator) EvaluateChebySpecial(ct *Ciphertext, n float64, cheby *ChebyshevInterpolation, evakey *EvaluationKey) (res *Ciphertext) {

	C := make(map[uint64]*Ciphertext)

	C[1] = ct.CopyNew().Ciphertext()

	eval.MultByConst(C[1], 2/((cheby.b-cheby.a)*n), C[1])
	eval.AddConst(C[1], (-cheby.a-cheby.b)/(cheby.b-cheby.a), C[1])
	eval.Rescale(C[1], eval.scale, C[1])

	return eval.evalCheby(cheby, C, evakey)
}

func (eval *evaluator) evalCheby(cheby *ChebyshevInterpolation, C map[uint64]*Ciphertext, evakey *EvaluationKey) (res *Ciphertext) {

	logDegree := uint64(bits.Len64(cheby.Degree()))
	logSplit := (logDegree >> 1) //optimalSplit(logDegree) //

	for i := uint64(2); i < (1 << logSplit); i++ {
		computePowerBasisCheby(i, C, eval, evakey)
	}

	for i := logSplit; i < logDegree; i++ {
		computePowerBasisCheby(1<<i, C, eval, evakey)
	}

	res = recurseCheby(eval.scale, logSplit, logDegree, &cheby.Poly, C, eval, evakey)

	C = nil

	return
}

func computePowerBasis(n uint64, C map[uint64]*Ciphertext, evaluator *evaluator, evakey *EvaluationKey) {

	if C[n] == nil {

		// Computes the index required to compute the asked ring evaluation
		a := uint64(math.Ceil(float64(n) / 2))
		b := n >> 1

		// Recurses on the given indexes
		computePowerBasis(a, C, evaluator, evakey)
		computePowerBasis(b, C, evaluator, evakey)

		// Computes C[n] = C[a]*C[b]
		C[n] = evaluator.MulRelinNew(C[a], C[b], evakey)

		evaluator.Rescale(C[n], evaluator.scale, C[n])
	}
}

func computePowerBasisCheby(n uint64, C map[uint64]*Ciphertext, evaluator *evaluator, evakey *EvaluationKey) {

	// Given a hash table with the first three evaluations of the Chebyshev ring at x in the interval a, b:
	// C0 = 1 (actually not stored in the hash table)
	// C1 = (2*x - a - b)/(b-a)
	// C2 = 2*C1*C1 - C0
	// Evaluates the nth degree Chebyshev ring in a recursive manner, storing intermediate results in the hashtable.
	// Consumes at most ceil(sqrt(n)) levels for an evaluation at Cn.
	// Uses the following property: for a given Chebyshev ring Cn = 2*Ca*Cb - Cc, n = a+b and c = abs(a-b)

	if C[n] == nil {

		// Computes the index required to compute the asked ring evaluation
		a := uint64(math.Ceil(float64(n) / 2))
		b := n >> 1
		c := uint64(math.Abs(float64(a) - float64(b)))

		// Recurses on the given indexes
		computePowerBasisCheby(a, C, evaluator, evakey)
		computePowerBasisCheby(b, C, evaluator, evakey)

		// Since C[0] is not stored (but rather seen as the constant 1), only recurses on c if c!= 0
		if c != 0 {
			computePowerBasisCheby(c, C, evaluator, evakey)
		}

		// Computes C[n] = C[a]*C[b]
		C[n] = evaluator.MulRelinNew(C[a], C[b], evakey)
		evaluator.Rescale(C[n], evaluator.scale, C[n])

		// Computes C[n] = 2*C[a]*C[b]
		evaluator.Add(C[n], C[n], C[n])

		// Computes C[n] = 2*C[a]*C[b] - C[c]
		if c == 0 {
			evaluator.AddConst(C[n], -1, C[n])
		} else {
			evaluator.Sub(C[n], C[c], C[n])
		}

	}
}

func splitCoeffs(coeffs *Poly, split uint64) (coeffsq, coeffsr *Poly) {

	// Splits a polynomial p such that p = q*C^degree + r.

	coeffsr = new(Poly)
	coeffsr.coeffs = make([]float64, split)
	if coeffs.maxDeg == coeffs.Degree() {
		coeffsr.maxDeg = split - 1
	} else {
		coeffsr.maxDeg = coeffs.maxDeg - (coeffs.Degree() - split + 1)
	}

	for i := uint64(0); i < split; i++ {
		coeffsr.coeffs[i] = coeffs.coeffs[i]
	}

	coeffsq = new(Poly)
	coeffsq.coeffs = make([]float64, coeffs.Degree()-split+1)
	coeffsq.maxDeg = coeffs.maxDeg

	coeffsq.coeffs[0] = coeffs.coeffs[split]
	for i := split + 1; i < coeffs.Degree()+1; i++ {
		coeffsq.coeffs[i-split] = coeffs.coeffs[i]
	}

	if coeffs.lead {
		coeffsq.lead = true
	}

	return coeffsq, coeffsr
}

func splitCoeffsCheby(coeffs *Poly, split uint64) (coeffsq, coeffsr *Poly) {

	// Splits a Chebyshev polynomial p such that p = q*C^degree + r, where q and r are a linear combination of a Chebyshev basis.
	coeffsr = new(Poly)
	coeffsr.coeffs = make([]float64, split)
	if coeffs.maxDeg == coeffs.Degree() {
		coeffsr.maxDeg = split - 1
	} else {
		coeffsr.maxDeg = coeffs.maxDeg - (coeffs.Degree() - split + 1)
	}

	for i := uint64(0); i < split; i++ {
		coeffsr.coeffs[i] = coeffs.coeffs[i]
	}

	coeffsq = new(Poly)
	coeffsq.coeffs = make([]float64, coeffs.Degree()-split+1)
	coeffsq.maxDeg = coeffs.maxDeg

	coeffsq.coeffs[0] = coeffs.coeffs[split]
	for i, j := split+1, uint64(1); i < coeffs.Degree()+1; i, j = i+1, j+1 {
		coeffsq.coeffs[i-split] = 2 * coeffs.coeffs[i]
		coeffsr.coeffs[split-j] -= coeffs.coeffs[i]
	}

	if coeffs.lead {
		coeffsq.lead = true
	}

	return coeffsq, coeffsr
}

func recurse(targetScale float64, logSplit, logDegree uint64, coeffs *Poly, C map[uint64]*Ciphertext, evaluator *evaluator, evakey *EvaluationKey) (res *Ciphertext) {
	// Recursively computes the evalution of the Chebyshev polynomial using a baby-set giant-step algorithm.
	if coeffs.Degree() < (1 << logSplit) {

		if coeffs.lead && coeffs.maxDeg > ((1<<logDegree)-(1<<(logSplit-1))) && logSplit > 1 {

			logDegree = uint64(bits.Len64(coeffs.Degree()))
			logSplit = logDegree >> 1

			return recurse(targetScale, logSplit, logDegree, coeffs, C, evaluator, evakey)
		}

		return evaluatePolyFromPowerBasis(targetScale, coeffs, C, evaluator)
	}

	var nextPower = uint64(1 << logSplit)
	for nextPower < (coeffs.Degree()>>1)+1 {
		nextPower <<= 1
	}

	coeffsq, coeffsr := splitCoeffs(coeffs, nextPower)

	level := C[nextPower].Level() - 1

	if coeffsq.maxDeg >= 1<<(logDegree-1) && coeffsq.lead {
		level++
	}

	currentQi := float64(evaluator.params.qi[level])

	res = recurse(targetScale*currentQi/C[nextPower].Scale(), logSplit, logDegree, coeffsq, C, evaluator, evakey)

	tmp := recurse(targetScale, logSplit, logDegree, coeffsr, C, evaluator, evakey)

	if res.Level() > tmp.Level() {
		for res.Level() != tmp.Level()+1 {
			evaluator.DropLevel(res, 1)
		}
	}

	evaluator.MulRelin(res, C[nextPower], evakey, res)

	if res.Level() > tmp.Level() {
		evaluator.Rescale(res, evaluator.scale, res)
		evaluator.Add(res, tmp, res)
	} else {
		evaluator.Add(res, tmp, res)
		evaluator.Rescale(res, evaluator.scale, res)
	}

	tmp = nil

	return
}

func recurseCheby(targetScale float64, logSplit, logDegree uint64, coeffs *Poly, C map[uint64]*Ciphertext, evaluator *evaluator, evakey *EvaluationKey) (res *Ciphertext) {
	// Recursively computes the evalution of the Chebyshev polynomial using a baby-set giant-step algorithm.
	if coeffs.Degree() < (1 << logSplit) {

		if coeffs.lead && coeffs.maxDeg > ((1<<logDegree)-(1<<(logSplit-1))) && logSplit > 1 {

			logDegree = uint64(bits.Len64(coeffs.Degree()))
			logSplit = logDegree >> 1

			return recurseCheby(targetScale, logSplit, logDegree, coeffs, C, evaluator, evakey)
		}

		return evaluatePolyFromPowerBasis(targetScale, coeffs, C, evaluator)
	}

	var nextPower = uint64(1 << logSplit)
	for nextPower < (coeffs.Degree()>>1)+1 {
		nextPower <<= 1
	}

	coeffsq, coeffsr := splitCoeffsCheby(coeffs, nextPower)

	level := C[nextPower].Level() - 1

	if coeffsq.maxDeg >= 1<<(logDegree-1) && coeffsq.lead {
		level++
	}

	currentQi := float64(evaluator.params.qi[level])

	res = recurseCheby(targetScale*currentQi/C[nextPower].Scale(), logSplit, logDegree, coeffsq, C, evaluator, evakey)

	var tmp *Ciphertext
	tmp = recurseCheby(targetScale, logSplit, logDegree, coeffsr, C, evaluator, evakey)

	if res.Level() > tmp.Level() {
		for res.Level() != tmp.Level()+1 {
			evaluator.DropLevel(res, 1)
		}
	}

	evaluator.MulRelin(res, C[nextPower], evakey, res)

	if res.Level() > tmp.Level() {
		evaluator.Rescale(res, evaluator.scale, res)
		evaluator.Add(res, tmp, res)
	} else {
		evaluator.Add(res, tmp, res)
		evaluator.Rescale(res, evaluator.scale, res)
	}

	tmp = nil

	return

}

func evaluatePolyFromPowerBasis(targetScale float64, coeffs *Poly, C map[uint64]*Ciphertext, evaluator *evaluator) (res *Ciphertext) {

	if coeffs.Degree() == 0 {

		res = NewCiphertext(evaluator.params, 1, C[1].Level(), targetScale)

		if math.Abs(coeffs.coeffs[0]) > 1e-14 {
			evaluator.AddConst(res, coeffs.coeffs[0], res)
		}

		return
	}

	currentQi := float64(evaluator.params.qi[C[coeffs.Degree()].Level()])

	ctScale := targetScale * currentQi

	res = NewCiphertext(evaluator.params, 1, C[coeffs.Degree()].Level(), ctScale)

	if math.Abs(coeffs.coeffs[0]) > 1e-14 {
		evaluator.AddConst(res, coeffs.coeffs[0], res)
	}

	for key := coeffs.Degree(); key > 0; key-- {

		if key != 0 && (math.Abs(coeffs.coeffs[key]) > 1e-14) {

			// Target scale * rescale-scale / power basis scale
			constScale := targetScale * currentQi / C[key].Scale()

			evaluator.multByIntegerAndAdd(C[key], int64(coeffs.coeffs[key]*constScale), res)
		}
	}

	evaluator.Rescale(res, evaluator.scale, res)

	return
}
