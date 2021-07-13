package ckks

import (
	"fmt"
	"math"
	"math/bits"
)

// Poly is a struct storing the coeffients of a polynomial
// that then can be evaluated on the ciphertext
type Poly struct {
	maxDeg int
	coeffs []complex128
	lead   bool
}

// NewPoly creates a new Poly from the input coefficients
func NewPoly(coeffs []complex128) (p *Poly) {

	p = new(Poly)
	p.coeffs = make([]complex128, len(coeffs))
	copy(p.coeffs, coeffs)
	p.maxDeg = len(coeffs) - 1
	p.lead = true

	return
}

// checkEnoughLevels checks that enough levels are available to evaluate the polynomial.
// Also checks if c is a gaussian integer or not. If not, then one more level is needed
// to evaluate the polynomial.
func checkEnoughLevels(levels int, pol *Poly, c complex128) (err error) {

	logDegree := int(math.Log2(float64(len(pol.coeffs))) + 0.5)

	if real(c) != float64(int64(real(c))) || imag(c) != float64(int64(imag(c))) {
		logDegree++
	}

	if levels < logDegree {
		return fmt.Errorf("%d levels < %d log(d) -> cannot evaluate", levels, logDegree)
	}

	return nil
}

// Degree returns the degree of the polynomial
func (p *Poly) Degree() int {
	return len(p.coeffs) - 1
}

// EvaluatePoly evaluates a polynomial in standard basis on the input Ciphertext in ceil(log2(deg+1)) levels.
// Returns an error if the input ciphertext does not have enough level to carry out the full polynomial evaluation.
// Returns an error if something is wrong with the scale.

func (eval *evaluator) EvaluatePoly(ct0 *Ciphertext, pol *Poly, targetScale float64) (opOut *Ciphertext, err error) {

	if err := checkEnoughLevels(ct0.Level(), pol, 1); err != nil {
		return ct0, err
	}

	C := make(map[int]*Ciphertext)

	C[1] = ct0.CopyNew()

	logDegree := bits.Len64(uint64(pol.Degree()))
	logSplit := (logDegree >> 1) //optimalSplit(logDegree) //

	for i := 2; i < (1 << logSplit); i++ {
		if err = computePowerBasis(i, C, eval); err != nil {
			return nil, err
		}
	}

	for i := logSplit; i < logDegree; i++ {
		if err = computePowerBasis(1<<i, C, eval); err != nil {
			return nil, err
		}
	}

	opOut, err = recurse(targetScale, logSplit, logDegree, pol, C, eval)

	C = nil
	return opOut, err
}

// EvaluateCheby evaluates a polynomial in Chebyshev basis on the input Ciphertext in ceil(log2(deg+1))+1 levels.
// Returns an error if the input ciphertext does not have enough level to carry out the full polynomial evaluation.
// Returns an error if something is wrong with the scale.
// A change of basis ct' = (2/(b-a)) * (ct + (-a-b)/(b-a)) is necessary before the polynomial evaluation to ensure correctness.

func (eval *evaluator) EvaluateCheby(op *Ciphertext, cheby *ChebyshevInterpolation, tartetScale float64) (opOut *Ciphertext, err error) {

	if err := checkEnoughLevels(op.Level(), &cheby.Poly, 1); err != nil {
		return op, err
	}

	C := make(map[int]*Ciphertext)

	C[1] = op.CopyNew()

	logDegree := int(bits.Len64(uint64(cheby.Degree())))
	logSplit := (logDegree >> 1) //optimalSplit(logDegree) //

	for i := 2; i < (1 << logSplit); i++ {
		if err = computePowerBasisCheby(i, C, eval); err != nil {
			return nil, err
		}
	}

	for i := logSplit; i < logDegree; i++ {
		if err = computePowerBasisCheby(1<<i, C, eval); err != nil {
			return nil, err
		}
	}

	opOut, err = recurseCheby(tartetScale, logSplit, logDegree, &cheby.Poly, C, eval)

	C = nil

	return opOut, err
}

func computePowerBasis(n int, C map[int]*Ciphertext, evaluator *evaluator) (err error) {

	if C[n] == nil {

		// Computes the index required to compute the asked ring evaluation
		a := int(math.Ceil(float64(n) / 2))
		b := n >> 1

		// Recurses on the given indexes
		if err = computePowerBasis(a, C, evaluator); err != nil {
			return err
		}
		if err = computePowerBasis(b, C, evaluator); err != nil {
			return err
		}

		// Computes C[n] = C[a]*C[b]
		C[n] = evaluator.MulRelinNew(C[a], C[b])

		if err = evaluator.Rescale(C[n], evaluator.scale, C[n]); err != nil {
			return err
		}
	}

	return nil
}

func computePowerBasisCheby(n int, C map[int]*Ciphertext, evaluator *evaluator) (err error) {

	// Given a hash table with the first three evaluations of the Chebyshev ring at x in the interval a, b:
	// C0 = 1 (actually not stored in the hash table)
	// C1 = (2*x - a - b)/(b-a)
	// C2 = 2*C1*C1 - C0
	// Evaluates the nth degree Chebyshev ring in a recursive manner, storing intermediate results in the hashtable.
	// Consumes at most ceil(sqrt(n)) levels for an evaluation at Cn.
	// Uses the following property: for a given Chebyshev ring Cn = 2*Ca*Cb - Cc, n = a+b and c = abs(a-b)

	if C[n] == nil {

		// Computes the index required to compute the asked ring evaluation
		a := int(math.Ceil(float64(n) / 2))
		b := n >> 1
		c := int(math.Abs(float64(a) - float64(b)))

		// Recurses on the given indexes
		if err = computePowerBasisCheby(a, C, evaluator); err != nil {
			return err
		}
		if err = computePowerBasisCheby(b, C, evaluator); err != nil {
			return err
		}

		// Since C[0] is not stored (but rather seen as the constant 1), only recurses on c if c!= 0
		if c != 0 {
			if err = computePowerBasisCheby(c, C, evaluator); err != nil {
				return err
			}
		}

		// Computes C[n] = C[a]*C[b]
		//fmt.Println("Mul", C[a].Level(), C[b].Level())
		C[n] = evaluator.MulRelinNew(C[a], C[b])
		if err = evaluator.Rescale(C[n], evaluator.scale, C[n]); err != nil {
			return err
		}

		// Computes C[n] = 2*C[a]*C[b]
		evaluator.Add(C[n], C[n], C[n])

		// Computes C[n] = 2*C[a]*C[b] - C[c]
		if c == 0 {
			evaluator.AddConst(C[n], -1, C[n])
		} else {
			evaluator.Sub(C[n], C[c], C[n])
		}

	}

	return nil
}

func splitCoeffs(coeffs *Poly, split int) (coeffsq, coeffsr *Poly) {

	// Splits a polynomial p such that p = q*C^degree + r.

	coeffsr = new(Poly)
	coeffsr.coeffs = make([]complex128, split)
	if coeffs.maxDeg == coeffs.Degree() {
		coeffsr.maxDeg = split - 1
	} else {
		coeffsr.maxDeg = coeffs.maxDeg - (coeffs.Degree() - split + 1)
	}

	for i := 0; i < split; i++ {
		coeffsr.coeffs[i] = coeffs.coeffs[i]
	}

	coeffsq = new(Poly)
	coeffsq.coeffs = make([]complex128, coeffs.Degree()-split+1)
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

func splitCoeffsCheby(coeffs *Poly, split int) (coeffsq, coeffsr *Poly) {

	// Splits a Chebyshev polynomial p such that p = q*C^degree + r, where q and r are a linear combination of a Chebyshev basis.
	coeffsr = new(Poly)
	coeffsr.coeffs = make([]complex128, split)
	if coeffs.maxDeg == coeffs.Degree() {
		coeffsr.maxDeg = split - 1
	} else {
		coeffsr.maxDeg = coeffs.maxDeg - (coeffs.Degree() - split + 1)
	}

	for i := 0; i < split; i++ {
		coeffsr.coeffs[i] = coeffs.coeffs[i]
	}

	coeffsq = new(Poly)
	coeffsq.coeffs = make([]complex128, coeffs.Degree()-split+1)
	coeffsq.maxDeg = coeffs.maxDeg

	coeffsq.coeffs[0] = coeffs.coeffs[split]
	for i, j := split+1, 1; i < coeffs.Degree()+1; i, j = i+1, j+1 {
		coeffsq.coeffs[i-split] = 2 * coeffs.coeffs[i]
		coeffsr.coeffs[split-j] -= coeffs.coeffs[i]
	}

	if coeffs.lead {
		coeffsq.lead = true
	}

	return coeffsq, coeffsr
}

func recurse(targetScale float64, logSplit, logDegree int, coeffs *Poly, C map[int]*Ciphertext, evaluator *evaluator) (res *Ciphertext, err error) {

	// Recursively computes the evalution of the Chebyshev polynomial using a baby-set giant-step algorithm.
	if coeffs.Degree() < (1 << logSplit) {

		if coeffs.lead && coeffs.maxDeg > ((1<<logDegree)-(1<<(logSplit-1))) && logSplit > 1 {

			logDegree = int(bits.Len64(uint64(coeffs.Degree())))
			logSplit = logDegree >> 1

			return recurse(targetScale, logSplit, logDegree, coeffs, C, evaluator)
		}

		return evaluatePolyFromPowerBasis(targetScale, coeffs, C, evaluator)
	}

	var nextPower = 1 << logSplit
	for nextPower < (coeffs.Degree()>>1)+1 {
		nextPower <<= 1
	}

	coeffsq, coeffsr := splitCoeffs(coeffs, nextPower)

	level := C[nextPower].Level() - 1

	if coeffsq.maxDeg >= 1<<(logDegree-1) && coeffsq.lead {
		level++
	}

	currentQi := float64(evaluator.params.Q()[level])

	//fmt.Printf("X^%2d: %d %d %t %d\n", nextPower, coeffsq.maxDeg, coeffsr.maxDeg, coeffsq.maxDeg >= 1<<(logDegree-1), level)
	//fmt.Printf("X^%2d: %f %f\n", nextPower, targetScale, targetScale* currentQi / C[nextPower].Scale())
	//fmt.Printf("X^%2d : qi %d %t %d %d\n", nextPower, level, coeffsq.lead, coeffsq.maxDeg, 1<<(logDegree-1))
	//fmt.Println()
	var tmp *Ciphertext
	if res, err = recurse(targetScale*currentQi/C[nextPower].Scale, logSplit, logDegree, coeffsq, C, evaluator); err != nil {
		return nil, err
	}

	if tmp, err = recurse(targetScale, logSplit, logDegree, coeffsr, C, evaluator); err != nil {
		return nil, err
	}

	if res.Level() > tmp.Level() {
		for res.Level() != tmp.Level()+1 {
			evaluator.DropLevel(res, 1)
		}
	}

	//fmt.Printf("X^%2d: (%d %f -> \n", nextPower, res.Level(), res.Scale())
	evaluator.MulRelin(res, C[nextPower], res)

	if res.Level() > tmp.Level() {
		if err = evaluator.Rescale(res, evaluator.scale, res); err != nil {
			return nil, err
		}
		//fmt.Printf("%f = %d) + (%d %f) = ", res.Scale(), res.Level(), tmp.Level(), tmp.Scale())
		evaluator.Add(res, tmp, res)
		//fmt.Printf("(%d %f) %f\n", res.Level(), res.Scale(), res.Scale()-tmp.Scale())
	} else {
		evaluator.Add(res, tmp, res)
		if err = evaluator.Rescale(res, evaluator.scale, res); err != nil {
			return nil, err
		}

	}

	tmp = nil

	return
}

func recurseCheby(targetScale float64, logSplit, logDegree int, coeffs *Poly, C map[int]*Ciphertext, evaluator *evaluator) (res *Ciphertext, err error) {

	// Recursively computes the evalution of the Chebyshev polynomial using a baby-set giant-step algorithm.
	if coeffs.Degree() < (1 << logSplit) {

		if coeffs.lead && coeffs.maxDeg > ((1<<logDegree)-(1<<(logSplit-1))) && logSplit > 1 {

			logDegree = int(bits.Len64(uint64(coeffs.Degree())))
			logSplit = logDegree >> 1

			return recurseCheby(targetScale, logSplit, logDegree, coeffs, C, evaluator)
		}

		return evaluatePolyFromPowerBasis(targetScale, coeffs, C, evaluator)
	}

	var nextPower = 1 << logSplit
	for nextPower < (coeffs.Degree()>>1)+1 {
		nextPower <<= 1
	}

	coeffsq, coeffsr := splitCoeffsCheby(coeffs, nextPower)

	level := C[nextPower].Level() - 1

	if coeffsq.maxDeg >= 1<<(logDegree-1) && coeffsq.lead {
		level++
	}

	currentQi := float64(evaluator.params.Q()[level])

	//fmt.Printf("X^%2d: %d %d %t %d\n", nextPower, coeffsq.maxDeg, coeffsr.maxDeg, coeffsq.maxDeg >= 1<<(logDegree-1), level)
	//fmt.Printf("X^%2d: %f %f\n", nextPower, targetScale, targetScale* currentQi / C[nextPower].Scale())
	//fmt.Printf("X^%2d : qi %d %t %d %d\n", nextPower, level, coeffsq.lead, coeffsq.maxDeg, 1<<(logDegree-1))
	//fmt.Println()

	if res, err = recurseCheby(targetScale*currentQi/C[nextPower].Scale, logSplit, logDegree, coeffsq, C, evaluator); err != nil {
		return nil, err
	}

	var tmp *Ciphertext
	if tmp, err = recurseCheby(targetScale, logSplit, logDegree, coeffsr, C, evaluator); err != nil {
		return nil, err
	}

	if res.Level() > tmp.Level() {
		for res.Level() != tmp.Level()+1 {
			evaluator.DropLevel(res, 1)
		}
	}

	//fmt.Printf("X^%2d: (%d %f -> \n", nextPower, res.Level(), res.Scale())
	evaluator.MulRelin(res, C[nextPower], res)

	if res.Level() > tmp.Level() {
		if err = evaluator.Rescale(res, evaluator.scale, res); err != nil {
			return nil, err
		}
		//fmt.Printf("%f = %d) + (%d %f) = ", res.Scale(), res.Level(), tmp.Level(), tmp.Scale())
		evaluator.Add(res, tmp, res)
		//fmt.Printf("(%d %f) %f\n", res.Level(), res.Scale(), res.Scale()-tmp.Scale())
	} else {
		evaluator.Add(res, tmp, res)
		if err = evaluator.Rescale(res, evaluator.scale, res); err != nil {
			return nil, err
		}
	}

	tmp = nil

	return

}

func evaluatePolyFromPowerBasis(targetScale float64, coeffs *Poly, C map[int]*Ciphertext, evaluator *evaluator) (res *Ciphertext, err error) {

	if coeffs.Degree() == 0 {

		res = NewCiphertext(evaluator.params, 1, C[1].Level(), targetScale)

		if math.Abs(real(coeffs.coeffs[0])) > 1e-14 || math.Abs(imag(coeffs.coeffs[0])) > 1e-14 {
			evaluator.AddConst(res, coeffs.coeffs[0], res)
		}

		return
	}

	currentQi := float64(evaluator.params.Q()[C[coeffs.Degree()].Level()])

	ctScale := targetScale * currentQi

	//fmt.Printf("%d %f\n", coeffs.maxDeg, targetScale)
	//fmt.Println("current Qi", evaluator.params.Qi[C[coeffs.Degree()].Level()])
	//fmt.Println(coeffs.Degree(), C[coeffs.Degree()].Level())

	res = NewCiphertext(evaluator.params, 1, C[coeffs.Degree()].Level(), ctScale)

	if math.Abs(real(coeffs.coeffs[0])) > 1e-14 || math.Abs(imag(coeffs.coeffs[0])) > 1e-14 {
		evaluator.AddConst(res, coeffs.coeffs[0], res)
	}

	for key := coeffs.Degree(); key > 0; key-- {

		if key != 0 && (math.Abs(real(coeffs.coeffs[key])) > 1e-14 || math.Abs(imag(coeffs.coeffs[key])) > 1e-14) {

			// Target scale * rescale-scale / power basis scale
			constScale := targetScale * currentQi / C[key].Scale

			cReal := int64(real(coeffs.coeffs[key]) * constScale)
			cImag := int64(imag(coeffs.coeffs[key]) * constScale)

			evaluator.MultByGaussianIntegerAndAdd(C[key], cReal, cImag, res)
		}
	}

	if err = evaluator.Rescale(res, evaluator.scale, res); err != nil {
		return nil, err
	}

	return
}
