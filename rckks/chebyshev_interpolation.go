package rckks

import (
	"math"
)

// ChebyshevInterpolation is a struct storing the coefficients, degree and range of a Chebyshev interpolation polynomial.
type ChebyshevInterpolation struct {
	Poly
	a float64
	b float64
}

// Approximate computes a Chebyshev approximation of the input function, for the range [-a, b] of degree degree.
// To be used in conjunction with the function EvaluateCheby.
func Approximate(function func(float64) float64, a, b float64, degree int) (cheby *ChebyshevInterpolation) {

	cheby = new(ChebyshevInterpolation)
	cheby.a = a
	cheby.b = b
	cheby.maxDeg = uint64(degree)
	cheby.lead = true

	nodes := chebyshevNodes(degree+1, a, b)

	fi := make([]float64, len(nodes))
	for i := range nodes {
		fi[i] = function(nodes[i])
	}

	cheby.coeffs = chebyCoeffs(nodes, fi, a, b)

	return
}

func chebyshevNodes(n int, a, b float64) (u []float64) {
	u = make([]float64, n)
	var x, y float64
	for k := 1; k < n+1; k++ {
		x = 0.5 * (a + b)
		y = 0.5 * (b - a)
		u[k-1] = x + y*math.Cos((float64(k)-0.5)*(3.141592653589793/float64(n)))
	}
	return
}

func evaluateChebyshevPolynomial(coeffs []float64, x float64, a, b float64) (y float64) {
	var u, Tprev, Tnext, T float64
	u = (2*x - a - b) / (b - a)
	Tprev = 1
	T = u
	y = coeffs[0]
	for i := 1; i < len(coeffs); i++ {
		y = y + T*coeffs[i]
		Tnext = 2*u*T - Tprev
		Tprev = T
		T = Tnext
	}
	return
}

func chebyCoeffs(nodes, fi []float64, a, b float64) (coeffs []float64) {

	var u, Tprev, T, Tnext float64

	n := len(nodes)

	coeffs = make([]float64, n)

	for i := 0; i < n; i++ {

		u = (2*nodes[i] - a - b) / (b - a)
		Tprev = 1
		T = u

		for j := 0; j < n; j++ {
			coeffs[j] += fi[i] * Tprev
			Tnext = 2*u*T - Tprev
			Tprev = T
			T = Tnext
		}
	}

	coeffs[0] /= float64(n)
	for i := 1; i < n; i++ {
		coeffs[i] *= (2.0 / float64(n))
	}

	return
}
