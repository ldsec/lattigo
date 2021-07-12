package dckks

import (
	"github.com/ldsec/lattigo/v2/ring"
	"math"
	"math/bits"
)

func extendBasisSmallNormAndCenter(ringQ, ringP *ring.Ring, polQ, polP *ring.Poly) {
	var coeff, Q, QHalf, sign uint64
	Q = ringQ.Modulus[0]
	QHalf = Q >> 1

	for j := 0; j < ringQ.N; j++ {

		coeff = polQ.Coeffs[0][j]

		sign = 1
		if coeff > QHalf {
			coeff = Q - coeff
			sign = 0
		}

		for i, pi := range ringP.Modulus {
			polP.Coeffs[i][j] = (coeff * sign) | (pi-coeff)*(sign^1)
		}
	}
}

// GetMinimumLevelForBootstrapping takes the security parameter lambda, the ciphertext scale, the number of parties and the moduli chain
// and returns the minimum level at which the collective refresh can be called with a security of at least 128-bits.
// It returns 3 parameters :
// minLevel : the minimum level at which the collective refresh must be called to ensure correctness
// logBound : the bit length of the masks to be sampled to mask the plaintext and ensure 128-bits of statistical indistinguishability
// ok 		: a boolean flag, which is set to false if no such instance exist
func GetMinimumLevelForBootstrapping(lambda int, scale float64, nParties int, moduli []uint64) (minLevel, logBound int, ok bool) {
	logBound = lambda + int(math.Ceil(math.Log2(scale)))
	maxBound := logBound + bits.Len64(uint64(nParties))
	minLevel = -1
	logQ := 0
	for i := 0; logQ < maxBound; i++ {
		if i >= len(moduli) {
			return 0, 0, false
		}

		logQ += bits.Len64(moduli[i])
		minLevel++
	}
	if len(moduli) < minLevel {
		return 0, 0, false
	}

	return minLevel, logBound, true
}
