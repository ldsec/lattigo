package ring

import (
	"math/bits"
	"unsafe"

	"github.com/ldsec/lattigo/v2/utils"
)

// GenGaloisParams generates the generators for the Galois endomorphisms.
func GenGaloisParams(n, gen uint64) (galElRotCol []uint64) {

	var m, mask uint64

	m = n << 1

	mask = m - 1

	galElRotCol = make([]uint64, n>>1)

	galElRotCol[0] = 1

	for i := uint64(1); i < n>>1; i++ {
		galElRotCol[i] = (galElRotCol[i-1] * gen) & mask
	}

	return
}

// PermuteNTTIndex computes the index table for PermuteNTT.
func PermuteNTTIndex(gen, power, N, NthRoot uint64) (index []uint64) {

	genPow := ModExp(gen, power, NthRoot)

	var mask, logNthRoot, tmp1, tmp2 uint64

	logNthRoot = uint64(bits.Len64(NthRoot) - 2)

	mask = NthRoot - 1

	index = make([]uint64, N)

	for i := uint64(0); i < N; i++ {
		tmp1 = 2*utils.BitReverse64(i, logNthRoot) + 1

		tmp2 = ((genPow * tmp1 & mask) - 1) >> 1

		index[i] = utils.BitReverse64(tmp2, logNthRoot)
	}

	return
}

// PermuteNTT applies the Galois transform on a polynomial in the NTT domain.
// It maps the coefficients x^i to x^(gen*i)
// It must be noted that the result cannot be in-place.
func PermuteNTT(polIn *Poly, gen, power, N, NthRoot uint64, polOut *Poly) {

	var tmp uint64

	index := PermuteNTTIndex(gen, power, N, NthRoot)

	for j := uint64(0); j < N; j++ {

		tmp = index[j]

		for i := 0; i < len(polIn.Coeffs); i++ {

			polOut.Coeffs[i][j] = polIn.Coeffs[i][tmp]
		}
	}
}

// PermuteNTTLvl applies the Galois transform on a polynomial in the NTT domain, up to a given level.
// It maps the coefficients x^i to x^(gen*i)
// It must be noted that the result cannot be in-place.
func PermuteNTTLvl(level uint64, polIn *Poly, gen uint64, polOut *Poly) {

	var N, tmp, mask, logN, tmp1, tmp2 uint64

	N = uint64(len(polIn.Coeffs[0]))

	logN = uint64(bits.Len64(N) - 1)

	mask = (N << 1) - 1

	index := make([]uint64, N)

	for i := uint64(0); i < N; i++ {
		tmp1 = 2*utils.BitReverse64(i, logN) + 1

		tmp2 = ((gen * tmp1 & mask) - 1) >> 1

		index[i] = utils.BitReverse64(tmp2, logN)
	}

	for j := uint64(0); j < N; j++ {

		tmp = index[j]

		for i := uint64(0); i < level+1; i++ {

			polOut.Coeffs[i][j] = polIn.Coeffs[i][tmp]
		}
	}
}

// PermuteNTTWithIndexLvl applies the Galois transform on a polynomial in the NTT domain, up to a given level.
// It maps the coefficients x^i to x^(gen*i) using the PermuteNTTIndex table.
// It must be noted that the result cannot be in-place.
func PermuteNTTWithIndexLvl(level uint64, polIn *Poly, index []uint64, polOut *Poly) {

	for j := uint64(0); j < uint64(len(polIn.Coeffs[0])); j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&index[j]))

		for i := uint64(0); i < level+1; i++ {

			z := (*[8]uint64)(unsafe.Pointer(&polOut.Coeffs[i][j]))
			y := polIn.Coeffs[i]

			z[0] = y[x[0]]
			z[1] = y[x[1]]
			z[2] = y[x[2]]
			z[3] = y[x[3]]
			z[4] = y[x[4]]
			z[5] = y[x[5]]
			z[6] = y[x[6]]
			z[7] = y[x[7]]
		}
	}
}

// PermuteNTTWithIndexAndAddNoModLvl applies the Galois transform on a polynomial in the NTT domain, up to a given level,
// and adds the result to the output polynomial without modular reduction.
// It maps the coefficients x^i to x^(gen*i) using the PermuteNTTIndex table.
// It must be noted that the result cannot be in-place.
func PermuteNTTWithIndexAndAddNoModLvl(level uint64, polIn *Poly, index []uint64, polOut *Poly) {

	for j := uint64(0); j < uint64(len(polIn.Coeffs[0])); j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&index[j]))

		for i := uint64(0); i < level+1; i++ {

			z := (*[8]uint64)(unsafe.Pointer(&polOut.Coeffs[i][j]))
			y := polIn.Coeffs[i]

			z[0] += y[x[0]]
			z[1] += y[x[1]]
			z[2] += y[x[2]]
			z[3] += y[x[3]]
			z[4] += y[x[4]]
			z[5] += y[x[5]]
			z[6] += y[x[6]]
			z[7] += y[x[7]]
		}
	}
}

// Permute applies the Galois transform on a polynomial outside of the NTT domain.
// It maps the coefficients x^i to x^(gen*i)
// It must be noted that the result cannot be in-place.
func (r *Ring) Permute(polIn *Poly, gen uint64, polOut *Poly) {

	var mask, index, indexRaw, logN, tmp uint64

	mask = r.N - 1

	logN = uint64(bits.Len64(mask))

	for i := uint64(0); i < r.N; i++ {

		indexRaw = i * gen

		index = indexRaw & mask

		tmp = (indexRaw >> logN) & 1

		for j, qi := range r.Modulus {

			polOut.Coeffs[j][index] = polIn.Coeffs[j][i]*(tmp^1) | (qi-polIn.Coeffs[j][i])*tmp
		}
	}
}
