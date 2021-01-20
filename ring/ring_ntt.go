package ring

import (
	"math/bits"
	"unsafe"
)

// NTT computes the NTT of p1 and returns the result on p2.
func (r *Ring) NTT(p1, p2 *Poly) {
	for x := range r.Modulus {
		NTT(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsi[x], r.Modulus[x], r.BredParams[x])
	}
}

// NTTLvl computes the NTT of p1 and returns the result on p2.
// The value level defines the number of moduli of the input polynomials.
func (r *Ring) NTTLvl(level uint64, p1, p2 *Poly) {
	for x := uint64(0); x < level+1; x++ {
		NTT(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsi[x], r.Modulus[x], r.BredParams[x])
	}
}

// NTTLazy computes the NTT of p1 and returns the result on p2.
func (r *Ring) NTTLazy(p1, p2 *Poly) {
	for x := range r.Modulus {
		NTTLazy(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsi[x], r.Modulus[x])
	}
}

// NTTLazyLvl computes the NTT of p1 and returns the result on p2.
// The value level defines the number of moduli of the input polynomials.
func (r *Ring) NTTLazyLvl(level uint64, p1, p2 *Poly) {
	for x := uint64(0); x < level+1; x++ {
		NTTLazy(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsi[x], r.Modulus[x])
	}
}

// butterfly computes X, Y = U + V*Psi, U - V*Psi mod Q.
func butterfly(U, V uint64, Psi FastBRedOperand, twoQ, fourQ, Q uint64) (uint64, uint64) {

	if U >= fourQ {
		U -= fourQ
	}

	V = FastBRed(V, Psi, Q)

	return U + V, U + twoQ - V
}

// NTT computes the NTT on the input coefficients using the input parameters.
func NTT(coeffsIn, coeffsOut []uint64, N uint64, nttPsi []FastBRedOperand, Q uint64, bredParams []uint64) {

	NTTLazy(coeffsIn, coeffsOut, N, nttPsi, Q)
	// Finish with an exact reduction
	for i := uint64(0); i < N; i = i + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&coeffsOut[i]))

		x[0] = BRedAdd(x[0], Q, bredParams)
		x[1] = BRedAdd(x[1], Q, bredParams)
		x[2] = BRedAdd(x[2], Q, bredParams)
		x[3] = BRedAdd(x[3], Q, bredParams)
		x[4] = BRedAdd(x[4], Q, bredParams)
		x[5] = BRedAdd(x[5], Q, bredParams)
		x[6] = BRedAdd(x[6], Q, bredParams)
		x[7] = BRedAdd(x[7], Q, bredParams)
	}
}

// NTTLazy computes the NTT on the input coefficients using the input parameters with output values in the range [0, 2q-1].
func NTTLazy(coeffsIn, coeffsOut []uint64, N uint64, nttPsi []FastBRedOperand, Q uint64) {
	var j1, j2, t uint64
	var F FastBRedOperand

	fourQ := 4 * Q
	twoQ := 2 * Q

	// Copy the result of the first round of butterflies on p2 with approximate reduction
	t = N >> 1
	F = nttPsi[1]

	var V uint64

	for j := uint64(0); j <= t-1; j = j + 8 {

		xin := (*[8]uint64)(unsafe.Pointer(&coeffsIn[j]))
		yin := (*[8]uint64)(unsafe.Pointer(&coeffsIn[j+t]))

		xout := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j]))
		yout := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j+t]))

		V = FastBRed(yin[0], F, Q)
		xout[0], yout[0] = xin[0]+V, xin[0]+twoQ-V

		V = FastBRed(yin[1], F, Q)
		xout[1], yout[1] = xin[1]+V, xin[1]+twoQ-V

		V = FastBRed(yin[2], F, Q)
		xout[2], yout[2] = xin[2]+V, xin[2]+twoQ-V

		V = FastBRed(yin[3], F, Q)
		xout[3], yout[3] = xin[3]+V, xin[3]+twoQ-V

		V = FastBRed(yin[4], F, Q)
		xout[4], yout[4] = xin[4]+V, xin[4]+twoQ-V

		V = FastBRed(yin[5], F, Q)
		xout[5], yout[5] = xin[5]+V, xin[5]+twoQ-V

		V = FastBRed(yin[6], F, Q)
		xout[6], yout[6] = xin[6]+V, xin[6]+twoQ-V

		V = FastBRed(yin[7], F, Q)
		xout[7], yout[7] = xin[7]+V, xin[7]+twoQ-V
	}

	// Continue the rest of the second to the n-1 butterflies on p2 with approximate reduction
	var reduce bool

	for m := uint64(2); m < N; m <<= 1 {

		reduce = (bits.Len64(m)&1 == 1)

		t >>= 1

		if t >= 8 {

			for i := uint64(0); i < m; i++ {

				j1 = (i * t) << 1

				j2 = j1 + t - 1

				F = nttPsi[m+i]

				if reduce {

					for j := j1; j <= j2; j = j + 8 {

						x := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j]))
						y := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j+t]))

						x[0], y[0] = butterfly(x[0], y[0], F, twoQ, fourQ, Q)
						x[1], y[1] = butterfly(x[1], y[1], F, twoQ, fourQ, Q)
						x[2], y[2] = butterfly(x[2], y[2], F, twoQ, fourQ, Q)
						x[3], y[3] = butterfly(x[3], y[3], F, twoQ, fourQ, Q)
						x[4], y[4] = butterfly(x[4], y[4], F, twoQ, fourQ, Q)
						x[5], y[5] = butterfly(x[5], y[5], F, twoQ, fourQ, Q)
						x[6], y[6] = butterfly(x[6], y[6], F, twoQ, fourQ, Q)
						x[7], y[7] = butterfly(x[7], y[7], F, twoQ, fourQ, Q)
					}

				} else {

					for j := j1; j <= j2; j = j + 8 {

						x := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j]))
						y := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j+t]))

						V = FastBRedConstant(y[0], F, Q)
						x[0], y[0] = x[0]+V, x[0]+twoQ-V

						V = FastBRedConstant(y[1], F, Q)
						x[1], y[1] = x[1]+V, x[1]+twoQ-V

						V = FastBRedConstant(y[2], F, Q)
						x[2], y[2] = x[2]+V, x[2]+twoQ-V

						V = FastBRedConstant(y[3], F, Q)
						x[3], y[3] = x[3]+V, x[3]+twoQ-V

						V = FastBRedConstant(y[4], F, Q)
						x[4], y[4] = x[4]+V, x[4]+twoQ-V

						V = FastBRedConstant(y[5], F, Q)
						x[5], y[5] = x[5]+V, x[5]+twoQ-V

						V = FastBRedConstant(y[6], F, Q)
						x[6], y[6] = x[6]+V, x[6]+twoQ-V

						V = FastBRedConstant(y[7], F, Q)
						x[7], y[7] = x[7]+V, x[7]+twoQ-V
					}
				}
			}

		} else if t == 4 {

			if reduce {

				for i := uint64(0); i < m; i = i + 2 {

					j1 = (i * t) << 1

					psi := (*[2]FastBRedOperand)(unsafe.Pointer(&nttPsi[m+i]))
					x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

					x[0], x[4] = butterfly(x[0], x[4], psi[0], twoQ, fourQ, Q)
					x[1], x[5] = butterfly(x[1], x[5], psi[0], twoQ, fourQ, Q)
					x[2], x[6] = butterfly(x[2], x[6], psi[0], twoQ, fourQ, Q)
					x[3], x[7] = butterfly(x[3], x[7], psi[0], twoQ, fourQ, Q)
					x[8], x[12] = butterfly(x[8], x[12], psi[1], twoQ, fourQ, Q)
					x[9], x[13] = butterfly(x[9], x[13], psi[1], twoQ, fourQ, Q)
					x[10], x[14] = butterfly(x[10], x[14], psi[1], twoQ, fourQ, Q)
					x[11], x[15] = butterfly(x[11], x[15], psi[1], twoQ, fourQ, Q)

				}
			} else {

				for i := uint64(0); i < m; i = i + 2 {

					j1 = (i * t) << 1

					psi := (*[2]FastBRedOperand)(unsafe.Pointer(&nttPsi[m+i]))
					x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

					V = FastBRedConstant(x[4], psi[0], Q)
					x[0], x[4] = x[0]+V, x[0]+twoQ-V

					V = FastBRedConstant(x[5], psi[0], Q)
					x[1], x[5] = x[1]+V, x[1]+twoQ-V

					V = FastBRedConstant(x[6], psi[0], Q)
					x[2], x[6] = x[2]+V, x[2]+twoQ-V

					V = FastBRedConstant(x[7], psi[0], Q)
					x[3], x[7] = x[3]+V, x[3]+twoQ-V

					V = FastBRedConstant(x[12], psi[1], Q)
					x[8], x[12] = x[8]+V, x[8]+twoQ-V

					V = FastBRedConstant(x[13], psi[1], Q)
					x[9], x[13] = x[9]+V, x[9]+twoQ-V

					V = FastBRedConstant(x[14], psi[1], Q)
					x[10], x[14] = x[10]+V, x[10]+twoQ-V

					V = FastBRedConstant(x[15], psi[1], Q)
					x[11], x[15] = x[11]+V, x[11]+twoQ-V

				}

			}

		} else if t == 2 {

			if reduce {

				for i := uint64(0); i < m; i = i + 4 {

					j1 = (i * t) << 1

					psi := (*[4]FastBRedOperand)(unsafe.Pointer(&nttPsi[m+i]))
					x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

					x[0], x[2] = butterfly(x[0], x[2], psi[0], twoQ, fourQ, Q)
					x[1], x[3] = butterfly(x[1], x[3], psi[0], twoQ, fourQ, Q)
					x[4], x[6] = butterfly(x[4], x[6], psi[1], twoQ, fourQ, Q)
					x[5], x[7] = butterfly(x[5], x[7], psi[1], twoQ, fourQ, Q)
					x[8], x[10] = butterfly(x[8], x[10], psi[2], twoQ, fourQ, Q)
					x[9], x[11] = butterfly(x[9], x[11], psi[2], twoQ, fourQ, Q)
					x[12], x[14] = butterfly(x[12], x[14], psi[3], twoQ, fourQ, Q)
					x[13], x[15] = butterfly(x[13], x[15], psi[3], twoQ, fourQ, Q)
				}
			} else {

				for i := uint64(0); i < m; i = i + 4 {

					j1 = (i * t) << 1

					psi := (*[4]FastBRedOperand)(unsafe.Pointer(&nttPsi[m+i]))
					x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

					V = FastBRedConstant(x[2], psi[0], Q)
					x[0], x[2] = x[0]+V, x[0]+twoQ-V

					V = FastBRedConstant(x[3], psi[0], Q)
					x[1], x[3] = x[1]+V, x[1]+twoQ-V

					V = FastBRedConstant(x[6], psi[1], Q)
					x[4], x[6] = x[4]+V, x[4]+twoQ-V

					V = FastBRedConstant(x[7], psi[1], Q)
					x[5], x[7] = x[5]+V, x[5]+twoQ-V

					V = FastBRedConstant(x[10], psi[2], Q)
					x[8], x[10] = x[8]+V, x[8]+twoQ-V

					V = FastBRedConstant(x[11], psi[2], Q)
					x[9], x[11] = x[9]+V, x[9]+twoQ-V

					V = FastBRedConstant(x[14], psi[3], Q)
					x[12], x[14] = x[12]+V, x[12]+twoQ-V

					V = FastBRedConstant(x[15], psi[3], Q)
					x[13], x[15] = x[13]+V, x[13]+twoQ-V
				}
			}

		} else {

			for i := uint64(0); i < m; i = i + 8 {

				psi := (*[8]FastBRedOperand)(unsafe.Pointer(&nttPsi[m+i]))
				x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[2*i]))

				x[0], x[1] = butterfly(x[0], x[1], psi[0], twoQ, fourQ, Q)
				x[2], x[3] = butterfly(x[2], x[3], psi[1], twoQ, fourQ, Q)
				x[4], x[5] = butterfly(x[4], x[5], psi[2], twoQ, fourQ, Q)
				x[6], x[7] = butterfly(x[6], x[7], psi[3], twoQ, fourQ, Q)
				x[8], x[9] = butterfly(x[8], x[9], psi[4], twoQ, fourQ, Q)
				x[10], x[11] = butterfly(x[10], x[11], psi[5], twoQ, fourQ, Q)
				x[12], x[13] = butterfly(x[12], x[13], psi[6], twoQ, fourQ, Q)
				x[14], x[15] = butterfly(x[14], x[15], psi[7], twoQ, fourQ, Q)
			}

			/*
				for i := uint64(0); i < m; i = i + 8 {

					psi := (*[8]FastBRedOperand)(unsafe.Pointer(&nttPsi[m+i]))
					x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[2*i]))

					V = FastBRedConstant(x[1], psi[0], Q)
					x[0], x[1] = x[0]+V, x[0]+twoQ-V

					V = FastBRedConstant(x[3], psi[1], Q)
					x[2], x[3] = x[2]+V, x[2]+twoQ-V

					V = FastBRedConstant(x[5], psi[2], Q)
					x[4], x[5] = x[4]+V, x[4]+twoQ-V

					V = FastBRedConstant(x[7], psi[3], Q)
					x[6], x[7] = x[6]+V, x[6]+twoQ-V

					V = FastBRedConstant(x[9], psi[4], Q)
					x[8], x[9] = x[8]+V, x[8]+twoQ-V

					V = FastBRedConstant(x[11], psi[5], Q)
					x[10], x[11] = x[10]+V, x[10]+twoQ-V

					V = FastBRedConstant(x[13], psi[6], Q)
					x[12], x[13] = x[12]+V, x[12]+twoQ-V

					V = FastBRedConstant(x[15], psi[7], Q)
					x[14], x[15] = x[14]+V, x[14]+twoQ-V
				}
			*/
		}

	}
}

// InvNTT computes the inverse-NTT of p1 and returns the result on p2.
func (r *Ring) InvNTT(p1, p2 *Poly) {
	for x := range r.Modulus {
		InvNTT(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsiInv[x], r.NttNInv[x], r.Modulus[x])
	}
}

// InvNTTLvl computes the inverse-NTT of p1 and returns the result on p2.
// The value level defines the number of moduli of the input polynomials.
func (r *Ring) InvNTTLvl(level uint64, p1, p2 *Poly) {
	for x := uint64(0); x < level+1; x++ {
		InvNTT(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsiInv[x], r.NttNInv[x], r.Modulus[x])
	}
}

// invbutterfly computes X, Y = U + V, (U - V) * Psi mod Q.
func invbutterfly(U, V uint64, Psi FastBRedOperand, twoQ, fourQ, Q uint64) (X, Y uint64) {
	X = U + V
	if X >= twoQ {
		X -= twoQ
	}
	Y = FastBRedConstant(U+fourQ-V, Psi, Q) // At the moment it is not possible to use FastBRedConstant if Q > 61 bits
	return
}

// InvNTT computes the InvNTT transformation on the input coefficients using the input parameters.
func InvNTT(coeffsIn, coeffsOut []uint64, N uint64, nttPsiInv []FastBRedOperand, nttNInv FastBRedOperand, Q uint64) {

	var j1, j2, h, t uint64
	var F FastBRedOperand

	// Copy the result of the first round of butterflies on p2 with approximate reduction
	t = 1
	h = N >> 1
	twoQ := Q << 1
	fourQ := Q << 2

	for i := uint64(0); i < h; i = i + 8 {

		psi := (*[8]FastBRedOperand)(unsafe.Pointer(&nttPsiInv[h+i]))
		xin := (*[16]uint64)(unsafe.Pointer(&coeffsIn[2*i]))
		xout := (*[16]uint64)(unsafe.Pointer(&coeffsOut[2*i]))

		xout[0], xout[1] = invbutterfly(xin[0], xin[1], psi[0], twoQ, fourQ, Q)
		xout[2], xout[3] = invbutterfly(xin[2], xin[3], psi[1], twoQ, fourQ, Q)
		xout[4], xout[5] = invbutterfly(xin[4], xin[5], psi[2], twoQ, fourQ, Q)
		xout[6], xout[7] = invbutterfly(xin[6], xin[7], psi[3], twoQ, fourQ, Q)
		xout[8], xout[9] = invbutterfly(xin[8], xin[9], psi[4], twoQ, fourQ, Q)
		xout[10], xout[11] = invbutterfly(xin[10], xin[11], psi[5], twoQ, fourQ, Q)
		xout[12], xout[13] = invbutterfly(xin[12], xin[13], psi[6], twoQ, fourQ, Q)
		xout[14], xout[15] = invbutterfly(xin[14], xin[15], psi[7], twoQ, fourQ, Q)
	}

	// Continue the rest of the second to the n-1 butterflies on p2 with approximate reduction
	t <<= 1
	for m := N >> 1; m > 1; m >>= 1 {

		j1 = 0
		h = m >> 1

		if t >= 8 {

			for i := uint64(0); i < h; i++ {

				j2 = j1 + t - 1

				F = nttPsiInv[h+i]

				for j := j1; j <= j2; j = j + 8 {

					x := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j]))
					y := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j+t]))

					x[0], y[0] = invbutterfly(x[0], y[0], F, twoQ, fourQ, Q)
					x[1], y[1] = invbutterfly(x[1], y[1], F, twoQ, fourQ, Q)
					x[2], y[2] = invbutterfly(x[2], y[2], F, twoQ, fourQ, Q)
					x[3], y[3] = invbutterfly(x[3], y[3], F, twoQ, fourQ, Q)
					x[4], y[4] = invbutterfly(x[4], y[4], F, twoQ, fourQ, Q)
					x[5], y[5] = invbutterfly(x[5], y[5], F, twoQ, fourQ, Q)
					x[6], y[6] = invbutterfly(x[6], y[6], F, twoQ, fourQ, Q)
					x[7], y[7] = invbutterfly(x[7], y[7], F, twoQ, fourQ, Q)
				}

				j1 = j1 + (t << 1)
			}

		} else if t == 4 {

			for i := uint64(0); i < h; i = i + 2 {

				psi := (*[2]FastBRedOperand)(unsafe.Pointer(&nttPsiInv[h+i]))
				x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

				x[0], x[4] = invbutterfly(x[0], x[4], psi[0], twoQ, fourQ, Q)
				x[1], x[5] = invbutterfly(x[1], x[5], psi[0], twoQ, fourQ, Q)
				x[2], x[6] = invbutterfly(x[2], x[6], psi[0], twoQ, fourQ, Q)
				x[3], x[7] = invbutterfly(x[3], x[7], psi[0], twoQ, fourQ, Q)
				x[8], x[12] = invbutterfly(x[8], x[12], psi[1], twoQ, fourQ, Q)
				x[9], x[13] = invbutterfly(x[9], x[13], psi[1], twoQ, fourQ, Q)
				x[10], x[14] = invbutterfly(x[10], x[14], psi[1], twoQ, fourQ, Q)
				x[11], x[15] = invbutterfly(x[11], x[15], psi[1], twoQ, fourQ, Q)

				j1 = j1 + (t << 2)
			}

		} else {

			for i := uint64(0); i < h; i = i + 4 {

				psi := (*[4]FastBRedOperand)(unsafe.Pointer(&nttPsiInv[h+i]))
				x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

				x[0], x[2] = invbutterfly(x[0], x[2], psi[0], twoQ, fourQ, Q)
				x[1], x[3] = invbutterfly(x[1], x[3], psi[0], twoQ, fourQ, Q)
				x[4], x[6] = invbutterfly(x[4], x[6], psi[1], twoQ, fourQ, Q)
				x[5], x[7] = invbutterfly(x[5], x[7], psi[1], twoQ, fourQ, Q)
				x[8], x[10] = invbutterfly(x[8], x[10], psi[2], twoQ, fourQ, Q)
				x[9], x[11] = invbutterfly(x[9], x[11], psi[2], twoQ, fourQ, Q)
				x[12], x[14] = invbutterfly(x[12], x[14], psi[3], twoQ, fourQ, Q)
				x[13], x[15] = invbutterfly(x[13], x[15], psi[3], twoQ, fourQ, Q)

				j1 = j1 + (t << 3)
			}
		}

		t <<= 1
	}

	// Finish with an exact reduction
	for i := uint64(0); i < N; i = i + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&coeffsOut[i]))

		x[0] = FastBRed(x[0], nttNInv, Q)
		x[1] = FastBRed(x[1], nttNInv, Q)
		x[2] = FastBRed(x[2], nttNInv, Q)
		x[3] = FastBRed(x[3], nttNInv, Q)
		x[4] = FastBRed(x[4], nttNInv, Q)
		x[5] = FastBRed(x[5], nttNInv, Q)
		x[6] = FastBRed(x[6], nttNInv, Q)
		x[7] = FastBRed(x[7], nttNInv, Q)
	}
}

// InvNTTLazy computes the inverse-NTT of p1 and returns the result on p2.
// Output values are in the range [0, 2q-1]
func (r *Ring) InvNTTLazy(p1, p2 *Poly) {
	for x := range r.Modulus {
		InvNTTLazy(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsiInv[x], r.NttNInv[x], r.Modulus[x])
	}
}

// InvNTTLazyLvl computes the inverse-NTT of p1 and returns the result on p2.
// The value level defines the number of moduli of the input polynomials.
// Output values are in the range [0, 2q-1]
func (r *Ring) InvNTTLazyLvl(level uint64, p1, p2 *Poly) {
	for x := uint64(0); x < level+1; x++ {
		InvNTTLazy(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsiInv[x], r.NttNInv[x], r.Modulus[x])
	}
}

// InvNTTLazy computes the InvNTT transformation on the input coefficients using the input parameters with output values in the range [0, 2q-1].
func InvNTTLazy(coeffsIn, coeffsOut []uint64, N uint64, nttPsiInv []FastBRedOperand, nttNInv FastBRedOperand, Q uint64) {

	var j1, j2, h, t uint64
	var F FastBRedOperand

	// Copy the result of the first round of butterflies on p2 with approximate reduction
	t = 1
	h = N >> 1

	twoQ := Q << 1
	fourQ := Q << 2

	for i := uint64(0); i < h; i = i + 8 {

		psi := (*[8]FastBRedOperand)(unsafe.Pointer(&nttPsiInv[h+i]))
		xin := (*[16]uint64)(unsafe.Pointer(&coeffsIn[2*i]))
		xout := (*[16]uint64)(unsafe.Pointer(&coeffsOut[2*i]))

		xout[0], xout[1] = invbutterfly(xin[0], xin[1], psi[0], twoQ, fourQ, Q)
		xout[2], xout[3] = invbutterfly(xin[2], xin[3], psi[1], twoQ, fourQ, Q)
		xout[4], xout[5] = invbutterfly(xin[4], xin[5], psi[2], twoQ, fourQ, Q)
		xout[6], xout[7] = invbutterfly(xin[6], xin[7], psi[3], twoQ, fourQ, Q)
		xout[8], xout[9] = invbutterfly(xin[8], xin[9], psi[4], twoQ, fourQ, Q)
		xout[10], xout[11] = invbutterfly(xin[10], xin[11], psi[5], twoQ, fourQ, Q)
		xout[12], xout[13] = invbutterfly(xin[12], xin[13], psi[6], twoQ, fourQ, Q)
		xout[14], xout[15] = invbutterfly(xin[14], xin[15], psi[7], twoQ, fourQ, Q)
	}

	// Continue the rest of the second to the n-1 butterflies on p2 with approximate reduction
	t <<= 1
	for m := N >> 1; m > 1; m >>= 1 {

		j1 = 0
		h = m >> 1

		if t >= 8 {

			for i := uint64(0); i < h; i++ {

				j2 = j1 + t - 1

				F = nttPsiInv[h+i]

				for j := j1; j <= j2; j = j + 8 {

					x := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j]))
					y := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j+t]))

					x[0], y[0] = invbutterfly(x[0], y[0], F, twoQ, fourQ, Q)
					x[1], y[1] = invbutterfly(x[1], y[1], F, twoQ, fourQ, Q)
					x[2], y[2] = invbutterfly(x[2], y[2], F, twoQ, fourQ, Q)
					x[3], y[3] = invbutterfly(x[3], y[3], F, twoQ, fourQ, Q)
					x[4], y[4] = invbutterfly(x[4], y[4], F, twoQ, fourQ, Q)
					x[5], y[5] = invbutterfly(x[5], y[5], F, twoQ, fourQ, Q)
					x[6], y[6] = invbutterfly(x[6], y[6], F, twoQ, fourQ, Q)
					x[7], y[7] = invbutterfly(x[7], y[7], F, twoQ, fourQ, Q)
				}

				j1 = j1 + (t << 1)
			}

		} else if t == 4 {

			for i := uint64(0); i < h; i = i + 2 {

				psi := (*[2]FastBRedOperand)(unsafe.Pointer(&nttPsiInv[h+i]))
				x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

				x[0], x[4] = invbutterfly(x[0], x[4], psi[0], twoQ, fourQ, Q)
				x[1], x[5] = invbutterfly(x[1], x[5], psi[0], twoQ, fourQ, Q)
				x[2], x[6] = invbutterfly(x[2], x[6], psi[0], twoQ, fourQ, Q)
				x[3], x[7] = invbutterfly(x[3], x[7], psi[0], twoQ, fourQ, Q)
				x[8], x[12] = invbutterfly(x[8], x[12], psi[1], twoQ, fourQ, Q)
				x[9], x[13] = invbutterfly(x[9], x[13], psi[1], twoQ, fourQ, Q)
				x[10], x[14] = invbutterfly(x[10], x[14], psi[1], twoQ, fourQ, Q)
				x[11], x[15] = invbutterfly(x[11], x[15], psi[1], twoQ, fourQ, Q)

				j1 = j1 + (t << 2)
			}

		} else {

			for i := uint64(0); i < h; i = i + 4 {

				psi := (*[4]FastBRedOperand)(unsafe.Pointer(&nttPsiInv[h+i]))
				x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

				x[0], x[2] = invbutterfly(x[0], x[2], psi[0], twoQ, fourQ, Q)
				x[1], x[3] = invbutterfly(x[1], x[3], psi[0], twoQ, fourQ, Q)
				x[4], x[6] = invbutterfly(x[4], x[6], psi[1], twoQ, fourQ, Q)
				x[5], x[7] = invbutterfly(x[5], x[7], psi[1], twoQ, fourQ, Q)
				x[8], x[10] = invbutterfly(x[8], x[10], psi[2], twoQ, fourQ, Q)
				x[9], x[11] = invbutterfly(x[9], x[11], psi[2], twoQ, fourQ, Q)
				x[12], x[14] = invbutterfly(x[12], x[14], psi[3], twoQ, fourQ, Q)
				x[13], x[15] = invbutterfly(x[13], x[15], psi[3], twoQ, fourQ, Q)

				j1 = j1 + (t << 3)
			}
		}

		t <<= 1
	}

	// Finish with an exact reduction
	for i := uint64(0); i < N; i = i + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&coeffsOut[i]))

		x[0] = FastBRedConstant(x[0], nttNInv, Q)
		x[1] = FastBRedConstant(x[1], nttNInv, Q)
		x[2] = FastBRedConstant(x[2], nttNInv, Q)
		x[3] = FastBRedConstant(x[3], nttNInv, Q)
		x[4] = FastBRedConstant(x[4], nttNInv, Q)
		x[5] = FastBRedConstant(x[5], nttNInv, Q)
		x[6] = FastBRedConstant(x[6], nttNInv, Q)
		x[7] = FastBRedConstant(x[7], nttNInv, Q)
	}
}
