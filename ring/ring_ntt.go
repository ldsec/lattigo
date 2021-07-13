package ring

import (
	"math/bits"
	"unsafe"
)

// NTT computes the NTT of p1 and returns the result on p2.
func (r *Ring) NTT(p1, p2 *Poly) {
	for x := range r.Modulus {
		NTT(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsi[x], r.Modulus[x], r.MredParams[x], r.BredParams[x])
	}
}

// NTTLvl computes the NTT of p1 and returns the result on p2.
// The value level defines the number of moduli of the input polynomials.
func (r *Ring) NTTLvl(level int, p1, p2 *Poly) {
	for x := 0; x < level+1; x++ {
		NTT(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsi[x], r.Modulus[x], r.MredParams[x], r.BredParams[x])
	}
}

// NTTLazy computes the NTT of p1 and returns the result on p2.
// Output values are in the range [0, 2q-1]
func (r *Ring) NTTLazy(p1, p2 *Poly) {
	for x := range r.Modulus {
		NTTLazy(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsi[x], r.Modulus[x], r.MredParams[x], r.BredParams[x])
	}
}

// NTTLazyLvl computes the NTT of p1 and returns the result on p2.
// The value level defines the number of moduli of the input polynomials.
// Output values are in the range [0, 2q-1]
func (r *Ring) NTTLazyLvl(level int, p1, p2 *Poly) {
	for x := 0; x < level+1; x++ {
		NTTLazy(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsi[x], r.Modulus[x], r.MredParams[x], r.BredParams[x])
	}
}

// InvNTT computes the inverse-NTT of p1 and returns the result on p2.
func (r *Ring) InvNTT(p1, p2 *Poly) {
	for x := range r.Modulus {
		InvNTT(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsiInv[x], r.NttNInv[x], r.Modulus[x], r.MredParams[x])
	}
}

// InvNTTLvl computes the inverse-NTT of p1 and returns the result on p2.
// The value level defines the number of moduli of the input polynomials.
func (r *Ring) InvNTTLvl(level int, p1, p2 *Poly) {
	for x := 0; x < level+1; x++ {
		InvNTT(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsiInv[x], r.NttNInv[x], r.Modulus[x], r.MredParams[x])
	}
}

// InvNTTLazy computes the inverse-NTT of p1 and returns the result on p2.
// Output values are in the range [0, 2q-1]
func (r *Ring) InvNTTLazy(p1, p2 *Poly) {
	for x := range r.Modulus {
		InvNTTLazy(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsiInv[x], r.NttNInv[x], r.Modulus[x], r.MredParams[x])
	}
}

// InvNTTLazyLvl computes the inverse-NTT of p1 and returns the result on p2.
// The value level defines the number of moduli of the input polynomials.
// Output values are in the range [0, 2q-1]
func (r *Ring) InvNTTLazyLvl(level int, p1, p2 *Poly) {
	for x := 0; x < level+1; x++ {
		InvNTTLazy(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsiInv[x], r.NttNInv[x], r.Modulus[x], r.MredParams[x])
	}
}

// butterfly computes X, Y = U + V*Psi, U - V*Psi mod Q.
func butterfly(U, V, Psi, twoQ, fourQ, Q, Qinv uint64) (uint64, uint64) {
	if U >= fourQ {
		U -= fourQ
	}
	V = MRedConstant(V, Psi, Q, Qinv)
	return U + V, U + twoQ - V
}

// NTT computes the NTT on the input coefficients using the input parameters.
func NTT(coeffsIn, coeffsOut []uint64, N int, nttPsi []uint64, Q, mredParams uint64, bredParams []uint64) {

	NTTLazy(coeffsIn, coeffsOut, N, nttPsi, Q, mredParams, bredParams)
	// Finish with an exact reduction
	for i := 0; i < N; i = i + 8 {

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
func NTTLazy(coeffsIn, coeffsOut []uint64, N int, nttPsi []uint64, Q, QInv uint64, bredParams []uint64) {
	var j1, j2, t int
	var F, V uint64

	fourQ := 4 * Q
	twoQ := 2 * Q

	// Copy the result of the first round of butterflies on p2 with approximate reduction
	t = N >> 1
	F = nttPsi[1]

	for j := 0; j <= t-1; j = j + 8 {

		xin := (*[8]uint64)(unsafe.Pointer(&coeffsIn[j]))
		yin := (*[8]uint64)(unsafe.Pointer(&coeffsIn[j+t]))

		xout := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j]))
		yout := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j+t]))

		V = MRedConstant(yin[0], F, Q, QInv)
		xout[0], yout[0] = xin[0]+V, xin[0]+twoQ-V

		V = MRedConstant(yin[1], F, Q, QInv)
		xout[1], yout[1] = xin[1]+V, xin[1]+twoQ-V

		V = MRedConstant(yin[2], F, Q, QInv)
		xout[2], yout[2] = xin[2]+V, xin[2]+twoQ-V

		V = MRedConstant(yin[3], F, Q, QInv)
		xout[3], yout[3] = xin[3]+V, xin[3]+twoQ-V

		V = MRedConstant(yin[4], F, Q, QInv)
		xout[4], yout[4] = xin[4]+V, xin[4]+twoQ-V

		V = MRedConstant(yin[5], F, Q, QInv)
		xout[5], yout[5] = xin[5]+V, xin[5]+twoQ-V

		V = MRedConstant(yin[6], F, Q, QInv)
		xout[6], yout[6] = xin[6]+V, xin[6]+twoQ-V

		V = MRedConstant(yin[7], F, Q, QInv)
		xout[7], yout[7] = xin[7]+V, xin[7]+twoQ-V
	}

	// Continue the rest of the second to the n-1 butterflies on p2 with approximate reduction
	var reduce bool

	for m := 2; m < N; m <<= 1 {

		reduce = (bits.Len64(uint64(m))&1 == 1)

		t >>= 1

		if t >= 8 {

			for i := 0; i < m; i++ {

				j1 = (i * t) << 1

				j2 = j1 + t - 1

				F = nttPsi[m+i]

				if reduce {

					for j := j1; j <= j2; j = j + 8 {

						x := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j]))
						y := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j+t]))

						x[0], y[0] = butterfly(x[0], y[0], F, twoQ, fourQ, Q, QInv)
						x[1], y[1] = butterfly(x[1], y[1], F, twoQ, fourQ, Q, QInv)
						x[2], y[2] = butterfly(x[2], y[2], F, twoQ, fourQ, Q, QInv)
						x[3], y[3] = butterfly(x[3], y[3], F, twoQ, fourQ, Q, QInv)
						x[4], y[4] = butterfly(x[4], y[4], F, twoQ, fourQ, Q, QInv)
						x[5], y[5] = butterfly(x[5], y[5], F, twoQ, fourQ, Q, QInv)
						x[6], y[6] = butterfly(x[6], y[6], F, twoQ, fourQ, Q, QInv)
						x[7], y[7] = butterfly(x[7], y[7], F, twoQ, fourQ, Q, QInv)
					}

				} else {

					for j := j1; j <= j2; j = j + 8 {

						x := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j]))
						y := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j+t]))

						V = MRedConstant(y[0], F, Q, QInv)
						x[0], y[0] = x[0]+V, x[0]+twoQ-V

						V = MRedConstant(y[1], F, Q, QInv)
						x[1], y[1] = x[1]+V, x[1]+twoQ-V

						V = MRedConstant(y[2], F, Q, QInv)
						x[2], y[2] = x[2]+V, x[2]+twoQ-V

						V = MRedConstant(y[3], F, Q, QInv)
						x[3], y[3] = x[3]+V, x[3]+twoQ-V

						V = MRedConstant(y[4], F, Q, QInv)
						x[4], y[4] = x[4]+V, x[4]+twoQ-V

						V = MRedConstant(y[5], F, Q, QInv)
						x[5], y[5] = x[5]+V, x[5]+twoQ-V

						V = MRedConstant(y[6], F, Q, QInv)
						x[6], y[6] = x[6]+V, x[6]+twoQ-V

						V = MRedConstant(y[7], F, Q, QInv)
						x[7], y[7] = x[7]+V, x[7]+twoQ-V
					}
				}
			}

		} else if t == 4 {

			if reduce {

				for i := 0; i < m; i = i + 2 {

					j1 = (i * t) << 1

					psi := (*[2]uint64)(unsafe.Pointer(&nttPsi[m+i]))
					x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

					x[0], x[4] = butterfly(x[0], x[4], psi[0], twoQ, fourQ, Q, QInv)
					x[1], x[5] = butterfly(x[1], x[5], psi[0], twoQ, fourQ, Q, QInv)
					x[2], x[6] = butterfly(x[2], x[6], psi[0], twoQ, fourQ, Q, QInv)
					x[3], x[7] = butterfly(x[3], x[7], psi[0], twoQ, fourQ, Q, QInv)
					x[8], x[12] = butterfly(x[8], x[12], psi[1], twoQ, fourQ, Q, QInv)
					x[9], x[13] = butterfly(x[9], x[13], psi[1], twoQ, fourQ, Q, QInv)
					x[10], x[14] = butterfly(x[10], x[14], psi[1], twoQ, fourQ, Q, QInv)
					x[11], x[15] = butterfly(x[11], x[15], psi[1], twoQ, fourQ, Q, QInv)

				}
			} else {

				for i := 0; i < m; i = i + 2 {

					j1 = (i * t) << 1

					psi := (*[2]uint64)(unsafe.Pointer(&nttPsi[m+i]))
					x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

					V = MRedConstant(x[4], psi[0], Q, QInv)
					x[0], x[4] = x[0]+V, x[0]+twoQ-V

					V = MRedConstant(x[5], psi[0], Q, QInv)
					x[1], x[5] = x[1]+V, x[1]+twoQ-V

					V = MRedConstant(x[6], psi[0], Q, QInv)
					x[2], x[6] = x[2]+V, x[2]+twoQ-V

					V = MRedConstant(x[7], psi[0], Q, QInv)
					x[3], x[7] = x[3]+V, x[3]+twoQ-V

					V = MRedConstant(x[12], psi[1], Q, QInv)
					x[8], x[12] = x[8]+V, x[8]+twoQ-V

					V = MRedConstant(x[13], psi[1], Q, QInv)
					x[9], x[13] = x[9]+V, x[9]+twoQ-V

					V = MRedConstant(x[14], psi[1], Q, QInv)
					x[10], x[14] = x[10]+V, x[10]+twoQ-V

					V = MRedConstant(x[15], psi[1], Q, QInv)
					x[11], x[15] = x[11]+V, x[11]+twoQ-V

				}

			}

		} else if t == 2 {

			if reduce {

				for i := 0; i < m; i = i + 4 {

					j1 = (i * t) << 1

					psi := (*[4]uint64)(unsafe.Pointer(&nttPsi[m+i]))
					x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

					x[0], x[2] = butterfly(x[0], x[2], psi[0], twoQ, fourQ, Q, QInv)
					x[1], x[3] = butterfly(x[1], x[3], psi[0], twoQ, fourQ, Q, QInv)
					x[4], x[6] = butterfly(x[4], x[6], psi[1], twoQ, fourQ, Q, QInv)
					x[5], x[7] = butterfly(x[5], x[7], psi[1], twoQ, fourQ, Q, QInv)
					x[8], x[10] = butterfly(x[8], x[10], psi[2], twoQ, fourQ, Q, QInv)
					x[9], x[11] = butterfly(x[9], x[11], psi[2], twoQ, fourQ, Q, QInv)
					x[12], x[14] = butterfly(x[12], x[14], psi[3], twoQ, fourQ, Q, QInv)
					x[13], x[15] = butterfly(x[13], x[15], psi[3], twoQ, fourQ, Q, QInv)
				}
			} else {

				for i := 0; i < m; i = i + 4 {

					j1 = (i * t) << 1

					psi := (*[4]uint64)(unsafe.Pointer(&nttPsi[m+i]))
					x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

					V = MRedConstant(x[2], psi[0], Q, QInv)
					x[0], x[2] = x[0]+V, x[0]+twoQ-V

					V = MRedConstant(x[3], psi[0], Q, QInv)
					x[1], x[3] = x[1]+V, x[1]+twoQ-V

					V = MRedConstant(x[6], psi[1], Q, QInv)
					x[4], x[6] = x[4]+V, x[4]+twoQ-V

					V = MRedConstant(x[7], psi[1], Q, QInv)
					x[5], x[7] = x[5]+V, x[5]+twoQ-V

					V = MRedConstant(x[10], psi[2], Q, QInv)
					x[8], x[10] = x[8]+V, x[8]+twoQ-V

					V = MRedConstant(x[11], psi[2], Q, QInv)
					x[9], x[11] = x[9]+V, x[9]+twoQ-V

					V = MRedConstant(x[14], psi[3], Q, QInv)
					x[12], x[14] = x[12]+V, x[12]+twoQ-V

					V = MRedConstant(x[15], psi[3], Q, QInv)
					x[13], x[15] = x[13]+V, x[13]+twoQ-V
				}
			}

		} else {

			for i := 0; i < m; i = i + 8 {

				psi := (*[8]uint64)(unsafe.Pointer(&nttPsi[m+i]))
				x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[2*i]))

				x[0], x[1] = butterfly(x[0], x[1], psi[0], twoQ, fourQ, Q, QInv)
				x[2], x[3] = butterfly(x[2], x[3], psi[1], twoQ, fourQ, Q, QInv)
				x[4], x[5] = butterfly(x[4], x[5], psi[2], twoQ, fourQ, Q, QInv)
				x[6], x[7] = butterfly(x[6], x[7], psi[3], twoQ, fourQ, Q, QInv)
				x[8], x[9] = butterfly(x[8], x[9], psi[4], twoQ, fourQ, Q, QInv)
				x[10], x[11] = butterfly(x[10], x[11], psi[5], twoQ, fourQ, Q, QInv)
				x[12], x[13] = butterfly(x[12], x[13], psi[6], twoQ, fourQ, Q, QInv)
				x[14], x[15] = butterfly(x[14], x[15], psi[7], twoQ, fourQ, Q, QInv)
			}

			/*
				for i := uint64(0); i < m; i = i + 8 {

					psi := (*[8]uint64)(unsafe.Pointer(&nttPsi[m+i]))
					x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[2*i]))

					V = MRedConstant(x[1], psi[0], Q, QInv)
					x[0], x[1] = x[0]+V, x[0]+twoQ-V

					V = MRedConstant(x[3], psi[1], Q, QInv)
					x[2], x[3] = x[2]+V, x[2]+twoQ-V

					V = MRedConstant(x[5], psi[2], Q, QInv)
					x[4], x[5] = x[4]+V, x[4]+twoQ-V

					V = MRedConstant(x[7], psi[3], Q, QInv)
					x[6], x[7] = x[6]+V, x[6]+twoQ-V

					V = MRedConstant(x[9], psi[4], Q, QInv)
					x[8], x[9] = x[8]+V, x[8]+twoQ-V

					V = MRedConstant(x[11], psi[5], Q, QInv)
					x[10], x[11] = x[10]+V, x[10]+twoQ-V

					V = MRedConstant(x[13], psi[6], Q, QInv)
					x[12], x[13] = x[12]+V, x[12]+twoQ-V

					V = MRedConstant(x[15], psi[7], Q, QInv)
					x[14], x[15] = x[14]+V, x[14]+twoQ-V
				}
			*/
		}
	}
}

// invbutterfly computes X, Y = U + V, (U - V) * Psi mod Q.
func invbutterfly(U, V, Psi, twoQ, fourQ, Q, Qinv uint64) (X, Y uint64) {
	X = U + V
	if X >= twoQ {
		X -= twoQ
	}
	Y = MRedConstant(U+fourQ-V, Psi, Q, Qinv) // At the moment it is not possible to use MRedConstant if Q > 61 bits
	return
}

// InvNTT computes the InvNTT transformation on the input coefficients using the input parameters.
func InvNTT(coeffsIn, coeffsOut []uint64, N int, nttPsiInv []uint64, nttNInv, Q, QInv uint64) {

	var j1, j2, h, t int
	var F uint64

	// Copy the result of the first round of butterflies on p2 with approximate reduction
	t = 1
	h = N >> 1
	twoQ := Q << 1
	fourQ := Q << 2

	for i := 0; i < h; i = i + 8 {

		psi := (*[8]uint64)(unsafe.Pointer(&nttPsiInv[h+i]))
		xin := (*[16]uint64)(unsafe.Pointer(&coeffsIn[2*i]))
		xout := (*[16]uint64)(unsafe.Pointer(&coeffsOut[2*i]))

		xout[0], xout[1] = invbutterfly(xin[0], xin[1], psi[0], twoQ, fourQ, Q, QInv)
		xout[2], xout[3] = invbutterfly(xin[2], xin[3], psi[1], twoQ, fourQ, Q, QInv)
		xout[4], xout[5] = invbutterfly(xin[4], xin[5], psi[2], twoQ, fourQ, Q, QInv)
		xout[6], xout[7] = invbutterfly(xin[6], xin[7], psi[3], twoQ, fourQ, Q, QInv)
		xout[8], xout[9] = invbutterfly(xin[8], xin[9], psi[4], twoQ, fourQ, Q, QInv)
		xout[10], xout[11] = invbutterfly(xin[10], xin[11], psi[5], twoQ, fourQ, Q, QInv)
		xout[12], xout[13] = invbutterfly(xin[12], xin[13], psi[6], twoQ, fourQ, Q, QInv)
		xout[14], xout[15] = invbutterfly(xin[14], xin[15], psi[7], twoQ, fourQ, Q, QInv)
	}

	// Continue the rest of the second to the n-1 butterflies on p2 with approximate reduction
	t <<= 1
	for m := N >> 1; m > 1; m >>= 1 {

		j1 = 0
		h = m >> 1

		if t >= 8 {

			for i := 0; i < h; i++ {

				j2 = j1 + t - 1

				F = nttPsiInv[h+i]

				for j := j1; j <= j2; j = j + 8 {

					x := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j]))
					y := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j+t]))

					x[0], y[0] = invbutterfly(x[0], y[0], F, twoQ, fourQ, Q, QInv)
					x[1], y[1] = invbutterfly(x[1], y[1], F, twoQ, fourQ, Q, QInv)
					x[2], y[2] = invbutterfly(x[2], y[2], F, twoQ, fourQ, Q, QInv)
					x[3], y[3] = invbutterfly(x[3], y[3], F, twoQ, fourQ, Q, QInv)
					x[4], y[4] = invbutterfly(x[4], y[4], F, twoQ, fourQ, Q, QInv)
					x[5], y[5] = invbutterfly(x[5], y[5], F, twoQ, fourQ, Q, QInv)
					x[6], y[6] = invbutterfly(x[6], y[6], F, twoQ, fourQ, Q, QInv)
					x[7], y[7] = invbutterfly(x[7], y[7], F, twoQ, fourQ, Q, QInv)
				}

				j1 = j1 + (t << 1)
			}

		} else if t == 4 {

			for i := 0; i < h; i = i + 2 {

				psi := (*[2]uint64)(unsafe.Pointer(&nttPsiInv[h+i]))
				x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

				x[0], x[4] = invbutterfly(x[0], x[4], psi[0], twoQ, fourQ, Q, QInv)
				x[1], x[5] = invbutterfly(x[1], x[5], psi[0], twoQ, fourQ, Q, QInv)
				x[2], x[6] = invbutterfly(x[2], x[6], psi[0], twoQ, fourQ, Q, QInv)
				x[3], x[7] = invbutterfly(x[3], x[7], psi[0], twoQ, fourQ, Q, QInv)
				x[8], x[12] = invbutterfly(x[8], x[12], psi[1], twoQ, fourQ, Q, QInv)
				x[9], x[13] = invbutterfly(x[9], x[13], psi[1], twoQ, fourQ, Q, QInv)
				x[10], x[14] = invbutterfly(x[10], x[14], psi[1], twoQ, fourQ, Q, QInv)
				x[11], x[15] = invbutterfly(x[11], x[15], psi[1], twoQ, fourQ, Q, QInv)

				j1 = j1 + (t << 2)
			}

		} else {

			for i := 0; i < h; i = i + 4 {

				psi := (*[4]uint64)(unsafe.Pointer(&nttPsiInv[h+i]))
				x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

				x[0], x[2] = invbutterfly(x[0], x[2], psi[0], twoQ, fourQ, Q, QInv)
				x[1], x[3] = invbutterfly(x[1], x[3], psi[0], twoQ, fourQ, Q, QInv)
				x[4], x[6] = invbutterfly(x[4], x[6], psi[1], twoQ, fourQ, Q, QInv)
				x[5], x[7] = invbutterfly(x[5], x[7], psi[1], twoQ, fourQ, Q, QInv)
				x[8], x[10] = invbutterfly(x[8], x[10], psi[2], twoQ, fourQ, Q, QInv)
				x[9], x[11] = invbutterfly(x[9], x[11], psi[2], twoQ, fourQ, Q, QInv)
				x[12], x[14] = invbutterfly(x[12], x[14], psi[3], twoQ, fourQ, Q, QInv)
				x[13], x[15] = invbutterfly(x[13], x[15], psi[3], twoQ, fourQ, Q, QInv)

				j1 = j1 + (t << 3)
			}
		}

		t <<= 1
	}

	// Finish with an exact reduction
	for i := 0; i < N; i = i + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&coeffsOut[i]))

		x[0] = MRed(x[0], nttNInv, Q, QInv)
		x[1] = MRed(x[1], nttNInv, Q, QInv)
		x[2] = MRed(x[2], nttNInv, Q, QInv)
		x[3] = MRed(x[3], nttNInv, Q, QInv)
		x[4] = MRed(x[4], nttNInv, Q, QInv)
		x[5] = MRed(x[5], nttNInv, Q, QInv)
		x[6] = MRed(x[6], nttNInv, Q, QInv)
		x[7] = MRed(x[7], nttNInv, Q, QInv)
	}
}

// InvNTTLazy computes the InvNTT transformation on the input coefficients using the input parameters with output values in the range [0, 2q-1].
func InvNTTLazy(coeffsIn, coeffsOut []uint64, N int, nttPsiInv []uint64, nttNInv, Q, mredParams uint64) {

	var j1, j2, h, t int
	var F uint64

	// Copy the result of the first round of butterflies on p2 with approximate reduction
	t = 1
	h = N >> 1

	twoQ := Q << 1
	fourQ := Q << 2

	for i := 0; i < h; i = i + 8 {

		psi := (*[8]uint64)(unsafe.Pointer(&nttPsiInv[h+i]))
		xin := (*[16]uint64)(unsafe.Pointer(&coeffsIn[2*i]))
		xout := (*[16]uint64)(unsafe.Pointer(&coeffsOut[2*i]))

		xout[0], xout[1] = invbutterfly(xin[0], xin[1], psi[0], twoQ, fourQ, Q, mredParams)
		xout[2], xout[3] = invbutterfly(xin[2], xin[3], psi[1], twoQ, fourQ, Q, mredParams)
		xout[4], xout[5] = invbutterfly(xin[4], xin[5], psi[2], twoQ, fourQ, Q, mredParams)
		xout[6], xout[7] = invbutterfly(xin[6], xin[7], psi[3], twoQ, fourQ, Q, mredParams)
		xout[8], xout[9] = invbutterfly(xin[8], xin[9], psi[4], twoQ, fourQ, Q, mredParams)
		xout[10], xout[11] = invbutterfly(xin[10], xin[11], psi[5], twoQ, fourQ, Q, mredParams)
		xout[12], xout[13] = invbutterfly(xin[12], xin[13], psi[6], twoQ, fourQ, Q, mredParams)
		xout[14], xout[15] = invbutterfly(xin[14], xin[15], psi[7], twoQ, fourQ, Q, mredParams)
	}

	// Continue the rest of the second to the n-1 butterflies on p2 with approximate reduction
	t <<= 1
	for m := N >> 1; m > 1; m >>= 1 {

		j1 = 0
		h = m >> 1

		if t >= 8 {

			for i := 0; i < h; i++ {

				j2 = j1 + t - 1

				F = nttPsiInv[h+i]

				for j := j1; j <= j2; j = j + 8 {

					x := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j]))
					y := (*[8]uint64)(unsafe.Pointer(&coeffsOut[j+t]))

					x[0], y[0] = invbutterfly(x[0], y[0], F, twoQ, fourQ, Q, mredParams)
					x[1], y[1] = invbutterfly(x[1], y[1], F, twoQ, fourQ, Q, mredParams)
					x[2], y[2] = invbutterfly(x[2], y[2], F, twoQ, fourQ, Q, mredParams)
					x[3], y[3] = invbutterfly(x[3], y[3], F, twoQ, fourQ, Q, mredParams)
					x[4], y[4] = invbutterfly(x[4], y[4], F, twoQ, fourQ, Q, mredParams)
					x[5], y[5] = invbutterfly(x[5], y[5], F, twoQ, fourQ, Q, mredParams)
					x[6], y[6] = invbutterfly(x[6], y[6], F, twoQ, fourQ, Q, mredParams)
					x[7], y[7] = invbutterfly(x[7], y[7], F, twoQ, fourQ, Q, mredParams)
				}

				j1 = j1 + (t << 1)
			}

		} else if t == 4 {

			for i := 0; i < h; i = i + 2 {

				psi := (*[2]uint64)(unsafe.Pointer(&nttPsiInv[h+i]))
				x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

				x[0], x[4] = invbutterfly(x[0], x[4], psi[0], twoQ, fourQ, Q, mredParams)
				x[1], x[5] = invbutterfly(x[1], x[5], psi[0], twoQ, fourQ, Q, mredParams)
				x[2], x[6] = invbutterfly(x[2], x[6], psi[0], twoQ, fourQ, Q, mredParams)
				x[3], x[7] = invbutterfly(x[3], x[7], psi[0], twoQ, fourQ, Q, mredParams)
				x[8], x[12] = invbutterfly(x[8], x[12], psi[1], twoQ, fourQ, Q, mredParams)
				x[9], x[13] = invbutterfly(x[9], x[13], psi[1], twoQ, fourQ, Q, mredParams)
				x[10], x[14] = invbutterfly(x[10], x[14], psi[1], twoQ, fourQ, Q, mredParams)
				x[11], x[15] = invbutterfly(x[11], x[15], psi[1], twoQ, fourQ, Q, mredParams)

				j1 = j1 + (t << 2)
			}

		} else {

			for i := 0; i < h; i = i + 4 {

				psi := (*[4]uint64)(unsafe.Pointer(&nttPsiInv[h+i]))
				x := (*[16]uint64)(unsafe.Pointer(&coeffsOut[j1]))

				x[0], x[2] = invbutterfly(x[0], x[2], psi[0], twoQ, fourQ, Q, mredParams)
				x[1], x[3] = invbutterfly(x[1], x[3], psi[0], twoQ, fourQ, Q, mredParams)
				x[4], x[6] = invbutterfly(x[4], x[6], psi[1], twoQ, fourQ, Q, mredParams)
				x[5], x[7] = invbutterfly(x[5], x[7], psi[1], twoQ, fourQ, Q, mredParams)
				x[8], x[10] = invbutterfly(x[8], x[10], psi[2], twoQ, fourQ, Q, mredParams)
				x[9], x[11] = invbutterfly(x[9], x[11], psi[2], twoQ, fourQ, Q, mredParams)
				x[12], x[14] = invbutterfly(x[12], x[14], psi[3], twoQ, fourQ, Q, mredParams)
				x[13], x[15] = invbutterfly(x[13], x[15], psi[3], twoQ, fourQ, Q, mredParams)

				j1 = j1 + (t << 3)
			}
		}

		t <<= 1
	}

	// Finish with an exact reduction
	for i := 0; i < N; i = i + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&coeffsOut[i]))

		x[0] = MRedConstant(x[0], nttNInv, Q, mredParams)
		x[1] = MRedConstant(x[1], nttNInv, Q, mredParams)
		x[2] = MRedConstant(x[2], nttNInv, Q, mredParams)
		x[3] = MRedConstant(x[3], nttNInv, Q, mredParams)
		x[4] = MRedConstant(x[4], nttNInv, Q, mredParams)
		x[5] = MRedConstant(x[5], nttNInv, Q, mredParams)
		x[6] = MRedConstant(x[6], nttNInv, Q, mredParams)
		x[7] = MRedConstant(x[7], nttNInv, Q, mredParams)
	}
}

///////////////////////////////////
/// For benchmark purposes only ///
///////////////////////////////////

// NTTBarrett performs the NTT operation using Barrett reduction.
// For benchmark purposes only.
func (r *Ring) NTTBarrett(p1, p2 *Poly) {
	for x := range r.Modulus {
		NTTBarrett(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsi[x], r.Modulus[x], r.BredParams[x])
	}
}

// InvNTTBarrett performs the inverse NTT operation using Barrett reduction.
// For benchmark purposes only.
func (r *Ring) InvNTTBarrett(p1, p2 *Poly) {
	for x := range r.Modulus {
		InvNTTBarrett(p1.Coeffs[x], p2.Coeffs[x], r.N, r.NttPsiInv[x], r.NttNInv[x], r.Modulus[x], r.BredParams[x])
	}
}

// butterflyBarrett computes X, Y = U + V*Psi, U - V*Psi mod Q using Barrett reduction.
// For benchmark purposes only.
func butterflyBarrett(U, V, Psi, Q uint64, bredParams []uint64) (X, Y uint64) {
	if U > 2*Q {
		U -= 2 * Q
	}
	V = BRedConstant(V, Psi, Q, bredParams)
	X = U + V
	Y = U + 2*Q - V
	return
}

// invbutterflyBarrett computes X, Y = U + V, (U - V) * Psi mod Q using Barrett reduction.
// For benchmark purposes only.
func invbutterflyBarrett(U, V, Psi, Q uint64, bredParams []uint64) (X, Y uint64) {
	X = U + V
	if X > 2*Q {
		X -= 2 * Q
	}
	Y = BRedConstant(U+2*Q-V, Psi, Q, bredParams) // At the moment it is not possible to use MRedConstant if Q > 61 bits
	return
}

// NTTBarrett computes the NTT using Barrett reduction.
// For benchmark purposes only.
func NTTBarrett(coeffsIn, coeffsOut []uint64, N int, nttPsi []uint64, Q uint64, bredParams []uint64) {
	var j1, j2, t int
	var F uint64

	t = N >> 1
	j2 = t - 1
	F = nttPsi[1]
	for j := 0; j <= j2; j++ {
		coeffsOut[j], coeffsOut[j+t] = butterflyBarrett(coeffsIn[j], coeffsIn[j+t], F, Q, bredParams)
	}

	for m := 2; m < N; m <<= 1 {
		t >>= 1
		for i := 0; i < m; i++ {

			j1 = (i * t) << 1

			j2 = j1 + t - 1

			F = nttPsi[m+i]

			for j := j1; j <= j2; j++ {
				coeffsOut[j], coeffsOut[j+t] = butterflyBarrett(coeffsOut[j], coeffsOut[j+t], F, Q, bredParams)
			}
		}
	}

	for i := 0; i < N; i++ {
		coeffsOut[i] = BRedAdd(coeffsOut[i], Q, bredParams)
	}
}

// InvNTTBarrett computes the Inverse NTT using Barrett reduction.
// For benchmark purposes only.
func InvNTTBarrett(coeffsIn, coeffsOut []uint64, N int, nttPsiInv []uint64, nttNInv, Q uint64, bredParams []uint64) {

	var j1, j2, h, t int
	var F uint64

	t = 1
	j1 = 0
	h = N >> 1

	for i := 0; i < h; i++ {

		j2 = j1

		F = nttPsiInv[h+i]

		for j := j1; j <= j2; j++ {
			coeffsOut[j], coeffsOut[j+t] = invbutterflyBarrett(coeffsIn[j], coeffsIn[j+t], F, Q, bredParams)
		}

		j1 = j1 + (t << 1)
	}

	t <<= 1
	for m := N >> 1; m > 1; m >>= 1 {

		j1 = 0
		h = m >> 1

		for i := 0; i < h; i++ {

			j2 = j1 + t - 1

			F = nttPsiInv[h+i]

			for j := j1; j <= j2; j++ {
				coeffsOut[j], coeffsOut[j+t] = invbutterflyBarrett(coeffsOut[j], coeffsOut[j+t], F, Q, bredParams)
			}

			j1 = j1 + (t << 1)
		}

		t <<= 1
	}

	for j := 0; j < N; j++ {
		coeffsOut[j] = BRed(coeffsOut[j], nttNInv, Q, bredParams)
	}
}
