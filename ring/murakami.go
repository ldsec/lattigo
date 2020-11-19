package ring

import (
	"errors"
	"unsafe"
)

// GenMurakamiParams generates the necessary params for the following map applies the following linear map :
// R[X + X^{-1}]/(X^{2N} + 1) <-> R[X]/(X^{N}-1) <-> R[X]/(X^{N} + 1)
// Papers :
// 1) Approximate Homomorphic Encryptionover the Conjugate-invariant Ring (https://eprint.iacr.org/2018/952)
// 2) Number Theoretic Transforms for Secure Signal Processing (https://arxiv.org/abs/1607.05229
func (r *Ring) GenMurakamiParams() error {

	if r.N == 0 || r.Modulus == nil {
		panic("error : invalid r parameters (missing)")
	}

	// Check if each qi is prime and if qi = 1 mod 2n
	for _, qi := range r.Modulus {
		if IsPrime(qi) == false || qi&((r.N<<2)-1) != 1 {
			r.allowsNTT = false
			return errors.New("warning : provided modulus does not allow Murakami (need 4nth root enabled)")
		}
	}

	r.Murakami = make([][]uint64, len(r.Modulus))
	r.MurakamiInv0 = make([][]uint64, len(r.Modulus))
	r.MurakamiInv1 = make([][]uint64, len(r.Modulus))
	r.MurakamiInv2 = make([]uint64, len(r.Modulus))

	for i, qi := range r.Modulus {

		mredParams := r.MredParams[i]
		bredParams := r.BredParams[i]

		r.Murakami[i] = make([]uint64, r.N)
		r.MurakamiInv0[i] = make([]uint64, r.N)
		r.MurakamiInv1[i] = make([]uint64, r.N)

		g := primitiveRoot(qi)

		_4n := uint64(r.N << 2)

		power := (qi - 1) / _4n
		powerInv := (qi - 1) - power

		w := MForm(ModExp(g, power, qi), qi, r.BredParams[i])
		wInv := MForm(ModExp(g, powerInv, qi), qi, r.BredParams[i])

		_2Inv := MForm(ModExp(2, qi-2, qi), qi, r.BredParams[i])

		for j := uint64(1); j < r.N>>1; j++ {

			w3i := modexpMontgomery(w, 3*j, qi, mredParams, bredParams)
			w2n3i := modexpMontgomery(w, 2*r.N-3*j, qi, mredParams, bredParams)

			wInv3i := modexpMontgomery(wInv, 3*j, qi, mredParams, bredParams)
			wInv2n3i := modexpMontgomery(wInv, 2*r.N-3*j, qi, mredParams, bredParams)

			r.Murakami[i][j] = w3i
			r.Murakami[i][r.N-j] = w2n3i
			r.Murakami[i][r.N-j] = MRed(r.Murakami[i][r.N-j], r.NttPsi[i][1], qi, mredParams)

			r.MurakamiInv0[i][j] = wInv2n3i
			r.MurakamiInv0[i][j] = MRed(r.MurakamiInv0[i][j], w3i, qi, mredParams)
			r.MurakamiInv0[i][r.N-j] = wInv3i
			r.MurakamiInv0[i][r.N-j] = MRed(r.MurakamiInv0[i][r.N-j], wInv3i, qi, mredParams)

			r.MurakamiInv1[i][j] = MRed(_2Inv, wInv3i, qi, mredParams)
			r.MurakamiInv1[i][r.N-j] = MRed(_2Inv, r.NttPsi[i][1], qi, mredParams)
			r.MurakamiInv1[i][r.N-j] = MRed(r.MurakamiInv1[i][r.N-j], w3i, qi, mredParams)
		}

		r.MurakamiInv2[i] = r.NttPsi[i][3] + r.NttPsiInv[i][3]
		r.MurakamiInv2[i] = MRed(r.MurakamiInv2[i], _2Inv, qi, mredParams)
		r.MurakamiInv2[i] = MRed(r.MurakamiInv2[i], r.NttPsi[i][1], qi, mredParams)
	}

	return nil

}

// MapXX2NToXNAndMurakami applies the followingliear map :
// R[X + X^{-1}]/(X^{2N} + 1) -> R[X]/(X^{N}-1) -> R[X]/(X^{N} + 1)
// Papers :
// 1) Approximate Homomorphic Encryptionover the Conjugate-invariant Ring (https://eprint.iacr.org/2018/952)
// 2) Number Theoretic Transforms for Secure Signal Processing (https://arxiv.org/abs/1607.05229)
func (r *Ring) MapXX2NToXNAndMurakami(level uint64, p1 *Poly) {
	// w is a 4nth root of unity
	//a[i], a[n-i] = (a[i]*w[i] + a[n-i]*w[i-n])*w[2i], (a[n-i]*w[n-i] + a[i]*w[m-i])*w[2*(n-i)]
	for i, qi := range r.Modulus[:level+1] {
		MapXX2NToXNAndMurakami(p1.Coeffs[i], r.Murakami[i], r.N, r.NttPsiInv[i][1], r.NttPsi[i][3]+r.NttPsi[i][2], qi, r.MredParams[i])
	}
}

// MapXNToXX2NAndMurakami applies the following linaer map :
// R[X]/(X^{N} + 1) -> R[X]/(X^{N}-1) -> R[X + X^{-1}]/(X^{2N} + 1)
// Papers :
// 1) Approximate Homomorphic Encryptionover the Conjugate-invariant Ring (https://eprint.iacr.org/2018/952)
// 2) Number Theoretic Transforms for Secure Signal Processing (https://arxiv.org/abs/1607.05229)
func (r *Ring) MapXNToXX2NAndMurakami(level uint64, p1 *Poly) {
	// w is a 4nth root of unity
	//a[i], a[n-i] = (a[i]*w[m-i]*w[m-2*i] + a[n-i]*w[i]*w[-2*(n-i)]), (a[n-i]*w[i-n]*w[-2*(n-i)] + a[i]*w[n-i]*w[-2i])
	for i, qi := range r.Modulus[:level+1] {
		MapXNToXX2NAndMurakami(p1.Coeffs[i], r.MurakamiInv0[i], r.MurakamiInv1[i], r.N, r.MurakamiInv2[i], qi, r.MredParams[i])
	}
}

// MapXX2NToXNAndMurakami applies the followingliear map :
// R[X + X^{-1}]/(X^{2N} + 1) -> R[X]/(X^{N}-1) -> R[X]/(X^{N} + 1)
// Papers :
// 1) Approximate Homomorphic Encryptionover the Conjugate-invariant Ring (https://eprint.iacr.org/2018/952)
// 2) Number Theoretic Transforms for Secure Signal Processing (https://arxiv.org/abs/1607.05229)
func MapXX2NToXNAndMurakami(p1tmp, m []uint64, N, psiInv, psi2And3, qi, qInv uint64) {

	for j := uint64(1); j < (N>>1)-7; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p1tmp[N-j-7]))
		m0 := (*[8]uint64)(unsafe.Pointer(&m[j]))
		m1 := (*[8]uint64)(unsafe.Pointer(&m[N-j-7]))

		x[0], y[7] = x[0]+MRed(y[7], psiInv, qi, qInv), y[7]+MRed(x[0], psiInv, qi, qInv)
		x[1], y[6] = x[1]+MRed(y[6], psiInv, qi, qInv), y[6]+MRed(x[1], psiInv, qi, qInv)
		x[2], y[5] = x[2]+MRed(y[5], psiInv, qi, qInv), y[5]+MRed(x[2], psiInv, qi, qInv)
		x[3], y[4] = x[3]+MRed(y[4], psiInv, qi, qInv), y[4]+MRed(x[3], psiInv, qi, qInv)
		x[4], y[3] = x[4]+MRed(y[3], psiInv, qi, qInv), y[3]+MRed(x[4], psiInv, qi, qInv)
		x[5], y[2] = x[5]+MRed(y[2], psiInv, qi, qInv), y[2]+MRed(x[5], psiInv, qi, qInv)
		x[6], y[1] = x[6]+MRed(y[1], psiInv, qi, qInv), y[1]+MRed(x[6], psiInv, qi, qInv)
		x[7], y[0] = x[7]+MRed(y[0], psiInv, qi, qInv), y[0]+MRed(x[7], psiInv, qi, qInv)

		x[0], y[7] = MRed(x[0], m0[0], qi, qInv), MRed(y[7], m1[7], qi, qInv)
		x[1], y[6] = MRed(x[1], m0[1], qi, qInv), MRed(y[6], m1[6], qi, qInv)
		x[2], y[5] = MRed(x[2], m0[2], qi, qInv), MRed(y[5], m1[5], qi, qInv)
		x[3], y[4] = MRed(x[3], m0[3], qi, qInv), MRed(y[4], m1[4], qi, qInv)
		x[4], y[3] = MRed(x[4], m0[4], qi, qInv), MRed(y[3], m1[3], qi, qInv)
		x[5], y[2] = MRed(x[5], m0[5], qi, qInv), MRed(y[2], m1[2], qi, qInv)
		x[6], y[1] = MRed(x[6], m0[6], qi, qInv), MRed(y[1], m1[1], qi, qInv)
		x[7], y[0] = MRed(x[7], m0[7], qi, qInv), MRed(y[0], m1[0], qi, qInv)
	}

	j := (N >> 1) - 7

	x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
	y := (*[8]uint64)(unsafe.Pointer(&p1tmp[N-j-6]))
	m0 := (*[8]uint64)(unsafe.Pointer(&m[j]))
	m1 := (*[8]uint64)(unsafe.Pointer(&m[N-j-6]))

	x[0], y[6] = x[0]+MRed(y[6], psiInv, qi, qInv), y[6]+MRed(x[0], psiInv, qi, qInv)
	x[1], y[5] = x[1]+MRed(y[5], psiInv, qi, qInv), y[5]+MRed(x[1], psiInv, qi, qInv)
	x[2], y[4] = x[2]+MRed(y[4], psiInv, qi, qInv), y[4]+MRed(x[2], psiInv, qi, qInv)
	x[3], y[3] = x[3]+MRed(y[3], psiInv, qi, qInv), y[3]+MRed(x[3], psiInv, qi, qInv)
	x[4], y[2] = x[4]+MRed(y[2], psiInv, qi, qInv), y[2]+MRed(x[4], psiInv, qi, qInv)
	x[5], y[1] = x[5]+MRed(y[1], psiInv, qi, qInv), y[1]+MRed(x[5], psiInv, qi, qInv)
	x[6], y[0] = x[6]+MRed(y[0], psiInv, qi, qInv), y[0]+MRed(x[6], psiInv, qi, qInv)

	x[0], y[6] = MRed(x[0], m0[0], qi, qInv), MRed(y[6], m1[6], qi, qInv)
	x[1], y[5] = MRed(x[1], m0[1], qi, qInv), MRed(y[5], m1[5], qi, qInv)
	x[2], y[4] = MRed(x[2], m0[2], qi, qInv), MRed(y[4], m1[4], qi, qInv)
	x[3], y[3] = MRed(x[3], m0[3], qi, qInv), MRed(y[3], m1[3], qi, qInv)
	x[4], y[2] = MRed(x[4], m0[4], qi, qInv), MRed(y[2], m1[2], qi, qInv)
	x[5], y[1] = MRed(x[5], m0[5], qi, qInv), MRed(y[1], m1[1], qi, qInv)
	x[6], y[0] = MRed(x[6], m0[6], qi, qInv), MRed(y[0], m1[0], qi, qInv)

	p1tmp[N>>1] = MRed(p1tmp[N>>1], psi2And3, qi, qInv)

}

// MapXNToXX2NAndMurakami applies the following linaer map :
// R[X]/(X^{N} + 1) -> R[X]/(X^{N}-1) -> R[X + X^{-1}]/(X^{2N} + 1)
// Papers :
// 1) Approximate Homomorphic Encryptionover the Conjugate-invariant Ring (https://eprint.iacr.org/2018/952)
// 2) Number Theoretic Transforms for Secure Signal Processing (https://arxiv.org/abs/1607.05229)
func MapXNToXX2NAndMurakami(p1tmp, m0, m1 []uint64, N, m2, qi, qInv uint64) {

	for j := uint64(1); j < (N>>1)-7; j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
		y := (*[8]uint64)(unsafe.Pointer(&p1tmp[N-j-7]))
		m00 := (*[8]uint64)(unsafe.Pointer(&m0[j]))
		m01 := (*[8]uint64)(unsafe.Pointer(&m0[N-j-7]))
		m10 := (*[8]uint64)(unsafe.Pointer(&m1[j]))
		m11 := (*[8]uint64)(unsafe.Pointer(&m1[N-j-7]))

		x[0], y[7] = x[0]+MRed(y[7], m00[0], qi, qInv), y[7]+MRed(x[0], m01[7], qi, qInv)
		x[1], y[6] = x[1]+MRed(y[6], m00[1], qi, qInv), y[6]+MRed(x[1], m01[6], qi, qInv)
		x[2], y[5] = x[2]+MRed(y[5], m00[2], qi, qInv), y[5]+MRed(x[2], m01[5], qi, qInv)
		x[3], y[4] = x[3]+MRed(y[4], m00[3], qi, qInv), y[4]+MRed(x[3], m01[4], qi, qInv)
		x[4], y[3] = x[4]+MRed(y[3], m00[4], qi, qInv), y[3]+MRed(x[4], m01[3], qi, qInv)
		x[5], y[2] = x[5]+MRed(y[2], m00[5], qi, qInv), y[2]+MRed(x[5], m01[2], qi, qInv)
		x[6], y[1] = x[6]+MRed(y[1], m00[6], qi, qInv), y[1]+MRed(x[6], m01[1], qi, qInv)
		x[7], y[0] = x[7]+MRed(y[0], m00[7], qi, qInv), y[0]+MRed(x[7], m01[0], qi, qInv)

		x[0], y[7] = MRed(x[0], m10[0], qi, qInv), MRed(y[7], m11[7], qi, qInv)
		x[1], y[6] = MRed(x[1], m10[1], qi, qInv), MRed(y[6], m11[6], qi, qInv)
		x[2], y[5] = MRed(x[2], m10[2], qi, qInv), MRed(y[5], m11[5], qi, qInv)
		x[3], y[4] = MRed(x[3], m10[3], qi, qInv), MRed(y[4], m11[4], qi, qInv)
		x[4], y[3] = MRed(x[4], m10[4], qi, qInv), MRed(y[3], m11[3], qi, qInv)
		x[5], y[2] = MRed(x[5], m10[5], qi, qInv), MRed(y[2], m11[2], qi, qInv)
		x[6], y[1] = MRed(x[6], m10[6], qi, qInv), MRed(y[1], m11[1], qi, qInv)
		x[7], y[0] = MRed(x[7], m10[7], qi, qInv), MRed(y[0], m11[0], qi, qInv)
	}

	j := (N >> 1) - 7
	x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
	y := (*[8]uint64)(unsafe.Pointer(&p1tmp[N-j-6]))
	m00 := (*[8]uint64)(unsafe.Pointer(&m0[j]))
	m01 := (*[8]uint64)(unsafe.Pointer(&m0[N-j-6]))
	m10 := (*[8]uint64)(unsafe.Pointer(&m1[j]))
	m11 := (*[8]uint64)(unsafe.Pointer(&m1[N-j-6]))

	x[0], y[6] = x[0]+MRed(y[6], m00[0], qi, qInv), y[6]+MRed(x[0], m01[6], qi, qInv)
	x[1], y[5] = x[1]+MRed(y[5], m00[1], qi, qInv), y[5]+MRed(x[1], m01[5], qi, qInv)
	x[2], y[4] = x[2]+MRed(y[4], m00[2], qi, qInv), y[4]+MRed(x[2], m01[4], qi, qInv)
	x[3], y[3] = x[3]+MRed(y[3], m00[3], qi, qInv), y[3]+MRed(x[3], m01[3], qi, qInv)
	x[4], y[2] = x[4]+MRed(y[2], m00[4], qi, qInv), y[2]+MRed(x[4], m01[2], qi, qInv)
	x[5], y[1] = x[5]+MRed(y[1], m00[5], qi, qInv), y[1]+MRed(x[5], m01[1], qi, qInv)
	x[6], y[0] = x[6]+MRed(y[0], m00[6], qi, qInv), y[0]+MRed(x[6], m01[0], qi, qInv)

	x[0], y[6] = MRed(x[0], m10[0], qi, qInv), MRed(y[6], m11[6], qi, qInv)
	x[1], y[5] = MRed(x[1], m10[1], qi, qInv), MRed(y[5], m11[5], qi, qInv)
	x[2], y[4] = MRed(x[2], m10[2], qi, qInv), MRed(y[4], m11[4], qi, qInv)
	x[3], y[3] = MRed(x[3], m10[3], qi, qInv), MRed(y[3], m11[3], qi, qInv)
	x[4], y[2] = MRed(x[4], m10[4], qi, qInv), MRed(y[2], m11[2], qi, qInv)
	x[5], y[1] = MRed(x[5], m10[5], qi, qInv), MRed(y[1], m11[1], qi, qInv)
	x[6], y[0] = MRed(x[6], m10[6], qi, qInv), MRed(y[0], m11[0], qi, qInv)

	p1tmp[N>>1] = MRed(p1tmp[N>>1], m2, qi, qInv)

}
