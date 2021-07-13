package ring

import (
	"math/big"
	"math/bits"
	"unsafe"

	"github.com/ldsec/lattigo/v2/utils"
)

func (r *Ring) minLevelTernary(p1, p2, p3 *Poly) int {
	return utils.MinInt(utils.MinInt(len(r.Modulus)-1, p1.Level()), utils.MinInt(p2.Level(), p3.Level()))
}

func (r *Ring) minLevelBinary(p1, p2 *Poly) int {
	return utils.MinInt(utils.MinInt(len(r.Modulus)-1, p1.Level()), p2.Level())
}

// Add adds p1 to p2 coefficient-wise and writes the result on p3.
func (r *Ring) Add(p1, p2, p3 *Poly) {
	r.AddLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// AddLvl adds p1 to p2 coefficient-wise for the moduli from
// q_0 up to q_level and writes the result on p3.
func (r *Ring) AddLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))

			z[0] = CRed(x[0]+y[0], qi)
			z[1] = CRed(x[1]+y[1], qi)
			z[2] = CRed(x[2]+y[2], qi)
			z[3] = CRed(x[3]+y[3], qi)
			z[4] = CRed(x[4]+y[4], qi)
			z[5] = CRed(x[5]+y[5], qi)
			z[6] = CRed(x[6]+y[6], qi)
			z[7] = CRed(x[7]+y[7], qi)
		}
	}
}

// AddNoMod adds p1 to p2 coefficient-wise without
// modular reduction and writes the result on p3.
func (r *Ring) AddNoMod(p1, p2, p3 *Poly) {
	r.AddNoModLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// AddNoModLvl adds p1 to p2 coefficient-wise without modular reduction
// for the moduli from q_0 up to q_level and writes the result on p3.
func (r *Ring) AddNoModLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))

			z[0] = x[0] + y[0]
			z[1] = x[1] + y[1]
			z[2] = x[2] + y[2]
			z[3] = x[3] + y[3]
			z[4] = x[4] + y[4]
			z[5] = x[5] + y[5]
			z[6] = x[6] + y[6]
			z[7] = x[7] + y[7]
		}
	}
}

// Sub subtracts p2 to p1 coefficient-wise and writes the result on p3.
func (r *Ring) Sub(p1, p2, p3 *Poly) {
	r.SubLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// SubLvl subtracts p2 to p1 coefficient-wise and writes the result on p3.
func (r *Ring) SubLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))

			z[0] = CRed((x[0]+qi)-y[0], qi)
			z[1] = CRed((x[1]+qi)-y[1], qi)
			z[2] = CRed((x[2]+qi)-y[2], qi)
			z[3] = CRed((x[3]+qi)-y[3], qi)
			z[4] = CRed((x[4]+qi)-y[4], qi)
			z[5] = CRed((x[5]+qi)-y[5], qi)
			z[6] = CRed((x[6]+qi)-y[6], qi)
			z[7] = CRed((x[7]+qi)-y[7], qi)
		}
	}
}

// SubNoMod subtracts p2 to p1 coefficient-wise without
// modular reduction and returns the result on p3.
func (r *Ring) SubNoMod(p1, p2, p3 *Poly) {
	r.SubNoModLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// SubNoModLvl subtracts p2 to p1 coefficient-wise without modular reduction
// for the moduli from q_0 up to q_level and writes the result on p3.
func (r *Ring) SubNoModLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))

			z[0] = x[0] + qi - y[0]
			z[1] = x[1] + qi - y[1]
			z[2] = x[2] + qi - y[2]
			z[3] = x[3] + qi - y[3]
			z[4] = x[4] + qi - y[4]
			z[5] = x[5] + qi - y[5]
			z[6] = x[6] + qi - y[6]
			z[7] = x[7] + qi - y[7]
		}
	}
}

// Neg sets all coefficients of p1 to their additive inverse and writes the result on p2.
func (r *Ring) Neg(p1, p2 *Poly) {
	r.NegLvl(r.minLevelBinary(p1, p2), p1, p2)
}

// NegLvl sets the coefficients of p1 to their additive inverse for
// the moduli from q_0 up to q_level and writes the result on p2.
func (r *Ring) NegLvl(level int, p1, p2 *Poly) {
	for i := 0; i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = qi - x[0]
			z[1] = qi - x[1]
			z[2] = qi - x[2]
			z[3] = qi - x[3]
			z[4] = qi - x[4]
			z[5] = qi - x[5]
			z[6] = qi - x[6]
			z[7] = qi - x[7]
		}
	}
}

// Reduce applies a modular reduction on the coefficients of p1 and writes the result on p2.
func (r *Ring) Reduce(p1, p2 *Poly) {
	r.ReduceLvl(r.minLevelBinary(p1, p2), p1, p2)
}

// ReduceLvl applies a modular reduction on the coefficients of p1
// for the moduli from q_0 up to q_level and writes the result on p2.
func (r *Ring) ReduceLvl(level int, p1, p2 *Poly) {
	for i := 0; i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		bredParams := r.BredParams[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = BRedAdd(x[0], qi, bredParams)
			z[1] = BRedAdd(x[1], qi, bredParams)
			z[2] = BRedAdd(x[2], qi, bredParams)
			z[3] = BRedAdd(x[3], qi, bredParams)
			z[4] = BRedAdd(x[4], qi, bredParams)
			z[5] = BRedAdd(x[5], qi, bredParams)
			z[6] = BRedAdd(x[6], qi, bredParams)
			z[7] = BRedAdd(x[7], qi, bredParams)
		}
	}
}

// ReduceConstant applies a modular reduction on the coefficients of p1 and writes the result on p2.
// Return values in [0, 2q-1]
func (r *Ring) ReduceConstant(p1, p2 *Poly) {
	r.ReduceConstantLvl(r.minLevelBinary(p1, p2), p1, p2)
}

// ReduceConstantLvl applies a modular reduction on the coefficients of p1
// for the moduli from q_0 up to q_level and writes the result on p2.
// Return values in [0, 2q-1]
func (r *Ring) ReduceConstantLvl(level int, p1, p2 *Poly) {
	for i := 0; i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		bredParams := r.BredParams[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = BRedAddConstant(x[0], qi, bredParams)
			z[1] = BRedAddConstant(x[1], qi, bredParams)
			z[2] = BRedAddConstant(x[2], qi, bredParams)
			z[3] = BRedAddConstant(x[3], qi, bredParams)
			z[4] = BRedAddConstant(x[4], qi, bredParams)
			z[5] = BRedAddConstant(x[5], qi, bredParams)
			z[6] = BRedAddConstant(x[6], qi, bredParams)
			z[7] = BRedAddConstant(x[7], qi, bredParams)
		}
	}
}

// Mod applies a modular reduction by m on the coefficients of p1 and writes the result on p2.
func (r *Ring) Mod(p1 *Poly, m uint64, p2 *Poly) {
	bredParams := BRedParams(m)
	for i := range r.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = BRedAdd(x[0], m, bredParams)
			z[1] = BRedAdd(x[1], m, bredParams)
			z[2] = BRedAdd(x[2], m, bredParams)
			z[3] = BRedAdd(x[3], m, bredParams)
			z[4] = BRedAdd(x[4], m, bredParams)
			z[5] = BRedAdd(x[5], m, bredParams)
			z[6] = BRedAdd(x[6], m, bredParams)
			z[7] = BRedAdd(x[7], m, bredParams)
		}
	}
}

// MulCoeffs multiplies p1 by p2 coefficient-wise, performs a
// Barrett modular reduction and writes the result on p3.
func (r *Ring) MulCoeffs(p1, p2, p3 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		bredParams := r.BredParams[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))

			z[0] = BRed(x[0], y[0], qi, bredParams)
			z[1] = BRed(x[1], y[1], qi, bredParams)
			z[2] = BRed(x[2], y[2], qi, bredParams)
			z[3] = BRed(x[3], y[3], qi, bredParams)
			z[4] = BRed(x[4], y[4], qi, bredParams)
			z[5] = BRed(x[5], y[5], qi, bredParams)
			z[6] = BRed(x[6], y[6], qi, bredParams)
			z[7] = BRed(x[7], y[7], qi, bredParams)
		}
	}
}

// MulCoeffsAndAdd multiplies p1 by p2 coefficient-wise with
// a Barret modular reduction and adds the result to p3.
func (r *Ring) MulCoeffsAndAdd(p1, p2, p3 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		bredParams := r.BredParams[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))

			z[0] = CRed(z[0]+BRed(x[0], y[0], qi, bredParams), qi)
			z[1] = CRed(z[1]+BRed(x[1], y[1], qi, bredParams), qi)
			z[2] = CRed(z[2]+BRed(x[2], y[2], qi, bredParams), qi)
			z[3] = CRed(z[3]+BRed(x[3], y[3], qi, bredParams), qi)
			z[4] = CRed(z[4]+BRed(x[4], y[4], qi, bredParams), qi)
			z[5] = CRed(z[5]+BRed(x[5], y[5], qi, bredParams), qi)
			z[6] = CRed(z[6]+BRed(x[6], y[6], qi, bredParams), qi)
			z[7] = CRed(z[7]+BRed(x[7], y[7], qi, bredParams), qi)
		}
	}
}

// MulCoeffsAndAddNoMod multiplies p1 by p2 coefficient-wise with a Barrett
// modular reduction and adds the result to p3 without modular reduction.
func (r *Ring) MulCoeffsAndAddNoMod(p1, p2, p3 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		bredParams := r.BredParams[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))

			z[0] += BRed(x[0], y[0], qi, bredParams)
			z[1] += BRed(x[1], y[1], qi, bredParams)
			z[2] += BRed(x[2], y[2], qi, bredParams)
			z[3] += BRed(x[3], y[3], qi, bredParams)
			z[4] += BRed(x[4], y[4], qi, bredParams)
			z[5] += BRed(x[5], y[5], qi, bredParams)
			z[6] += BRed(x[6], y[6], qi, bredParams)
			z[7] += BRed(x[7], y[7], qi, bredParams)
		}
	}
}

// MulCoeffsMontgomery multiplies p1 by p2 coefficient-wise with a
// Montgomery modular reduction and returns the result on p3.
func (r *Ring) MulCoeffsMontgomery(p1, p2, p3 *Poly) {
	r.MulCoeffsMontgomeryLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// MulCoeffsMontgomeryLvl multiplies p1 by p2 coefficient-wise with a Montgomery
// modular reduction for the moduli from q_0 up to q_level and returns the result on p3.
func (r *Ring) MulCoeffsMontgomeryLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))

			z[0] = MRed(x[0], y[0], qi, mredParams)
			z[1] = MRed(x[1], y[1], qi, mredParams)
			z[2] = MRed(x[2], y[2], qi, mredParams)
			z[3] = MRed(x[3], y[3], qi, mredParams)
			z[4] = MRed(x[4], y[4], qi, mredParams)
			z[5] = MRed(x[5], y[5], qi, mredParams)
			z[6] = MRed(x[6], y[6], qi, mredParams)
			z[7] = MRed(x[7], y[7], qi, mredParams)
		}
	}
}

// MulCoeffsMontgomeryConstant multiplies p1 by p2 coefficient-wise with a
// constant-time Montgomery modular reduction and writes the result on p3.
func (r *Ring) MulCoeffsMontgomeryConstant(p1, p2, p3 *Poly) {
	r.MulCoeffsMontgomeryConstantLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// MulCoeffsMontgomeryConstantLvl multiplies p1 by p2 coefficient-wise with a Montgomery
// modular reduction for the moduli from q_0 up to q_level and returns the result on p3.
func (r *Ring) MulCoeffsMontgomeryConstantLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))

			z[0] = MRedConstant(x[0], y[0], qi, mredParams)
			z[1] = MRedConstant(x[1], y[1], qi, mredParams)
			z[2] = MRedConstant(x[2], y[2], qi, mredParams)
			z[3] = MRedConstant(x[3], y[3], qi, mredParams)
			z[4] = MRedConstant(x[4], y[4], qi, mredParams)
			z[5] = MRedConstant(x[5], y[5], qi, mredParams)
			z[6] = MRedConstant(x[6], y[6], qi, mredParams)
			z[7] = MRedConstant(x[7], y[7], qi, mredParams)
		}
	}
}

// MulCoeffsMontgomeryAndAdd multiplies p1 by p2 coefficient-wise with a
// Montgomery modular reduction and adds the result to p3.
func (r *Ring) MulCoeffsMontgomeryAndAdd(p1, p2, p3 *Poly) {
	r.MulCoeffsMontgomeryAndAddLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// MulCoeffsMontgomeryAndAddLvl multiplies p1 by p2 coefficient-wise with a Montgomery
// modular reduction for the moduli from q_0 up to q_level and adds the result to p3.
func (r *Ring) MulCoeffsMontgomeryAndAddLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))

			z[0] = CRed(z[0]+MRed(x[0], y[0], qi, mredParams), qi)
			z[1] = CRed(z[1]+MRed(x[1], y[1], qi, mredParams), qi)
			z[2] = CRed(z[2]+MRed(x[2], y[2], qi, mredParams), qi)
			z[3] = CRed(z[3]+MRed(x[3], y[3], qi, mredParams), qi)
			z[4] = CRed(z[4]+MRed(x[4], y[4], qi, mredParams), qi)
			z[5] = CRed(z[5]+MRed(x[5], y[5], qi, mredParams), qi)
			z[6] = CRed(z[6]+MRed(x[6], y[6], qi, mredParams), qi)
			z[7] = CRed(z[7]+MRed(x[7], y[7], qi, mredParams), qi)
		}
	}
}

// MulCoeffsMontgomeryAndAddNoMod multiplies p1 by p2 coefficient-wise with a
// Montgomery modular reduction and adds the result to p3 without modular reduction.
func (r *Ring) MulCoeffsMontgomeryAndAddNoMod(p1, p2, p3 *Poly) {
	r.MulCoeffsMontgomeryAndAddNoModLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// MulCoeffsMontgomeryAndAddNoModLvl multiplies p1 by p2 coefficient-wise with a Montgomery modular
// reduction for the moduli from q_0 up to q_level and adds the result to p3 without modular reduction.
func (r *Ring) MulCoeffsMontgomeryAndAddNoModLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))

			z[0] += MRed(x[0], y[0], qi, mredParams)
			z[1] += MRed(x[1], y[1], qi, mredParams)
			z[2] += MRed(x[2], y[2], qi, mredParams)
			z[3] += MRed(x[3], y[3], qi, mredParams)
			z[4] += MRed(x[4], y[4], qi, mredParams)
			z[5] += MRed(x[5], y[5], qi, mredParams)
			z[6] += MRed(x[6], y[6], qi, mredParams)
			z[7] += MRed(x[7], y[7], qi, mredParams)
		}
	}
}

// MulCoeffsMontgomeryConstantAndAddNoMod multiplies p1 by p2 coefficient-wise with a
// Montgomery modular reduction and adds the result to p3 without modular reduction.
// Return values in [0, 3q-1]
func (r *Ring) MulCoeffsMontgomeryConstantAndAddNoMod(p1, p2, p3 *Poly) {
	r.MulCoeffsMontgomeryConstantAndAddNoModLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// MulCoeffsMontgomeryConstantAndAddNoModLvl multiplies p1 by p2 coefficient-wise with a constant-time Montgomery
// modular reduction for the moduli from q_0 up to q_level and adds the result to p3 without modular reduction.
// Return values in [0, 3q-1]
func (r *Ring) MulCoeffsMontgomeryConstantAndAddNoModLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))

			z[0] += MRedConstant(x[0], y[0], qi, mredParams)
			z[1] += MRedConstant(x[1], y[1], qi, mredParams)
			z[2] += MRedConstant(x[2], y[2], qi, mredParams)
			z[3] += MRedConstant(x[3], y[3], qi, mredParams)
			z[4] += MRedConstant(x[4], y[4], qi, mredParams)
			z[5] += MRedConstant(x[5], y[5], qi, mredParams)
			z[6] += MRedConstant(x[6], y[6], qi, mredParams)
			z[7] += MRedConstant(x[7], y[7], qi, mredParams)
		}
	}
}

// MulCoeffsMontgomeryAndSub multiplies p1 by p2 coefficient-wise with
// a Montgomery modular reduction and subtracts the result from p3.
func (r *Ring) MulCoeffsMontgomeryAndSub(p1, p2, p3 *Poly) {
	r.MulCoeffsMontgomeryAndSubLvl(r.minLevelTernary(p1, p2, p3), p1, p2, p3)
}

// MulCoeffsMontgomeryAndSubLvl multiplies p1 by p2 coefficient-wise with
// a Montgomery modular reduction and subtracts the result from p3.
func (r *Ring) MulCoeffsMontgomeryAndSubLvl(level int, p1, p2, p3 *Poly) {
	for i := 0; i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))

			z[0] = CRed(z[0]+(qi-MRed(x[0], y[0], qi, mredParams)), qi)
			z[1] = CRed(z[1]+(qi-MRed(x[1], y[1], qi, mredParams)), qi)
			z[2] = CRed(z[2]+(qi-MRed(x[2], y[2], qi, mredParams)), qi)
			z[3] = CRed(z[3]+(qi-MRed(x[3], y[3], qi, mredParams)), qi)
			z[4] = CRed(z[4]+(qi-MRed(x[4], y[4], qi, mredParams)), qi)
			z[5] = CRed(z[5]+(qi-MRed(x[5], y[5], qi, mredParams)), qi)
			z[6] = CRed(z[6]+(qi-MRed(x[6], y[6], qi, mredParams)), qi)
			z[7] = CRed(z[7]+(qi-MRed(x[7], y[7], qi, mredParams)), qi)
		}
	}
}

// MulCoeffsMontgomeryAndSubNoMod multiplies p1 by p2 coefficient-wise with a Montgomery
// modular reduction and subtracts the result from p3 without modular reduction.
func (r *Ring) MulCoeffsMontgomeryAndSubNoMod(p1, p2, p3 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))

			z[0] += (qi - MRed(x[0], y[0], qi, mredParams))
			z[1] += (qi - MRed(x[1], y[1], qi, mredParams))
			z[2] += (qi - MRed(x[2], y[2], qi, mredParams))
			z[3] += (qi - MRed(x[3], y[3], qi, mredParams))
			z[4] += (qi - MRed(x[4], y[4], qi, mredParams))
			z[5] += (qi - MRed(x[5], y[5], qi, mredParams))
			z[6] += (qi - MRed(x[6], y[6], qi, mredParams))
			z[7] += (qi - MRed(x[7], y[7], qi, mredParams))
		}
	}
}

// MulCoeffsConstant multiplies p1 by p2 coefficient-wise with a constant-time
// Barrett modular reduction and writes the result on p3.
func (r *Ring) MulCoeffsConstant(p1, p2, p3 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp, p3tmp := p1.Coeffs[i], p2.Coeffs[i], p3.Coeffs[i]
		bredParams := r.BredParams[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))

			z[0] = BRedConstant(x[0], y[0], qi, bredParams)
			z[1] = BRedConstant(x[1], y[1], qi, bredParams)
			z[2] = BRedConstant(x[2], y[2], qi, bredParams)
			z[3] = BRedConstant(x[3], y[3], qi, bredParams)
			z[4] = BRedConstant(x[4], y[4], qi, bredParams)
			z[5] = BRedConstant(x[5], y[5], qi, bredParams)
			z[6] = BRedConstant(x[6], y[6], qi, bredParams)
			z[7] = BRedConstant(x[7], y[7], qi, bredParams)
		}
	}
}

// AddScalar adds a scalar to each coefficient of p1 and writes the result on p2.
func (r *Ring) AddScalar(p1 *Poly, scalar uint64, p2 *Poly) {
	for i, Qi := range r.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p1.Coeffs[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = CRed(x[0]+scalar, Qi)
			z[1] = CRed(x[1]+scalar, Qi)
			z[2] = CRed(x[2]+scalar, Qi)
			z[3] = CRed(x[3]+scalar, Qi)
			z[4] = CRed(x[4]+scalar, Qi)
			z[5] = CRed(x[5]+scalar, Qi)
			z[6] = CRed(x[6]+scalar, Qi)
			z[7] = CRed(x[7]+scalar, Qi)
		}
	}
}

// AddScalarBigint adds a big.Int scalar to each coefficient of p1 and writes the result on p2.
func (r *Ring) AddScalarBigint(p1 *Poly, scalar *big.Int, p2 *Poly) {
	tmp := new(big.Int)
	for i, Qi := range r.Modulus {
		scalarQi := tmp.Mod(scalar, NewUint(Qi)).Uint64()
		p1tmp, p2tmp := p1.Coeffs[i], p1.Coeffs[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = CRed(x[0]+scalarQi, Qi)
			z[1] = CRed(x[1]+scalarQi, Qi)
			z[2] = CRed(x[2]+scalarQi, Qi)
			z[3] = CRed(x[3]+scalarQi, Qi)
			z[4] = CRed(x[4]+scalarQi, Qi)
			z[5] = CRed(x[5]+scalarQi, Qi)
			z[6] = CRed(x[6]+scalarQi, Qi)
			z[7] = CRed(x[7]+scalarQi, Qi)
		}
	}
}

// SubScalar subtracts a scalar from each coefficient of p1 and writes the result on p2.
func (r *Ring) SubScalar(p1 *Poly, scalar uint64, p2 *Poly) {
	for i, Qi := range r.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p1.Coeffs[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = CRed(x[0]+Qi-scalar, Qi)
			z[1] = CRed(x[1]+Qi-scalar, Qi)
			z[2] = CRed(x[2]+Qi-scalar, Qi)
			z[3] = CRed(x[3]+Qi-scalar, Qi)
			z[4] = CRed(x[4]+Qi-scalar, Qi)
			z[5] = CRed(x[5]+Qi-scalar, Qi)
			z[6] = CRed(x[6]+Qi-scalar, Qi)
			z[7] = CRed(x[7]+Qi-scalar, Qi)
		}
	}
}

// SubScalarBigint subtracts a big.Int scalar from each coefficient of p1 and writes the result on p2.
func (r *Ring) SubScalarBigint(p1 *Poly, scalar *big.Int, p2 *Poly) {
	tmp := new(big.Int)
	for i, Qi := range r.Modulus {
		scalarQi := tmp.Mod(scalar, NewUint(Qi)).Uint64()
		p1tmp, p2tmp := p1.Coeffs[i], p1.Coeffs[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = CRed(x[0]+Qi-scalarQi, Qi)
			z[1] = CRed(x[1]+Qi-scalarQi, Qi)
			z[2] = CRed(x[2]+Qi-scalarQi, Qi)
			z[3] = CRed(x[3]+Qi-scalarQi, Qi)
			z[4] = CRed(x[4]+Qi-scalarQi, Qi)
			z[5] = CRed(x[5]+Qi-scalarQi, Qi)
			z[6] = CRed(x[6]+Qi-scalarQi, Qi)
			z[7] = CRed(x[7]+Qi-scalarQi, Qi)
		}
	}
}

// MulScalar multiplies each coefficient of p1 by a scalar and writes the result on p2.
func (r *Ring) MulScalar(p1 *Poly, scalar uint64, p2 *Poly) {
	r.MulScalarLvl(r.minLevelBinary(p1, p2), p1, scalar, p2)
}

// MulScalarLvl multiplies each coefficient of p1 by a scalar for the moduli from q_0 up to q_level and writes the result on p2.
func (r *Ring) MulScalarLvl(level int, p1 *Poly, scalar uint64, p2 *Poly) {
	for i := 0; i < level+1; i++ {
		Qi := r.Modulus[i]
		scalarMont := MForm(BRedAdd(scalar, Qi, r.BredParams[i]), Qi, r.BredParams[i])
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = MRed(x[0], scalarMont, Qi, mredParams)
			z[1] = MRed(x[1], scalarMont, Qi, mredParams)
			z[2] = MRed(x[2], scalarMont, Qi, mredParams)
			z[3] = MRed(x[3], scalarMont, Qi, mredParams)
			z[4] = MRed(x[4], scalarMont, Qi, mredParams)
			z[5] = MRed(x[5], scalarMont, Qi, mredParams)
			z[6] = MRed(x[6], scalarMont, Qi, mredParams)
			z[7] = MRed(x[7], scalarMont, Qi, mredParams)
		}
	}
}

// MulScalarBigint multiplies each coefficient of p1 by a big.Int scalar and writes the result on p2.
func (r *Ring) MulScalarBigint(p1 *Poly, scalar *big.Int, p2 *Poly) {
	r.MulScalarBigintLvl(r.minLevelBinary(p1, p2), p1, scalar, p2)
}

// MulScalarBigintLvl multiplies each coefficient of p1 by a big.Int scalar
//for the moduli from q_0 up to q_level and writes the result on p2.
func (r *Ring) MulScalarBigintLvl(level int, p1 *Poly, scalar *big.Int, p2 *Poly) {
	scalarQi := new(big.Int)
	for i := 0; i < level+1; i++ {
		Qi := r.Modulus[i]
		scalarQi.Mod(scalar, NewUint(Qi))
		scalarMont := MForm(BRedAdd(scalarQi.Uint64(), Qi, r.BredParams[i]), Qi, r.BredParams[i])
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = MRed(x[0], scalarMont, Qi, mredParams)
			z[1] = MRed(x[1], scalarMont, Qi, mredParams)
			z[2] = MRed(x[2], scalarMont, Qi, mredParams)
			z[3] = MRed(x[3], scalarMont, Qi, mredParams)
			z[4] = MRed(x[4], scalarMont, Qi, mredParams)
			z[5] = MRed(x[5], scalarMont, Qi, mredParams)
			z[6] = MRed(x[6], scalarMont, Qi, mredParams)
			z[7] = MRed(x[7], scalarMont, Qi, mredParams)
		}
	}
}

// Shift circulary shifts the coefficients of the polynomial p1 by n positions to the left and writes the result on p2.
func (r *Ring) Shift(p1 *Poly, n int, p2 *Poly) {
	mask := (1 << r.N) - 1
	for i := range r.Modulus {
		p2.Coeffs[i] = append(p1.Coeffs[i][(n&mask):], p1.Coeffs[i][:(n&mask)]...)
	}
}

// MForm switches p1 to the Montgomery domain and writes the result on p2.
func (r *Ring) MForm(p1, p2 *Poly) {
	r.MFormLvl(r.minLevelBinary(p1, p2), p1, p2)
}

// MFormLvl switches p1 to the Montgomery domain for the moduli from q_0 up to q_level and writes the result on p2.
func (r *Ring) MFormLvl(level int, p1, p2 *Poly) {
	for i := 0; i < level+1; i++ {
		qi := r.Modulus[i]
		bredParams := r.BredParams[i]
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = MForm(x[0], qi, bredParams)
			z[1] = MForm(x[1], qi, bredParams)
			z[2] = MForm(x[2], qi, bredParams)
			z[3] = MForm(x[3], qi, bredParams)
			z[4] = MForm(x[4], qi, bredParams)
			z[5] = MForm(x[5], qi, bredParams)
			z[6] = MForm(x[6], qi, bredParams)
			z[7] = MForm(x[7], qi, bredParams)
		}
	}
}

// InvMForm switches back p1 from the Montgomery domain to the conventional domain and writes the result on p2.
func (r *Ring) InvMForm(p1, p2 *Poly) {
	r.InvMFormLvl(r.minLevelBinary(p1, p2), p1, p2)
}

// InvMFormLvl switches back p1 from the Montgomery domain to the conventional domain and writes the result on p2.
func (r *Ring) InvMFormLvl(level int, p1, p2 *Poly) {
	for i, qi := range r.Modulus[:level+1] {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = InvMForm(x[0], qi, mredParams)
			z[1] = InvMForm(x[1], qi, mredParams)
			z[2] = InvMForm(x[2], qi, mredParams)
			z[3] = InvMForm(x[3], qi, mredParams)
			z[4] = InvMForm(x[4], qi, mredParams)
			z[5] = InvMForm(x[5], qi, mredParams)
			z[6] = InvMForm(x[6], qi, mredParams)
			z[7] = InvMForm(x[7], qi, mredParams)
		}
	}
}

// MulByPow2New multiplies p1 by 2^pow2 and returns the result in a new polynomial p2.
func (r *Ring) MulByPow2New(p1 *Poly, pow2 int) (p2 *Poly) {
	p2 = r.NewPoly()
	r.MulByPow2(p1, pow2, p2)
	return
}

// MulByPow2 multiplies p1 by 2^pow2 and writes the result on p2.
func (r *Ring) MulByPow2(p1 *Poly, pow2 int, p2 *Poly) {
	r.MulByPow2Lvl(r.minLevelBinary(p1, p2), p1, pow2, p2)
}

// MulByPow2Lvl multiplies p1 by 2^pow2 for the moduli from q_0 up to q_level and writes the result on p2.
func (r *Ring) MulByPow2Lvl(level int, p1 *Poly, pow2 int, p2 *Poly) {
	r.MFormLvl(level, p1, p2)
	for i := 0; i < level+1; i++ {
		qi := r.Modulus[i]
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = PowerOf2(x[0], pow2, qi, mredParams)
			z[1] = PowerOf2(x[1], pow2, qi, mredParams)
			z[2] = PowerOf2(x[2], pow2, qi, mredParams)
			z[3] = PowerOf2(x[3], pow2, qi, mredParams)
			z[4] = PowerOf2(x[4], pow2, qi, mredParams)
			z[5] = PowerOf2(x[5], pow2, qi, mredParams)
			z[6] = PowerOf2(x[6], pow2, qi, mredParams)
			z[7] = PowerOf2(x[7], pow2, qi, mredParams)
		}
	}
}

// MultByMonomialNew multiplies p1 by x^monomialDeg and writes the result on a new polynomial p2.
func (r *Ring) MultByMonomialNew(p1 *Poly, monomialDeg int) (p2 *Poly) {
	p2 = r.NewPoly()
	r.MultByMonomial(p1, monomialDeg, p2)
	return
}

// MultByMonomial multiplies p1 by x^monomialDeg and writes the result on p2.
func (r *Ring) MultByMonomial(p1 *Poly, monomialDeg int, p2 *Poly) {

	shift := monomialDeg % (r.N << 1)

	if shift == 0 {

		for i := range r.Modulus {
			p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
			for j := 0; j < r.N; j++ {
				p2tmp[j] = p1tmp[j]
			}
		}

	} else {

		tmpx := r.NewPoly()

		if shift < r.N {

			for i := range r.Modulus {
				p1tmp, tmpxT := p1.Coeffs[i], tmpx.Coeffs[i]
				for j := 0; j < r.N; j++ {
					tmpxT[j] = p1tmp[j]
				}
			}

		} else {

			for i, qi := range r.Modulus {
				p1tmp, tmpxT := p1.Coeffs[i], tmpx.Coeffs[i]
				for j := 0; j < r.N; j++ {
					tmpxT[j] = qi - p1tmp[j]
				}
			}
		}

		shift %= r.N

		for i, qi := range r.Modulus {
			p2tmp, tmpxT := p2.Coeffs[i], tmpx.Coeffs[i]
			for j := 0; j < shift; j++ {
				p2tmp[j] = qi - tmpxT[r.N-shift+j]
			}
		}

		for i := range r.Modulus {
			p2tmp, tmpxT := p2.Coeffs[i], tmpx.Coeffs[i]
			for j := shift; j < r.N; j++ {
				p2tmp[j] = tmpxT[j-shift]

			}
		}
	}
}

// MulByVectorMontgomery multiplies p1 by a vector of uint64 coefficients and writes the result on p2.
func (r *Ring) MulByVectorMontgomery(p1 *Poly, vector []uint64, p2 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&vector[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = MRed(x[0], y[0], qi, mredParams)
			z[1] = MRed(x[1], y[1], qi, mredParams)
			z[2] = MRed(x[2], y[2], qi, mredParams)
			z[3] = MRed(x[3], y[3], qi, mredParams)
			z[4] = MRed(x[4], y[4], qi, mredParams)
			z[5] = MRed(x[5], y[5], qi, mredParams)
			z[6] = MRed(x[6], y[6], qi, mredParams)
			z[7] = MRed(x[7], y[7], qi, mredParams)
		}
	}
}

// MulByVectorMontgomeryAndAddNoMod multiplies p1 by a vector of uint64 coefficients and adds the result on p2 without modular reduction.
func (r *Ring) MulByVectorMontgomeryAndAddNoMod(p1 *Poly, vector []uint64, p2 *Poly) {
	for i, qi := range r.Modulus {
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		mredParams := r.MredParams[i]
		for j := 0; j < r.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&vector[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] += MRed(x[0], y[0], qi, mredParams)
			z[1] += MRed(x[1], y[1], qi, mredParams)
			z[2] += MRed(x[2], y[2], qi, mredParams)
			z[3] += MRed(x[3], y[3], qi, mredParams)
			z[4] += MRed(x[4], y[4], qi, mredParams)
			z[5] += MRed(x[5], y[5], qi, mredParams)
			z[6] += MRed(x[6], y[6], qi, mredParams)
			z[7] += MRed(x[7], y[7], qi, mredParams)
		}
	}
}

// BitReverse applies a bit reverse permutation on the coefficients of p1 and writes the result on p2.
// In can safely be used for in-place permutation.
func (r *Ring) BitReverse(p1, p2 *Poly) {
	bitLenOfN := uint64(bits.Len64(uint64(r.N)) - 1)

	if p1 != p2 {
		for i := range r.Modulus {
			p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
			for j := 0; j < r.N; j++ {
				p2tmp[utils.BitReverse64(uint64(j), bitLenOfN)] = p1tmp[j]
			}
		}
	} else { // In place in case p1 = p2
		for x := range r.Modulus {
			p2tmp := p2.Coeffs[x]
			for i := 0; i < r.N; i++ {
				j := utils.BitReverse64(uint64(i), bitLenOfN)
				if i < int(j) {
					p2tmp[i], p2tmp[j] = p2tmp[j], p2tmp[i]
				}
			}
		}
	}
}

// Rotate applies a Galois automorphism on p1 in NTT form,
// rotating the coefficients to the right by n positions, and writes the result on p2.
// It requires the data to be permuted in bit-reversal order before applying the NTT.
func (r *Ring) Rotate(p1 *Poly, n int, p2 *Poly) {

	var root, gal uint64

	n &= (1 << r.N) - 1

	for i, qi := range r.Modulus {

		mredParams := r.MredParams[i]

		root = MRed(r.PsiMont[i], r.PsiMont[i], qi, mredParams)

		root = modexpMontgomery(root, n, qi, mredParams, r.BredParams[i])

		gal = MForm(1, qi, r.BredParams[i])

		p1tmp, p2tmp := p1.Coeffs[i], p1.Coeffs[i]

		for j := 1; j < r.N; j++ {

			gal = MRed(gal, root, qi, mredParams)

			p2tmp[j] = MRed(p1tmp[j], gal, qi, mredParams)

		}
	}
}
