// Package ring implements RNS-accelerated modular arithmetic operations for polynomials, including:
// RNS basis extension; RNS rescaling; number theoretic transform (NTT); uniform, Gaussian and ternary sampling.
package ring

import (
	"bytes"
	"encoding/gob"
	"errors"
	"fmt"
	"math/big"
	"math/bits"

	"github.com/ldsec/lattigo/v2/utils"
)

// Ring is a structure that keeps all the variables required to operate on a polynomial represented in this ring.
type Ring struct {

	// Polynomial nb.Coefficients
	N int

	// Moduli
	Modulus []uint64

	// 2^bit_length(Qi) - 1
	Mask []uint64

	// Indicates whether NTT can be used with the current ring.
	AllowsNTT bool

	// Product of the Moduli
	ModulusBigint *big.Int

	// Fast reduction parameters
	BredParams [][]uint64
	MredParams []uint64

	RescaleParams [][]uint64

	//NTT Parameters
	PsiMont    []uint64 //2N-th primitive root in Montgomery form
	PsiInvMont []uint64 //2N-th inverse primitive root in Montgomery form

	NttPsi    [][]uint64 //powers of the inverse of the 2N-th primitive root in Montgomery form (in bit-reversed order)
	NttPsiInv [][]uint64 //powers of the inverse of the 2N-th primitive root in Montgomery form (in bit-reversed order)
	NttNInv   []uint64   //[N^-1] mod Qi in Montgomery form
}

// NewRing creates a new RNS Ring with degree N and coefficient moduli Moduli. N must be a power of two larger than 8. Moduli should be
// a non-empty []uint64 with distinct prime elements. For the Ring instance to support NTT operation, these elements must also be equal
// to 1 modulo 2*N. Non-nil r and error are returned in the case of non NTT-enabling parameters.
func NewRing(N int, Moduli []uint64) (r *Ring, err error) {
	r = new(Ring)
	err = r.setParameters(N, Moduli)
	if err != nil {
		return nil, err
	}
	return r, r.genNTTParams()
}

// setParameters initializes a *Ring by setting the required precomputed values (except for the NTT-related values, which are set by the
// genNTTParams function).
func (r *Ring) setParameters(N int, Modulus []uint64) error {

	// Checks if N is a power of 2
	if (N < 16) || (N&(N-1)) != 0 && N != 0 {
		return errors.New("invalid ring degree (must be a power of 2 >= 8)")
	}

	if len(Modulus) == 0 {
		return errors.New("invalid modulus (must be a non-empty []uint64)")
	}

	if !utils.AllDistinct(Modulus) {
		return errors.New("invalid modulus (moduli are not distinct)")
	}

	r.AllowsNTT = false

	r.N = N

	r.Modulus = make([]uint64, len(Modulus))
	r.Mask = make([]uint64, len(Modulus))

	for i, qi := range Modulus {
		r.Modulus[i] = qi
		r.Mask[i] = (1 << uint64(bits.Len64(qi))) - 1
	}

	// Compute the bigQ
	r.ModulusBigint = NewInt(1)
	for _, qi := range r.Modulus {
		r.ModulusBigint.Mul(r.ModulusBigint, NewUint(qi))
	}

	// Compute the fast reduction parameters
	r.BredParams = make([][]uint64, len(r.Modulus))
	r.MredParams = make([]uint64, len(r.Modulus))

	for i, qi := range r.Modulus {

		// Compute the fast modular reduction parameters for the Ring
		r.BredParams[i] = BRedParams(qi)

		// If qi is not a power of 2, we can compute the MRedParams (otherwise, it
		// would return an error as there is no valid Montgomery form mod a power of 2)
		if (qi&(qi-1)) != 0 && qi != 0 {
			r.MredParams[i] = MRedParams(qi)
		}
	}

	return nil
}

// genNTTParams checks that N has been correctly initialized, and checks that each modulus is a prime congruent to 1 mod 2N (i.e. NTT-friendly).
// Then, it computes the variables required for the NTT. The purpose of ValidateParameters is to validate that the moduli allow the NTT, and to compute the
// NTT parameters.
func (r *Ring) genNTTParams() error {

	if r.AllowsNTT {
		return nil
	}

	if r.N == 0 || r.Modulus == nil {
		return errors.New("invalid r parameters (missing)")
	}

	// Check if each qi is prime and equal to 1 mod NthRoot
	for i, qi := range r.Modulus {
		if !IsPrime(qi) {
			return fmt.Errorf("invalid modulus (Modulus[%d] is not prime)", i)
		}

		if int(qi)&((r.N<<1)-1) != 1 {
			r.AllowsNTT = false
			return fmt.Errorf("invalid modulus (Modulus[%d] != 1 mod 2N)", i)
		}
	}

	r.RescaleParams = make([][]uint64, len(r.Modulus)-1)

	for j := len(r.Modulus) - 1; j > 0; j-- {

		r.RescaleParams[j-1] = make([]uint64, j)

		for i := 0; i < j; i++ {

			r.RescaleParams[j-1][i] = MForm(r.Modulus[i]-ModExp(r.Modulus[j], int(r.Modulus[i]-2), r.Modulus[i]), r.Modulus[i], r.BredParams[i])
		}
	}

	r.PsiMont = make([]uint64, len(r.Modulus))
	r.PsiInvMont = make([]uint64, len(r.Modulus))
	r.NttPsi = make([][]uint64, len(r.Modulus))
	r.NttPsiInv = make([][]uint64, len(r.Modulus))
	r.NttNInv = make([]uint64, len(r.Modulus))

	bitLenofN := bits.Len64(uint64(r.N)) - 1

	for i, qi := range r.Modulus {

		// 1.1 Compute N^(-1) mod Q in Montgomery form
		r.NttNInv[i] = MForm(ModExp(uint64(r.N), int(qi-2), qi), qi, r.BredParams[i])

		// 1.2 Compute Psi and PsiInv in Montgomery form
		r.NttPsi[i] = make([]uint64, r.N)
		r.NttPsiInv[i] = make([]uint64, r.N)

		// Finds a 2N-th primitive Root
		g := primitiveRoot(qi)

		_2n := r.N << 1

		power := (int(qi) - 1) / _2n
		powerInv := (int(qi) - 1) - power

		// Computes Psi and PsiInv in Montgomery form
		PsiMont := MForm(ModExp(g, power, qi), qi, r.BredParams[i])
		PsiInvMont := MForm(ModExp(g, powerInv, qi), qi, r.BredParams[i])

		r.PsiMont[i] = PsiMont
		r.PsiInvMont[i] = PsiInvMont

		r.NttPsi[i][0] = MForm(1, qi, r.BredParams[i])
		r.NttPsiInv[i][0] = MForm(1, qi, r.BredParams[i])

		// Compute nttPsi[j] = nttPsi[j-1]*Psi and nttPsiInv[j] = nttPsiInv[j-1]*PsiInv
		for j := 1; j < r.N; j++ {

			indexReversePrev := utils.BitReverse64(uint64(j-1), uint64(bitLenofN))
			indexReverseNext := utils.BitReverse64(uint64(j), uint64(bitLenofN))

			r.NttPsi[i][indexReverseNext] = MRed(r.NttPsi[i][indexReversePrev], PsiMont, qi, r.MredParams[i])
			r.NttPsiInv[i][indexReverseNext] = MRed(r.NttPsiInv[i][indexReversePrev], PsiInvMont, qi, r.MredParams[i])
		}
	}

	r.AllowsNTT = true

	return nil
}

// Minimal required information to recover the full ring. Used to import and export the ring.
type ringParams struct {
	N       int
	Modulus []uint64
}

// MarshalBinary encodes the target ring on a slice of bytes.
func (r *Ring) MarshalBinary() ([]byte, error) {

	parameters := ringParams{r.N, r.Modulus}

	var buf bytes.Buffer
	enc := gob.NewEncoder(&buf)
	if err := enc.Encode(parameters); err != nil {
		return nil, err
	}
	return buf.Bytes(), nil
}

// UnmarshalBinary decodes a slice of bytes on the target Ring.
func (r *Ring) UnmarshalBinary(data []byte) error {

	parameters := ringParams{}

	reader := bytes.NewReader(data)
	dec := gob.NewDecoder(reader)
	if err := dec.Decode(&parameters); err != nil {
		return err
	}

	if err := r.setParameters(parameters.N, parameters.Modulus); err != nil {
		return err
	}
	if err := r.genNTTParams(); err != nil {
		return err
	}

	return nil
}

// NewPoly creates a new polynomial with all coefficients set to 0.
func (r *Ring) NewPoly() *Poly {
	p := new(Poly)

	p.Coeffs = make([][]uint64, len(r.Modulus))
	for i := 0; i < len(r.Modulus); i++ {
		p.Coeffs[i] = make([]uint64, r.N)
	}

	return p
}

// NewPolyLvl creates a new polynomial with all coefficients set to 0.
func (r *Ring) NewPolyLvl(level int) *Poly {
	p := new(Poly)

	p.Coeffs = make([][]uint64, level+1)
	for i := 0; i < level+1; i++ {
		p.Coeffs[i] = make([]uint64, r.N)
	}

	return p
}

// SetCoefficientsInt64 sets the coefficients of p1 from an int64 array.
func (r *Ring) SetCoefficientsInt64(coeffs []int64, p1 *Poly) {
	for i, coeff := range coeffs {
		for j, Qi := range r.Modulus {
			p1.Coeffs[j][i] = CRed(uint64((coeff%int64(Qi) + int64(Qi))), Qi)
		}
	}
}

// SetCoefficientsUint64 sets the coefficients of p1 from an uint64 array.
func (r *Ring) SetCoefficientsUint64(coeffs []uint64, p1 *Poly) {
	for i, coeff := range coeffs {
		for j, Qi := range r.Modulus {
			p1.Coeffs[j][i] = coeff % Qi
		}
	}
}

// SetCoefficientsString parses an array of string as Int variables, and sets the
// coefficients of p1 with these Int variables.
func (r *Ring) SetCoefficientsString(coeffs []string, p1 *Poly) {
	QiBigint := new(big.Int)
	coeffTmp := new(big.Int)
	for i, Qi := range r.Modulus {
		QiBigint.SetUint64(Qi)
		for j, coeff := range coeffs {
			p1.Coeffs[i][j] = coeffTmp.Mod(NewIntFromString(coeff), QiBigint).Uint64()
		}
	}
}

// SetCoefficientsBigint sets the coefficients of p1 from an array of Int variables.
func (r *Ring) SetCoefficientsBigint(coeffs []*big.Int, p1 *Poly) {
	QiBigint := new(big.Int)
	coeffTmp := new(big.Int)
	for i, Qi := range r.Modulus {
		QiBigint.SetUint64(Qi)
		for j, coeff := range coeffs {
			p1.Coeffs[i][j] = coeffTmp.Mod(coeff, QiBigint).Uint64()

		}
	}
}

// SetCoefficientsBigintLvl sets the coefficients of p1 from an array of Int variables.
func (r *Ring) SetCoefficientsBigintLvl(level int, coeffs []*big.Int, p1 *Poly) {

	QiBigint := new(big.Int)
	coeffTmp := new(big.Int)
	for i := 0; i < level+1; i++ {
		QiBigint.SetUint64(r.Modulus[i])
		for j, coeff := range coeffs {
			p1.Coeffs[i][j] = coeffTmp.Mod(coeff, QiBigint).Uint64()

		}
	}
}

// PolyToString reconstructs p1 and returns the result in an array of string.
func (r *Ring) PolyToString(p1 *Poly) []string {

	coeffsBigint := make([]*big.Int, r.N)
	r.PolyToBigint(p1, coeffsBigint)
	coeffsString := make([]string, len(coeffsBigint))

	for i := range coeffsBigint {
		coeffsString[i] = coeffsBigint[i].String()
	}

	return coeffsString
}

// PolyToBigint reconstructs p1 and returns the result in an array of Int.
func (r *Ring) PolyToBigint(p1 *Poly, coeffsBigint []*big.Int) {
	r.PolyToBigintLvl(p1.Level(), p1, coeffsBigint)
}

// PolyToBigintLvl reconstructs p1 and returns the result in an array of Int.
func (r *Ring) PolyToBigintLvl(level int, p1 *Poly, coeffsBigint []*big.Int) {
	var qi uint64

	crtReconstruction := make([]*big.Int, level+1)

	QiB := new(big.Int)
	tmp := new(big.Int)
	modulusBigint := NewUint(1)

	for i := 0; i < level+1; i++ {

		qi = r.Modulus[i]
		QiB.SetUint64(qi)

		modulusBigint.Mul(modulusBigint, QiB)

		crtReconstruction[i] = new(big.Int)
		crtReconstruction[i].Quo(r.ModulusBigint, QiB)
		tmp.ModInverse(crtReconstruction[i], QiB)
		tmp.Mod(tmp, QiB)
		crtReconstruction[i].Mul(crtReconstruction[i], tmp)
	}

	for x := 0; x < r.N; x++ {

		tmp.SetUint64(0)
		coeffsBigint[x] = new(big.Int)

		for i := 0; i < level+1; i++ {
			coeffsBigint[x].Add(coeffsBigint[x], tmp.Mul(NewUint(p1.Coeffs[i][x]), crtReconstruction[i]))
		}

		coeffsBigint[x].Mod(coeffsBigint[x], modulusBigint)
	}
}

// PolyToBigintCenteredLvl reconstructs p1 and returns the result in an array of Int.
// Coefficients are centered around Q/2
func (r *Ring) PolyToBigintCenteredLvl(level int, p1 *Poly, coeffsBigint []*big.Int) {
	var qi uint64

	crtReconstruction := make([]*big.Int, level+1)

	QiB := new(big.Int)
	tmp := new(big.Int)
	modulusBigint := NewUint(1)

	for i := 0; i < level+1; i++ {

		qi = r.Modulus[i]
		QiB.SetUint64(qi)

		modulusBigint.Mul(modulusBigint, QiB)

		crtReconstruction[i] = new(big.Int)
		crtReconstruction[i].Quo(r.ModulusBigint, QiB)
		tmp.ModInverse(crtReconstruction[i], QiB)
		tmp.Mod(tmp, QiB)
		crtReconstruction[i].Mul(crtReconstruction[i], tmp)
	}

	modulusBigintHalf := new(big.Int)
	modulusBigintHalf.Rsh(modulusBigint, 1)

	var sign int
	for x := 0; x < r.N; x++ {

		tmp.SetUint64(0)
		coeffsBigint[x].SetUint64(0)

		for i := 0; i < level+1; i++ {
			coeffsBigint[x].Add(coeffsBigint[x], tmp.Mul(NewUint(p1.Coeffs[i][x]), crtReconstruction[i]))
		}

		coeffsBigint[x].Mod(coeffsBigint[x], modulusBigint)

		// Centers the coefficients
		sign = coeffsBigint[x].Cmp(modulusBigintHalf)

		if sign == 1 || sign == 0 {
			coeffsBigint[x].Sub(coeffsBigint[x], modulusBigint)
		}
	}
}

// Equal checks if p1 = p2 in the given Ring.
func (r *Ring) Equal(p1, p2 *Poly) bool {

	for i := 0; i < len(r.Modulus); i++ {
		if len(p1.Coeffs[i]) != len(p2.Coeffs[i]) {
			return false
		}
	}

	r.Reduce(p1, p1)
	r.Reduce(p2, p2)

	for i := 0; i < len(r.Modulus); i++ {
		for j := 0; j < r.N; j++ {
			if p1.Coeffs[i][j] != p2.Coeffs[i][j] {
				return false
			}
		}
	}

	return true
}

// EqualLvl checks if p1 = p2 in the given Ring, up to a given level.
func (r *Ring) EqualLvl(level int, p1, p2 *Poly) bool {

	for i := 0; i < level+1; i++ {
		if len(p1.Coeffs[i]) != len(p2.Coeffs[i]) {
			return false
		}
	}

	r.ReduceLvl(level, p1, p1)
	r.ReduceLvl(level, p2, p2)

	for i := 0; i < level+1; i++ {
		for j := 0; j < r.N; j++ {
			if p1.Coeffs[i][j] != p2.Coeffs[i][j] {
				return false
			}
		}
	}

	return true
}
