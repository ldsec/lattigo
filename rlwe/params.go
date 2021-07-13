package rlwe

import (
	"encoding/json"
	"fmt"
	"math"
	"math/big"
	"math/bits"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// MaxLogN is the log2 of the largest supported polynomial modulus degree.
const MaxLogN = 16

// MinLogN is the log2 of the smallest supported polynomial modulus degree (needed to ensure the NTT correctness).
const MinLogN = 4

// MaxModuliCount is the largest supported number of moduli in the RNS representation.
const MaxModuliCount = 34

// MaxModuliSize is the largest bit-length supported for the moduli in the RNS representation.
const MaxModuliSize = 60

// DefaultSigma is the default error distribution standard deviation
const DefaultSigma = 3.2

// GaloisGen is an integer of order N=2^d modulo M=2N and that spans Z_M with the integer -1.
// The j-th ring automorphism takes the root zeta to zeta^(5j).
const GaloisGen uint64 = 5

// ParametersLiteral is a literal representation of BFV parameters.  It has public
// fields and is used to express unchecked user-defined parameters literally into
// Go programs. The NewParametersFromLiteral function is used to generate the actual
// checked parameters from the literal representation.
type ParametersLiteral struct {
	LogN  int
	Q     []uint64
	P     []uint64
	LogQ  []int `json:",omitempty"`
	LogP  []int `json:",omitempty"`
	Sigma float64
}

// Parameters represents a set of generic RLWE parameters. Its fields are private and
// immutable. See ParametersLiteral for user-specified parameters.
type Parameters struct {
	logN   int
	qi     []uint64
	pi     []uint64
	sigma  float64
	ringQ  *ring.Ring
	ringP  *ring.Ring
	ringQP *ring.Ring
}

var (
	// TestPN12QP109 is a set of default parameters with logN=12 and logQP=109
	TestPN12QP109 = ParametersLiteral{
		LogN:  12,
		Q:     []uint64{0x7ffffec001, 0x40002001}, // 39 + 39 bits
		P:     []uint64{0x8000016001},             // 30 bits
		Sigma: DefaultSigma,
	}
	// TestPN13QP218 is a set of default parameters with logN=13 and logQP=218
	TestPN13QP218 = ParametersLiteral{
		LogN:  13,
		Q:     []uint64{0x3fffffffef8001, 0x4000000011c001, 0x40000000120001}, // 54 + 54 + 54 bits
		P:     []uint64{0x7ffffffffb4001},                                     // 55 bits
		Sigma: DefaultSigma,
	}

	// TestPN14QP438 is a set of default parameters with logN=14 and logQP=438
	TestPN14QP438 = ParametersLiteral{
		LogN: 14,
		Q: []uint64{0x100000000060001, 0x80000000068001, 0x80000000080001,
			0x3fffffffef8001, 0x40000000120001, 0x3fffffffeb8001}, // 56 + 55 + 55 + 54 + 54 + 54 bits
		P:     []uint64{0x80000000130001, 0x7fffffffe90001}, // 55 + 55 bits
		Sigma: DefaultSigma,
	}

	// TestPN15QP880 is a set of default parameters with logN=15 and logQP=880
	TestPN15QP880 = ParametersLiteral{
		LogN: 15,
		Q: []uint64{0x7ffffffffe70001, 0x7ffffffffe10001, 0x7ffffffffcc0001, // 59 + 59 + 59 bits
			0x400000000270001, 0x400000000350001, 0x400000000360001, // 58 + 58 + 58 bits
			0x3ffffffffc10001, 0x3ffffffffbe0001, 0x3ffffffffbd0001, // 58 + 58 + 58 bits
			0x4000000004d0001, 0x400000000570001, 0x400000000660001}, // 58 + 58 + 58 bits
		P:     []uint64{0xffffffffffc0001, 0x10000000001d0001, 0x10000000006e0001}, // 60 + 60 + 60 bits
		Sigma: DefaultSigma,
	}
)

// NewParameters returns a new set of generic RLWE parameters from the given ring degree logn, moduli q and p, and
// error distribution parameter sigma. It returns the empty parameters Parameters{} and a non-nil error if the
// specified parameters are invalid.
func NewParameters(logn int, q, p []uint64, sigma float64) (Parameters, error) {
	var err error
	if err = checkSizeParams(logn, len(q), len(p)); err != nil {
		return Parameters{}, err
	}

	// Checks if moduli are valid
	if err = CheckModuli(q, p, logn); err != nil {
		return Parameters{}, err
	}

	params := Parameters{
		logN:  logn,
		pi:    make([]uint64, len(p)),
		qi:    make([]uint64, len(q)),
		sigma: sigma,
	}

	if params.ringQ, err = ring.NewRing(1<<logn, q); err != nil {
		return Parameters{}, err
	}

	if len(p) != 0 {
		if params.ringP, err = ring.NewRing(1<<logn, p); err != nil {
			return Parameters{}, err
		}
	}

	if params.ringQP, err = ring.NewRing(1<<logn, append(q, p...)); err != nil {
		return Parameters{}, err
	}

	copy(params.qi, q)
	copy(params.pi, p)
	return params, nil
}

// NewParametersFromLiteral instantiate a set of generic RLWE parameters from a ParametersLiteral specification.
// It returns the empty parameters Parameters{} and a non-nil error if the specified parameters are invalid.
func NewParametersFromLiteral(paramDef ParametersLiteral) (Parameters, error) {
	switch {
	case paramDef.Q != nil && paramDef.LogQ == nil && paramDef.P != nil && paramDef.LogP == nil:
		return NewParameters(paramDef.LogN, paramDef.Q, paramDef.P, paramDef.Sigma)
	case paramDef.LogQ != nil && paramDef.Q == nil && paramDef.LogP != nil && paramDef.P == nil:
		q, p, err := GenModuli(paramDef.LogN, paramDef.LogQ, paramDef.LogP)
		if err != nil {
			return Parameters{}, err
		}
		return NewParameters(paramDef.LogN, q, p, paramDef.Sigma)
	default:
		return Parameters{}, fmt.Errorf("invalid parameter literal")
	}
}

// N returns the ring degree
func (p Parameters) N() int {
	return 1 << p.logN
}

// LogN returns the log of the degree of the polynomial ring
func (p Parameters) LogN() int {
	return p.logN
}

// RingQ returns a pointer to ringQ
func (p Parameters) RingQ() *ring.Ring {
	return p.ringQ
}

// RingP returns a pointer to ringP
func (p Parameters) RingP() *ring.Ring {
	return p.ringP
}

// RingQP returns a pointer to ringQP
func (p Parameters) RingQP() *ring.Ring {
	return p.ringQP
}

// Sigma returns standard deviation of the noise distribution
func (p Parameters) Sigma() float64 {
	return p.sigma
}

// MaxLevel returns the maximum level of a ciphertext
func (p Parameters) MaxLevel() int {
	return p.QCount() - 1
}

// Q returns a new slice with the factors of the ciphertext modulus q
func (p Parameters) Q() []uint64 {
	qi := make([]uint64, len(p.qi))
	copy(qi, p.qi)
	return qi
}

// QCount returns the number of factors of the ciphertext modulus Q
func (p Parameters) QCount() int {
	return len(p.qi)
}

// QBigInt return the ciphertext-space modulus Q in big.Integer, reconstructed, representation.
func (p Parameters) QBigInt() *big.Int {
	q := big.NewInt(1)
	for _, qi := range p.qi {
		q.Mul(q, new(big.Int).SetUint64(qi))
	}
	return q
}

// P returns a new slice with the factors of the ciphertext modulus extension P
func (p Parameters) P() []uint64 {
	pi := make([]uint64, len(p.pi))
	copy(pi, p.pi)
	return pi
}

// PCount returns the number of factors of the ciphertext modulus extension P
func (p Parameters) PCount() int {
	return len(p.pi)
}

// PBigInt return the ciphertext-space extention modulus P in big.Integer, reconstructed, representation.
func (p Parameters) PBigInt() *big.Int {
	pInt := big.NewInt(1)
	for _, pi := range p.pi {
		pInt.Mul(pInt, new(big.Int).SetUint64(pi))
	}
	return pInt
}

// QP return the extended ciphertext-space modulus QP in RNS representation.
func (p Parameters) QP() []uint64 {
	qp := make([]uint64, len(p.qi)+len(p.pi))
	copy(qp, p.qi)
	copy(qp[len(p.qi):], p.pi)
	return qp
}

// QPCount returns the number of factors of the ciphertext modulus + the modulus extension P
func (p Parameters) QPCount() int {
	return len(p.qi) + len(p.pi)
}

// QPBigInt return the extended ciphertext-space modulus QP in big.Integer, reconstructed, representation.
func (p Parameters) QPBigInt() *big.Int {
	pqInt := p.QBigInt()
	pqInt.Mul(pqInt, p.PBigInt())
	return pqInt
}

// LogQ returns the size of the extended modulus Q in bits
func (p Parameters) LogQ() int {
	tmp := ring.NewUint(1)
	for _, qi := range p.qi {
		tmp.Mul(tmp, ring.NewUint(qi))
	}
	return tmp.BitLen()
}

// LogP returns the size of the extended modulus P in bits
func (p Parameters) LogP() int {
	tmp := ring.NewUint(1)
	for _, pi := range p.pi {
		tmp.Mul(tmp, ring.NewUint(pi))
	}
	return tmp.BitLen()
}

// LogQP returns the size of the extended modulus QP in bits
func (p Parameters) LogQP() int {
	tmp := ring.NewUint(1)
	for _, qi := range p.qi {
		tmp.Mul(tmp, ring.NewUint(qi))
	}
	for _, pi := range p.pi {
		tmp.Mul(tmp, ring.NewUint(pi))
	}
	return tmp.BitLen()
}

// Alpha returns the number of moduli in in P
func (p Parameters) Alpha() int {
	return p.PCount()
}

// Beta returns the number of element in the RNS decomposition basis: Ceil(lenQi / lenPi)
func (p Parameters) Beta() int {
	if p.Alpha() != 0 {
		return int(math.Ceil(float64(p.QCount()) / float64(p.Alpha())))
	}

	return 1
}

// QiOverflowMargin returns floor(2^64 / max(Qi)), i.e. the number of times elements of Z_max{Qi} can
// be added together before overflowing 2^64.
func (p *Parameters) QiOverflowMargin(level int) int {
	return int(math.Exp2(64) / float64(utils.MaxSliceUint64(p.qi[:level+1])))
}

// PiOverflowMargin returns floor(2^64 / max(Pi)), i.e. the number of times elements of Z_max{Pi} can
// be added together before overflowing 2^64.
func (p *Parameters) PiOverflowMargin() int {
	return int(math.Exp2(64) / float64(utils.MaxSliceUint64(p.pi)))
}

// GaloisElementForColumnRotationBy returns the galois element for plaintext
// column rotations by k position to the left. Providing a negative k is
// equivalent to a right rotation.
func (p Parameters) GaloisElementForColumnRotationBy(k int) uint64 {
	twoN := 1 << (p.logN + 1)
	mask := twoN - 1
	kRed := k & mask
	return ring.ModExp(GaloisGen, kRed, uint64(twoN))
}

// GaloisElementForRowRotation returns the galois element for generating the row
// rotation automorphism
func (p Parameters) GaloisElementForRowRotation() uint64 {
	return (1 << (p.logN + 1)) - 1
}

// GaloisElementsForRowInnerSum returns a list of all galois elements required to
// perform an InnerSum operation. This corresponds to all the left rotations by
// k-positions where k is a power of two and the row-rotation element.
func (p Parameters) GaloisElementsForRowInnerSum() (galEls []uint64) {
	galEls = make([]uint64, p.logN+1)
	galEls[0] = p.GaloisElementForRowRotation()
	for i := 0; i < int(p.logN)-1; i++ {
		galEls[i+1] = p.GaloisElementForColumnRotationBy(1 << i)
	}
	return galEls
}

// InverseGaloisElement takes a galois element and returns the galois element
//  corresponding to the inverse automorphism
func (p Parameters) InverseGaloisElement(galEl uint64) uint64 {
	twoN := 1 << (p.logN + 1)
	return ring.ModExp(galEl, twoN-1, uint64(twoN))
}

// Equals checks two Parameter structs for equality.
func (p Parameters) Equals(other Parameters) bool {
	res := p.logN == other.logN
	res = res && utils.EqualSliceUint64(p.qi, other.qi)
	res = res && utils.EqualSliceUint64(p.pi, other.pi)
	res = res && (p.sigma == other.sigma)
	return res
}

// CopyNew makes a deep copy of the receiver and returns it.
func (p Parameters) CopyNew() Parameters {
	qi, pi := p.qi, p.pi
	p.qi, p.pi = make([]uint64, len(p.qi)), make([]uint64, len(p.pi))
	copy(p.qi, qi)
	copy(p.pi, pi)
	return p
}

// MarshalBinary returns a []byte representation of the parameter set.
func (p Parameters) MarshalBinary() ([]byte, error) {
	if p.LogN() == 0 { // if N is 0, then p is the zero value
		return []byte{}, nil
	}

	// 1 byte : logN
	// 1 byte : #Q
	// 1 byte : #P
	// 8 byte : sigma
	// 8 * (#Q) : Q
	// 8 * (#P) : P
	b := utils.NewBuffer(make([]byte, 0, p.MarshalBinarySize()))
	b.WriteUint8(uint8(p.logN))
	b.WriteUint8(uint8(len(p.qi)))
	b.WriteUint8(uint8(len(p.pi)))
	b.WriteUint64(math.Float64bits(p.sigma))
	b.WriteUint64Slice(p.qi)
	b.WriteUint64Slice(p.pi)
	return b.Bytes(), nil
}

// UnmarshalBinary decodes a []byte into a parameter set struct.
func (p *Parameters) UnmarshalBinary(data []byte) error {
	if len(data) < 11 {
		return fmt.Errorf("invalid rlwe.Parameter serialization")
	}
	b := utils.NewBuffer(data)
	logN := int(b.ReadUint8())
	lenQ := int(b.ReadUint8())
	lenP := int(b.ReadUint8())
	sigma := math.Float64frombits(b.ReadUint64())

	if err := checkSizeParams(logN, lenQ, lenP); err != nil {
		return err
	}

	qi := make([]uint64, lenQ)
	pi := make([]uint64, lenP)
	b.ReadUint64Slice(qi)
	b.ReadUint64Slice(pi)

	var err error
	*p, err = NewParameters(logN, qi, pi, sigma)
	return err
}

// MarshalBinarySize returns the length of the []byte encoding of the reciever.
func (p Parameters) MarshalBinarySize() int {
	return 11 + (len(p.qi)+len(p.pi))<<3
}

// MarshalJSON returns a JSON representation of this parameter set. See `Marshal` from the `encoding/json` package.
func (p Parameters) MarshalJSON() ([]byte, error) {
	return json.Marshal(&ParametersLiteral{LogN: p.logN, Q: p.qi, P: p.pi, Sigma: p.sigma})
}

// UnmarshalJSON reads a JSON representation of a parameter set into the receiver Parameter. See `Unmarshal` from the `encoding/json` package.
func (p *Parameters) UnmarshalJSON(data []byte) (err error) {
	var params ParametersLiteral
	if err = json.Unmarshal(data, &params); err != nil {
		return err
	}
	*p, err = NewParametersFromLiteral(params)
	return
}

// CheckModuli checks that the provided q and p correspond to a valid moduli chain
func CheckModuli(q, p []uint64, logN int) error {

	if len(q) > MaxModuliCount {
		return fmt.Errorf("#Qi is larger than %d", MaxModuliCount)
	}

	if len(p) > MaxModuliCount {
		return fmt.Errorf("#Pi is larger than %d", MaxModuliCount)
	}

	for i, qi := range q {
		if uint64(bits.Len64(qi)-1) > MaxModuliSize+1 {
			return fmt.Errorf("Qi bit-size (i=%d) is larger than %d", i, MaxModuliSize)
		}
	}

	for i, pi := range p {
		if uint64(bits.Len64(pi)-1) > MaxModuliSize+2 {
			return fmt.Errorf("Pi bit-size (i=%d) is larger than %d", i, MaxModuliSize)
		}
	}

	N := uint64(1 << logN)

	for i, qi := range q {
		if !ring.IsPrime(qi) || qi&((N<<1)-1) != 1 {
			return fmt.Errorf("Qi (i=%d) is not an NTT prime", i)
		}
	}

	for i, pi := range p {
		if !ring.IsPrime(pi) || pi&((N<<1)-1) != 1 {
			return fmt.Errorf("Pi (i=%d) is not an NTT prime", i)
		}
	}

	return nil
}

func checkSizeParams(logN int, lenQ, lenP int) error {
	if logN > MaxLogN {
		return fmt.Errorf("logN=%d is larger than MaxLogN=%d", logN, MaxLogN)
	}
	if logN < MinLogN {
		return fmt.Errorf("logN=%d is smaller than MinLogN=%d", logN, MinLogN)
	}
	if lenQ > MaxModuliCount {
		return fmt.Errorf("lenQ=%d is larger than MaxModuliCount=%d", lenQ, MaxModuliCount)
	}
	if lenP > MaxModuliCount {
		return fmt.Errorf("lenP=%d is larger than MaxModuliCount=%d", lenP, MaxModuliCount)
	}
	return nil
}

func checkModuliLogSize(logQ, logP []int) error {

	for i, qi := range logQ {
		if qi <= 0 || qi > MaxModuliSize {
			return fmt.Errorf("logQ[%d]=%d is not in ]0, %d]", i, qi, MaxModuliSize)
		}
	}

	for i, pi := range logP {
		if pi <= 0 || pi > MaxModuliSize+1 {
			return fmt.Errorf("logP[%d]=%d is not in ]0,%d]", i, pi, MaxModuliSize+1)
		}
	}

	return nil
}

// GenModuli generates a valid moduli chain from the provided moduli sizes.
func GenModuli(logN int, logQ, logP []int) (q, p []uint64, err error) {

	if err = checkSizeParams(logN, len(logQ), len(logP)); err != nil {
		return
	}

	if err = checkModuliLogSize(logQ, logP); err != nil {
		return
	}

	// Extracts all the different primes bit size and maps their number
	primesbitlen := make(map[int]int)
	for _, qi := range logQ {
		primesbitlen[qi]++
	}

	for _, pj := range logP {
		primesbitlen[pj]++
	}

	// For each bit-size, finds that many primes
	primes := make(map[int][]uint64)
	for key, value := range primesbitlen {
		primes[key] = ring.GenerateNTTPrimes(int(key), 2<<logN, int(value))
	}

	// Assigns the primes to the moduli chain
	for _, qi := range logQ {
		q = append(q, primes[qi][0])
		primes[qi] = primes[qi][1:]
	}

	// Assigns the primes to the special primes list for the extended ring
	for _, pj := range logP {
		p = append(p, primes[pj][0])
		primes[pj] = primes[pj][1:]
	}

	return
}
