package rckks

import (
	"errors"
	"fmt"
	"math"
	"math/big"
	"math/bits"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// MaxLogN is the log2 of the largest supported polynomial modulus degree.
const MaxLogN = 16

// MaxModuliCount is the largest supported number of moduli in the RNS representation.
const MaxModuliCount = 34

// MaxModuliSize is the largest bit-length supported for the moduli in the RNS representation.
const MaxModuliSize = 60

// DefaultSigma is the default error distribution standard deviation
const DefaultSigma = 3.2

// Name of the different default parameter sets
const (
	// PN12QP109 is the index in DefaultParams for logQP = 109
	PN12QP109 = iota
	// PN13QP218 is the index in DefaultParams for logQP = 218
	PN13QP218
	// PN14QP438 is the index in DefaultParams for logQP = 438
	PN14QP438
	// PN15QP880 is the index in DefaultParams for logQP = 880
	PN15QP880
	// PN16QP1761 is the index in DefaultParams for logQP = 1761
	PN16QP1761

	// PN12QP101pq is the index in DefaultParams for logQP = 101 (post quantum)
	PN12QP101pq
	// PN13QP202pq is the index in DefaultParams for logQP = 202 (post quantum)
	PN13QP202pq
	// PN14QP411pq is the index in DefaultParams for logQP = 411 (post quantum)
	PN14QP411pq
	// PN15QP827pq is the index in DefaultParams for logQP = 827 (post quantum)
	PN15QP827pq
	// PN16QP1654pq is the index in DefaultParams for logQP = 1654 (post quantum)
	PN16QP1654pq
)

// DefaultParams is a set of default CKKS parameters ensuring 128 bit security.
var DefaultParams = []*Parameters{

	//LogQi = 109
	{logN: 12,
		logSlots: 12,
		qi: []uint64{0x1ffffe0001, // 37 + 32
			0x100014001},
		pi:    []uint64{0x4000038001}, // 38
		scale: 1 << 32,
		sigma: DefaultSigma,
	},

	//LogQi = 218
	{logN: 13,
		logSlots: 13,
		qi: []uint64{0x200038001, // 33 + 5 x 30
			0x3ffe8001,
			0x40020001,
			0x40038001,
			0x3ffc0001,
			0x40080001,
		},
		pi:    []uint64{0x800008001}, // 35
		scale: 1 << 30,
		sigma: DefaultSigma,
	},

	//LogQiP = 438
	{logN: 14,
		logSlots: 14, // 45 + 9 x 34
		qi: []uint64{0x2000000a0001, 0x3fffd0001,
			0x400060001, 0x3fff90001, 0x400080001,
			0x400180001, 0x3ffd20001, 0x400300001,
			0x400360001, 0x4003e0001},
		pi:    []uint64{0x80000050001, 0x7ffffdb0001}, // 43, 43
		scale: 1 << 34,
		sigma: DefaultSigma,
	},

	//LogQi = 880
	{logN: 15,
		logSlots: 15, // 50 + 17 x 40
		qi: []uint64{0x4000000120001, 0x10000140001, 0xffffe80001, 0xffffc40001, 0x100003e0001,
			0xffffb20001, 0x10000500001, 0xffff940001, 0xffff8a0001, 0xffff820001,
			0xffff780001, 0x10000960001, 0x10000a40001, 0xffff580001, 0x10000b60001,
			0xffff480001, 0xffff420001, 0xffff340001},
		pi:    []uint64{0x3ffffffd20001, 0x4000000420001, 0x3ffffffb80001}, // 50, 50, 50
		scale: 1 << 40,
		sigma: DefaultSigma,
	},

	//LogQi = 1761
	{logN: 16,
		logSlots: 16, // 55 + 33 x 45
		qi: []uint64{0x80000000080001, 0x200000440001, 0x200000500001, 0x1fffff980001, 0x200000c80001,
			0x1ffffeb40001, 0x1ffffe640001, 0x200001a00001, 0x200001e80001, 0x1ffffe0c0001,
			0x200002480001, 0x200002800001, 0x1ffffd800001, 0x200002900001, 0x1ffffd700001,
			0x2000029c0001, 0x1ffffcf00001, 0x200003140001, 0x1ffffcc80001, 0x1ffffcb40001,
			0x1ffffc980001, 0x200003740001, 0x200003800001, 0x200003d40001, 0x1ffffc200001,
			0x1ffffc140001, 0x200004100001, 0x200004180001, 0x1ffffbc40001, 0x200004700001,
			0x1ffffb900001, 0x200004cc0001, 0x1ffffb240001, 0x200004e80001},
		pi:    []uint64{0x80000000440001, 0x80000000500001, 0x7fffffff380001, 0x80000000e00001}, // 4 x 55
		scale: 1 << 45,
		sigma: DefaultSigma,
	},

	//LogQi = 100.99992948699881
	{logN: 12,
		logSlots: 12, // 35 + 30
		qi:       []uint64{0x800004001, 0x3fff4001},
		pi:       []uint64{0xffffc4001}, // 36
		scale:    1 << 30,
		sigma:    DefaultSigma,
	},

	//LogQi = 201.99146237773812
	{logN: 13,
		logSlots: 13, // 33 + 5 x 27
		qi:       []uint64{0x200038001, 0x8008001, 0x8020001, 0x7fa8001, 0x80e0001, 0x7e90001},
		pi:       []uint64{0x400018001}, // 34
		scale:    1 << 27,
		sigma:    DefaultSigma,
	},

	//LogQiP = 411.0002505571375
	{logN: 14,
		logSlots: 14, // 40 + 9 x 33
		qi: []uint64{0x10000140001, 0x1fff90001, 0x200080001, 0x1fff60001, 0x200100001,
			0x1fff00001, 0x1ffef0001, 0x1ffe60001, 0x2001d0001, 0x2002e0001},

		pi:    []uint64{0x1ffffe0001, 0x1ffffc0001}, // 37, 37
		scale: 1 << 33,
		sigma: DefaultSigma,
	},

	//LogQi = 827.0000336017748
	{logN: 15,
		logSlots: 15, // 46 + 17 x 38
		qi: []uint64{0x400000060001, 0x3fffe80001, 0x4000300001, 0x3fffb80001, 0x40004a0001,
			0x3fffb20001, 0x4000540001, 0x4000560001, 0x3fff900001, 0x4000720001,
			0x3fff8e0001, 0x4000800001, 0x40008a0001, 0x3fff6c0001, 0x40009e0001,
			0x3fff300001, 0x3fff1c0001, 0x4000fc0001},
		pi:    []uint64{0x2000000a0001, 0x2000000e0001, 0x1fffffc20001}, // 3 x 45
		scale: 1 << 38,
		sigma: DefaultSigma,
	},

	//LogQi = 1654.0000043447205
	{logN: 16,
		logSlots: 16, // 55 + 31 x 45
		qi: []uint64{0x80000000080001, 0x200000440001, 0x200000500001, 0x1fffff980001, 0x200000c80001,
			0x1ffffeb40001, 0x1ffffe640001, 0x200001a00001, 0x200001e80001, 0x1ffffe0c0001,
			0x200002480001, 0x200002800001, 0x1ffffd800001, 0x200002900001, 0x1ffffd700001,
			0x2000029c0001, 0x1ffffcf00001, 0x200003140001, 0x1ffffcc80001, 0x1ffffcb40001,
			0x1ffffc980001, 0x200003740001, 0x200003800001, 0x200003d40001, 0x1ffffc200001,
			0x1ffffc140001, 0x200004100001, 0x200004180001, 0x1ffffbc40001, 0x200004700001,
			0x1ffffb900001, 0x200004cc0001},
		pi:    []uint64{0x80000001c0001, 0x80000002c0001, 0x8000000500001, 0x7ffffff9c0001}, // 4 x 51
		scale: 1 << 45,
		sigma: DefaultSigma,
	},
}

// Moduli stores the NTT primes of the RNS representation.
type Moduli struct {
	Qi []uint64 // Ciphertext prime moduli
	Pi []uint64 // Keys additional prime moduli
}

// Print prints the moduli in hexadimal
func (m *Moduli) Print() {
	for _, qi := range m.Qi {
		fmt.Printf("0x%x,\n", qi)
	}
	fmt.Println()

	for _, pj := range m.Pi {
		fmt.Printf("0x%x,\n", pj)
	}
	fmt.Println()
}

// Copy creates a copy of the target Moduli.
func (m *Moduli) Copy() Moduli {

	qi := make([]uint64, len(m.Qi))
	copy(qi, m.Qi)

	pi := make([]uint64, len(m.Pi))
	copy(pi, m.Pi)

	return Moduli{qi, pi}
}

// LogModuli stores the bit-length of the NTT primes of the RNS representation.
type LogModuli struct {
	LogQi []uint64 // Ciphertext prime moduli bit-size
	LogPi []uint64 // Keys additional prime moduli bit-size
}

// Copy creates a copy of the target Moduli.
func (m *LogModuli) Copy() LogModuli {

	LogQi := make([]uint64, len(m.LogQi))
	copy(LogQi, m.LogQi)

	LogPi := make([]uint64, len(m.LogPi))
	copy(LogPi, m.LogPi)

	return LogModuli{LogQi, LogPi}
}

// Parameters represents a given parameter set for the BFV cryptosystem.
type Parameters struct {
	qi       []uint64
	pi       []uint64
	logN     uint64 // Ring degree (power of 2)
	logSlots uint64
	scale    float64
	sigma    float64 // Gaussian sampling variance
}

// NewParametersFromModuli creates a new Parameters struct and returns a pointer to it.
func NewParametersFromModuli(logN uint64, m *Moduli) (p *Parameters, err error) {
	p = new(Parameters)

	if (logN < 3) || (logN > MaxLogN) {
		return nil, fmt.Errorf("invalid polynomial ring log degree: %d", logN)
	}

	p.logN = logN

	if err = checkModuli(m, logN); err != nil {
		return nil, err
	}

	p.qi = make([]uint64, len(m.Qi), len(m.Qi))
	copy(p.qi, m.Qi)
	p.pi = make([]uint64, len(m.Pi), len(m.Pi))
	copy(p.pi, m.Pi)

	p.sigma = DefaultSigma

	return p, nil

}

// NewParametersFromLogModuli creates a new Parameters struct and returns a pointer to it.
func NewParametersFromLogModuli(logN uint64, lm *LogModuli) (p *Parameters, err error) {

	if err = checkLogModuli(lm); err != nil {
		return nil, err
	}

	// If LogModuli is valid and then generates the moduli
	return NewParametersFromModuli(logN, genModuli(lm, logN+1))
}

// NewPolyQ returns a new empty polynomial of degree 2^LogN in basis Qi.
func (p *Parameters) NewPolyQ() *ring.Poly {
	return ring.NewPoly(p.N(), p.QiCount())
}

// NewPolyP returns a new empty polynomial of degree 2^LogN in basis Pi.
func (p *Parameters) NewPolyP() *ring.Poly {
	return ring.NewPoly(p.N(), p.PiCount())
}

// NewPolyQP returns a new empty polynomial of degree 2^LogN in basis Qi + Pi.
func (p *Parameters) NewPolyQP() *ring.Poly {
	return ring.NewPoly(p.N(), p.QPiCount())
}

// N returns the ring degree
func (p *Parameters) N() uint64 {
	return 1 << p.logN
}

// LogN returns the log of the degree of the polynomial ring
func (p *Parameters) LogN() uint64 {
	return p.logN
}

// LogSlots returns the log of the number of slots
func (p *Parameters) LogSlots() uint64 {
	return p.logSlots
}

// MaxLevel returns the maximum ciphertext level
func (p *Parameters) MaxLevel() uint64 {
	return p.QiCount() - 1
}

// Levels returns then number of total levels enabled by the parameters
func (p *Parameters) Levels() uint64 {
	return p.QiCount()
}

// Slots returns number of avaliable plaintext slots
func (p *Parameters) Slots() uint64 {
	return 1 << p.logSlots
}

// MaxSlots returns the theoretical maximum of plaintext slots allowed by the ring degree
func (p *Parameters) MaxSlots() uint64 {
	return p.N()
}

// MaxLogSlots returns the log of the maximum number of slots enabled by the parameters
func (p *Parameters) MaxLogSlots() uint64 {
	return p.logN
}

// Sigma returns standard deviation of the noise distribution
func (p *Parameters) Sigma() float64 {
	return p.sigma
}

// Scale returns the default plaintext/ciphertext scale
func (p *Parameters) Scale() float64 {
	return p.scale
}

// SetScale sets the default plaintext/ciphertext scale
func (p *Parameters) SetScale(scale float64) {
	p.scale = scale
}

// SetLogSlots sets the value logSlots of the parameters.
func (p *Parameters) SetLogSlots(logSlots uint64) (err error) {
	if (logSlots == 0) || (logSlots > p.MaxLogSlots()) {
		return fmt.Errorf("slots cannot be greater than LogN-1")
	}

	p.logSlots = logSlots

	return nil
}

// SetSigma sets the value sigma of the parameters
func (p *Parameters) SetSigma(sigma float64) {
	p.sigma = sigma
}

// LogModuli generates a LogModuli struct from the parameters' Moduli struct and returns it.
func (p *Parameters) LogModuli() (lm *LogModuli) {

	lm = new(LogModuli)

	lm.LogQi = make([]uint64, len(p.qi), len(p.qi))
	for i := range p.qi {
		lm.LogQi[i] = uint64(math.Round(math.Log2(float64(p.qi[i]))))
	}

	lm.LogPi = make([]uint64, len(p.pi), len(p.pi))
	for i := range p.pi {
		lm.LogPi[i] = uint64(math.Round(math.Log2(float64(p.pi[i]))))
	}

	return
}

// Moduli returns a struct Moduli with the moduli of the parameters
func (p *Parameters) Moduli() (m *Moduli) {
	m = new(Moduli)
	m.Qi = make([]uint64, p.QiCount(), p.QiCount())
	m.Pi = make([]uint64, p.PiCount(), p.PiCount())
	copy(m.Qi, p.qi)
	copy(m.Pi, p.pi)
	return
}

// Qi returns a new slice with the factors of the ciphertext modulus q
func (p *Parameters) Qi() []uint64 {
	qi := make([]uint64, len(p.qi))
	copy(qi, p.qi)
	return qi
}

// QiCount returns the number of factors of the ciphertext modulus Q
func (p *Parameters) QiCount() uint64 {
	return uint64(len(p.qi))
}

// Pi returns a new slice with the factors of the ciphertext modulus extention P
func (p *Parameters) Pi() []uint64 {
	pi := make([]uint64, len(p.pi))
	copy(pi, p.pi)
	return pi
}

// PiCount returns the number of factors of the ciphertext modulus extention P
func (p *Parameters) PiCount() uint64 {
	return uint64(len(p.pi))
}

// QPiCount returns the number of factors of the ciphertext modulus + the extention modulus P
func (p *Parameters) QPiCount() uint64 {
	return uint64(len(p.qi) + len(p.pi))
}

// LogQP returns the size of the extended modulus QP in bits
func (p *Parameters) LogQP() uint64 {
	tmp := ring.NewUint(1)
	for _, qi := range p.qi {
		tmp.Mul(tmp, ring.NewUint(qi))
	}
	for _, pi := range p.pi {
		tmp.Mul(tmp, ring.NewUint(pi))
	}
	return uint64(tmp.BitLen())
}

// LogQLvl returns the size of the modulus Q in bits at a specific level
func (p *Parameters) LogQLvl(level uint64) uint64 {
	tmp := p.QLvl(level)
	return uint64(tmp.BitLen())
}

// QLvl returns the product of the moduli at the given level as a big.Int
func (p *Parameters) QLvl(level uint64) *big.Int {
	tmp := ring.NewUint(1)
	for _, qi := range p.qi[:level+1] {
		tmp.Mul(tmp, ring.NewUint(qi))
	}
	return tmp
}

// LogQ returns the size of the modulus Q in bits
func (p *Parameters) LogQ() uint64 {
	return p.LogQLvl(p.QiCount() - 1)
}

// Q returns the product of all the moduli as a big.Int
func (p *Parameters) Q() *big.Int {
	return p.QLvl(p.QiCount() - 1)
}

// LogP returns the size of the modulus P in bits
func (p *Parameters) LogP() uint64 {
	tmp := ring.NewUint(1)
	for _, pi := range p.pi {
		tmp.Mul(tmp, ring.NewUint(pi))
	}
	return uint64(tmp.BitLen())
}

// LogQAlpha returns the size in bits of the sum of the norm of
// each element of the special RNS decomposition basis for the
// key-switching.
// LogQAlpha is the size of the element that is multipled by the
// error during the keyswitching and then divided by P.
// LogQAlpha should be smaller than P or the error added during
// the key-switching wont be negligible.
func (p *Parameters) LogQAlpha() uint64 {

	alpha := p.PiCount()

	if alpha == 0 {
		return 0
	}

	res := ring.NewUint(0)
	var j uint64
	for i := uint64(0); i < p.QiCount(); i = i + alpha {

		j = i + alpha
		if j > p.QiCount() {
			j = p.QiCount()
		}

		tmp := ring.NewUint(1)
		for _, qi := range p.pi[i:j] {
			tmp.Mul(tmp, ring.NewUint(qi))
		}

		res.Add(res, tmp)
	}

	return uint64(res.BitLen())
}

// Alpha returns the number of moduli in in P
func (p *Parameters) Alpha() uint64 {
	return p.PiCount()
}

// Beta returns the number of element in the RNS decomposition basis: Ceil(lenQi / lenPi)
func (p *Parameters) Beta() uint64 {
	if p.Alpha() != 0 {
		return uint64(math.Ceil(float64(p.QiCount()) / float64(p.Alpha())))
	}

	return 0
}

// Copy creates a copy of the target parameters.
func (p *Parameters) Copy() (paramsCopy *Parameters) {

	paramsCopy = new(Parameters)
	paramsCopy.logN = p.logN
	paramsCopy.logSlots = p.logSlots
	paramsCopy.scale = p.scale
	paramsCopy.sigma = p.sigma
	paramsCopy.qi = make([]uint64, len(p.qi), len(p.qi))
	copy(paramsCopy.qi, p.qi)
	paramsCopy.pi = make([]uint64, len(p.pi), len(p.pi))
	copy(paramsCopy.pi, p.pi)
	return
}

// Equals compares two sets of parameters for equality.
func (p *Parameters) Equals(other *Parameters) (res bool) {

	if p == other {
		return true
	}

	res = p.logN == other.logN
	res = res && (p.logSlots == other.logSlots)
	res = res && (p.scale == other.scale)
	res = res && (p.sigma == other.sigma)
	res = res && utils.EqualSliceUint64(p.qi, other.qi)
	res = res && utils.EqualSliceUint64(p.pi, other.pi)
	return
}

// MarshalBinary returns a []byte representation of the parameter set.
func (p *Parameters) MarshalBinary() ([]byte, error) {
	if p.logN == 0 { // if N is 0, then p is the zero value
		return []byte{}, nil
	}

	b := utils.NewBuffer(make([]byte, 0, 21+(p.QPiCount())<<3))

	b.WriteUint8(uint8(p.logN))
	b.WriteUint8(uint8(p.logSlots))
	b.WriteUint64(math.Float64bits(p.scale))
	b.WriteUint64(math.Float64bits(p.sigma))
	b.WriteUint8(uint8(p.QiCount()))
	b.WriteUint8(uint8(p.PiCount()))
	b.WriteUint64Slice(p.qi)
	b.WriteUint64Slice(p.pi)

	return b.Bytes(), nil
}

// UnmarshalBinary decodes a []byte into a parameter set struct
func (p *Parameters) UnmarshalBinary(data []byte) (err error) {

	if len(data) < 3 {
		return errors.New("invalid parameters encoding")
	}

	b := utils.NewBuffer(data)

	p.logN = uint64(b.ReadUint8())

	if p.logN > MaxLogN {
		return fmt.Errorf("LogN larger than %d", MaxLogN)
	}

	p.logSlots = uint64(b.ReadUint8())

	if p.logSlots > p.logN {
		return fmt.Errorf("LogSlots larger than %d", MaxLogN-1)
	}

	p.scale = math.Float64frombits(b.ReadUint64())
	p.sigma = math.Float64frombits(b.ReadUint64())

	lenQi := b.ReadUint8()
	lenPi := b.ReadUint8()

	p.qi = make([]uint64, lenQi, lenQi)
	p.pi = make([]uint64, lenPi, lenPi)

	b.ReadUint64Slice(p.qi)
	b.ReadUint64Slice(p.pi)

	if err = checkModuli(p.Moduli(), p.logN); err != nil {
		return err
	}

	return nil
}

func checkModuli(m *Moduli, logN uint64) error {

	if len(m.Qi) > MaxModuliCount {
		return fmt.Errorf("#Qi is larger than %d", MaxModuliCount)
	}

	if len(m.Pi) > MaxModuliCount {
		return fmt.Errorf("#Pi is larger than %d", MaxModuliCount)
	}

	for i, qi := range m.Qi {
		if uint64(bits.Len64(qi)-1) > MaxModuliSize+1 {
			return fmt.Errorf("Qi bit-size (i=%d) is larger than %d", i, MaxModuliSize)
		}
	}

	for i, pi := range m.Pi {
		if uint64(bits.Len64(pi)-1) > MaxModuliSize+2 {
			return fmt.Errorf("Pi bit-size (i=%d) is larger than %d", i, MaxModuliSize)
		}
	}

	N := uint64(1 << logN)

	for i, qi := range m.Qi {
		if !ring.IsPrime(qi) || qi&((N<<1)-1) != 1 {
			return fmt.Errorf("Qi (i=%d) is not an NTT prime", i)
		}
	}

	for i, pi := range m.Pi {
		if !ring.IsPrime(pi) || pi&((N<<1)-1) != 1 {
			return fmt.Errorf("Pi (i=%d) is not an NTT prime", i)
		}
	}

	return nil
}

func checkLogModuli(m *LogModuli) error {

	if len(m.LogQi) > MaxModuliCount {
		return fmt.Errorf("#LogQi is larger than %d", MaxModuliCount)
	}

	if len(m.LogPi) > MaxModuliCount {
		return fmt.Errorf("#LogPi is larger than %d", MaxModuliCount)
	}

	for i, qi := range m.LogQi {
		if qi > MaxModuliSize {
			return fmt.Errorf("LogQi (i=%d) is larger than %d", i, MaxModuliSize)
		}
	}

	for i, pi := range m.LogPi {
		if pi > MaxModuliSize+1 {
			return fmt.Errorf("LogPi (i=%d) is larger than %d", i, MaxModuliSize)
		}
	}

	return nil
}

func genModuli(lm *LogModuli, logN uint64) (m *Moduli) {

	m = new(Moduli)

	// Extracts all the different primes bit size and maps their number
	primesbitlen := make(map[uint64]uint64)
	for _, qi := range lm.LogQi {
		primesbitlen[qi]++
	}

	for _, pj := range lm.LogPi {
		primesbitlen[pj]++
	}

	// For each bit-size, finds that many primes
	primes := make(map[uint64][]uint64)
	for key, value := range primesbitlen {
		primes[key] = ring.GenerateNTTPrimes(key, logN, value)
	}

	// Assigns the primes to the ckks moduli chain
	m.Qi = make([]uint64, len(lm.LogQi))
	for i, qi := range lm.LogQi {
		m.Qi[i] = primes[qi][0]
		primes[qi] = primes[qi][1:]
	}

	// Assigns the primes to the special primes list for the extended ring
	m.Pi = make([]uint64, len(lm.LogPi))
	for i, pj := range lm.LogPi {
		m.Pi[i] = primes[pj][0]
		primes[pj] = primes[pj][1:]
	}

	return m
}
