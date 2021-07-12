package ring

import (
	"encoding/binary"

	"github.com/ldsec/lattigo/v2/utils"
)

// UniformSampler wraps a util.PRNG and represents the state of a sampler of uniform polynomials.
type UniformSampler struct {
	baseSampler
	randomBufferN []byte
}

// NewUniformSampler creates a new instance of UniformSampler from a PRNG and ring definition.
func NewUniformSampler(prng utils.PRNG, baseRing *Ring) *UniformSampler {
	uniformSampler := new(UniformSampler)
	uniformSampler.baseRing = baseRing
	uniformSampler.prng = prng
	uniformSampler.randomBufferN = make([]byte, baseRing.N)
	return uniformSampler
}

// Read generates a new polynomial with coefficients following a uniform distribution over [0, Qi-1].
func (uniformSampler *UniformSampler) Read(Pol *Poly) {

	var randomUint, mask, qi uint64
	var ptr int

	uniformSampler.prng.Clock(uniformSampler.randomBufferN)

	for j := range uniformSampler.baseRing.Modulus {

		qi = uniformSampler.baseRing.Modulus[j]

		// Start by computing the mask
		mask = uniformSampler.baseRing.Mask[j]

		ptmp := Pol.Coeffs[j]

		// Iterate for each modulus over each coefficient
		for i := 0; i < uniformSampler.baseRing.N; i++ {

			// Sample an integer between [0, qi-1]
			for {

				// Refill the pool if it runs empty
				if ptr == uniformSampler.baseRing.N {
					uniformSampler.prng.Clock(uniformSampler.randomBufferN)
					ptr = 0
				}

				// Read bytes from the pool
				randomUint = binary.BigEndian.Uint64(uniformSampler.randomBufferN[ptr:ptr+8]) & mask
				ptr += 8

				// If the integer is between [0, qi-1], break the loop
				if randomUint < qi {
					break
				}
			}

			ptmp[i] = randomUint
		}
	}
}

// ReadLvl generates a new polynomial with coefficients following a uniform distribution over [0, Qi-1].
func (uniformSampler *UniformSampler) ReadLvl(level int, Pol *Poly) {

	var randomUint, mask, qi uint64
	var ptr int

	uniformSampler.prng.Clock(uniformSampler.randomBufferN)

	for j := 0; j < level+1; j++ {

		qi = uniformSampler.baseRing.Modulus[j]

		// Start by computing the mask
		mask = uniformSampler.baseRing.Mask[j]

		ptmp := Pol.Coeffs[j]

		// Iterate for each modulus over each coefficient
		for i := 0; i < uniformSampler.baseRing.N; i++ {

			// Sample an integer between [0, qi-1]
			for {

				// Refill the pool if it runs empty
				if ptr == uniformSampler.baseRing.N {
					uniformSampler.prng.Clock(uniformSampler.randomBufferN)
					ptr = 0
				}

				// Read bytes from the pool
				randomUint = binary.BigEndian.Uint64(uniformSampler.randomBufferN[ptr:ptr+8]) & mask
				ptr += 8

				// If the integer is between [0, qi-1], break the loop
				if randomUint < qi {
					break
				}
			}

			ptmp[i] = randomUint
		}
	}
}

// ReadNew generates a new polynomial with coefficients following a uniform distribution over [0, Qi-1].
// Polynomial is created at the max level.
func (uniformSampler *UniformSampler) ReadNew() (Pol *Poly) {
	Pol = uniformSampler.baseRing.NewPoly()
	uniformSampler.Read(Pol)
	return
}

// ReadLvlNew generates a new polynomial with coefficients following a uniform distribution over [0, Qi-1].
// Polynomial is created at the specified level.
func (uniformSampler *UniformSampler) ReadLvlNew(level int) (Pol *Poly) {
	Pol = uniformSampler.baseRing.NewPolyLvl(level)
	uniformSampler.ReadLvl(level, Pol)
	return
}

// RandUniform samples a uniform randomInt variable in the range [0, mask] until randomInt is in the range [0, v-1].
// mask needs to be of the form 2^n -1.
func RandUniform(prng utils.PRNG, v uint64, mask uint64) (randomInt uint64) {
	for {
		randomInt = randInt64(prng, mask)
		if randomInt < v {
			return randomInt
		}
	}
}

// randInt32 samples a uniform variable in the range [0, mask], where mask is of the form 2^n-1, with n in [0, 32].
func randInt32(prng utils.PRNG, mask uint64) uint64 {

	// generate random 4 bytes
	randomBytes := make([]byte, 4)
	prng.Clock(randomBytes)

	// convert 4 bytes to a uint32
	randomUint32 := uint64(binary.BigEndian.Uint32(randomBytes))

	// return required bits
	return mask & randomUint32
}

// randInt64 samples a uniform variable in the range [0, mask], where mask is of the form 2^n-1, with n in [0, 64].
func randInt64(prng utils.PRNG, mask uint64) uint64 {

	// generate random 8 bytes
	randomBytes := make([]byte, 8)
	prng.Clock(randomBytes)

	// convert 8 bytes to a uint64
	randomUint64 := binary.BigEndian.Uint64(randomBytes)

	// return required bits
	return mask & randomUint64
}
