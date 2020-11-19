package drckks

import (
	"github.com/ldsec/lattigo/v2/rckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
	"math/big"
)

// PermuteProtocol is a struct storing the parameters for the PermuteProtocol protocol.
type PermuteProtocol struct {
	drckksContext   *drckksContext
	encoder         rckks.EncoderBigComplex
	tmp             *ring.Poly
	maskBigint      []*big.Int
	maskFloat       []*big.Float
	maskComplex     []*ring.Complex
	gaussianSampler *ring.GaussianSampler
}

// NewPermuteProtocol creates a new instance of the PermuteProtocol.
func NewPermuteProtocol(params *rckks.Parameters) (pp *PermuteProtocol) {

	prec := uint64(256)

	pp = new(PermuteProtocol)
	pp.encoder = rckks.NewEncoderBigComplex(params, prec)
	drckksContext := newDrckksContext(params)
	pp.drckksContext = drckksContext
	pp.tmp = drckksContext.ringQ.NewPoly()
	pp.maskBigint = make([]*big.Int, drckksContext.n)
	pp.maskComplex = make([]*ring.Complex, drckksContext.n)

	for i := uint64(0); i < drckksContext.n; i++ {
		pp.maskComplex[i] = new(ring.Complex)
		pp.maskComplex[i][0] = new(big.Float)
		pp.maskComplex[i][1] = new(big.Float)
		pp.maskComplex[i][0].SetPrec(uint(prec))
		pp.maskComplex[i][1].SetPrec(uint(prec))
	}

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	pp.gaussianSampler = ring.NewGaussianSampler(prng, drckksContext.ringQ, params.Sigma(), uint64(6*params.Sigma()))

	return
}

// AllocateShares allocates the shares of the Refresh protocol.
func (pp *PermuteProtocol) AllocateShares(levelStart uint64) (RefreshShareDecrypt, RefreshShareRecrypt) {
	return pp.drckksContext.ringQ.NewPolyLvl(levelStart), pp.drckksContext.ringQ.NewPoly()
}

func (pp *PermuteProtocol) permuteWithIndex(permutation []uint64, values []*ring.Complex) {
	tmp := make([]*ring.Complex, len(values))

	for i := range values {
		tmp[i] = values[permutation[i]].Copy()
	}

	for i := range values {
		values[i] = tmp[i]
	}
}

// GenShares generates the decryption and recryption shares of the Refresh protocol.
func (pp *PermuteProtocol) GenShares(sk *ring.Poly, levelStart, nParties uint64, ciphertext *rckks.Ciphertext, crs *ring.Poly, slots uint64, permutation []uint64, shareDecrypt RefreshShareDecrypt, shareRecrypt RefreshShareRecrypt) {

	ringQ := pp.drckksContext.ringQ

	bound := ring.NewUint(ringQ.Modulus[0])
	for i := uint64(1); i < levelStart+1; i++ {
		bound.Mul(bound, ring.NewUint(ringQ.Modulus[i]))
	}

	bound.Quo(bound, ring.NewUint(2*nParties))
	boundHalf := new(big.Int).Rsh(bound, 1)

	maxSlots := pp.drckksContext.n

	// Samples the whole N coefficients for h0
	// Samples in the Slot domain
	var sign int
	for i := uint64(0); i < maxSlots; i++ {

		pp.maskBigint[i] = ring.RandInt(bound)

		sign = pp.maskBigint[i].Cmp(boundHalf)
		if sign == 1 || sign == 0 {
			pp.maskBigint[i].Sub(pp.maskBigint[i], bound)
		}

		if i == 0 {
			pp.maskComplex[i][0].SetInt(pp.maskBigint[i])
			pp.maskComplex[i][1].SetFloat64(0.0)
		} else {
			pp.maskComplex[i][0].SetInt(pp.maskBigint[i])
			pp.maskComplex[maxSlots-i][1].SetInt(pp.maskBigint[i])
			pp.maskComplex[maxSlots-i][1].Neg(pp.maskComplex[maxSlots-i][1])
		}
	}

	// h0 = mask (at level min)
	ringQ.SetCoefficientsBigintLvl(levelStart, pp.maskBigint, shareDecrypt)
	rckks.NTTRCKKSLvl(ringQ, levelStart, shareDecrypt, shareDecrypt)
	// h0 = sk*c1 + mask
	ringQ.MulCoeffsMontgomeryAndAddLvl(levelStart, sk, ciphertext.Value()[1], shareDecrypt)
	// h0 = sk*c1 + mask + e0
	pp.gaussianSampler.Read(pp.tmp)
	rckks.NTTRCKKS(ringQ, pp.tmp, pp.tmp)
	ringQ.AddLvl(levelStart, shareDecrypt, pp.tmp, shareDecrypt)

	// h1 = pi(mask) (at level max)
	pp.encoder.FFT(pp.maskComplex, slots)            // Switches to the coeff domain
	pp.permuteWithIndex(permutation, pp.maskComplex) // Applies the permutation
	pp.encoder.InvFFT(pp.maskComplex, slots)         // Switches back to the

	for i := uint64(0); i < maxSlots; i++ {
		pp.maskComplex[i].Real().Int(pp.maskBigint[i])
	}

	ringQ.SetCoefficientsBigint(pp.maskBigint, shareRecrypt)

	rckks.NTTRCKKS(ringQ, shareRecrypt, shareRecrypt)

	// h1 = sk*a + mask
	ringQ.MulCoeffsMontgomeryAndAdd(sk, crs, shareRecrypt)

	// h1 = sk*a + mask + e1
	pp.gaussianSampler.Read(pp.tmp)
	rckks.NTTRCKKS(ringQ, pp.tmp, pp.tmp)
	ringQ.Add(shareRecrypt, pp.tmp, shareRecrypt)

	// h1 = -sk*c1 - mask - e1
	ringQ.Neg(shareRecrypt, shareRecrypt)

	pp.tmp.Zero()
}

// Aggregate adds share1 with share2 on shareOut.
func (pp *PermuteProtocol) Aggregate(share1, share2, shareOut *ring.Poly) {
	pp.drckksContext.ringQ.AddLvl(uint64(len(share1.Coeffs)-1), share1, share2, shareOut)
}

// Decrypt operates a masked decryption on the ciphertext with the given decryption share.
func (pp *PermuteProtocol) Decrypt(ciphertext *rckks.Ciphertext, shareDecrypt RefreshShareDecrypt) {
	pp.drckksContext.ringQ.AddLvl(ciphertext.Level(), ciphertext.Value()[0], shareDecrypt, ciphertext.Value()[0])
}

// Permute takes a masked decrypted ciphertext at modulus Q_0 and returns the same masked decrypted ciphertext at modulus Q_L, with Q_0 << Q_L.
// Operates a permutation of the plaintext slots.
func (pp *PermuteProtocol) Permute(ciphertext *rckks.Ciphertext, permutation []uint64, slots uint64) {
	drckksContext := pp.drckksContext
	ringQ := pp.drckksContext.ringQ

	rckks.InvNTTRCKKSLvl(ringQ, ciphertext.Level(), ciphertext.Value()[0], ciphertext.Value()[0])

	ringQ.PolyToBigint(ciphertext.Value()[0], pp.maskBigint)

	QStart := ring.NewUint(ringQ.Modulus[0])
	for i := uint64(1); i < ciphertext.Level()+1; i++ {
		QStart.Mul(QStart, ring.NewUint(ringQ.Modulus[i]))
	}
	QHalf := new(big.Int).Rsh(QStart, 1)

	maxSlots := pp.drckksContext.n
	//gap := maxSlots / slots

	var sign int
	for i := uint64(0); i < maxSlots; i++ {

		// Centers the value around the current modulus
		sign = pp.maskBigint[i].Cmp(QHalf)
		if sign == 1 || sign == 0 {
			pp.maskBigint[i].Sub(pp.maskBigint[i], QStart)
		}

		if i == 0 {
			pp.maskComplex[i][0].SetInt(pp.maskBigint[i])
			pp.maskComplex[i][1].SetFloat64(0.0)
		} else {
			pp.maskComplex[i][0].SetInt(pp.maskBigint[i])
			pp.maskComplex[maxSlots-i][1].SetInt(pp.maskBigint[i])
			pp.maskComplex[maxSlots-i][1].Neg(pp.maskComplex[maxSlots-i][1])
		}
	}

	pp.encoder.FFT(pp.maskComplex, slots)

	pp.permuteWithIndex(permutation, pp.maskComplex)

	pp.encoder.InvFFT(pp.maskComplex, slots)

	for i := uint64(0); i < maxSlots; i++ {
		pp.maskComplex[i].Real().Int(pp.maskBigint[i])
	}

	for ciphertext.Level() != drckksContext.params.MaxLevel() {
		ciphertext.Value()[0].Coeffs = append(ciphertext.Value()[0].Coeffs, make([][]uint64, 1)...)
		ciphertext.Value()[0].Coeffs[ciphertext.Level()] = make([]uint64, drckksContext.n)
	}

	ringQ.SetCoefficientsBigintLvl(ciphertext.Level(), pp.maskBigint, ciphertext.Value()[0])

	rckks.NTTRCKKSLvl(ringQ, ciphertext.Level(), ciphertext.Value()[0], ciphertext.Value()[0])
}

// Recrypt operates a masked recryption on the masked decrypted ciphertext.
func (pp *PermuteProtocol) Recrypt(ciphertext *rckks.Ciphertext, crs *ring.Poly, shareRecrypt RefreshShareRecrypt) {

	pp.drckksContext.ringQ.Add(ciphertext.Value()[0], shareRecrypt, ciphertext.Value()[0])

	ciphertext.Value()[1] = crs.CopyNew()
}
