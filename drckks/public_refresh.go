package drckks

import (
	"math/big"

	"github.com/ldsec/lattigo/v2/rckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// RefreshProtocol is a struct storing the parameters for the Refresh protocol.
type RefreshProtocol struct {
	drckksContext   *drckksContext
	tmp             *ring.Poly
	maskBigint      []*big.Int
	gaussianSampler *ring.GaussianSampler
}

// RefreshShareDecrypt is a struct storing the masked decryption share.
type RefreshShareDecrypt *ring.Poly

// RefreshShareRecrypt is a struct storing the masked recryption share.
type RefreshShareRecrypt *ring.Poly

// NewRefreshProtocol creates a new instance of the Refresh protocol.
func NewRefreshProtocol(params *rckks.Parameters) (refreshProtocol *RefreshProtocol) {

	refreshProtocol = new(RefreshProtocol)
	drckksContext := newDrckksContext(params)
	refreshProtocol.drckksContext = drckksContext
	refreshProtocol.tmp = drckksContext.ringQ.NewPoly()
	refreshProtocol.maskBigint = make([]*big.Int, drckksContext.n)
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	refreshProtocol.gaussianSampler = ring.NewGaussianSampler(prng, drckksContext.ringQ, params.Sigma(), uint64(6*params.Sigma()))

	return
}

// AllocateShares allocates the shares of the Refresh protocol.
func (refreshProtocol *RefreshProtocol) AllocateShares(levelStart uint64) (RefreshShareDecrypt, RefreshShareRecrypt) {
	return refreshProtocol.drckksContext.ringQ.NewPolyLvl(levelStart), refreshProtocol.drckksContext.ringQ.NewPoly()
}

// GenShares generates the decryption and recryption shares of the Refresh protocol.
func (refreshProtocol *RefreshProtocol) GenShares(sk *ring.Poly, levelStart, nParties uint64, ciphertext *rckks.Ciphertext, crs *ring.Poly, shareDecrypt RefreshShareDecrypt, shareRecrypt RefreshShareRecrypt) {

	ringQ := refreshProtocol.drckksContext.ringQ

	bound := ring.NewUint(ringQ.Modulus[0])
	for i := uint64(1); i < levelStart+1; i++ {
		bound.Mul(bound, ring.NewUint(ringQ.Modulus[i]))
	}

	bound.Quo(bound, ring.NewUint(2*nParties))
	boundHalf := new(big.Int).Rsh(bound, 1)

	var sign int
	for i := range refreshProtocol.maskBigint {
		refreshProtocol.maskBigint[i] = ring.RandInt(bound)
		sign = refreshProtocol.maskBigint[i].Cmp(boundHalf)
		if sign == 1 || sign == 0 {
			refreshProtocol.maskBigint[i].Sub(refreshProtocol.maskBigint[i], bound)
		}
	}

	// h0 = mask (at level min)
	ringQ.SetCoefficientsBigintLvl(levelStart, refreshProtocol.maskBigint, shareDecrypt)
	// h1 = mask (at level max)
	ringQ.SetCoefficientsBigint(refreshProtocol.maskBigint, shareRecrypt)

	for i := range refreshProtocol.maskBigint {
		refreshProtocol.maskBigint[i].SetUint64(0)
	}

	rckks.NTTRCKKSLvl(ringQ, levelStart, shareDecrypt, shareDecrypt)
	rckks.NTTRCKKS(ringQ, shareRecrypt, shareRecrypt)

	// h0 = sk*c1 + mask
	ringQ.MulCoeffsMontgomeryAndAddLvl(levelStart, sk, ciphertext.Value()[1], shareDecrypt)

	// h1 = sk*a + mask
	ringQ.MulCoeffsMontgomeryAndAdd(sk, crs, shareRecrypt)

	// h0 = sk*c1 + mask + e0
	refreshProtocol.gaussianSampler.ReadLvl(levelStart, refreshProtocol.tmp)
	rckks.NTTRCKKSLvl(ringQ, levelStart, refreshProtocol.tmp, refreshProtocol.tmp)
	ringQ.AddLvl(levelStart, shareDecrypt, refreshProtocol.tmp, shareDecrypt)

	// h1 = sk*a + mask + e1
	refreshProtocol.gaussianSampler.Read(refreshProtocol.tmp)
	rckks.NTTRCKKS(ringQ, refreshProtocol.tmp, refreshProtocol.tmp)
	ringQ.Add(shareRecrypt, refreshProtocol.tmp, shareRecrypt)

	// h1 = -sk*c1 - mask - e0
	ringQ.Neg(shareRecrypt, shareRecrypt)

	refreshProtocol.tmp.Zero()
}

// Aggregate adds share1 with share2 on shareOut.
func (refreshProtocol *RefreshProtocol) Aggregate(share1, share2, shareOut *ring.Poly) {
	refreshProtocol.drckksContext.ringQ.AddLvl(uint64(len(share1.Coeffs)-1), share1, share2, shareOut)
}

// Decrypt operates a masked decryption on the ciphertext with the given decryption share.
func (refreshProtocol *RefreshProtocol) Decrypt(ciphertext *rckks.Ciphertext, shareDecrypt RefreshShareDecrypt) {
	refreshProtocol.drckksContext.ringQ.AddLvl(ciphertext.Level(), ciphertext.Value()[0], shareDecrypt, ciphertext.Value()[0])
}

// Recode takes a masked decrypted ciphertext at modulus Q_0 and returns the same masked decrypted ciphertext at modulus Q_L, with Q_0 << Q_L.
func (refreshProtocol *RefreshProtocol) Recode(ciphertext *rckks.Ciphertext) {
	drckksContext := refreshProtocol.drckksContext
	ringQ := refreshProtocol.drckksContext.ringQ

	rckks.InvNTTRCKKSLvl(ringQ, ciphertext.Level(), ciphertext.Value()[0], ciphertext.Value()[0])

	ringQ.PolyToBigint(ciphertext.Value()[0], refreshProtocol.maskBigint)

	QStart := ring.NewUint(ringQ.Modulus[0])
	for i := uint64(1); i < ciphertext.Level()+1; i++ {
		QStart.Mul(QStart, ring.NewUint(ringQ.Modulus[i]))
	}

	QHalf := new(big.Int).Rsh(QStart, 1)

	for ciphertext.Level() != drckksContext.params.MaxLevel() {
		ciphertext.Value()[0].Coeffs = append(ciphertext.Value()[0].Coeffs, make([][]uint64, 1)...)
		ciphertext.Value()[0].Coeffs[ciphertext.Level()] = make([]uint64, drckksContext.n)
	}

	var sign int
	for i := uint64(0); i < drckksContext.n; i++ {
		sign = refreshProtocol.maskBigint[i].Cmp(QHalf)
		if sign == 1 || sign == 0 {
			refreshProtocol.maskBigint[i].Sub(refreshProtocol.maskBigint[i], QStart)
		}
	}

	ringQ.SetCoefficientsBigintLvl(ciphertext.Level(), refreshProtocol.maskBigint, ciphertext.Value()[0])

	rckks.NTTRCKKSLvl(ringQ, ciphertext.Level(), ciphertext.Value()[0], ciphertext.Value()[0])
}

// Recrypt operates a masked recryption on the masked decrypted ciphertext.
func (refreshProtocol *RefreshProtocol) Recrypt(ciphertext *rckks.Ciphertext, crs *ring.Poly, shareRecrypt RefreshShareRecrypt) {

	refreshProtocol.drckksContext.ringQ.Add(ciphertext.Value()[0], shareRecrypt, ciphertext.Value()[0])
	crs.Coeffs = crs.Coeffs[:ciphertext.Level()+1]
	ciphertext.Value()[1] = crs.CopyNew()
}
