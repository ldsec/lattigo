package dckks

import (
	"math/big"

	"encoding/binary"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// MaskedTransformProtocol is a struct storing the parameters for the MaskedTransformProtocol protocol.
type MaskedTransformProtocol struct {
	e2s E2SProtocol
	s2e S2EProtocol

	defaultScale *big.Int

	ringQ         *ring.Ring
	tmpMask       []*big.Int
	tmpBigComplex []*ring.Complex
	encoder       ckks.EncoderBigComplex
}

// MaskedTransformFunc is a method template for linear transforms that can be
// evaluated on a ciphertext during its collective refresh
type MaskedTransformFunc func(coeffsIn, coeffsOut []*ring.Complex)

// MaskedTransformShare is a struct storing the decryption and recryption shares.
type MaskedTransformShare struct {
	e2sShare drlwe.CKSShare
	s2eShare drlwe.CKSShare
}

// MarshalBinary encodes a RefreshShare on a slice of bytes.
func (share *MaskedTransformShare) MarshalBinary() (data []byte, err error) {
	var e2sData, s2eData []byte
	if e2sData, err = share.e2sShare.MarshalBinary(); err != nil {
		return nil, err
	}
	if s2eData, err = share.s2eShare.MarshalBinary(); err != nil {
		return nil, err
	}
	data = make([]byte, 8)
	binary.LittleEndian.PutUint64(data, uint64(len(e2sData)))
	data = append(data, e2sData...)
	data = append(data, s2eData...)
	return data, nil
}

// UnmarshalBinary decodes a marshaled RefreshShare on the target RefreshShare.
func (share *MaskedTransformShare) UnmarshalBinary(data []byte) error {

	e2sDataLen := binary.LittleEndian.Uint64(data[:8])

	if err := share.e2sShare.UnmarshalBinary(data[8 : e2sDataLen+8]); err != nil {
		return err
	}
	if err := share.s2eShare.UnmarshalBinary(data[8+e2sDataLen:]); err != nil {
		return err
	}
	return nil
}

// NewMaskedTransformProtocol creates a new instance of the PermuteProtocol.
func NewMaskedTransformProtocol(params ckks.Parameters, precision int, sigmaSmudging float64) (rfp *MaskedTransformProtocol) {

	rfp = new(MaskedTransformProtocol)
	rfp.e2s = *NewE2SProtocol(params, sigmaSmudging)
	rfp.s2e = *NewS2EProtocol(params, sigmaSmudging)

	rfp.defaultScale = new(big.Int)
	ring.NewFloat(params.Scale(), precision).Int(rfp.defaultScale)

	rfp.ringQ = rfp.e2s.ringQ

	rfp.tmpMask = make([]*big.Int, rfp.ringQ.N)
	rfp.tmpBigComplex = make([]*ring.Complex, rfp.ringQ.N>>1)
	for i := range rfp.tmpBigComplex {
		rfp.tmpMask[i*2] = new(big.Int)
		rfp.tmpMask[i*2+1] = new(big.Int)
		rfp.tmpBigComplex[i] = ring.NewComplex(ring.NewFloat(0, precision), ring.NewFloat(0, precision))
	}
	rfp.encoder = ckks.NewEncoderBigComplex(params, precision)
	return
}

// AllocateShare allocates the shares of the PermuteProtocol
func (rfp *MaskedTransformProtocol) AllocateShare(levelDecrypt, levelRecrypt int) *MaskedTransformShare {
	return &MaskedTransformShare{*rfp.e2s.AllocateShare(levelDecrypt), *rfp.s2e.AllocateShare(levelRecrypt)}
}

// GenShares generates the shares of the PermuteProtocol
// This protocol requires additional inputs which are :
// logBound : the bit length of the masks
// logSlots : the bit length of the number of slots
//
// The method "GetMinimumLevelForBootstrapping" should be used to get the minimum level at which the masked transform can be called while still ensure 128-bits of security, as well as the
// value for logBound.
func (rfp *MaskedTransformProtocol) GenShares(sk *rlwe.SecretKey, logBound, logSlots int, ct *ckks.Ciphertext, crs *ring.Poly, transform MaskedTransformFunc, shareOut *MaskedTransformShare) {

	if ct.Level() != shareOut.e2sShare.Value.Level() {
		panic("ciphertext level must be equal to e2sShare")
	}

	if crs.Level() != shareOut.s2eShare.Value.Level() {
		panic("crs level must be equal to s2eShare")
	}

	// Generates the decryption share
	// Returns [M_i] on rfp.tmpMask and [a*s_i -M_i + e] on e2sShare
	rfp.e2s.GenShare(sk, logBound, logSlots, ct, &rlwe.AdditiveShareBigint{Value: rfp.tmpMask}, &shareOut.e2sShare)

	slots := 1 << logSlots
	gap := rfp.ringQ.N / (2 * slots)

	// Applies LT(M_i)
	if transform != nil {
		// Extracts sparse coefficients
		for i, jdx, idx := 0, rfp.ringQ.N>>1, 0; i < slots; i, jdx, idx = i+1, jdx+gap, idx+gap {
			rfp.tmpBigComplex[idx][0].SetInt(rfp.tmpMask[idx])
			rfp.tmpBigComplex[idx][1].SetInt(rfp.tmpMask[jdx])
		}

		// Decodes
		rfp.encoder.FFT(rfp.tmpBigComplex, 1<<logSlots)

		// Applies the linear transform
		transform(rfp.tmpBigComplex, rfp.tmpBigComplex)

		// Recodes
		rfp.encoder.InvFFT(rfp.tmpBigComplex, 1<<logSlots)

		// Puts the coefficient back
		for i, jdx, idx := 0, rfp.ringQ.N>>1, 0; i < slots; i, jdx, idx = i+1, jdx+gap, idx+gap {
			rfp.tmpBigComplex[i].Real().Int(rfp.tmpMask[idx])
			rfp.tmpBigComplex[i].Imag().Int(rfp.tmpMask[jdx])
		}
	}

	// Applies LT(M_i) * diffscale
	inputScaleInt := new(big.Int)
	ring.NewFloat(ct.Scale, 256).Int(inputScaleInt)

	// Scales the mask by the ratio between the two scales
	for i := 0; i < rfp.ringQ.N; i++ {
		if i%gap == 0 {
			rfp.tmpMask[i].Mul(rfp.tmpMask[i], rfp.defaultScale)
			rfp.tmpMask[i].Quo(rfp.tmpMask[i], inputScaleInt)
		}
	}

	// Returns [-a*s_i + LT(M_i) * diffscale + e] on s2eShare
	rfp.s2e.GenShare(sk, crs, &rlwe.AdditiveShareBigint{Value: rfp.tmpMask}, &shareOut.s2eShare)
}

// Aggregate sums share1 and share2 on shareOut.
func (rfp *MaskedTransformProtocol) Aggregate(share1, share2, shareOut *MaskedTransformShare) {

	if share1.e2sShare.Value.Level() != share2.e2sShare.Value.Level() || share1.e2sShare.Value.Level() != shareOut.e2sShare.Value.Level() {
		panic("all e2s shares must be at the same level")
	}

	if share1.s2eShare.Value.Level() != share2.s2eShare.Value.Level() || share1.s2eShare.Value.Level() != shareOut.s2eShare.Value.Level() {
		panic("all s2e shares must be at the same level")
	}

	rfp.ringQ.AddLvl(share1.e2sShare.Value.Level(), share1.e2sShare.Value, share2.e2sShare.Value, shareOut.e2sShare.Value)
	rfp.ringQ.AddLvl(share1.s2eShare.Value.Level(), share1.s2eShare.Value, share2.s2eShare.Value, shareOut.s2eShare.Value)
}

// Transform applies Decrypt, Recode and Recrypt on the input ciphertext.
func (rfp *MaskedTransformProtocol) Transform(ct *ckks.Ciphertext, logSlots int, transform MaskedTransformFunc, crs *ring.Poly, share *MaskedTransformShare, ciphertextOut *ckks.Ciphertext) {

	if ct.Level() != share.e2sShare.Value.Level() {
		panic("ciphertext level and e2s level must be the same")
	}

	if crs.Level() != share.s2eShare.Value.Level() {
		panic("crs level and s2e level must be the same")
	}

	maxLevel := crs.Level()

	// Returns -sum(M_i) + x (outside of the NTT domain)
	rfp.e2s.GetShare(nil, &share.e2sShare, ct, &rlwe.AdditiveShareBigint{Value: rfp.tmpMask})

	slots := 1 << logSlots
	gap := rfp.ringQ.N / (2 * slots)

	// Returns LT(-sum(M_i) + x)
	if transform != nil {
		// Extracts sparse coefficients
		for i, jdx, idx := 0, rfp.ringQ.N>>1, 0; i < slots; i, jdx, idx = i+1, jdx+gap, idx+gap {
			rfp.tmpBigComplex[idx][0].SetInt(rfp.tmpMask[idx])
			rfp.tmpBigComplex[idx][1].SetInt(rfp.tmpMask[jdx])
		}

		// Decodes
		rfp.encoder.FFT(rfp.tmpBigComplex, 1<<logSlots)

		// Applies the linear transform
		transform(rfp.tmpBigComplex, rfp.tmpBigComplex)

		// Recodes
		rfp.encoder.InvFFT(rfp.tmpBigComplex, 1<<logSlots)

		// Puts the coefficient back
		for i, jdx, idx := 0, rfp.ringQ.N>>1, 0; i < slots; i, jdx, idx = i+1, jdx+gap, idx+gap {
			rfp.tmpBigComplex[i].Real().Int(rfp.tmpMask[idx])
			rfp.tmpBigComplex[i].Imag().Int(rfp.tmpMask[jdx])
		}
	}

	// Returns LT(-sum(M_i) + x) * diffscale
	inputScaleInt := new(big.Int)
	ring.NewFloat(ct.Scale, 256).Int(inputScaleInt)

	// Scales the mask by the ratio between the two scales
	for i := range rfp.tmpMask {
		rfp.tmpMask[i].Mul(rfp.tmpMask[i], rfp.defaultScale)
		rfp.tmpMask[i].Quo(rfp.tmpMask[i], inputScaleInt)
	}

	// Extend the levels of the ciphertext for future allocation
	for ciphertextOut.Level() != crs.Level() {
		level := ciphertextOut.Level() + 1

		ciphertextOut.Value[0].Coeffs = append(ciphertextOut.Value[0].Coeffs, make([][]uint64, 1)...)
		ciphertextOut.Value[0].Coeffs[level] = make([]uint64, rfp.ringQ.N)

		ciphertextOut.Value[1].Coeffs = append(ciphertextOut.Value[1].Coeffs, make([][]uint64, 1)...)
		ciphertextOut.Value[1].Coeffs[level] = make([]uint64, rfp.ringQ.N)
	}

	// Sets LT(-sum(M_i) + x) * diffscale in the RNS domain
	rfp.ringQ.SetCoefficientsBigintLvl(maxLevel, rfp.tmpMask, ciphertextOut.Value[0])

	// Sets LT(-sum(M_i) + x) * diffscale in the NTT domain
	rfp.ringQ.NTTLvl(maxLevel, ciphertextOut.Value[0], ciphertextOut.Value[0])

	// LT(-sum(M_i) + x) * diffscale + [-a*s + LT(M_i) * diffscale + e] = [-a*s + LT(x) * diffscale + e]
	rfp.ringQ.AddLvl(maxLevel, ciphertextOut.Value[0], share.s2eShare.Value, ciphertextOut.Value[0])

	// Copies the result on the out ciphertext
	rfp.s2e.GetEncryption(&drlwe.CKSShare{Value: ciphertextOut.Value[0]}, crs, ciphertextOut)
}
