package dbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// MaskedTransformProtocol is a struct storing the parameters for the MaskedTransformProtocol protocol.
type MaskedTransformProtocol struct {
	e2s E2SProtocol
	s2e S2EProtocol

	encoder bfv.Encoder
	ringQ   *ring.Ring
	ringT   *ring.Ring

	tmpPt bfv.Plaintext
	//	tmpPtRt     bfv.PlaintextRingT
	tmpMask     *ring.Poly
	tmpMaskPerm *ring.Poly
}

// MaskedTransformFunc is type of functions that can be operated on masked BFV plaintexts as a part of the
// Masked Transform Protocol.
type MaskedTransformFunc func(ptIn, ptOut bfv.PlaintextRingT)

// MaskedTransformShare is a struct storing the decryption and recryption shares.
type MaskedTransformShare struct {
	e2sShare drlwe.CKSShare
	s2eShare drlwe.CKSShare
}

// MarshalBinary encodes a RefreshShare on a slice of bytes.
func (share *MaskedTransformShare) MarshalBinary() ([]byte, error) {
	e2sData, err := share.e2sShare.MarshalBinary()
	if err != nil {
		return nil, err
	}
	s2eData, err := share.s2eShare.MarshalBinary()
	if err != nil {
		return nil, err
	}
	return append(e2sData, s2eData...), nil
}

// UnmarshalBinary decodes a marshaled RefreshShare on the target RefreshShare.
func (share *MaskedTransformShare) UnmarshalBinary(data []byte) error {
	shareLen := len(data) >> 1
	share.e2sShare.UnmarshalBinary(data[:shareLen])
	share.s2eShare.UnmarshalBinary(data[shareLen:])
	return nil
}

// NewMaskedTransformProtocol creates a new instance of the PermuteProtocol.
func NewMaskedTransformProtocol(params bfv.Parameters, sigmaSmudging float64) (rfp *MaskedTransformProtocol) {

	rfp = new(MaskedTransformProtocol)
	rfp.e2s = *NewE2SProtocol(params, sigmaSmudging)
	rfp.s2e = *NewS2EProtocol(params, sigmaSmudging)

	rfp.encoder = rfp.e2s.encoder

	rfp.ringQ = rfp.e2s.ringQ
	rfp.ringT = rfp.e2s.ringT
	rfp.tmpPt = *bfv.NewPlaintext(params)
	//rfp.tmpPtRt = *bfv.NewPlaintextRingT(params)
	rfp.tmpMask = rfp.ringT.NewPoly()
	rfp.tmpMaskPerm = rfp.ringT.NewPoly()
	return
}

// AllocateShare allocates the shares of the PermuteProtocol
func (rfp *MaskedTransformProtocol) AllocateShare() *MaskedTransformShare {
	level := len(rfp.ringQ.Modulus) - 1
	return &MaskedTransformShare{*rfp.e2s.AllocateShare(level), *rfp.s2e.AllocateShare(level)}
}

// GenShares generates the shares of the PermuteProtocol
func (rfp *MaskedTransformProtocol) GenShares(sk *rlwe.SecretKey, ciphertext *bfv.Ciphertext, crs *ring.Poly, transform MaskedTransformFunc, shareOut *MaskedTransformShare) {
	rfp.e2s.GenShare(sk, ciphertext, &rlwe.AdditiveShare{Value: *rfp.tmpMask}, &shareOut.e2sShare)
	mask := rfp.tmpMask
	if transform != nil {
		transform(bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: mask}}, bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: rfp.tmpMaskPerm}})
		mask = rfp.tmpMaskPerm
	}
	rfp.s2e.GenShare(sk, crs, &rlwe.AdditiveShare{Value: *mask}, &shareOut.s2eShare)
}

// Aggregate sums share1 and share2 on shareOut.
func (rfp *MaskedTransformProtocol) Aggregate(share1, share2, shareOut *MaskedTransformShare) {
	rfp.ringQ.Add(share1.e2sShare.Value, share2.e2sShare.Value, shareOut.e2sShare.Value)
	rfp.ringQ.Add(share1.s2eShare.Value, share2.s2eShare.Value, shareOut.s2eShare.Value)
}

// Transform applies Decrypt, Recode and Recrypt on the input ciphertext.
func (rfp *MaskedTransformProtocol) Transform(ciphertext *bfv.Ciphertext, transform MaskedTransformFunc, crs *ring.Poly, share *MaskedTransformShare, ciphertextOut *bfv.Ciphertext) {
	rfp.e2s.GetShare(nil, &share.e2sShare, ciphertext, &rlwe.AdditiveShare{Value: *rfp.tmpMask}) // tmpMask RingT(m - sum M_i)
	mask := rfp.tmpMask
	if transform != nil {
		transform(bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: mask}}, bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: rfp.tmpMaskPerm}})
		mask = rfp.tmpMaskPerm
	}
	//rfp.tmpPtRt.Value.Copy(mask)
	rfp.encoder.ScaleUp(&bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: mask}}, &rfp.tmpPt)
	//rfp.encoder.EncodeUint(mask.Coeffs[0], &rfp.tmpPt) // tmpMask RingQ( Delta*(m - sum M_i))
	rfp.ringQ.Add(rfp.tmpPt.Value, share.s2eShare.Value, ciphertextOut.Value[0])
	rfp.s2e.GetEncryption(&drlwe.CKSShare{Value: ciphertextOut.Value[0]}, crs, ciphertextOut)
}
