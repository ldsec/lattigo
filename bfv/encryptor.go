package bfv

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// Encryptor an encryption interface for the BFV scheme.
type Encryptor interface {
	Encrypt(plaintext *Plaintext, ciphertext *Ciphertext)
	EncryptNew(plaintext *Plaintext) *Ciphertext
	EncryptFromCRP(plaintext *Plaintext, crp *ring.Poly, ctOut *Ciphertext)
	EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext
}

type encryptor struct {
	rlwe.Encryptor
	params Parameters
}

// NewEncryptor instantiates a new Encryptor for the BFV scheme. The key argument can
// be either a *rlwe.PublicKey or a *rlwe.SecretKey.
func NewEncryptor(params Parameters, key interface{}) Encryptor {
	return &encryptor{rlwe.NewEncryptor(params.Parameters, key), params}
}

// NewFastEncryptor instantiates a new Encryptor for the BFV scheme.
// This encryptor's Encrypt method first encrypts zero in Q and then adds the plaintext.
// This method is faster than the normal encryptor but result in a noisier ciphertext.
func NewFastEncryptor(params Parameters, key *rlwe.PublicKey) Encryptor {
	return &encryptor{rlwe.NewFastEncryptor(params.Parameters, key), params}
}

// Encrypt encrypts the input plaintext and write the result on ctOut.
// The encryption algorithm depends on how the receiver encryptor was initialized (see
// NewEncryptor and NewFastEncryptor).
func (enc *encryptor) Encrypt(plaintext *Plaintext, ctOut *Ciphertext) {
	enc.Encryptor.Encrypt(&rlwe.Plaintext{Value: plaintext.Value}, &rlwe.Ciphertext{Value: ctOut.Value})
}

// EncryptNew encrypts the input plaintext returns the result as a newly allocated ciphertext.
// The encryption algorithm depends on how the receiver encryptor was initialized (see
// NewEncryptor and NewFastEncryptor).
func (enc *encryptor) EncryptNew(plaintext *Plaintext) *Ciphertext {
	ct := NewCiphertext(enc.params, 1)
	enc.Encryptor.Encrypt(plaintext.Plaintext, ct.Ciphertext)
	return ct
}

// EncryptFromCRP encrypts the input plaintext and writes the result in ctOut.
// The passed crp is always treated as being in the NTT domain.
func (enc *encryptor) EncryptFromCRP(plaintext *Plaintext, crp *ring.Poly, ctOut *Ciphertext) {
	enc.Encryptor.EncryptFromCRP(&rlwe.Plaintext{Value: plaintext.Value}, crp, &rlwe.Ciphertext{Value: ctOut.Value})
}

// EncryptFromCRPNew encrypts the input plaintext and returns the result as a newly allocated ciphertext.
// The passed crp is always treated as being in the NTT domain.
func (enc *encryptor) EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext {
	ct := NewCiphertext(enc.params, 1)
	enc.Encryptor.EncryptFromCRP(&rlwe.Plaintext{Value: plaintext.Value}, crp, ct.Ciphertext)
	return ct
}
