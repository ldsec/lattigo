package ckks

import (
	"github.com/ldsec/lattigo/v2/rlwe"
)

// Decryptor is an interface wrapping a rlwe.Decryptor.
type Decryptor interface {
	DecryptNew(ciphertext *Ciphertext) (plaintext *Plaintext)
	Decrypt(ciphertext *Ciphertext, plaintext *Plaintext)
}
type decryptor struct {
	rlwe.Decryptor
	params Parameters
}

// NewDecryptor instantiates a Decryptor for the CKKS scheme.
func NewDecryptor(params Parameters, sk *rlwe.SecretKey) Decryptor {
	return &decryptor{rlwe.NewDecryptor(params.Parameters, sk), params}
}

// Decrypt decrypts the ciphertext and write the result in ptOut.
func (dec *decryptor) DecryptNew(ciphertext *Ciphertext) (plaintext *Plaintext) {
	pt := NewPlaintext(dec.params, ciphertext.Level(), ciphertext.Scale)
	dec.Decryptor.Decrypt(ciphertext.Ciphertext, pt.Plaintext)
	return pt
}

// DecryptNew decrypts the ciphertext and returns the result in a newly allocated Plaintext.
func (dec *decryptor) Decrypt(ciphertext *Ciphertext, plaintext *Plaintext) {
	dec.Decryptor.Decrypt(&rlwe.Ciphertext{Value: ciphertext.Value}, &rlwe.Plaintext{Value: plaintext.Value})
	plaintext.Scale = ciphertext.Scale
}
