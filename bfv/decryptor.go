package bfv

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

// NewDecryptor instantiates a Decryptor for the BFV scheme.
func NewDecryptor(params Parameters, sk *rlwe.SecretKey) Decryptor {
	return &decryptor{rlwe.NewDecryptor(params.Parameters, sk), params}
}

// Decrypt decrypts the ciphertext and write the result in ptOut.
func (dec *decryptor) Decrypt(ct *Ciphertext, ptOut *Plaintext) {
	dec.Decryptor.Decrypt(&rlwe.Ciphertext{Value: ct.Value}, &rlwe.Plaintext{Value: ptOut.Value})
}

// DecryptNew decrypts the ciphertext and returns the result in a newly allocated Plaintext.
func (dec *decryptor) DecryptNew(ct *Ciphertext) (ptOut *Plaintext) {
	pt := NewPlaintext(dec.params)
	dec.Decryptor.Decrypt(ct.Ciphertext, pt.Plaintext)
	return pt
}
