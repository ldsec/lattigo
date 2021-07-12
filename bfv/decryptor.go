package bfv

import (
	"github.com/ldsec/lattigo/v2/rlwe"
)

// Decryptor is a struct wrapping a rlwe.Decryptor.
type Decryptor struct {
	rlwe.Decryptor
}

// NewDecryptor instantiates a new rlwe.Decryptor wrapped for BFV ciphertexts and plaintexts.
func NewDecryptor(params Parameters, sk *rlwe.SecretKey) *Decryptor {
	return &Decryptor{Decryptor: rlwe.NewDecryptor(params.Parameters, sk)}
}

// DecryptNew calls rlwe.Decryptor.DecryptNew.
func (decryptor *Decryptor) DecryptNew(ciphertext *Ciphertext) (plaintext *Plaintext) {
	return &Plaintext{Plaintext: decryptor.Decryptor.DecryptNew(&rlwe.Element{Value: ciphertext.Value})}
}

// Decrypt calls rlwe.Decryptor.Decrypt.
func (decryptor *Decryptor) Decrypt(ciphertext *Ciphertext, plaintext *Plaintext) {
	decryptor.Decryptor.Decrypt(&rlwe.Element{Value: ciphertext.Value}, &rlwe.Plaintext{Value: plaintext.Value})
}
