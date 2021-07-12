package ckks

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// Encryptor a struct wrapping an rlwe.Encryptor.
type Encryptor struct {
	rlwe.Encryptor
}

// NewEncryptor instanciate a new wrapped rlwe.Encryptor for CKKS ciphertexts and plaintexts.
func NewEncryptor(params Parameters, key interface{}) *Encryptor {
	switch key := key.(type) {
	case *rlwe.PublicKey:
		return &Encryptor{rlwe.NewEncryptorFromPk(params.Parameters, key)}
	case *rlwe.SecretKey:
		return &Encryptor{rlwe.NewEncryptorFromSk(params.Parameters, key)}
	default:
		panic("key must be either rlwe.PublicKey or rlwe.SecretKey")
	}
}

// EncryptNew calls rlwe.Encryptor.EncryptNew.
func (encryptor *Encryptor) EncryptNew(plaintext *Plaintext) *Ciphertext {
	return &Ciphertext{Element: encryptor.Encryptor.EncryptNTTNew(&rlwe.Plaintext{Value: plaintext.Value}), Scale: plaintext.Scale}
}

// Encrypt calls rlwe.Encryptor.Encrypt.
func (encryptor *Encryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {
	encryptor.Encryptor.Encrypt(&rlwe.Plaintext{Value: plaintext.Value}, &rlwe.Element{Value: ciphertext.Value})
	ciphertext.Scale = plaintext.Scale
}

// EncryptFastNew calls rlwe.Encryptor.EncryptFastNew.
func (encryptor *Encryptor) EncryptFastNew(plaintext *Plaintext) *Ciphertext {
	return &Ciphertext{Element: encryptor.Encryptor.EncryptFastNTTNew(&rlwe.Plaintext{Value: plaintext.Value}), Scale: plaintext.Scale}
}

// EncryptFast calls rlwe.Encryptor.EncryptFast.
func (encryptor *Encryptor) EncryptFast(plaintext *Plaintext, ciphertext *Ciphertext) {
	encryptor.Encryptor.EncryptFast(&rlwe.Plaintext{Value: plaintext.Value}, &rlwe.Element{Value: ciphertext.Value})
	ciphertext.Scale = plaintext.Scale
}

// EncryptFromCRP calls rlwe.Encryptor.EncryptFromCRP.
func (encryptor *Encryptor) EncryptFromCRP(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly) {
	encryptor.Encryptor.EncryptFromCRP(&rlwe.Plaintext{Value: plaintext.Value}, &rlwe.Element{Value: ciphertext.Value}, crp)
	ciphertext.Scale = plaintext.Scale
}

// EncryptFromCRPNew calls rlwe.Encryptor.EncryptFromCRPNew.
func (encryptor *Encryptor) EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext {
	return &Ciphertext{Element: encryptor.Encryptor.EncryptFromCRPNTTNew(&rlwe.Plaintext{Value: plaintext.Value}, crp), Scale: plaintext.Scale}
}
