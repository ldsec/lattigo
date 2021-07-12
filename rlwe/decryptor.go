package rlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// Decryptor is an interface for decrypting Ciphertexts. A Decryptor stores the secret-key.
type Decryptor interface {
	// DecryptNew decrypts the ciphertext and returns a newly created
	// plaintext. A Horner method is used for evaluating the decryption.
	// The level of the output plaintext is ciphertext.Level().
	DecryptNew(ciphertext *Element) (plaintext *Plaintext)

	// DecryptNew decrypts the ciphertext and returns a newly created
	// plaintext. A Horner method is used for evaluating the decryption.
	// The level of the output plaintext is ciphertext.Level().
	// Plaintext output is in the NTT domain.
	DecryptNTTNew(ciphertext *Element) (plaintext *Plaintext)

	// Decrypt decrypts the ciphertext and returns the result on the provided
	// receiver plaintext. A Horner method is used for evaluating the
	// decryption.
	// The level of the output plaintext is min(ciphertext.Level(), plaintext.Level())
	// Output domain will match plaintext.Value.IsNTT value.
	Decrypt(ciphertext *Element, plaintext *Plaintext)
}

// decryptor is a structure used to decrypt ciphertext. It stores the secret-key.
type decryptor struct {
	params Parameters
	ringQ  *ring.Ring
	pool   *ring.Poly
	sk     *SecretKey
}

// NewDecryptor instantiates a new Decryptor that will be able to decrypt ciphertexts
// encrypted under the provided secret-key.
func NewDecryptor(params Parameters, sk *SecretKey) Decryptor {

	if sk.Value[0].Degree() != params.N() {
		panic("secret_key is invalid for the provided parameters")
	}

	return &decryptor{
		params: params,
		ringQ:  params.RingQ(),
		pool:   params.RingQ().NewPoly(),
		sk:     sk,
	}
}

func (decryptor *decryptor) DecryptNew(ciphertext *Element) (plaintext *Plaintext) {

	plaintext = NewPlaintext(decryptor.params, ciphertext.Level())

	decryptor.Decrypt(ciphertext, plaintext)

	return plaintext
}

func (decryptor *decryptor) DecryptNTTNew(ciphertext *Element) (plaintext *Plaintext) {

	plaintext = NewPlaintext(decryptor.params, ciphertext.Level())
	plaintext.Value.IsNTT = true

	decryptor.Decrypt(ciphertext, plaintext)

	return plaintext
}

func (decryptor *decryptor) Decrypt(ciphertext *Element, plaintext *Plaintext) {

	ringQ := decryptor.ringQ

	level := utils.MinInt(ciphertext.Level(), plaintext.Level())

	plaintext.Value.Coeffs = plaintext.Value.Coeffs[:level+1]

	if ciphertext.Value[0].IsNTT {
		ring.CopyValuesLvl(level, ciphertext.Value[ciphertext.Degree()], plaintext.Value)
	} else {
		ringQ.NTTLazyLvl(level, ciphertext.Value[ciphertext.Degree()], plaintext.Value)
	}

	for i := ciphertext.Degree(); i > 0; i-- {

		ringQ.MulCoeffsMontgomeryLvl(level, plaintext.Value, decryptor.sk.Value[0], plaintext.Value)

		if !ciphertext.Value[0].IsNTT {
			ringQ.NTTLazyLvl(level, ciphertext.Value[i-1], decryptor.pool)
			ringQ.AddLvl(level, plaintext.Value, decryptor.pool, plaintext.Value)
		} else {
			ringQ.AddLvl(level, plaintext.Value, ciphertext.Value[i-1], plaintext.Value)
		}

		if i&7 == 7 {
			ringQ.ReduceLvl(level, plaintext.Value, plaintext.Value)
		}
	}

	if (ciphertext.Degree())&7 != 7 {
		ringQ.ReduceLvl(level, plaintext.Value, plaintext.Value)
	}

	if !plaintext.Value.IsNTT {
		ringQ.InvNTTLvl(level, plaintext.Value, plaintext.Value)
	}
}
