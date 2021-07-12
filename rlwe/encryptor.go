package rlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// Encryptor in an interface for encryptors
//
// encrypt with pk : ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk : ciphertext = [-a*sk + m + e, a]
type Encryptor interface {
	// EncryptNew encrypts the input plaintext using the stored key and returns
	// the result on a newly created ciphertext. The encryption is done by first
	// encrypting zero in QP, dividing by P and then adding the plaintext.
	// The level of the output ciphertext is plaintext.Level().
	EncryptNew(plaintext *Plaintext) *Element

	// EncryptNTTNew encrypts the input plaintext using the stored key and returns
	// the result on a newly created ciphertext. The encryption is done by first
	// encrypting zero in QP, dividing by P and then adding the plaintext.
	// The level of the output ciphertext is plaintext.Level().
	// Output will be in the NTT domain.
	EncryptNTTNew(plaintext *Plaintext) *Element

	// Encrypt encrypts the input plaintext using the stored key, and returns
	// the result on the receiver ciphertext. The encryption is done by first
	// encrypting zero in QP, dividing by P and then adding the plaintext.
	// The level of the output ciphertext is min(plaintext.Level(), ciphertext.Level()).
	// Output domain will match ciphertext.Value[0].IsNTT value.
	Encrypt(plaintext *Plaintext, ciphertext *Element)

	// EncryptFastNew encrypts the input plaintext using the stored key and returns
	// the result on a newly created ciphertext. The encryption is done by first
	// encrypting zero in Q and then adding the plaintext.
	// The level of the output ciphertext is plaintext.Level().
	EncryptFastNew(plaintext *Plaintext) *Element

	// EncryptFastNTTNew encrypts the input plaintext using the stored key and returns
	// the result on a newly created ciphertext. The encryption is done by first
	// encrypting zero in Q and then adding the plaintext.
	// The level of the output ciphertext is plaintext.Level().
	// Output will be in the NTT domain.
	EncryptFastNTTNew(plaintext *Plaintext) *Element

	// EncryptFast encrypts the input plaintext using the stored-key, and returns
	// the result on the receiver ciphertext. The encryption is done by first
	// encrypting zero in Q and then adding the plaintext.
	// The level of the output ciphertext is min(plaintext.Level(), ciphertext.Level()).
	// Output domain will match ciphertext.Value[0].IsNTT value.
	EncryptFast(plaintext *Plaintext, ciphertext *Element)

	// EncryptFromCRPNew encrypts the input plaintext using the stored key and returns
	// the result on a newly created ciphertext. The encryption is done by first encrypting
	// zero in QP, using the provided polynomial as the uniform polynomial, dividing by P and
	// then adding the plaintext.
	// CRP is always treated as being in the NTT domain.
	// The level of the output ciphertext is min(plaintext.Level(), len(CRP.Coeffs)-1).
	EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Element

	// EncryptFromCRPNTTNew encrypts the input plaintext using the stored key and returns
	// the result on a newly created ciphertext. The encryption is done by first encrypting
	// zero in QP, using the provided polynomial as the uniform polynomial, dividing by P and
	// then adding the plaintext.
	// The level of the output ciphertext is min(plaintext.Level(), len(CRP.Coeffs)-1).
	// CRP is always treated as being in the NTT domain.
	// Output will be in the NTT domain.
	EncryptFromCRPNTTNew(plaintext *Plaintext, crp *ring.Poly) *Element

	// EncryptFromCRP encrypts the input plaintext using the stored key and returns
	// the result tge receiver ciphertext. The encryption is done by first encrypting
	// zero in QP, using the provided polynomial as the uniform polynomial, dividing by P and
	// then adding the plaintext.
	// The level of the output ciphertext is min(plaintext.Level(), ciphertext.Level(), len(CRP.Coeffs)-1).
	// CRP is always treated as being in the NTT domain.
	// Output domain will match ciphertext.Value[0].IsNTT value.
	EncryptFromCRP(plaintext *Plaintext, ciphertext *Element, crp *ring.Poly)
}

// encryptor is a struct used to encrypt Plaintexts. It stores the public-key and/or secret-key.
type encryptor struct {
	params Parameters

	ringQ *ring.Ring
	ringP *ring.Ring

	poolQ [1]*ring.Poly
	poolP [3]*ring.Poly

	baseconverter   *ring.FastBasisExtender
	gaussianSampler *ring.GaussianSampler
	ternarySampler  *ring.TernarySampler
	uniformSampler  *ring.UniformSampler
}

type pkEncryptor struct {
	encryptor
	pk *PublicKey
}

type skEncryptor struct {
	encryptor
	sk *SecretKey
}

// NewEncryptorFromPk creates a new Encryptor with the provided public-key.
// This Encryptor can be used to encrypt Plaintexts, using the stored key.
func NewEncryptorFromPk(params Parameters, pk *PublicKey) Encryptor {
	enc := newEncryptor(params)

	if pk.Value[0][0].Degree() != params.N() || pk.Value[1][0].Degree() != params.N() {
		panic("cannot newEncryptor: pk ring degree does not match params ring degree")
	}

	return &pkEncryptor{enc, pk}
}

// NewEncryptorFromSk creates a new Encryptor with the provided secret-key.
// This Encryptor can be used to encrypt Plaintexts, using the stored key.
func NewEncryptorFromSk(params Parameters, sk *SecretKey) Encryptor {
	enc := newEncryptor(params)

	if sk.Value[0].Degree() != params.N() {
		panic("cannot newEncryptor: sk ring degree does not match params ring degree")
	}

	return &skEncryptor{enc, sk}
}

func newEncryptor(params Parameters) encryptor {

	ringQ := params.RingQ()
	ringP := params.RingP()

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	var baseconverter *ring.FastBasisExtender
	var poolP [3]*ring.Poly
	if params.PCount() != 0 {
		baseconverter = ring.NewFastBasisExtender(ringQ, ringP)
		poolP = [3]*ring.Poly{ringP.NewPoly(), ringP.NewPoly(), ringP.NewPoly()}
	}

	return encryptor{
		params:          params,
		ringQ:           ringQ,
		ringP:           ringP,
		poolQ:           [1]*ring.Poly{ringQ.NewPoly()},
		poolP:           poolP,
		baseconverter:   baseconverter,
		gaussianSampler: ring.NewGaussianSampler(prng, ringQ, params.Sigma(), int(6*params.Sigma())),
		ternarySampler:  ring.NewTernarySampler(prng, ringQ, 0.5, false),
		uniformSampler:  ring.NewUniformSampler(prng, ringQ),
	}
}

// EncryptNew encrypts the input Plaintext using the stored key and returns
// the result on a newly created Ciphertext.
//
// encrypt with pk: ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk: ciphertext = [-a*sk + m + e, a]
func (encryptor *pkEncryptor) EncryptNew(plaintext *Plaintext) *Element {

	if encryptor.baseconverter == nil {
		panic("Cannot EncryptNew : modulus P is empty -> use instead EncryptFastNew")
	}

	ciphertext := NewElement(encryptor.params, 1, plaintext.Level())
	encryptor.encrypt(plaintext, ciphertext, false)

	return ciphertext
}

// EncryptNTTNew encrypts the input Plaintext using the stored key and returns
// the result on a newly created Ciphertext.
//
// encrypt with pk: ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk: ciphertext = [-a*sk + m + e, a]
func (encryptor *pkEncryptor) EncryptNTTNew(plaintext *Plaintext) *Element {

	if encryptor.baseconverter == nil {
		panic("Cannot EncryptNew : modulus P is empty -> use instead EncryptFastNew")
	}

	ciphertext := NewElement(encryptor.params, 1, plaintext.Level())
	ciphertext.Value[0].IsNTT = true
	ciphertext.Value[1].IsNTT = true
	encryptor.encrypt(plaintext, ciphertext, false)

	return ciphertext
}

func (encryptor *pkEncryptor) Encrypt(plaintext *Plaintext, ciphertext *Element) {

	if encryptor.baseconverter == nil {
		panic("Cannot Encrypt : modulus P is empty -> use instead EncryptFast")
	}

	encryptor.encrypt(plaintext, ciphertext, false)
}

func (encryptor *pkEncryptor) EncryptFastNew(plaintext *Plaintext) *Element {
	ciphertext := NewElement(encryptor.params, 1, plaintext.Level())
	encryptor.encrypt(plaintext, ciphertext, true)

	return ciphertext
}

func (encryptor *pkEncryptor) EncryptFastNTTNew(plaintext *Plaintext) *Element {
	ciphertext := NewElement(encryptor.params, 1, plaintext.Level())
	ciphertext.Value[0].IsNTT = true
	ciphertext.Value[1].IsNTT = true
	encryptor.encrypt(plaintext, ciphertext, true)

	return ciphertext
}

func (encryptor *pkEncryptor) EncryptFast(plaintext *Plaintext, ciphertext *Element) {
	encryptor.encrypt(plaintext, ciphertext, true)
}

func (encryptor *pkEncryptor) EncryptFromCRP(plaintext *Plaintext, ciphertext *Element, crp *ring.Poly) {
	panic("Cannot encrypt with CRP using an encryptor created with the public-key")
}

func (encryptor *pkEncryptor) EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Element {
	panic("Cannot encrypt with CRP using an encryptor created with the public-key")
}

func (encryptor *pkEncryptor) EncryptFromCRPNTTNew(plaintext *Plaintext, crp *ring.Poly) *Element {
	panic("Cannot encrypt with CRP using an encryptor created with the public-key")
}

// Encrypt encrypts the input Plaintext using the stored key, and returns the result
// on the receiver Element.
//
// encrypt with pk: ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk: ciphertext = [-a*sk + m + e, a]
func (encryptor *pkEncryptor) encrypt(plaintext *Plaintext, ciphertext *Element, fast bool) {

	lvl := utils.MinInt(plaintext.Level(), ciphertext.Level())

	poolQ0 := encryptor.poolQ[0]
	poolP0 := encryptor.poolP[0]
	poolP1 := encryptor.poolP[1]
	poolP2 := encryptor.poolP[2]

	// We sample a R-WLE instance (encryption of zero) over the extended ring (ciphertext ring + special prime)

	ringQ := encryptor.ringQ

	ciphertextNTT := ciphertext.Value[0].IsNTT

	if fast {

		encryptor.ternarySampler.ReadLvl(lvl, poolQ0)
		ringQ.NTTLvl(lvl, poolQ0, poolQ0)
		ringQ.MFormLvl(lvl, poolQ0, poolQ0)

		// ct0 = u*pk0
		ringQ.MulCoeffsMontgomeryLvl(lvl, poolQ0, encryptor.pk.Value[0][0], ciphertext.Value[0])
		// ct1 = u*pk1
		ringQ.MulCoeffsMontgomeryLvl(lvl, poolQ0, encryptor.pk.Value[1][0], ciphertext.Value[1])

		if ciphertextNTT {

			// ct1 = u*pk1 + e1
			encryptor.gaussianSampler.ReadLvl(lvl, poolQ0)
			ringQ.NTTLvl(lvl, poolQ0, poolQ0)
			ringQ.AddLvl(lvl, ciphertext.Value[1], poolQ0, ciphertext.Value[1])

			// ct0 = u*pk0 + e0
			encryptor.gaussianSampler.ReadLvl(lvl, poolQ0)

			if !plaintext.Value.IsNTT {
				ringQ.AddLvl(lvl, poolQ0, plaintext.Value, poolQ0)
				ringQ.NTTLvl(lvl, poolQ0, poolQ0)
				ringQ.AddLvl(lvl, ciphertext.Value[0], poolQ0, ciphertext.Value[0])
			} else {
				ringQ.NTTLvl(lvl, poolQ0, poolQ0)
				ringQ.AddLvl(lvl, ciphertext.Value[0], poolQ0, ciphertext.Value[0])
				ringQ.AddLvl(lvl, ciphertext.Value[0], plaintext.Value, ciphertext.Value[0])
			}

		} else {

			ringQ.InvNTTLvl(lvl, ciphertext.Value[0], ciphertext.Value[0])
			ringQ.InvNTTLvl(lvl, ciphertext.Value[1], ciphertext.Value[1])

			// ct[0] = pk[0]*u + e0
			encryptor.gaussianSampler.ReadAndAddLvl(ciphertext.Level(), ciphertext.Value[0])

			// ct[1] = pk[1]*u + e1
			encryptor.gaussianSampler.ReadAndAddLvl(ciphertext.Level(), ciphertext.Value[1])

			if !plaintext.Value.IsNTT {
				ringQ.AddLvl(lvl, ciphertext.Value[0], plaintext.Value, ciphertext.Value[0])
			} else {
				ringQ.InvNTTLvl(lvl, plaintext.Value, poolQ0)
				ringQ.AddLvl(lvl, ciphertext.Value[0], poolQ0, ciphertext.Value[0])
			}
		}
	} else {

		ringP := encryptor.ringP

		encryptor.ternarySampler.ReadLvl(lvl, poolQ0)
		ExtendBasisSmallNormAndCenter(ringQ, ringP, poolQ0, poolP0)

		// (#Q + #P) NTT
		ringQ.NTTLvl(lvl, poolQ0, poolQ0)
		ringP.NTT(poolP0, poolP0)

		ringQ.MFormLvl(lvl, poolQ0, poolQ0)
		ringP.MForm(poolP0, poolP0)

		// ct0 = u*pk0
		// ct1 = u*pk1
		ringQ.MulCoeffsMontgomeryLvl(lvl, poolQ0, encryptor.pk.Value[0][0], ciphertext.Value[0])
		ringQ.MulCoeffsMontgomeryLvl(lvl, poolQ0, encryptor.pk.Value[1][0], ciphertext.Value[1])
		ringP.MulCoeffsMontgomery(poolP0, encryptor.pk.Value[1][1], poolP1)
		ringP.MulCoeffsMontgomery(poolP0, encryptor.pk.Value[0][1], poolP0)

		// 2*(#Q + #P) NTT
		ringQ.InvNTTLvl(lvl, ciphertext.Value[0], ciphertext.Value[0])
		ringQ.InvNTTLvl(lvl, ciphertext.Value[1], ciphertext.Value[1])
		ringP.InvNTT(poolP0, poolP0)
		ringP.InvNTT(poolP1, poolP1)

		encryptor.gaussianSampler.ReadLvl(lvl, poolQ0)
		ExtendBasisSmallNormAndCenter(ringQ, ringP, poolQ0, poolP2)
		ringQ.AddLvl(lvl, ciphertext.Value[0], poolQ0, ciphertext.Value[0])
		ringP.Add(poolP0, poolP2, poolP0)

		encryptor.gaussianSampler.ReadLvl(lvl, poolQ0)
		ExtendBasisSmallNormAndCenter(ringQ, ringP, poolQ0, poolP2)
		ringQ.AddLvl(lvl, ciphertext.Value[1], poolQ0, ciphertext.Value[1])
		ringP.Add(poolP1, poolP2, poolP1)

		// ct0 = (u*pk0 + e0)/P
		encryptor.baseconverter.ModDownSplitPQ(lvl, ciphertext.Value[0], poolP0, ciphertext.Value[0])

		// ct1 = (u*pk1 + e1)/P
		encryptor.baseconverter.ModDownSplitPQ(lvl, ciphertext.Value[1], poolP1, ciphertext.Value[1])

		if ciphertextNTT {

			if !plaintext.Value.IsNTT {
				ringQ.AddLvl(lvl, ciphertext.Value[0], plaintext.Value, ciphertext.Value[0])
			}

			// 2*#Q NTT
			ringQ.NTTLvl(lvl, ciphertext.Value[0], ciphertext.Value[0])
			ringQ.NTTLvl(lvl, ciphertext.Value[1], ciphertext.Value[1])

			if plaintext.Value.IsNTT {
				// ct0 = (u*pk0 + e0)/P + m
				ringQ.AddLvl(lvl, ciphertext.Value[0], plaintext.Value, ciphertext.Value[0])
			}

		} else {

			if !plaintext.Value.IsNTT {
				ringQ.AddLvl(lvl, ciphertext.Value[0], plaintext.Value, ciphertext.Value[0])
			} else {
				ringQ.InvNTTLvl(lvl, plaintext.Value, poolQ0)
				ringQ.AddLvl(lvl, ciphertext.Value[0], poolQ0, ciphertext.Value[0])
			}
		}
	}

	ciphertext.Value[1].IsNTT = ciphertext.Value[0].IsNTT

	ciphertext.Value[0].Coeffs = ciphertext.Value[0].Coeffs[:lvl+1]
	ciphertext.Value[1].Coeffs = ciphertext.Value[1].Coeffs[:lvl+1]
}

func (encryptor *skEncryptor) EncryptNew(plaintext *Plaintext) *Element {
	ciphertext := NewElement(encryptor.params, 1, plaintext.Level())
	encryptor.Encrypt(plaintext, ciphertext)
	return ciphertext
}

func (encryptor *skEncryptor) EncryptNTTNew(plaintext *Plaintext) *Element {
	ciphertext := NewElement(encryptor.params, 1, plaintext.Level())
	ciphertext.Value[0].IsNTT = true
	ciphertext.Value[1].IsNTT = true
	encryptor.Encrypt(plaintext, ciphertext)
	return ciphertext
}

func (encryptor *skEncryptor) Encrypt(plaintext *Plaintext, ciphertext *Element) {
	encryptor.uniformSampler.ReadLvl(utils.MinInt(plaintext.Level(), ciphertext.Level()), ciphertext.Value[1])
	encryptor.encrypt(plaintext, ciphertext)
}

func (encryptor *skEncryptor) EncryptFastNew(plaintext *Plaintext) *Element {
	panic("Cannot Encrypt : SkEncryptor doesn't support EncryptFastNew() -> use instead EncryptNew()")
}

func (encryptor *skEncryptor) EncryptFastNTTNew(plaintext *Plaintext) *Element {
	panic("Cannot Encrypt : SkEncryptor doesn't support EncryptFastNew() -> use instead EncryptNew()")
}

func (encryptor *skEncryptor) EncryptFast(plaintext *Plaintext, ciphertext *Element) {
	panic("Cannot Encrypt : SkEncryptor doesn't support EncryptFast() -> use instead Encrypt()")
}

func (encryptor *skEncryptor) EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Element {
	ciphertext := NewElement(encryptor.params, 1, plaintext.Level())
	encryptor.EncryptFromCRP(plaintext, ciphertext, crp)
	return ciphertext
}

func (encryptor *skEncryptor) EncryptFromCRPNTTNew(plaintext *Plaintext, crp *ring.Poly) *Element {
	ciphertext := NewElement(encryptor.params, 1, plaintext.Level())
	ciphertext.Value[0].IsNTT = true
	ciphertext.Value[1].IsNTT = true
	encryptor.EncryptFromCRP(plaintext, ciphertext, crp)
	return ciphertext
}

func (encryptor *skEncryptor) EncryptFromCRP(plaintext *Plaintext, ciphertext *Element, crp *ring.Poly) {
	ring.CopyValues(crp, ciphertext.Value[1])
	encryptor.encrypt(plaintext, ciphertext)
}

func (encryptor *skEncryptor) encrypt(plaintext *Plaintext, ciphertext *Element) {

	ringQ := encryptor.ringQ

	lvl := utils.MinInt(plaintext.Level(), ciphertext.Level())

	poolQ0 := encryptor.poolQ[0]

	ciphertextNTT := ciphertext.Value[0].IsNTT

	ringQ.MulCoeffsMontgomeryLvl(lvl, ciphertext.Value[1], encryptor.sk.Value[0], ciphertext.Value[0])
	ringQ.NegLvl(lvl, ciphertext.Value[0], ciphertext.Value[0])

	if ciphertextNTT {

		encryptor.gaussianSampler.ReadLvl(lvl, poolQ0)

		if plaintext.Value.IsNTT {
			ringQ.NTTLvl(lvl, poolQ0, poolQ0)
			ringQ.AddLvl(lvl, ciphertext.Value[0], poolQ0, ciphertext.Value[0])
			ringQ.AddLvl(lvl, ciphertext.Value[0], plaintext.Value, ciphertext.Value[0])
		} else {
			ringQ.AddLvl(lvl, poolQ0, plaintext.Value, poolQ0)
			ringQ.NTTLvl(lvl, poolQ0, poolQ0)
			ringQ.AddLvl(lvl, ciphertext.Value[0], poolQ0, ciphertext.Value[0])
		}

		ciphertext.Value[0].IsNTT = true
		ciphertext.Value[1].IsNTT = true

	} else {

		if plaintext.Value.IsNTT {
			ringQ.AddLvl(lvl, ciphertext.Value[0], plaintext.Value, ciphertext.Value[0])
			ringQ.InvNTTLvl(lvl, ciphertext.Value[0], ciphertext.Value[0])

		} else {
			ringQ.InvNTTLvl(lvl, ciphertext.Value[0], ciphertext.Value[0])
			ringQ.AddLvl(lvl, ciphertext.Value[0], plaintext.Value, ciphertext.Value[0])
		}

		encryptor.gaussianSampler.ReadAndAddLvl(ciphertext.Level(), ciphertext.Value[0])

		ringQ.InvNTTLvl(lvl, ciphertext.Value[1], ciphertext.Value[1])

		ciphertext.Value[0].IsNTT = false
		ciphertext.Value[1].IsNTT = false

	}

	ciphertext.Value[0].Coeffs = ciphertext.Value[0].Coeffs[:lvl+1]
	ciphertext.Value[1].Coeffs = ciphertext.Value[1].Coeffs[:lvl+1]
}
