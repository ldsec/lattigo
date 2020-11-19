package rckks

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
	EncryptNew(plaintext *Plaintext) *Ciphertext

	// Encrypt encrypts the input plaintext using the stored key, and returns
	// the result on the receiver ciphertext. The encryption is done by first
	// encrypting zero in QP, dividing by P and then adding the plaintext.
	Encrypt(plaintext *Plaintext, ciphertext *Ciphertext)

	// EncryptFastNew encrypts the input plaintext using the stored key and returns
	// the result on a newly created ciphertext. The encryption is done by first
	// encrypting zero in Q and then adding the plaintext.
	EncryptFastNew(plaintext *Plaintext) *Ciphertext

	// EncryptFsat encrypts the input plaintext using the stored-key, and returns
	// the result onthe receiver ciphertext. The encryption is done by first
	// encrypting zero in Q and then adding the plaintext.
	EncryptFast(plaintext *Plaintext, ciphertext *Ciphertext)

	// EncryptFromCRPNew encrypts the input plaintext using the stored key and returns
	// the result on a newly created ciphertext. The encryption is done by first encrypting
	// zero in QP, using the provided polynomial as the uniform polynomial, dividing by P and
	// then adding the plaintext.
	EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext

	// EncryptFromCRP encrypts the input plaintext using the stored key and returns
	// the result tge receiver ciphertext. The encryption is done by first encrypting
	// zero in QP, using the provided polynomial as the uniform polynomial, dividing by P and
	// then adding the plaintext.
	EncryptFromCRP(plaintext *Plaintext, ciphertetx *Ciphertext, crp *ring.Poly)

	// EncryptFromCRPNew encrypts the input plaintext using the stored key and returns
	// the result on a newly created ciphertext. The encryption is done by first encrypting
	// zero in Q, using the provided polynomial as the uniform polynomial, and
	// then adding the plaintext.
	EncryptFromCRPFastNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext

	// EncryptFromCRP encrypts the input plaintext using the stored key and returns
	// the result tge receiver ciphertext. The encryption is done by first encrypting
	// zero in Q, using the provided polynomial as the uniform polynomial, and
	// then adding the plaintext.
	EncryptFromCRPFast(plaintext *Plaintext, ciphertetx *Ciphertext, crp *ring.Poly)
}

// encryptor is a struct used to encrypt Plaintexts. It stores the public-key and/or secret-key.
type encryptor struct {
	params *Parameters

	ringQ  *ring.Ring
	ringQP *ring.Ring

	polypool [3]*ring.Poly

	baseconverter              *ring.FastBasisExtender
	gaussianSamplerQ           *ring.GaussianSampler
	uniformSamplerQ            *ring.UniformSampler
	ternarySamplerMontgomeryQ  *ring.TernarySampler
	gaussianSamplerQP          *ring.GaussianSampler
	uniformSamplerQP           *ring.UniformSampler
	ternarySamplerMontgomeryQP *ring.TernarySampler
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
func NewEncryptorFromPk(params *Parameters, pk *PublicKey) Encryptor {
	enc := newEncryptor(params)

	if uint64(pk.pk[0].GetDegree()) != params.N() || uint64(pk.pk[1].GetDegree()) != params.N() {
		panic("cannot newEncrpytor: pk ring degree does not match params ring degree")
	}

	return &pkEncryptor{enc, pk}
}

// NewEncryptorFromSk creates a new Encryptor with the provided secret-key.
// This Encryptor can be used to encrypt Plaintexts, using the stored key.
func NewEncryptorFromSk(params *Parameters, sk *SecretKey) Encryptor {
	enc := newEncryptor(params)

	if uint64(sk.sk.GetDegree()) != params.N() {
		panic("cannot newEncryptor: sk ring degree does not match params ring degree")
	}

	return &skEncryptor{enc, sk}
}

func newEncryptor(params *Parameters) encryptor {

	var q, qp *ring.Ring
	var err error
	if q, err = ring.NewRing(params.N(), params.qi); err != nil {
		panic(err)
	}

	if q, err = ring.NewRingWithNthRoot(params.N(), params.N()<<2, params.qi); err != nil {
		panic(err)
	}

	var baseconverter *ring.FastBasisExtender
	if params.PiCount() != 0 {

		if qp, err = ring.NewRingWithNthRoot(params.N(), params.N()<<2, append(params.qi, params.pi...)); err != nil {
			panic(err)
		}

		p, err := ring.NewRingWithNthRoot(params.N(), params.N()<<2, params.pi)
		if err != nil {
			panic(err)
		}

		baseconverter = ring.NewFastBasisExtender(q, p)
	}

	prng, err := utils.NewKeyedPRNG(nil)
	if err != nil {
		panic(err)
	}

	return encryptor{
		params:                     params.Copy(),
		ringQ:                      q,
		ringQP:                     qp,
		polypool:                   [3]*ring.Poly{qp.NewPoly(), qp.NewPoly(), qp.NewPoly()},
		baseconverter:              baseconverter,
		gaussianSamplerQ:           ring.NewGaussianSampler(prng, q, params.sigma, uint64(6*params.sigma)),
		uniformSamplerQ:            ring.NewUniformSampler(prng, q),
		ternarySamplerMontgomeryQ:  ring.NewTernarySampler(prng, q, 0.5, true),
		gaussianSamplerQP:          ring.NewGaussianSampler(prng, qp, params.sigma, uint64(6*params.sigma)),
		uniformSamplerQP:           ring.NewUniformSampler(prng, qp),
		ternarySamplerMontgomeryQP: ring.NewTernarySampler(prng, qp, 0.5, true),
	}
}

// EncryptNew encrypts the input Plaintext using the stored key and returns
// the result on a newly created Ciphertext.
//
// encrypt with pk: ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk: ciphertext = [-a*sk + m + e, a]
func (encryptor *pkEncryptor) EncryptNew(plaintext *Plaintext) *Ciphertext {

	if encryptor.baseconverter == nil {
		panic("Cannot EncryptNew : modulus P is empty -> use instead EncryptFastNew")
	}

	ciphertext := NewCiphertext(encryptor.params, 1, plaintext.Level(), plaintext.Scale())
	encryptor.encrypt(plaintext, ciphertext, false)

	return ciphertext
}

func (encryptor *pkEncryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {

	if encryptor.baseconverter == nil {
		panic("Cannot Encrypt : modulus P is empty -> use instead EncryptFast")
	}

	encryptor.encrypt(plaintext, ciphertext, false)
}

func (encryptor *pkEncryptor) EncryptFastNew(plaintext *Plaintext) *Ciphertext {
	ciphertext := NewCiphertext(encryptor.params, 1, plaintext.Level(), plaintext.Scale())
	encryptor.encrypt(plaintext, ciphertext, true)

	return ciphertext
}

func (encryptor *pkEncryptor) EncryptFast(plaintext *Plaintext, ciphertext *Ciphertext) {
	encryptor.encrypt(plaintext, ciphertext, true)
}

func (encryptor *pkEncryptor) EncryptFromCRP(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly) {
	panic("Cannot encrypt with CRP using an encryptor created with the public-key")
}

func (encryptor *pkEncryptor) EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext {
	panic("Cannot encrypt with CRP using an encryptor created with the public-key")
}

func (encryptor *pkEncryptor) EncryptFromCRPFast(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly) {
	panic("Cannot encrypt with CRP using an encryptor created with the public-key")
}

func (encryptor *pkEncryptor) EncryptFromCRPFastNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext {
	panic("Cannot encrypt with CRP using an encryptor created with the public-key")
}

// Encrypt encrypts the input Plaintext using the stored key, and returns the result
// on the receiver Ciphertext.
//
// encrypt with pk: ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk: ciphertext = [-a*sk + m + e, a]
func (encryptor *pkEncryptor) encrypt(plaintext *Plaintext, ciphertext *Ciphertext, fast bool) {

	// We sample a R-WLE instance (encryption of zero) over the extended ring (ciphertext ring + special prime)

	ringQ := encryptor.ringQ

	if fast {

		level := encryptor.params.QiCount() - 1

		encryptor.ternarySamplerMontgomeryQ.Read(encryptor.polypool[2])

		NTTRCKKS(ringQ, encryptor.polypool[2], encryptor.polypool[2])

		// ct0 = u*pk0
		ringQ.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[0], ciphertext.value[0])
		// ct1 = u*pk1
		ringQ.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[1], ciphertext.value[1])

		// ct1 = u*pk1 + e1
		encryptor.gaussianSamplerQ.ReadLvl(level, encryptor.polypool[0])

		NTTRCKKS(ringQ, encryptor.polypool[0], encryptor.polypool[0])

		ringQ.Add(ciphertext.value[1], encryptor.polypool[0], ciphertext.value[1])

		if !plaintext.isNTT {

			// ct0 = u*pk0 + e0
			encryptor.gaussianSamplerQ.ReadLvl(level, encryptor.polypool[0])

			// ct0 = (u*pk0 + e0)/P + m
			ringQ.Add(encryptor.polypool[0], plaintext.value, encryptor.polypool[0])

			NTTRCKKS(ringQ, encryptor.polypool[0], encryptor.polypool[0])

			ringQ.Add(ciphertext.value[0], encryptor.polypool[0], ciphertext.value[0])

		} else {
			// ct0 = u*pk0 + e0
			encryptor.gaussianSamplerQ.ReadLvl(level, encryptor.polypool[0])

			NTTRCKKS(ringQ, encryptor.polypool[0], encryptor.polypool[0])

			ringQ.Add(ciphertext.value[0], encryptor.polypool[0], ciphertext.value[0])
			ringQ.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])
		}

	} else {

		ringQP := encryptor.ringQP

		level := uint64(len(ringQP.Modulus) - 1)

		encryptor.ternarySamplerMontgomeryQP.Read(encryptor.polypool[2])

		NTTRCKKS(ringQP, encryptor.polypool[2], encryptor.polypool[2])

		// ct0 = u*pk0
		ringQP.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[0], encryptor.polypool[0])
		// ct1 = u*pk1
		ringQP.MulCoeffsMontgomery(encryptor.polypool[2], encryptor.pk.pk[1], encryptor.polypool[1])

		// 2*(#Q + #P) NTT
		InvNTTRCKKS(ringQP, encryptor.polypool[0], encryptor.polypool[0])

		InvNTTRCKKS(ringQP, encryptor.polypool[1], encryptor.polypool[1])

		// ct0 = u*pk0 + e0
		encryptor.gaussianSamplerQP.ReadAndAddLvl(level, encryptor.polypool[0])
		// ct1 = u*pk1 + e1
		encryptor.gaussianSamplerQP.ReadAndAddLvl(level, encryptor.polypool[1])

		// ct0 = (u*pk0 + e0)/P
		encryptor.baseconverter.ModDownPQ(plaintext.Level(), encryptor.polypool[0], ciphertext.value[0])

		// ct1 = (u*pk1 + e1)/P
		encryptor.baseconverter.ModDownPQ(plaintext.Level(), encryptor.polypool[1], ciphertext.value[1])

		if !plaintext.isNTT {
			ringQ.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])
		}

		level = uint64(len(ringQ.Modulus) - 1)

		// 2*#Q NTT
		NTTRCKKS(ringQ, ciphertext.value[0], ciphertext.value[0])
		NTTRCKKS(ringQ, ciphertext.value[1], ciphertext.value[1])

		if plaintext.isNTT {
			// ct0 = (u*pk0 + e0)/P + m
			ringQ.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])
		}
	}

	ciphertext.isNTT = true
}

func (encryptor *skEncryptor) EncryptNew(plaintext *Plaintext) *Ciphertext {

	if encryptor.baseconverter == nil {
		panic("Cannot EncryptNew : modulus P is empty -> use instead EncryptFastNew")
	}

	ciphertext := NewCiphertext(encryptor.params, 1, plaintext.Level(), plaintext.Scale())
	encryptor.Encrypt(plaintext, ciphertext)
	return ciphertext
}

func (encryptor *skEncryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {
	encryptor.encryptSample(plaintext, ciphertext)
}

func (encryptor *skEncryptor) EncryptFastNew(plaintext *Plaintext) *Ciphertext {
	panic("Cannot Encrypt : SkEncryptor doesn't support EncryptFastNew() -> use instead EncryptNew()")
}

func (encryptor *skEncryptor) EncryptFast(plaintext *Plaintext, ciphertext *Ciphertext) {
	panic("Cannot Encrypt : SkEncryptor doesn't support EncryptFast() -> use instead Encrypt()")
}

func (encryptor *skEncryptor) EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext {
	ciphertext := NewCiphertext(encryptor.params, 1, plaintext.Level(), plaintext.Scale())
	encryptor.EncryptFromCRP(plaintext, ciphertext, crp)
	return ciphertext
}

func (encryptor *skEncryptor) EncryptFromCRP(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly) {
	encryptor.encryptFromCRP(plaintext, ciphertext, crp)
}

func (encryptor *skEncryptor) EncryptFromCRPFastNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext {
	panic("Cannot Encrypt : SkEncryptor doesn't support EncryptFromCRPFastNew() -> use instead EncryptFromCRPNew()")
}

func (encryptor *skEncryptor) EncryptFromCRPFast(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly) {
	panic("Cannot Encrypt : SkEncryptor doesn't support EncryptFromCRPFast() -> use instead EncryptFromCRP()")

}

func (encryptor *skEncryptor) encryptSample(plaintext *Plaintext, ciphertext *Ciphertext) {
	encryptor.uniformSamplerQ.Read(ciphertext.value[1])
	encryptor.encrypt(plaintext, ciphertext, ciphertext.value[1])
}

func (encryptor *skEncryptor) encryptFromCRP(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly) {
	encryptor.ringQ.Copy(crp, ciphertext.value[1])
	encryptor.encrypt(plaintext, ciphertext, ciphertext.value[1])
}

func (encryptor *skEncryptor) encrypt(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly) {

	ringQ := encryptor.ringQ

	level := encryptor.params.QiCount() - 1

	ringQ.MulCoeffsMontgomery(ciphertext.value[1], encryptor.sk.sk, ciphertext.value[0])
	ringQ.Neg(ciphertext.value[0], ciphertext.value[0])

	if plaintext.isNTT {
		encryptor.gaussianSamplerQ.ReadLvl(level, encryptor.polypool[0])

		NTTRCKKS(ringQ, encryptor.polypool[0], encryptor.polypool[0])

		ringQ.Add(ciphertext.value[0], encryptor.polypool[0], ciphertext.value[0])
		ringQ.Add(ciphertext.value[0], plaintext.value, ciphertext.value[0])
	} else {
		encryptor.gaussianSamplerQ.ReadLvl(level, encryptor.polypool[0])
		ringQ.Add(encryptor.polypool[0], plaintext.value, encryptor.polypool[0])

		NTTRCKKS(ringQ, encryptor.polypool[0], encryptor.polypool[0])

		ringQ.Add(ciphertext.value[0], encryptor.polypool[0], ciphertext.value[0])
	}

	ciphertext.isNTT = true
}
