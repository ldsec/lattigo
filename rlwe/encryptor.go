package rlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// Encryptor a generic RLWE encryption interface.
type Encryptor interface {
	// Encrypt encrypts the input plaintext and write the result on ctOut.
	// The encryption algorithm depends on the implementor.
	Encrypt(pt *Plaintext, ctOut *Ciphertext)

	// EncryptFromCRP encrypts the input plaintext and writes the result in ctOut.
	// The encryption algorithm depends on the implementor.
	EncryptFromCRP(pt *Plaintext, crp *ring.Poly, ctOut *Ciphertext)
}

// encryptorBase is a struct used to encrypt Plaintexts. It stores the public-key and/or secret-key.
type encryptorBase struct {
	params Parameters

	ringQ *ring.Ring
	ringP *ring.Ring

	poolQ [1]*ring.Poly
	poolP [3]*ring.Poly

	gaussianSampler *ring.GaussianSampler
	ternarySampler  *ring.TernarySampler
	uniformSampler  *ring.UniformSampler
}

type pkEncryptor struct {
	encryptorBase
	pk            *PublicKey
	baseconverter *ring.FastBasisExtender
}

type pkFastEncryptor struct {
	encryptorBase
	pk *PublicKey
}

type skEncryptor struct {
	encryptorBase
	sk *SecretKey
}

// NewEncryptor instatiates a new generic RLWE Encryptor. The key argument can
// be either a *rlwe.PublicKey or a *rlwe.SecretKey.
func NewEncryptor(params Parameters, key interface{}) Encryptor {
	switch key := key.(type) {
	case *PublicKey:
		if key.Value[0].Degree() != params.N() || key.Value[1].Degree() != params.N() {
			panic("cannot newEncryptor: pk ring degree does not match params ring degree")
		}
		encryptorBase := newEncryptorBase(params)
		if params.PCount() > 0 {
			baseconverter := ring.NewFastBasisExtender(params.ringQ, params.ringP)
			return &pkEncryptor{encryptorBase, key, baseconverter}
		}
		return &pkFastEncryptor{encryptorBase, key}
	case *SecretKey:
		if key.Value.Degree() != params.N() {
			panic("cannot newEncryptor: sk ring degree does not match params ring degree")
		}
		return &skEncryptor{newEncryptorBase(params), key}
	default:
		panic("key must be either rlwe.PublicKey or rlwe.SecretKey")
	}
}

// NewFastEncryptor instantiates a new  generic RLWE Encryptor.
// This encryptor's Encrypt method first encrypts zero in Q and then adds the plaintext.
// This method is faster than the normal encryptor but result in a noisier ciphertext.
func NewFastEncryptor(params Parameters, key *PublicKey) Encryptor {
	return &pkFastEncryptor{newEncryptorBase(params), key}
}

// Encrypt encrypts the input Plaintext and write the result in ctOut.
func (encryptor *pkEncryptor) Encrypt(plaintext *Plaintext, ctOut *Ciphertext) {

	lvl := utils.MinInt(plaintext.Level(), ctOut.Level())

	poolQ0 := encryptor.poolQ[0]
	poolP0 := encryptor.poolP[0]
	poolP1 := encryptor.poolP[1]
	poolP2 := encryptor.poolP[2]

	// We sample a R-WLE instance (encryption of zero) over the extended ring (ciphertext ring + special prime)

	ringQ := encryptor.ringQ
	ringP := encryptor.ringP

	ciphertextNTT := ctOut.Value[0].IsNTT

	encryptor.ternarySampler.ReadLvl(lvl, poolQ0)
	extendBasisSmallNormAndCenter(ringQ, ringP, poolQ0, poolP0)

	// (#Q + #P) NTT
	ringQ.NTTLvl(lvl, poolQ0, poolQ0)
	ringP.NTT(poolP0, poolP0)

	ringQ.MFormLvl(lvl, poolQ0, poolQ0)
	ringP.MForm(poolP0, poolP0)

	pk0P := new(ring.Poly)
	pk1P := new(ring.Poly)
	pk0P.Coeffs = encryptor.pk.Value[0].Coeffs[len(ringQ.Modulus):]
	pk1P.Coeffs = encryptor.pk.Value[1].Coeffs[len(ringQ.Modulus):]

	// ct0 = u*pk0
	// ct1 = u*pk1
	ringQ.MulCoeffsMontgomeryLvl(lvl, poolQ0, encryptor.pk.Value[0], ctOut.Value[0])
	ringQ.MulCoeffsMontgomeryLvl(lvl, poolQ0, encryptor.pk.Value[1], ctOut.Value[1])
	ringP.MulCoeffsMontgomery(poolP0, pk1P, poolP1)
	ringP.MulCoeffsMontgomery(poolP0, pk0P, poolP0)

	// 2*(#Q + #P) NTT
	ringQ.InvNTTLvl(lvl, ctOut.Value[0], ctOut.Value[0])
	ringQ.InvNTTLvl(lvl, ctOut.Value[1], ctOut.Value[1])
	ringP.InvNTT(poolP0, poolP0)
	ringP.InvNTT(poolP1, poolP1)

	encryptor.gaussianSampler.ReadLvl(lvl, poolQ0)
	extendBasisSmallNormAndCenter(ringQ, ringP, poolQ0, poolP2)
	ringQ.AddLvl(lvl, ctOut.Value[0], poolQ0, ctOut.Value[0])
	ringP.Add(poolP0, poolP2, poolP0)

	encryptor.gaussianSampler.ReadLvl(lvl, poolQ0)
	extendBasisSmallNormAndCenter(ringQ, ringP, poolQ0, poolP2)
	ringQ.AddLvl(lvl, ctOut.Value[1], poolQ0, ctOut.Value[1])
	ringP.Add(poolP1, poolP2, poolP1)

	// ct0 = (u*pk0 + e0)/P
	encryptor.baseconverter.ModDownSplitPQ(lvl, ctOut.Value[0], poolP0, ctOut.Value[0])

	// ct1 = (u*pk1 + e1)/P
	encryptor.baseconverter.ModDownSplitPQ(lvl, ctOut.Value[1], poolP1, ctOut.Value[1])

	if ciphertextNTT {

		if !plaintext.Value.IsNTT {
			ringQ.AddLvl(lvl, ctOut.Value[0], plaintext.Value, ctOut.Value[0])
		}

		// 2*#Q NTT
		ringQ.NTTLvl(lvl, ctOut.Value[0], ctOut.Value[0])
		ringQ.NTTLvl(lvl, ctOut.Value[1], ctOut.Value[1])

		if plaintext.Value.IsNTT {
			// ct0 = (u*pk0 + e0)/P + m
			ringQ.AddLvl(lvl, ctOut.Value[0], plaintext.Value, ctOut.Value[0])
		}

	} else {

		if !plaintext.Value.IsNTT {
			ringQ.AddLvl(lvl, ctOut.Value[0], plaintext.Value, ctOut.Value[0])
		} else {
			ringQ.InvNTTLvl(lvl, plaintext.Value, poolQ0)
			ringQ.AddLvl(lvl, ctOut.Value[0], poolQ0, ctOut.Value[0])
		}
	}

	ctOut.Value[1].IsNTT = ctOut.Value[0].IsNTT
	ctOut.Value[0].Coeffs = ctOut.Value[0].Coeffs[:lvl+1]
	ctOut.Value[1].Coeffs = ctOut.Value[1].Coeffs[:lvl+1]
}

// Encrypt encrypts the input Plaintext and write the result in ctOut.
// It first encrypts zero in Q and then adds the plaintext.
// This method is faster than the normal encryptor but result in a noisier ciphertext.
func (encryptor *pkFastEncryptor) Encrypt(plaintext *Plaintext, ctOut *Ciphertext) {
	lvl := utils.MinInt(plaintext.Level(), ctOut.Level())

	poolQ0 := encryptor.poolQ[0]

	ringQ := encryptor.ringQ

	ciphertextNTT := ctOut.Value[0].IsNTT

	encryptor.ternarySampler.ReadLvl(lvl, poolQ0)
	ringQ.NTTLvl(lvl, poolQ0, poolQ0)
	ringQ.MFormLvl(lvl, poolQ0, poolQ0)

	// ct0 = u*pk0
	ringQ.MulCoeffsMontgomeryLvl(lvl, poolQ0, encryptor.pk.Value[0], ctOut.Value[0])
	// ct1 = u*pk1
	ringQ.MulCoeffsMontgomeryLvl(lvl, poolQ0, encryptor.pk.Value[1], ctOut.Value[1])

	if ciphertextNTT {

		// ct1 = u*pk1 + e1
		encryptor.gaussianSampler.ReadLvl(lvl, poolQ0)
		ringQ.NTTLvl(lvl, poolQ0, poolQ0)
		ringQ.AddLvl(lvl, ctOut.Value[1], poolQ0, ctOut.Value[1])

		// ct0 = u*pk0 + e0
		encryptor.gaussianSampler.ReadLvl(lvl, poolQ0)

		if !plaintext.Value.IsNTT {
			ringQ.AddLvl(lvl, poolQ0, plaintext.Value, poolQ0)
			ringQ.NTTLvl(lvl, poolQ0, poolQ0)
			ringQ.AddLvl(lvl, ctOut.Value[0], poolQ0, ctOut.Value[0])
		} else {
			ringQ.NTTLvl(lvl, poolQ0, poolQ0)
			ringQ.AddLvl(lvl, ctOut.Value[0], poolQ0, ctOut.Value[0])
			ringQ.AddLvl(lvl, ctOut.Value[0], plaintext.Value, ctOut.Value[0])
		}

	} else {

		ringQ.InvNTTLvl(lvl, ctOut.Value[0], ctOut.Value[0])
		ringQ.InvNTTLvl(lvl, ctOut.Value[1], ctOut.Value[1])

		// ct[0] = pk[0]*u + e0
		encryptor.gaussianSampler.ReadAndAddLvl(ctOut.Level(), ctOut.Value[0])

		// ct[1] = pk[1]*u + e1
		encryptor.gaussianSampler.ReadAndAddLvl(ctOut.Level(), ctOut.Value[1])

		if !plaintext.Value.IsNTT {
			ringQ.AddLvl(lvl, ctOut.Value[0], plaintext.Value, ctOut.Value[0])
		} else {
			ringQ.InvNTTLvl(lvl, plaintext.Value, poolQ0)
			ringQ.AddLvl(lvl, ctOut.Value[0], poolQ0, ctOut.Value[0])
		}
	}

	ctOut.Value[1].IsNTT = ctOut.Value[0].IsNTT

	ctOut.Value[0].Coeffs = ctOut.Value[0].Coeffs[:lvl+1]
	ctOut.Value[1].Coeffs = ctOut.Value[1].Coeffs[:lvl+1]
}

// Encrypt encrypts the input Plaintext and write the result in ctOut.
func (encryptor *skEncryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {
	encryptor.uniformSampler.ReadLvl(utils.MinInt(plaintext.Level(), ciphertext.Level()), ciphertext.Value[1])
	encryptor.encrypt(plaintext, ciphertext)
}

// EncryptFromCRP encrypts the input Plaintext given the uniformly random element c1 and write the result in ctOut.
func (encryptor *skEncryptor) EncryptFromCRP(plaintext *Plaintext, crp *ring.Poly, ctOut *Ciphertext) {
	ring.CopyValues(crp, ctOut.Value[1])
	encryptor.encrypt(plaintext, ctOut)
}

func (encryptor *skEncryptor) encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {

	ringQ := encryptor.ringQ

	lvl := utils.MinInt(plaintext.Level(), ciphertext.Level())

	poolQ0 := encryptor.poolQ[0]

	ciphertextNTT := ciphertext.Value[0].IsNTT

	ringQ.MulCoeffsMontgomeryLvl(lvl, ciphertext.Value[1], encryptor.sk.Value, ciphertext.Value[0])
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

func newEncryptorBase(params Parameters) encryptorBase {

	ringQ := params.RingQ()
	ringP := params.RingP()

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	var poolP [3]*ring.Poly
	if params.PCount() != 0 {
		poolP = [3]*ring.Poly{ringP.NewPoly(), ringP.NewPoly(), ringP.NewPoly()}
	}

	return encryptorBase{
		params:          params,
		ringQ:           ringQ,
		ringP:           ringP,
		poolQ:           [1]*ring.Poly{ringQ.NewPoly()},
		poolP:           poolP,
		gaussianSampler: ring.NewGaussianSampler(prng, ringQ, params.Sigma(), int(6*params.Sigma())),
		ternarySampler:  ring.NewTernarySampler(prng, ringQ, 0.5, false),
		uniformSampler:  ring.NewUniformSampler(prng, ringQ),
	}
}

func (encryptor *encryptorBase) EncryptFromCRP(plaintext *Plaintext, crp *ring.Poly, ctOut *Ciphertext) {
	panic("Cannot encrypt with CRP using an encryptor created with the public-key")
}

func extendBasisSmallNormAndCenter(ringQ, ringP *ring.Ring, polQ, polP *ring.Poly) {
	var coeff, Q, QHalf, sign uint64
	Q = ringQ.Modulus[0]
	QHalf = Q >> 1

	for j := 0; j < ringQ.N; j++ {

		coeff = polQ.Coeffs[0][j]

		sign = 1
		if coeff > QHalf {
			coeff = Q - coeff
			sign = 0
		}

		for i, pi := range ringP.Modulus {
			polP.Coeffs[i][j] = (coeff * sign) | (pi-coeff)*(sign^1)
		}
	}
}
