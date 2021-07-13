package dbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

// E2SProtocol is the structure storing the parameters and temporary buffers
// required by the encryption-to-shares protocol.
type E2SProtocol struct {
	CKSProtocol

	ringQ       *ring.Ring
	ringT       *ring.Ring
	maskSampler *ring.UniformSampler
	encoder     bfv.Encoder

	zero              *rlwe.SecretKey
	tmpPlaintextRingT *bfv.PlaintextRingT
	tmpPlaintext      *bfv.Plaintext
}

// NewE2SProtocol creates a new E2SProtocol struct from the passed BFV parameters.
func NewE2SProtocol(params bfv.Parameters, sigmaSmudging float64) *E2SProtocol {
	e2s := new(E2SProtocol)
	e2s.CKSProtocol = *NewCKSProtocol(params, sigmaSmudging)
	e2s.ringQ = params.RingQ()
	e2s.ringT = params.RingT()
	e2s.encoder = bfv.NewEncoder(params)
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	e2s.maskSampler = ring.NewUniformSampler(prng, e2s.ringT)
	e2s.zero = rlwe.NewSecretKey(params.Parameters)
	e2s.tmpPlaintext = bfv.NewPlaintext(params)
	e2s.tmpPlaintextRingT = bfv.NewPlaintextRingT(params)
	return e2s
}

// GenShare generates a party's share in the encryption-to-shares protocol. This share consist in the additive secret-share of the party
// which is written in secretShareOut and in the public masked-decryption share written in publicShareOut.
func (e2s *E2SProtocol) GenShare(sk *rlwe.SecretKey, ct *bfv.Ciphertext, secretShareOut *rlwe.AdditiveShare, publicShareOut *drlwe.CKSShare) {
	e2s.CKSProtocol.GenShare(sk, e2s.zero, ct.Ciphertext, publicShareOut)
	e2s.maskSampler.Read(&secretShareOut.Value)
	e2s.encoder.ScaleUp(&bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: &secretShareOut.Value}}, e2s.tmpPlaintext)
	e2s.ringQ.Sub(publicShareOut.Value, e2s.tmpPlaintext.Value, publicShareOut.Value)
}

// GetShare is the final step of the encryption-to-share protocol. It performs the masked decryption of the target ciphertext followed by a
// the removal of the caller's secretShare as generated in the GenShare method.
// If the caller is not secret-key-share holder (i.e., didn't generate a decryption share), `secretShare` can be set to nil.
// Therefore, in order to obtain an additive sharing of the message, only one party should call this method, and the other parties should use
// the secretShareOut output of the GenShare method.
func (e2s *E2SProtocol) GetShare(secretShare *rlwe.AdditiveShare, aggregatePublicShare *drlwe.CKSShare, ct *bfv.Ciphertext, secretShareOut *rlwe.AdditiveShare) {
	e2s.ringQ.Add(aggregatePublicShare.Value, ct.Value[0], e2s.tmpPlaintext.Value)
	e2s.encoder.ScaleDown(e2s.tmpPlaintext, e2s.tmpPlaintextRingT)
	if secretShare != nil {
		e2s.ringT.Add(&secretShare.Value, e2s.tmpPlaintextRingT.Value, &secretShareOut.Value)
	} else {
		secretShareOut.Value.Copy(e2s.tmpPlaintextRingT.Value)
	}
}

// S2EProtocol is the structure storing the parameters and temporary buffers
// required by the shares-to-encryption protocol.
type S2EProtocol struct {
	CKSProtocol

	ringQ   *ring.Ring
	encoder bfv.Encoder

	zero         *rlwe.SecretKey
	tmpPlaintext *bfv.Plaintext
}

// NewS2EProtocol creates a new S2EProtocol struct from the passed BFV parameters.
func NewS2EProtocol(params bfv.Parameters, sigmaSmudging float64) *S2EProtocol {
	s2e := new(S2EProtocol)
	s2e.CKSProtocol = *NewCKSProtocol(params, sigmaSmudging)
	s2e.ringQ = params.RingQ()
	s2e.encoder = bfv.NewEncoder(params)
	s2e.zero = rlwe.NewSecretKey(params.Parameters)
	s2e.tmpPlaintext = bfv.NewPlaintext(params)
	return s2e
}

// GenShare generates a party's in the shares-to-encryption protocol given the party's secret-key share `sk`, a common
// polynomial sampled from the CRS `c1` and the party's secret share of the message.
func (s2e *S2EProtocol) GenShare(sk *rlwe.SecretKey, c1 *ring.Poly, secretShare *rlwe.AdditiveShare, c0ShareOut *drlwe.CKSShare) {
	s2e.encoder.ScaleUp(&bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: &secretShare.Value}}, s2e.tmpPlaintext)
	s2e.CKSProtocol.GenShare(s2e.zero, sk, &rlwe.Ciphertext{Value: []*ring.Poly{c0ShareOut.Value, c1}}, c0ShareOut)
	s2e.ringQ.Add(c0ShareOut.Value, s2e.tmpPlaintext.Value, c0ShareOut.Value)
}

// GetEncryption computes the final encryption of the secret-shared message when provided with the aggregation `c0Agg` of the parties'
// share in the protocol and with the common, CRS-sampled polynomial `c1`.
func (s2e *S2EProtocol) GetEncryption(c0Agg *drlwe.CKSShare, c1 *ring.Poly, ctOut *bfv.Ciphertext) {
	if ctOut.Degree() != 1 {
		panic("ctOut must have degree 1.")
	}
	ctOut.Value[0].Copy(c0Agg.Value)
	ctOut.Value[1].Copy(c1)
}
