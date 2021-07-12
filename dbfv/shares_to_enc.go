package dbfv

import (
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

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
	e2s := new(S2EProtocol)
	e2s.CKSProtocol = *NewCKSProtocol(params, sigmaSmudging)
	e2s.ringQ = params.RingQ()
	e2s.encoder = bfv.NewEncoder(params)
	e2s.zero = rlwe.NewSecretKey(params.Parameters)
	e2s.tmpPlaintext = bfv.NewPlaintext(params)
	return e2s
}

// GenShare generates a party's in the shares-to-encryption protocol given the party's secret-key share `sk`, a common
// polynomial sampled from the CRS `c1` and the party's secret share of the message.
func (s2e *S2EProtocol) GenShare(sk *rlwe.SecretKey, c1 *ring.Poly, secretShare rlwe.AdditiveShare, c0ShareOut *drlwe.CKSShare) {
	s2e.encoder.ScaleUp(&bfv.PlaintextRingT{Plaintext: &rlwe.Plaintext{Value: &secretShare.Value}}, s2e.tmpPlaintext)
	s2e.CKSProtocol.GenShare(s2e.zero, sk, &rlwe.Element{Value: []*ring.Poly{c0ShareOut.Value, c1}}, c0ShareOut)
	s2e.ringQ.Add(c0ShareOut.Value, s2e.tmpPlaintext.Value, c0ShareOut.Value)
}

// GetEncryption computes the final encryption of the secret-shared message when provided with the aggregation `c0Agg` of the parties'
// share in the protocol and with the common, CRS-sampled polynomial `c1`.
func (s2e *S2EProtocol) GetEncryption(c0Agg *drlwe.CKSShare, c1 *ring.Poly, ctOut bfv.Ciphertext) {
	if ctOut.Degree() != 1 {
		panic("ctOut must have degree 1.")
	}
	ctOut.Value[0].Copy(c0Agg.Value)
	ctOut.Value[1].Copy(c1)
}
