package dckks

import (
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"math/big"
)

// S2EProtocol is the structure storing the parameters and temporary buffers
// required by the shares-to-encryption protocol.
type S2EProtocol struct {
	CKSProtocol
	ringQ    *ring.Ring
	tmp      *ring.Poly
	ssBigint []*big.Int
	zero     *rlwe.SecretKey
}

// NewS2EProtocol creates a new S2EProtocol struct from the passed CKKS parameters.
func NewS2EProtocol(params ckks.Parameters, sigmaSmudging float64) *S2EProtocol {
	s2e := new(S2EProtocol)
	s2e.CKSProtocol = *NewCKSProtocol(params, sigmaSmudging)
	s2e.ringQ = params.RingQ()
	s2e.tmp = s2e.ringQ.NewPoly()
	s2e.ssBigint = make([]*big.Int, s2e.ringQ.N)
	s2e.zero = rlwe.NewSecretKey(params.Parameters)
	return s2e
}

// AllocateShare allocates a share of the S2E protocol
func (s2e S2EProtocol) AllocateShare(level int) (share *drlwe.CKSShare) {
	share = s2e.CKSProtocol.AllocateShare(level)
	share.Value.IsNTT = true
	return
}

// GenShare generates a party's in the shares-to-encryption protocol given the party's secret-key share `sk`, a common
// polynomial sampled from the CRS `c1` and the party's secret share of the message.
func (s2e *S2EProtocol) GenShare(sk *rlwe.SecretKey, c1 *ring.Poly, secretShare rlwe.AdditiveShareBigint, c0ShareOut *drlwe.CKSShare) {

	if c1.Level() != c0ShareOut.Value.Level() {
		panic("c1 and c0ShareOut level must be equal")
	}

	// Generates an encryption share
	c1.IsNTT = true
	s2e.CKSProtocol.GenShare(s2e.zero, sk, &rlwe.Element{Value: []*ring.Poly{nil, c1}}, c0ShareOut)

	s2e.ringQ.SetCoefficientsBigintLvl(c1.Level(), secretShare.Value, s2e.tmp)
	s2e.ringQ.NTTLvl(c1.Level(), s2e.tmp, s2e.tmp)
	s2e.ringQ.AddLvl(c1.Level(), c0ShareOut.Value, s2e.tmp, c0ShareOut.Value)
}

// GetEncryption computes the final encryption of the secret-shared message when provided with the aggregation `c0Agg` of the parties'
// share in the protocol and with the common, CRS-sampled polynomial `c1`.
func (s2e *S2EProtocol) GetEncryption(c0Agg *drlwe.CKSShare, c1 *ring.Poly, ctOut ckks.Ciphertext) {

	if ctOut.Degree() != 1 {
		panic("ctOut must have degree 1.")
	}

	if c0Agg.Value.Level() != c1.Level() {
		panic("c0Agg level must be equal to c1 level")
	}

	if ctOut.Level() != c1.Level() {
		panic("ctOut level must be equal to c1 level")
	}

	ctOut.Value[0].Copy(c0Agg.Value)
	ctOut.Value[1].Copy(c1)
}
