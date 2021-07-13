package dckks

import (
	"math/big"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/drlwe"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// E2SProtocol is the structure storing the parameters and temporary buffers
// required by the encryption-to-shares protocol.
type E2SProtocol struct {
	CKSProtocol

	ringQ      *ring.Ring
	zero       *rlwe.SecretKey
	maskBigint []*big.Int
	pool       *ring.Poly
}

// NewE2SProtocol creates a new E2SProtocol struct from the passed CKKS parameters.
func NewE2SProtocol(params ckks.Parameters, sigmaSmudging float64) *E2SProtocol {
	e2s := new(E2SProtocol)
	e2s.CKSProtocol = *NewCKSProtocol(params, sigmaSmudging)
	e2s.ringQ = params.RingQ()
	e2s.zero = rlwe.NewSecretKey(params.Parameters)
	e2s.maskBigint = make([]*big.Int, params.N())
	for i := range e2s.maskBigint {
		e2s.maskBigint[i] = new(big.Int)
	}
	e2s.pool = e2s.ringQ.NewPoly()
	return e2s
}

// AllocateShare allocates a share of the E2S protocol
func (e2s E2SProtocol) AllocateShare(level int) (share *drlwe.CKSShare) {
	share = e2s.CKSProtocol.AllocateShare(level)
	share.Value.IsNTT = true
	return
}

// GenShare generates a party's share in the encryption-to-shares protocol. This share consist in the additive secret-share of the party
// which is written in secretShareOut and in the public masked-decryption share written in publicShareOut.
// This protocol requires additional inputs which are :
// logBound : the bit length of the masks
// logSlots : the bit length of the number of slots
//
// The method "GetMinimumLevelForBootstrapping" should be used to get the minimum level at which E2S can be called while still ensure 128-bits of security, as well as the
// value for logBound.
func (e2s *E2SProtocol) GenShare(sk *rlwe.SecretKey, logBound, logSlots int, ct *ckks.Ciphertext, secretShareOut *rlwe.AdditiveShareBigint, publicShareOut *drlwe.CKSShare) {

	ringQ := e2s.ringQ

	// Get the upperbound on the norm
	// Ensures that bound >= 2^{128+logbound}
	bound := ring.NewUint(1)
	bound.Lsh(bound, uint(logBound))

	boundMax := ring.NewUint(ringQ.Modulus[0])
	for i := 1; i < ct.Level()+1; i++ {
		boundMax.Mul(boundMax, ring.NewUint(ringQ.Modulus[i]))
	}

	var sign int

	sign = bound.Cmp(boundMax)

	if sign == 1 || bound.Cmp(boundMax) == 1 {
		panic("ciphertext level is not large enough for refresh correctness")
	}

	gap := ringQ.N / (2 << logSlots)

	boundHalf := new(big.Int).Rsh(bound, 1)

	// Generate the mask in Z[Y] for Y = X^{N/(2*slots)}
	for i := 0; i < ringQ.N; i++ {

		if i%gap == 0 {
			e2s.maskBigint[i] = ring.RandInt(bound)
			sign = e2s.maskBigint[i].Cmp(boundHalf)
			if sign == 1 || sign == 0 {
				e2s.maskBigint[i].Sub(e2s.maskBigint[i], bound)
			}
		} else {
			e2s.maskBigint[i].SetUint64(0)
		}

	}

	// Set the mask on the out secret-share
	for i := range e2s.maskBigint {
		secretShareOut.Value[i].Set(e2s.maskBigint[i])
	}

	// Encrypt the mask
	// Generates an encryption of zero and subtracts the mask
	e2s.CKSProtocol.GenShare(sk, e2s.zero, ct.Ciphertext, publicShareOut)
	// Puts the mask in a poly
	e2s.ringQ.SetCoefficientsBigintLvl(ct.Level(), secretShareOut.Value, e2s.pool)
	// NTT the poly
	e2s.ringQ.NTTLvl(ct.Level(), e2s.pool, e2s.pool)
	// Substracts the mask to the encryption of zero
	e2s.ringQ.SubLvl(ct.Level(), publicShareOut.Value, e2s.pool, publicShareOut.Value)
}

// GetShare is the final step of the encryption-to-share protocol. It performs the masked decryption of the target ciphertext followed by a
// the removal of the caller's secretShare as generated in the GenShare method.
// If the caller is not secret-key-share holder (i.e., didn't generate a decryption share), `secretShare` can be set to nil.
// Therefore, in order to obtain an additive sharing of the message, only one party should call this method, and the other parties should use
// the secretShareOut output of the GenShare method.
func (e2s *E2SProtocol) GetShare(secretShare *rlwe.AdditiveShareBigint, aggregatePublicShare *drlwe.CKSShare, ct *ckks.Ciphertext, secretShareOut *rlwe.AdditiveShareBigint) {

	e2s.pool.Zero()

	// Adds the decryption share on the ciphertext and stores the result in a pool
	e2s.ringQ.AddLvl(ct.Level(), aggregatePublicShare.Value, ct.Value[0], e2s.pool)

	// Switches the LSSS RNS NTT ciphertext outside of the NTT domain
	e2s.ringQ.InvNTTLvl(ct.Level(), e2s.pool, e2s.pool)

	// Switches the LSSS RNS ciphertext outside of the RNS domain
	e2s.ringQ.PolyToBigintCenteredLvl(ct.Level(), e2s.pool, e2s.maskBigint)

	// Substracts the last mask
	if secretShare != nil {
		a := secretShareOut.Value
		b := e2s.maskBigint
		c := secretShare.Value
		for i := range secretShareOut.Value {
			a[i].Add(c[i], b[i])
		}
	} else {
		a := secretShareOut.Value
		b := e2s.maskBigint
		for i := range secretShareOut.Value {
			a[i].Set(b[i])
		}
	}
}

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
func (s2e *S2EProtocol) GenShare(sk *rlwe.SecretKey, c1 *ring.Poly, secretShare *rlwe.AdditiveShareBigint, c0ShareOut *drlwe.CKSShare) {

	if c1.Level() != c0ShareOut.Value.Level() {
		panic("c1 and c0ShareOut level must be equal")
	}

	// Generates an encryption share
	c1.IsNTT = true
	s2e.CKSProtocol.GenShare(s2e.zero, sk, &rlwe.Ciphertext{Value: []*ring.Poly{nil, c1}}, c0ShareOut)

	s2e.ringQ.SetCoefficientsBigintLvl(c1.Level(), secretShare.Value, s2e.tmp)
	s2e.ringQ.NTTLvl(c1.Level(), s2e.tmp, s2e.tmp)
	s2e.ringQ.AddLvl(c1.Level(), c0ShareOut.Value, s2e.tmp, c0ShareOut.Value)
}

// GetEncryption computes the final encryption of the secret-shared message when provided with the aggregation `c0Agg` of the parties'
// share in the protocol and with the common, CRS-sampled polynomial `c1`.
func (s2e *S2EProtocol) GetEncryption(c0Agg *drlwe.CKSShare, c1 *ring.Poly, ctOut *ckks.Ciphertext) {

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
