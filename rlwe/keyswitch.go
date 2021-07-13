package rlwe

import (
	"github.com/ldsec/lattigo/v2/ring"
	"math"
)

// KeySwitcher is a struct for RLWE key-switching.
type KeySwitcher struct {
	*Parameters
	*keySwitcherBuffer
	Baseconverter *ring.FastBasisExtender
	Decomposer    *ring.Decomposer
}

type keySwitcherBuffer struct {
	// PoolQ[0]/PoolP[0] : on the fly decomp(c2)
	// PoolQ[1-5]/PoolP[1-5] : available
	PoolQ       [6]*ring.Poly
	PoolP       [6]*ring.Poly
	PoolInvNTT  *ring.Poly
	PoolDecompQ []*ring.Poly // Memory pool for the basis extension in hoisting
	PoolDecompP []*ring.Poly // Memory pool for the basis extension in hoisting
}

func newKeySwitcherBuffer(params Parameters) *keySwitcherBuffer {

	buff := new(keySwitcherBuffer)
	beta := params.Beta()
	ringQ := params.RingQ()
	ringP := params.RingP()

	buff.PoolQ = [6]*ring.Poly{ringQ.NewPoly(), ringQ.NewPoly(), ringQ.NewPoly(), ringQ.NewPoly(), ringQ.NewPoly(), ringQ.NewPoly()}
	buff.PoolP = [6]*ring.Poly{ringP.NewPoly(), ringP.NewPoly(), ringP.NewPoly(), ringP.NewPoly(), ringP.NewPoly(), ringP.NewPoly()}

	buff.PoolInvNTT = ringQ.NewPoly()

	buff.PoolDecompQ = make([]*ring.Poly, beta)
	buff.PoolDecompP = make([]*ring.Poly, beta)

	for i := 0; i < beta; i++ {
		buff.PoolDecompQ[i] = ringQ.NewPoly()
		buff.PoolDecompP[i] = ringP.NewPoly()
	}

	return buff
}

// NewKeySwitcher creates a new KeySwitcher.
func NewKeySwitcher(params Parameters) *KeySwitcher {
	ks := new(KeySwitcher)
	ks.Parameters = &params
	ks.Baseconverter = ring.NewFastBasisExtender(params.RingQ(), params.RingP())
	ks.Decomposer = ring.NewDecomposer(params.RingQ().Modulus, params.RingP().Modulus)
	ks.keySwitcherBuffer = newKeySwitcherBuffer(params)
	return ks
}

// ShallowCopy creates a copy of a KeySwitcher, only reallocating the memory pool.
func (ks *KeySwitcher) ShallowCopy() *KeySwitcher {
	return &KeySwitcher{
		Parameters:        ks.Parameters,
		Decomposer:        ks.Decomposer,
		keySwitcherBuffer: newKeySwitcherBuffer(*ks.Parameters),
		Baseconverter:     ks.Baseconverter.ShallowCopy(),
	}
}

// SwitchKeysInPlace applies the general key-switching procedure of the form [c0 + cx*evakey[0], c1 + cx*evakey[1]]
// Will return the result in the same NTT domain as the input cx.
func (ks *KeySwitcher) SwitchKeysInPlace(level int, cx *ring.Poly, evakey *SwitchingKey, p0, p1 *ring.Poly) {
	ks.SwitchKeysInPlaceNoModDown(level, cx, evakey, p0, ks.PoolP[1], p1, ks.PoolP[2])

	if cx.IsNTT {
		ks.Baseconverter.ModDownSplitNTTPQ(level, p0, ks.PoolP[1], p0)
		ks.Baseconverter.ModDownSplitNTTPQ(level, p1, ks.PoolP[2], p1)
	} else {

		ks.ringQ.InvNTTLazyLvl(level, p0, p0)
		ks.ringQ.InvNTTLazyLvl(level, p1, p1)
		ks.ringP.InvNTTLazy(ks.PoolP[1], ks.PoolP[1])
		ks.ringP.InvNTTLazy(ks.PoolP[2], ks.PoolP[2])

		ks.Baseconverter.ModDownSplitPQ(level, p0, ks.PoolP[1], p0)
		ks.Baseconverter.ModDownSplitPQ(level, p1, ks.PoolP[2], p1)
	}
}

// DecomposeNTT applies the full RNS basis decomposition for all q_alpha_i on c2.
// Expects the IsNTT flag of c2 to correctly reflect the domain of c2.
// PoolDecompQ and PoolDecompQ are vectors of polynomials (mod Q and mod P) that store the
// special RNS decomposition of c2 (in the NTT domain)
func (ks *KeySwitcher) DecomposeNTT(levelQ int, c2 *ring.Poly, PoolDecompQ, PoolDecompP []*ring.Poly) {

	ringQ := ks.RingQ()

	var polyNTT, polyInvNTT *ring.Poly

	if c2.IsNTT {
		polyNTT = c2
		polyInvNTT = ks.PoolInvNTT
		ringQ.InvNTTLvl(levelQ, polyNTT, polyInvNTT)
	} else {
		polyNTT = ks.PoolInvNTT
		polyInvNTT = c2
		ringQ.NTTLvl(levelQ, polyInvNTT, polyNTT)
	}

	alpha := ks.Parameters.PCount()
	beta := int(math.Ceil(float64(levelQ+1) / float64(alpha)))

	for i := 0; i < beta; i++ {
		ks.DecomposeSingleNTT(levelQ, i, polyNTT, polyInvNTT, PoolDecompQ[i], PoolDecompP[i])
	}
}

// DecomposeSingleNTT takes the input polynomial c2 (c2NTT and c2InvNTT, respectively in the NTT and out of the NTT domain)
// modulo q_alpha_beta, and returns the result on c2QiQ are c2QiP the receiver polynomials
// respectively mod Q and mod P (in the NTT domain)
func (ks *KeySwitcher) DecomposeSingleNTT(level, beta int, c2NTT, c2InvNTT, c2QiQ, c2QiP *ring.Poly) {

	ringQ := ks.RingQ()
	ringP := ks.RingP()

	ks.Decomposer.DecomposeAndSplit(level, beta, c2InvNTT, c2QiQ, c2QiP)

	p0idxst := beta * len(ringP.Modulus)
	p0idxed := p0idxst + ks.Decomposer.Xalpha()[beta]

	// c2_qi = cx mod qi mod qi
	for x := 0; x < level+1; x++ {
		if p0idxst <= x && x < p0idxed {
			copy(c2QiQ.Coeffs[x], c2NTT.Coeffs[x])
		} else {
			ring.NTTLazy(c2QiQ.Coeffs[x], c2QiQ.Coeffs[x], ringQ.N, ringQ.NttPsi[x], ringQ.Modulus[x], ringQ.MredParams[x], ringQ.BredParams[x])
		}
	}
	// c2QiP = c2 mod qi mod pj
	ringP.NTTLazy(c2QiP, c2QiP)
}

// SwitchKeysInPlaceNoModDown applies the key-switch to the polynomial cx :
//
// pool2 = dot(decomp(cx) * evakey[0]) mod QP (encrypted input is multiplied by P factor)
// pool3 = dot(decomp(cx) * evakey[1]) mod QP (encrypted input is multiplied by P factor)
//
// Expects the flag IsNTT of cx to correctly reflect the domain of cx.
func (ks *KeySwitcher) SwitchKeysInPlaceNoModDown(level int, cx *ring.Poly, evakey *SwitchingKey, pool2Q, pool2P, pool3Q, pool3P *ring.Poly) {

	var reduce int

	ringQ := ks.ringQ
	ringP := ks.ringP

	c2QiQ := ks.PoolQ[0]
	c2QiP := ks.PoolP[0]

	var cxNTT, cxInvNTT *ring.Poly
	if cx.IsNTT {
		cxNTT = cx
		cxInvNTT = ks.PoolInvNTT
		ringQ.InvNTTLvl(level, cxNTT, cxInvNTT)
	} else {
		cxNTT = ks.PoolInvNTT
		cxInvNTT = cx
		ringQ.NTTLvl(level, cxInvNTT, cxNTT)
	}

	evakey0Q := new(ring.Poly)
	evakey1Q := new(ring.Poly)
	evakey0P := new(ring.Poly)
	evakey1P := new(ring.Poly)

	reduce = 0

	alpha := len(ringP.Modulus)
	beta := int(math.Ceil(float64(level+1) / float64(alpha)))

	QiOverF := ks.Parameters.QiOverflowMargin(level) >> 1
	PiOverF := ks.Parameters.PiOverflowMargin() >> 1

	// Key switching with CRT decomposition for the Qi
	for i := 0; i < beta; i++ {

		ks.DecomposeSingleNTT(level, i, cxNTT, cxInvNTT, c2QiQ, c2QiP)

		evakey0Q.Coeffs = evakey.Value[i][0].Coeffs[:level+1]
		evakey1Q.Coeffs = evakey.Value[i][1].Coeffs[:level+1]
		evakey0P.Coeffs = evakey.Value[i][0].Coeffs[len(ringQ.Modulus):]
		evakey1P.Coeffs = evakey.Value[i][1].Coeffs[len(ringQ.Modulus):]

		if i == 0 {
			ringQ.MulCoeffsMontgomeryConstantLvl(level, evakey0Q, c2QiQ, pool2Q)
			ringQ.MulCoeffsMontgomeryConstantLvl(level, evakey1Q, c2QiQ, pool3Q)
			ringP.MulCoeffsMontgomeryConstant(evakey0P, c2QiP, pool2P)
			ringP.MulCoeffsMontgomeryConstant(evakey1P, c2QiP, pool3P)
		} else {
			ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(level, evakey0Q, c2QiQ, pool2Q)
			ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(level, evakey1Q, c2QiQ, pool3Q)
			ringP.MulCoeffsMontgomeryConstantAndAddNoMod(evakey0P, c2QiP, pool2P)
			ringP.MulCoeffsMontgomeryConstantAndAddNoMod(evakey1P, c2QiP, pool3P)
		}

		if reduce%QiOverF == QiOverF-1 {
			ringQ.ReduceLvl(level, pool2Q, pool2Q)
			ringQ.ReduceLvl(level, pool3Q, pool3Q)
		}

		if reduce%PiOverF == PiOverF-1 {
			ringP.Reduce(pool2P, pool2P)
			ringP.Reduce(pool3P, pool3P)
		}

		reduce++
	}

	if reduce%QiOverF != 0 {
		ringQ.ReduceLvl(level, pool2Q, pool2Q)
		ringQ.ReduceLvl(level, pool3Q, pool3Q)
	}

	if reduce%PiOverF != 0 {
		ringP.Reduce(pool2P, pool2P)
		ringP.Reduce(pool3P, pool3P)
	}
}

// KeyswitchHoisted applies the key-switch to the decomposed polynomial c2 mod QP (PoolDecompQ and PoolDecompP)
// and divides the result by P, reducing the basis from QP to Q.
//
// pool2 = dot(PoolDecompQ||PoolDecompP * evakey[0]) mod Q
// pool3 = dot(PoolDecompQ||PoolDecompP * evakey[1]) mod Q
func (ks *KeySwitcher) KeyswitchHoisted(level int, PoolDecompQ, PoolDecompP []*ring.Poly, evakey *SwitchingKey, pool2Q, pool3Q, pool2P, pool3P *ring.Poly) {

	ks.KeyswitchHoistedNoModDown(level, PoolDecompQ, PoolDecompP, evakey, pool2Q, pool3Q, pool2P, pool3P)

	// Computes pool2Q = pool2Q/pool2P and pool3Q = pool3Q/pool3P
	ks.Baseconverter.ModDownSplitNTTPQ(level, pool2Q, pool2P, pool2Q)
	ks.Baseconverter.ModDownSplitNTTPQ(level, pool3Q, pool3P, pool3Q)
}

// KeyswitchHoistedNoModDown applies the key-switch to the decomposed polynomial c2 mod QP (PoolDecompQ and PoolDecompP)
//
// pool2 = dot(PoolDecompQ||PoolDecompP * evakey[0]) mod QP
// pool3 = dot(PoolDecompQ||PoolDecompP * evakey[1]) mod QP
func (ks *KeySwitcher) KeyswitchHoistedNoModDown(level int, PoolDecompQ, PoolDecompP []*ring.Poly, evakey *SwitchingKey, pool2Q, pool3Q, pool2P, pool3P *ring.Poly) {

	ringQ := ks.ringQ
	ringP := ks.ringP

	alpha := len(ringP.Modulus)
	beta := int(math.Ceil(float64(level+1) / float64(alpha)))

	evakey0Q := new(ring.Poly)
	evakey1Q := new(ring.Poly)
	evakey0P := new(ring.Poly)
	evakey1P := new(ring.Poly)

	QiOverF := ks.Parameters.QiOverflowMargin(level) >> 1
	PiOverF := ks.Parameters.PiOverflowMargin() >> 1

	// Key switching with CRT decomposition for the Qi
	var reduce int
	for i := 0; i < beta; i++ {

		evakey0Q.Coeffs = evakey.Value[i][0].Coeffs[:level+1]
		evakey1Q.Coeffs = evakey.Value[i][1].Coeffs[:level+1]
		evakey0P.Coeffs = evakey.Value[i][0].Coeffs[len(ringQ.Modulus):]
		evakey1P.Coeffs = evakey.Value[i][1].Coeffs[len(ringQ.Modulus):]

		if i == 0 {
			ringQ.MulCoeffsMontgomeryConstantLvl(level, evakey0Q, PoolDecompQ[i], pool2Q)
			ringQ.MulCoeffsMontgomeryConstantLvl(level, evakey1Q, PoolDecompQ[i], pool3Q)
			ringP.MulCoeffsMontgomeryConstant(evakey0P, PoolDecompP[i], pool2P)
			ringP.MulCoeffsMontgomeryConstant(evakey1P, PoolDecompP[i], pool3P)
		} else {
			ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(level, evakey0Q, PoolDecompQ[i], pool2Q)
			ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(level, evakey1Q, PoolDecompQ[i], pool3Q)
			ringP.MulCoeffsMontgomeryConstantAndAddNoMod(evakey0P, PoolDecompP[i], pool2P)
			ringP.MulCoeffsMontgomeryConstantAndAddNoMod(evakey1P, PoolDecompP[i], pool3P)
		}

		if reduce%QiOverF == QiOverF-1 {
			ringQ.ReduceLvl(level, pool2Q, pool2Q)
			ringQ.ReduceLvl(level, pool3Q, pool3Q)
		}

		if reduce%PiOverF == PiOverF-1 {
			ringP.Reduce(pool2P, pool2P)
			ringP.Reduce(pool3P, pool3P)
		}

		reduce++
	}

	if reduce%QiOverF != 0 {
		ringQ.ReduceLvl(level, pool2Q, pool2Q)
		ringQ.ReduceLvl(level, pool3Q, pool3Q)
	}

	if reduce%PiOverF != 0 {
		ringP.Reduce(pool2P, pool2P)
		ringP.Reduce(pool3P, pool3P)
	}
}
