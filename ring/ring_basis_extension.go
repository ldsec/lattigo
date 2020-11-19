package ring

import (
	"math"
	"math/big"
	"math/bits"
	"unsafe"
)

// FastBasisExtender stores the necessary parameters for RNS basis extension.
// The used algorithm is from https://eprint.iacr.org/2018/117.pdf.
type FastBasisExtender struct {
	RingQ           *Ring
	RingP           *Ring
	ParamsQP        *ModUpParams
	ParamsPQ        *ModUpParams
	ModDownParamsPQ []uint64
	ModDownParamsQP []uint64
	PolyPoolQ       *Poly
	PolyPoolP       *Poly
}

// ModUpParams stores the necessary parameters to apply ModUp Q->P or P->Q on a polynomial
type ModUpParams struct {
	Q []uint64
	P []uint64

	//Parameters for basis extension from Q to P
	// (Q/Qi)^-1) (mod each Qi) (in Montgomery form)
	qibMont []uint64
	// Q/qi (mod each Pj) (in Montgomery form)
	qispjMont [][]uint64
	// Q*v (mod each Pj) for v in [1,...,k] where k is the number of Pj moduli
	qpjInv [][]uint64

	bredParamsQ [][]uint64
	mredParamsQ []uint64

	bredParamsP [][]uint64
	mredParamsP []uint64
}

func genModDownParams(ringP, ringQ *Ring) (params []uint64) {

	params = make([]uint64, len(ringP.Modulus))

	bredParams := ringP.GetBredParams()
	tmp := new(big.Int)
	for i, Qi := range ringP.Modulus {

		params[i] = tmp.Mod(ringQ.ModulusBigint, NewUint(Qi)).Uint64()
		params[i] = ModExp(params[i], Qi-2, Qi)
		params[i] = MForm(params[i], Qi, bredParams[i])
	}

	return
}

// NewFastBasisExtender creates a new FastBasisExtender, enabling RNS basis extension from Q to P and P to Q.
func NewFastBasisExtender(ringQ, ringP *Ring) *FastBasisExtender {

	newParams := new(FastBasisExtender)

	newParams.RingQ = ringQ
	newParams.RingP = ringP

	newParams.ParamsQP = basisextenderparameters(ringQ.Modulus, ringP.Modulus)
	newParams.ParamsPQ = basisextenderparameters(ringP.Modulus, ringQ.Modulus)

	newParams.ModDownParamsPQ = genModDownParams(ringQ, ringP)
	newParams.ModDownParamsQP = genModDownParams(ringP, ringQ)

	newParams.PolyPoolQ = ringQ.NewPoly()
	newParams.PolyPoolP = ringP.NewPoly()

	return newParams
}

func basisextenderparameters(Q, P []uint64) (params *ModUpParams) {

	params = new(ModUpParams)

	params.Q = make([]uint64, len(Q))
	params.bredParamsQ = make([][]uint64, len(Q))
	params.mredParamsQ = make([]uint64, len(Q))
	for i, qi := range Q {
		params.Q[i] = Q[i]
		params.bredParamsQ[i] = BRedParams(qi)
		params.mredParamsQ[i] = MRedParams(qi)
	}

	params.P = make([]uint64, len(P))
	params.bredParamsP = make([][]uint64, len(P))
	params.mredParamsP = make([]uint64, len(P))
	for i, pj := range P {
		params.P[i] = P[i]
		params.bredParamsP[i] = BRedParams(pj)
		params.mredParamsP[i] = MRedParams(pj)
	}

	tmp := new(big.Int)
	QiB := new(big.Int)
	QiStar := new(big.Int)
	QiBarre := new(big.Int)

	modulusbigint := NewUint(1)
	for _, qi := range Q {
		modulusbigint.Mul(modulusbigint, NewUint(qi))
	}

	params.qibMont = make([]uint64, len(Q))
	params.qispjMont = make([][]uint64, len(P))

	for i := range P {
		params.qispjMont[i] = make([]uint64, len(Q))
	}

	for i, qi := range Q {

		QiB.SetUint64(qi)
		QiStar.Quo(modulusbigint, QiB)
		QiBarre.ModInverse(QiStar, QiB)
		QiBarre.Mod(QiBarre, QiB)

		// (Q/Qi)^-1) * r (mod Qi) (in Montgomery form)
		params.qibMont[i] = MForm(QiBarre.Uint64(), qi, params.bredParamsQ[i])

		for j, pj := range P {
			// (Q/qi * r) (mod Pj) (in Montgomery form)
			params.qispjMont[j][i] = MForm(tmp.Mod(QiStar, NewUint(pj)).Uint64(), pj, params.bredParamsP[j])
		}
	}

	params.qpjInv = make([][]uint64, len(P))
	for j, pj := range P {
		params.qpjInv[j] = make([]uint64, len(Q)+1)
		// Correction Term (v*Q) mod each Pj
		v := pj - tmp.Mod(modulusbigint, NewUint(pj)).Uint64()
		params.qpjInv[j][0] = 0
		for i := 1; i < len(Q)+1; i++ {
			params.qpjInv[j][i] = CRed(params.qpjInv[j][i-1]+v, pj)
		}
	}

	return
}

// ModUpSplitQP extends the RNS basis of a polynomial from Q to QP.
// Given a polynomial with coefficients in basis {Q0,Q1....Qlevel},
// it extends its basis from {Q0,Q1....Qlevel} to {Q0,Q1....Qlevel,P0,P1...Pj}
func (basisextender *FastBasisExtender) ModUpSplitQP(level uint64, p1, p2 *Poly) {
	ModUpExact(p1.Coeffs[:level+1], p2.Coeffs[:uint64(len(basisextender.ParamsQP.P))], basisextender.ParamsQP)
}

// ModUpSplitPQ extends the RNS basis of a polynomial from P to PQ.
// Given a polynomial with coefficients in basis {P0,P1....Plevel},
// it extends its basis from {P0,P1....Plevel} to {Q0,Q1...Qj}
func (basisextender *FastBasisExtender) ModUpSplitPQ(level uint64, p1, p2 *Poly) {
	ModUpExact(p1.Coeffs[:level+1], p2.Coeffs[:uint64(len(basisextender.ParamsPQ.P))], basisextender.ParamsPQ)
}

// ModDownNTTPQ reduces the basis RNS of a polynomial in the NTT domain
// from QP to Q and divides its coefficients by P.
// Given a polynomial with coefficients in basis {Q0,Q1....Qlevel,P0,P1...Pj},
// it reduces its basis from {Q0,Q1....Qlevel,P0,P1...Pj} to {Q0,Q1....Qlevel}
// and performs a rounded integer division of the result by P.
// Inputs must be in the NTT domain.
func (basisextender *FastBasisExtender) ModDownNTTPQ(level uint64, p1, p2 *Poly) {

	ringQ := basisextender.RingQ
	ringP := basisextender.RingP
	modDownParams := basisextender.ModDownParamsPQ
	polypool := basisextender.PolyPoolQ
	nQi := len(ringQ.Modulus)
	nPj := len(ringP.Modulus)

	// First we get the P basis part of p1 out of the NTT domain
	for j := 0; j < nPj; j++ {
		InvNTTLazy(p1.Coeffs[nQi+j], p1.Coeffs[nQi+j], ringP.N, ringP.GetNttPsiInv()[j], ringP.GetNttNInv()[j], ringP.Modulus[j], ringP.GetMredParams()[j])
	}

	// Then we target this P basis of p1 and convert it to a Q basis (at the "level" of p1) and copy it on polypool
	// polypool is now the representation of the P basis of p1 but in basis Q (at the "level" of p1)
	ModUpExact(p1.Coeffs[nQi:nQi+nPj], polypool.Coeffs[:level+1], basisextender.ParamsPQ)

	// Finally, for each level of p1 (and polypool since they now share the same basis) we compute p2 = (P^-1) * (p1 - polypool) mod Q
	for i := uint64(0); i < level+1; i++ {

		qi := ringQ.Modulus[i]
		twoqi := qi << 1
		p1tmp := p1.Coeffs[i]
		p2tmp := p2.Coeffs[i]
		p3tmp := polypool.Coeffs[i]
		params := qi - modDownParams[i]
		mredParams := ringQ.MredParams[i]
		bredParams := ringQ.BredParams[i]

		// First we switch back the relevant polypool CRT array back to the NTT domain
		NTTLazy(p3tmp, p3tmp, ringQ.N, ringQ.GetNttPsi()[i], qi, mredParams, bredParams)

		// Then for each coefficient we compute (P^-1) * (p1[i][j] - polypool[i][j]) mod qi
		for j := uint64(0); j < ringQ.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = MRed(y[0]+twoqi-x[0], params, qi, mredParams)
			z[1] = MRed(y[1]+twoqi-x[1], params, qi, mredParams)
			z[2] = MRed(y[2]+twoqi-x[2], params, qi, mredParams)
			z[3] = MRed(y[3]+twoqi-x[3], params, qi, mredParams)
			z[4] = MRed(y[4]+twoqi-x[4], params, qi, mredParams)
			z[5] = MRed(y[5]+twoqi-x[5], params, qi, mredParams)
			z[6] = MRed(y[6]+twoqi-x[6], params, qi, mredParams)
			z[7] = MRed(y[7]+twoqi-x[7], params, qi, mredParams)
		}
	}

	// In total we do len(P) + len(Q) NTT, which is optimal (linear in the number of moduli of P and Q)
}

// ModDownSplitNTTPQ reduces the basis of a polynomial.
// Given a polynomial with coefficients in basis {Q0,Q1....Qi} and {P0,P1...Pj},
// it reduces its basis from {Q0,Q1....Qi} and {P0,P1...Pj} to {Q0,Q1....Qi}
// and does a rounded integer division of the result by P.
// Inputs must be in the NTT domain.
func (basisextender *FastBasisExtender) ModDownSplitNTTPQ(level uint64, p1Q, p1P, p2 *Poly) {

	ringQ := basisextender.RingQ
	ringP := basisextender.RingP
	modDownParams := basisextender.ModDownParamsPQ
	polypool := basisextender.PolyPoolQ

	// First we get the P basis part of p1 out of the NTT domain
	ringP.InvNTTLazy(p1P, p1P)

	// Then we target this P basis of p1 and convert it to a Q basis (at the "level" of p1) and copy it on polypool
	// polypool is now the representation of the P basis of p1 but in basis Q (at the "level" of p1)
	ModUpExact(p1P.Coeffs, polypool.Coeffs[:level+1], basisextender.ParamsPQ)

	// Finally, for each level of p1 (and polypool since they now share the same basis) we compute p2 = (P^-1) * (p1 - polypool) mod Q
	for i := uint64(0); i < level+1; i++ {

		qi := ringQ.Modulus[i]
		twoqi := qi << 1
		p1tmp := p1Q.Coeffs[i]
		p2tmp := p2.Coeffs[i]
		p3tmp := polypool.Coeffs[i]
		params := qi - modDownParams[i]
		mredParams := ringQ.MredParams[i]
		bredParams := ringQ.BredParams[i]

		// First we switch back the relevant polypool CRT array back to the NTT domain
		NTTLazy(p3tmp, p3tmp, ringQ.N, ringQ.GetNttPsi()[i], ringQ.Modulus[i], mredParams, bredParams)

		// Then for each coefficient we compute (P^-1) * (p1[i][j] - polypool[i][j]) mod qi
		for j := uint64(0); j < ringQ.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = MRed(y[0]+twoqi-x[0], params, qi, mredParams)
			z[1] = MRed(y[1]+twoqi-x[1], params, qi, mredParams)
			z[2] = MRed(y[2]+twoqi-x[2], params, qi, mredParams)
			z[3] = MRed(y[3]+twoqi-x[3], params, qi, mredParams)
			z[4] = MRed(y[4]+twoqi-x[4], params, qi, mredParams)
			z[5] = MRed(y[5]+twoqi-x[5], params, qi, mredParams)
			z[6] = MRed(y[6]+twoqi-x[6], params, qi, mredParams)
			z[7] = MRed(y[7]+twoqi-x[7], params, qi, mredParams)
		}
	}

	// In total we do len(P) + len(Q) NTT, which is optimal (linear in the number of moduli of P and Q)
}

// ModDownPQ reduces the basis of a polynomial.
// Given a polynomial with coefficients in basis {Q0,Q1....Qlevel,P0,P1...Pj},
// it reduces its basis from {Q0,Q1....Qlevel,P0,P1...Pj} to {Q0,Q1....Qlevel}
// and does a rounded integer division of the result by P.
func (basisextender *FastBasisExtender) ModDownPQ(level uint64, p1, p2 *Poly) {

	ringQ := basisextender.RingQ
	modDownParams := basisextender.ModDownParamsPQ
	polypool := basisextender.PolyPoolQ
	nPi := uint64(len(basisextender.ParamsQP.P))

	// We target this P basis of p1 and convert it to a Q basis (at the "level" of p1) and copy it on polypool
	// polypool is now the representation of the P basis of p1 but in basis Q (at the "level" of p1)
	ModUpExact(p1.Coeffs[level+1:level+1+nPi], polypool.Coeffs[:level+1], basisextender.ParamsPQ)

	// Finally, for each level of p1 (and polypool since they now share the same basis) we compute p2 = (P^-1) * (p1 - polypool) mod Q
	for i := uint64(0); i < level+1; i++ {

		qi := ringQ.Modulus[i]
		twoqi := qi << 1
		p1tmp := p1.Coeffs[i]
		p2tmp := p2.Coeffs[i]
		p3tmp := polypool.Coeffs[i]
		params := qi - modDownParams[i]
		mredParams := ringQ.MredParams[i]

		// Then for each coefficient we compute (P^-1) * (p1[i][j] - polypool[i][j]) mod qi
		for j := uint64(0); j < ringQ.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = MRed(y[0]+twoqi-x[0], params, qi, mredParams)
			z[1] = MRed(y[1]+twoqi-x[1], params, qi, mredParams)
			z[2] = MRed(y[2]+twoqi-x[2], params, qi, mredParams)
			z[3] = MRed(y[3]+twoqi-x[3], params, qi, mredParams)
			z[4] = MRed(y[4]+twoqi-x[4], params, qi, mredParams)
			z[5] = MRed(y[5]+twoqi-x[5], params, qi, mredParams)
			z[6] = MRed(y[6]+twoqi-x[6], params, qi, mredParams)
			z[7] = MRed(y[7]+twoqi-x[7], params, qi, mredParams)
		}
	}

	// In total we do len(P) + len(Q) NTT, which is optimal (linear in the number of moduli of P and Q)
}

// ModDownSplitPQ reduces the basis of a polynomial.
// Given a polynomial with coefficients in basis {Q0,Q1....Qlevel} and {P0,P1...Pj},
// it reduces its basis from {Q0,Q1....Qlevel} and {P0,P1...Pj} to {Q0,Q1....Qlevel}
// and does a rounded integer division of the result by P.
func (basisextender *FastBasisExtender) ModDownSplitPQ(level uint64, p1Q, p1P, p2 *Poly) {

	ringQ := basisextender.RingQ
	modDownParams := basisextender.ModDownParamsPQ
	polypool := basisextender.PolyPoolQ

	// Then we target this P basis of p1 and convert it to a Q basis (at the "level" of p1) and copy it on polypool
	// polypool is now the representation of the P basis of p1 but in basis Q (at the "level" of p1)
	ModUpExact(p1P.Coeffs, polypool.Coeffs[:level+1], basisextender.ParamsPQ)

	// Finally, for each level of p1 (and polypool since they now share the same basis) we compute p2 = (P^-1) * (p1 - polypool) mod Q
	for i := uint64(0); i < level+1; i++ {

		qi := ringQ.Modulus[i]
		twoqi := qi << 1
		p1tmp := p1Q.Coeffs[i]
		p2tmp := p2.Coeffs[i]
		p3tmp := polypool.Coeffs[i]
		params := qi - modDownParams[i]
		mredParams := ringQ.MredParams[i]

		// Then for each coefficient we compute (P^-1) * (p1[i][j] - polypool[i][j]) mod qi
		for j := uint64(0); j < ringQ.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = MRed(y[0]+twoqi-x[0], params, qi, mredParams)
			z[1] = MRed(y[1]+twoqi-x[1], params, qi, mredParams)
			z[2] = MRed(y[2]+twoqi-x[2], params, qi, mredParams)
			z[3] = MRed(y[3]+twoqi-x[3], params, qi, mredParams)
			z[4] = MRed(y[4]+twoqi-x[4], params, qi, mredParams)
			z[5] = MRed(y[5]+twoqi-x[5], params, qi, mredParams)
			z[6] = MRed(y[6]+twoqi-x[6], params, qi, mredParams)
			z[7] = MRed(y[7]+twoqi-x[7], params, qi, mredParams)
		}
	}

	// In total we do len(P) + len(Q) NTT, which is optimal (linear in the number of moduli of P and Q)
}

// ModDownSplitQP reduces the basis of a polynomial.
// Given a polynomial with coefficients in basis {Q0,Q1....QlevelQ} and {P0,P1...PlevelP},
// it reduces its basis from {Q0,Q1....QlevelQ} and {P0,P1...PlevelP} to {P0,P1...PlevelP}
// and does a floored integer division of the result by Q.
func (basisextender *FastBasisExtender) ModDownSplitQP(levelQ, levelP uint64, p1Q, p1P, p2 *Poly) {

	ringP := basisextender.RingP
	modDownParams := basisextender.ModDownParamsQP
	polypool := basisextender.PolyPoolP

	// Then we target this P basis of p1 and convert it to a Q basis (at the "level" of p1) and copy it on polypool
	// polypool is now the representation of the P basis of p1 but in basis Q (at the "level" of p1)
	basisextender.ModUpSplitQP(levelQ, p1Q, polypool)

	// Finally, for each level of p1 (and polypool since they now share the same basis) we compute p2 = (P^-1) * (p1 - polypool) mod Q
	for i := uint64(0); i < levelP+1; i++ {

		qi := ringP.Modulus[i]
		twoqi := qi << 1
		p1tmp := p1P.Coeffs[i]
		p2tmp := p2.Coeffs[i]
		p3tmp := polypool.Coeffs[i]
		params := qi - modDownParams[i]
		mredParams := ringP.MredParams[i]

		// Then for each coefficient we compute (P^-1) * (p1[i][j] - polypool[i][j]) mod qi
		for j := uint64(0); j < ringP.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			y := (*[8]uint64)(unsafe.Pointer(&p3tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = MRed(y[0]+twoqi-x[0], params, qi, mredParams)
			z[1] = MRed(y[1]+twoqi-x[1], params, qi, mredParams)
			z[2] = MRed(y[2]+twoqi-x[2], params, qi, mredParams)
			z[3] = MRed(y[3]+twoqi-x[3], params, qi, mredParams)
			z[4] = MRed(y[4]+twoqi-x[4], params, qi, mredParams)
			z[5] = MRed(y[5]+twoqi-x[5], params, qi, mredParams)
			z[6] = MRed(y[6]+twoqi-x[6], params, qi, mredParams)
			z[7] = MRed(y[7]+twoqi-x[7], params, qi, mredParams)
		}
	}

	// In total we do len(P) + len(Q) NTT, which is optimal (linear in the number of moduli of P and Q)
}

// ModUpExact applies a basis extension from Q to P to the polynomial
// and approximates the number of overflow of Q using floating point arithmetic.
// Caution, returns the values in [0, 2q-1]
func ModUpExact(p1, p2 [][]uint64, params *ModUpParams) {

	var v [8]uint64
	var y0, y1, y2, y3, y4, y5, y6, y7 [32]uint64

	// We loop over each coefficient and apply the basis extension
	for x := uint64(0); x < uint64(len(p1[0])); x = x + 8 {

		reconstructRNS(uint64(len(p1)), x, p1, &v, &y0, &y1, &y2, &y3, &y4, &y5, &y6, &y7, params.Q, params.mredParamsQ, params.qibMont)

		for j := 0; j < len(p2); j++ {

			pj := params.P[j]
			qInv := params.mredParamsP[j]
			qpjInv := params.qpjInv[j]
			qispjMont := params.qispjMont[j]

			res := (*[8]uint64)(unsafe.Pointer(&p2[j][x]))

			multSum(res, &v, &y0, &y1, &y2, &y3, &y4, &y5, &y6, &y7, uint64(len(p1)), pj, qInv, qpjInv, qispjMont)
		}
	}
}

// Decomposer is a structure that stores the parameters of the arbitrary decomposer.
// This decomposer takes a p(x)_Q (in basis Q) and returns p(x) mod qi in basis QP, where
// qi = prod(Q_i) for 0<=i<=L, where L is the number of factors in P.
type Decomposer struct {
	nQprimes    uint64
	nPprimes    uint64
	alpha       uint64
	beta        uint64
	xalpha      []uint64
	ModUpParams [][]*ModUpParams
	QInt        *big.Int
	PInt        *big.Int
}

// Xalpha returns a slice that contains all the values of #Qi/#Pi.
func (decomposer *Decomposer) Xalpha() (xalpha []uint64) {
	return decomposer.xalpha
}

// NewDecomposer creates a new Decomposer.
func NewDecomposer(Q, P []uint64) (decomposer *Decomposer) {
	decomposer = new(Decomposer)

	decomposer.nQprimes = uint64(len(Q))
	decomposer.nPprimes = uint64(len(P))

	decomposer.QInt = NewUint(1)
	for i := range Q {
		decomposer.QInt.Mul(decomposer.QInt, NewUint(Q[i]))
	}

	decomposer.PInt = NewUint(1)
	for i := range P {
		decomposer.PInt.Mul(decomposer.PInt, NewUint(P[i]))
	}

	decomposer.alpha = uint64(len(P))
	decomposer.beta = uint64(math.Ceil(float64(len(Q)) / float64(decomposer.alpha)))

	decomposer.xalpha = make([]uint64, decomposer.beta)
	for i := range decomposer.xalpha {
		decomposer.xalpha[i] = decomposer.alpha
	}

	if uint64(len(Q))%decomposer.alpha != 0 {
		decomposer.xalpha[decomposer.beta-1] = uint64(len(Q)) % decomposer.alpha
	}

	decomposer.ModUpParams = make([][]*ModUpParams, decomposer.beta)

	// Create a basis extension for each possible combination of [Qi,Pj] according to xalpha
	for i := uint64(0); i < decomposer.beta; i++ {

		decomposer.ModUpParams[i] = make([]*ModUpParams, decomposer.xalpha[i]-1)

		for j := uint64(0); j < decomposer.xalpha[i]-1; j++ {

			Qi := make([]uint64, j+2)
			Pi := make([]uint64, len(Q)+len(P))

			for k := uint64(0); k < j+2; k++ {
				Qi[k] = Q[i*decomposer.alpha+k]
			}

			for k := 0; k < len(Q); k++ {
				Pi[k] = Q[k]
			}

			for k := len(Q); k < len(Q)+len(P); k++ {
				Pi[k] = P[k-len(Q)]
			}

			decomposer.ModUpParams[i][j] = basisextenderparameters(Qi, Pi)
		}
	}

	return
}

// DecomposeAndSplit decomposes a polynomial p(x) in basis Q, reduces it modulo qi, and returns
// the result in basis QP separately.
func (decomposer *Decomposer) DecomposeAndSplit(level, crtDecompLevel uint64, p0, p1Q, p1P *Poly) {

	alphai := decomposer.xalpha[crtDecompLevel]

	p0idxst := crtDecompLevel * decomposer.alpha
	p0idxed := p0idxst + alphai

	// First we check if the vector can simply by coping and rearranging elements (the case where no reconstruction is needed)
	if (p0idxed > level+1 && (level+1)%decomposer.nPprimes == 1) || alphai == 1 {

		for j := uint64(0); j < level+1; j++ {
			copy(p1Q.Coeffs[j], p0.Coeffs[p0idxst])
		}

		for j := uint64(0); j < decomposer.nPprimes; j++ {
			copy(p1P.Coeffs[j], p0.Coeffs[p0idxst])
		}

		// Otherwise, we apply a fast exact base conversion for the reconstruction
	} else {

		var index uint64
		if level >= alphai+crtDecompLevel*decomposer.alpha {
			index = decomposer.xalpha[crtDecompLevel] - 2
		} else {
			index = (level - 1) % decomposer.alpha
		}

		params := decomposer.ModUpParams[crtDecompLevel][index]

		var v [8]uint64
		var vi [8]float64
		var y0, y1, y2, y3, y4, y5, y6, y7 [32]uint64
		var qibMont, qi, pj, mredParams uint64
		var qif float64

		// We loop over each coefficient and apply the basis extension
		for x := uint64(0); x < uint64(len(p0.Coeffs[0])); x = x + 8 {

			vi[0], vi[1], vi[2], vi[3], vi[4], vi[5], vi[6], vi[7] = 0, 0, 0, 0, 0, 0, 0, 0

			// Coefficients to be decomposed
			for i, j := uint64(0), p0idxst; i < index+2; i, j = i+1, j+1 {

				qibMont = params.qibMont[i]
				qi = params.Q[i]
				mredParams = params.mredParamsQ[i]
				qif = float64(qi)

				px := (*[8]uint64)(unsafe.Pointer(&p0.Coeffs[j][x]))
				py := (*[8]uint64)(unsafe.Pointer(&p1Q.Coeffs[j][x]))

				// For the coefficients to be decomposed, we can simply copy them
				py[0], py[1], py[2], py[3], py[4], py[5], py[6], py[7] = px[0], px[1], px[2], px[3], px[4], px[5], px[6], px[7]

				y0[i] = MRed(px[0], qibMont, qi, mredParams)
				y1[i] = MRed(px[1], qibMont, qi, mredParams)
				y2[i] = MRed(px[2], qibMont, qi, mredParams)
				y3[i] = MRed(px[3], qibMont, qi, mredParams)
				y4[i] = MRed(px[4], qibMont, qi, mredParams)
				y5[i] = MRed(px[5], qibMont, qi, mredParams)
				y6[i] = MRed(px[6], qibMont, qi, mredParams)
				y7[i] = MRed(px[7], qibMont, qi, mredParams)

				// Computation of the correction term v * Q%pi
				vi[0] += float64(y0[i]) / qif
				vi[1] += float64(y1[i]) / qif
				vi[2] += float64(y2[i]) / qif
				vi[3] += float64(y3[i]) / qif
				vi[4] += float64(y4[i]) / qif
				vi[5] += float64(y5[i]) / qif
				vi[6] += float64(y6[i]) / qif
				vi[7] += float64(y7[i]) / qif
			}

			// Index of the correction term
			v[0] = uint64(vi[0])
			v[1] = uint64(vi[1])
			v[2] = uint64(vi[2])
			v[3] = uint64(vi[3])
			v[4] = uint64(vi[4])
			v[5] = uint64(vi[5])
			v[6] = uint64(vi[6])
			v[7] = uint64(vi[7])

			// Coefficients of index smaller than the ones to be decomposed
			for j := uint64(0); j < p0idxst; j++ {

				pj = params.P[j]
				qInv := params.mredParamsP[j]
				qpjInv := params.qpjInv[j]
				qispjMont := params.qispjMont[j]

				res := (*[8]uint64)(unsafe.Pointer(&p1Q.Coeffs[j][x]))

				multSum(res, &v, &y0, &y1, &y2, &y3, &y4, &y5, &y6, &y7, index+2, pj, qInv, qpjInv, qispjMont)
			}

			// Coefficients of index greater than the ones to be decomposed
			for j := decomposer.alpha * crtDecompLevel; j < level+1; j = j + 1 {

				pj = params.P[j]
				qInv := params.mredParamsP[j]
				qpjInv := params.qpjInv[j]
				qispjMont := params.qispjMont[j]

				res := (*[8]uint64)(unsafe.Pointer(&p1Q.Coeffs[j][x]))

				multSum(res, &v, &y0, &y1, &y2, &y3, &y4, &y5, &y6, &y7, index+2, pj, qInv, qpjInv, qispjMont)
			}

			// Coefficients of the special primes Pi
			for j, u := uint64(0), decomposer.nQprimes; j < decomposer.nPprimes; j, u = j+1, u+1 {

				pj = params.P[u]
				qInv := params.mredParamsP[u]
				qpjInv := params.qpjInv[u]
				qispjMont := params.qispjMont[u]

				res := (*[8]uint64)(unsafe.Pointer(&p1P.Coeffs[j][x]))

				multSum(res, &v, &y0, &y1, &y2, &y3, &y4, &y5, &y6, &y7, index+2, pj, qInv, qpjInv, qispjMont)
			}
		}
	}
}

func reconstructRNS(index, x uint64, p [][]uint64, v *[8]uint64, y0, y1, y2, y3, y4, y5, y6, y7 *[32]uint64, Q, QInv, QbMont []uint64) {

	var vi [8]float64
	var qi, qiInv, qibMont uint64
	var qif float64

	for i := uint64(0); i < index; i++ {

		qibMont = QbMont[i]
		qi = Q[i]
		qiInv = QInv[i]
		qif = float64(qi)

		y0[i] = MRed(p[i][x+0], qibMont, qi, qiInv)
		y1[i] = MRed(p[i][x+1], qibMont, qi, qiInv)
		y2[i] = MRed(p[i][x+2], qibMont, qi, qiInv)
		y3[i] = MRed(p[i][x+3], qibMont, qi, qiInv)
		y4[i] = MRed(p[i][x+4], qibMont, qi, qiInv)
		y5[i] = MRed(p[i][x+5], qibMont, qi, qiInv)
		y6[i] = MRed(p[i][x+6], qibMont, qi, qiInv)
		y7[i] = MRed(p[i][x+7], qibMont, qi, qiInv)

		// Computation of the correction term v * Q%pi
		vi[0] += float64(y0[i]) / qif
		vi[1] += float64(y1[i]) / qif
		vi[2] += float64(y2[i]) / qif
		vi[3] += float64(y3[i]) / qif
		vi[4] += float64(y4[i]) / qif
		vi[5] += float64(y5[i]) / qif
		vi[6] += float64(y6[i]) / qif
		vi[7] += float64(y7[i]) / qif
	}

	v[0] = uint64(vi[0])
	v[1] = uint64(vi[1])
	v[2] = uint64(vi[2])
	v[3] = uint64(vi[3])
	v[4] = uint64(vi[4])
	v[5] = uint64(vi[5])
	v[6] = uint64(vi[6])
	v[7] = uint64(vi[7])
}

// Caution, returns the values in [0, 2q-1]
func multSum(res, v *[8]uint64, y0, y1, y2, y3, y4, y5, y6, y7 *[32]uint64, index, pj, qInv uint64, qpjInv, qispjMont []uint64) {

	var rlo, rhi [8]uint64
	var mhi, mlo, c, hhi uint64

	// Accumulates the sum on uint128 and does a lazy montgomery reduction at the end
	for i := uint64(0); i < index; i++ {

		mhi, mlo = bits.Mul64(y0[i], qispjMont[i])
		rlo[0], c = bits.Add64(rlo[0], mlo, 0)
		rhi[0] += mhi + c

		mhi, mlo = bits.Mul64(y1[i], qispjMont[i])
		rlo[1], c = bits.Add64(rlo[1], mlo, 0)
		rhi[1] += mhi + c

		mhi, mlo = bits.Mul64(y2[i], qispjMont[i])
		rlo[2], c = bits.Add64(rlo[2], mlo, 0)
		rhi[2] += mhi + c

		mhi, mlo = bits.Mul64(y3[i], qispjMont[i])
		rlo[3], c = bits.Add64(rlo[3], mlo, 0)
		rhi[3] += mhi + c

		mhi, mlo = bits.Mul64(y4[i], qispjMont[i])
		rlo[4], c = bits.Add64(rlo[4], mlo, 0)
		rhi[4] += mhi + c

		mhi, mlo = bits.Mul64(y5[i], qispjMont[i])
		rlo[5], c = bits.Add64(rlo[5], mlo, 0)
		rhi[5] += mhi + c

		mhi, mlo = bits.Mul64(y6[i], qispjMont[i])
		rlo[6], c = bits.Add64(rlo[6], mlo, 0)
		rhi[6] += mhi + c

		mhi, mlo = bits.Mul64(y7[i], qispjMont[i])
		rlo[7], c = bits.Add64(rlo[7], mlo, 0)
		rhi[7] += mhi + c
	}

	hhi, _ = bits.Mul64(rlo[0]*qInv, pj)
	res[0] = rhi[0] - hhi + pj + qpjInv[v[0]]

	hhi, _ = bits.Mul64(rlo[1]*qInv, pj)
	res[1] = rhi[1] - hhi + pj + qpjInv[v[1]]

	hhi, _ = bits.Mul64(rlo[2]*qInv, pj)
	res[2] = rhi[2] - hhi + pj + qpjInv[v[2]]

	hhi, _ = bits.Mul64(rlo[3]*qInv, pj)
	res[3] = rhi[3] - hhi + pj + qpjInv[v[3]]

	hhi, _ = bits.Mul64(rlo[4]*qInv, pj)
	res[4] = rhi[4] - hhi + pj + qpjInv[v[4]]

	hhi, _ = bits.Mul64(rlo[5]*qInv, pj)
	res[5] = rhi[5] - hhi + pj + qpjInv[v[5]]

	hhi, _ = bits.Mul64(rlo[6]*qInv, pj)
	res[6] = rhi[6] - hhi + pj + qpjInv[v[6]]

	hhi, _ = bits.Mul64(rlo[7]*qInv, pj)
	res[7] = rhi[7] - hhi + pj + qpjInv[v[7]]
}
