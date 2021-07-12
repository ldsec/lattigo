package main

import (
	"fmt"
	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
	"time"
)

func main() {
	var err error
	ParametersLiteral := ckks.ParametersLiteral{
		LogN:     12,
		LogSlots: 10,
		Scale:    1 << 30,
		Sigma:    rlwe.DefaultSigma,
		Q: []uint64{
			0xffff820001,       // 40 Q0
			0x2000000a0001,     // 45 CtS
			0x2000000e0001,     // 45 CtS
			0x10000140001,      // 40
			0xffffe80001,       // 40
			0xffffc40001,       // 40
			0x100003e0001,      // 40
			0xffffb20001,       // 40
			0x10000500001,      // 40
			0xffff940001,       // 40
			0xffff8a0001,       // 40
			0xfffffffff840001,  // 60 Sine (double angle)
			0x1000000000860001, // 60 Sine (double angle)
			0xfffffffff6a0001,  // 60 Sine
			0x1000000000980001, // 60 Sine
			0xfffffffff5a0001,  // 60 Sine
			0x1000000000b00001, // 60 Sine
			0x1000000000ce0001, // 60 Sine
			0xfffffffff2a0001,  // 60 Sine
			0x100000000060001,  // 58 Repack
		},
		P: []uint64{
			0x1fffffffffe00001, // Pi 61
		},
	}

	SlotsToCoeffsParameters := ckks.SlotsToCoeffsParameters{
		LevelStart:  2,
		BSGSRatio:   16.0,
		BitReversed: false,
		ScalingFactor: [][]float64{
			{0x2000000a0001},
			{0x2000000e0001},
		},
	}

	H := 256

	var params ckks.Parameters
	if params, err = ckks.NewParametersFromLiteral(ParametersLiteral); err != nil {
		panic(err)
	}

	ringQ := params.RingQ()
	Q := ringQ.Modulus[0]

	encoder := ckks.NewEncoder(params)
	kgen := ckks.NewKeyGenerator(params)
	sk := kgen.GenSecretKeySparse(H)
	encryptor := ckks.NewEncryptor(params, sk)
	decryptor := ckks.NewDecryptor(params, sk)

	SlotsToCoeffsMatrix := SlotsToCoeffsParameters.GenSlotsToCoeffsMatrix(&params, params.LogSlots(), 1.0, encoder)

	// Gets the rotations indexes for SlotsToCoeffs
	rotations := SlotsToCoeffsParameters.RotationsForSlotsToCoeffs(&params, params.LogSlots())

	rotKey := kgen.GenRotationKeysForRotations(rotations, true, sk)
	rlk := kgen.GenRelinearizationKey(sk)

	eval := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: rotKey})

	values := make([]complex128, params.Slots())
	for i := range values {
		values[i] = complex(float64(i+1)/float64(params.Slots()), 1+float64(i+1)/float64(params.Slots()))
	}

	plaintext := ckks.NewPlaintext(params, params.MaxLevel(), params.Scale())
	// Must encode with 2*Slots because a real vector is returned
	encoder.Encode(plaintext, values, utils.MinInt(params.LogSlots()+1, params.LogN()-1))
	ct := encryptor.EncryptNew(plaintext)

	// Homomorphic Decoding
	ct = ckks.SlotsToCoeffs(ct, nil, SlotsToCoeffsMatrix, eval)

	// Decrypt and print coefficient domain
	coeffsFloat := encoder.DecodeCoeffsPublic(decryptor.DecryptNew(ct), 0)
	valuesFloat := make([]complex128, params.Slots())
	gap := params.N() / (2 * params.Slots())
	for i, idx := 0, 0; i < params.Slots(); i, idx = i+1, idx+gap {
		valuesFloat[i] = complex(coeffsFloat[idx], coeffsFloat[idx+(params.N()>>1)])
	}
	ckks.SliceBitReverseInPlaceComplex128(valuesFloat, params.Slots())

	// RLWE -> LWE Extraction
	lweReal, lweImag := ExtractLWESamplesBitReversed(ct, params)

	// Encode the LWE samples
	lweEncoded := make([]complex128, params.Slots())
	for i := 0; i < params.Slots(); i++ {
		lweEncoded[i] = complex(float64(lweReal[i].b[0]), float64(lweImag[i].b[0]))
	}

	ptLWE := ckks.NewPlaintext(params, params.MaxLevel(), 1.0)

	encoder.EncodeNTT(ptLWE, lweEncoded, params.LogSlots())

	// Encode the secret-key
	skInvNTT := ringQ.NewPolyLvl(0)
	ring.CopyValues(sk.Value, skInvNTT)
	ringQ.InvNTTLvl(0, skInvNTT, skInvNTT)

	// Visual of some values
	fmt.Println(valuesFloat[0], valuesFloat[params.Slots()-1])
	fmt.Println(complex(DecryptLWE(ringQ, lweReal[0], params.Scale(), skInvNTT), DecryptLWE(ringQ, lweImag[0], params.Scale(), skInvNTT)),
		complex(DecryptLWE(ringQ, lweReal[params.Slots()-1], params.Scale(), skInvNTT), DecryptLWE(ringQ, lweImag[params.Slots()-1], params.Scale(), skInvNTT)))

	ringQ.InvMFormLvl(0, skInvNTT, skInvNTT)

	fmt.Println("Encrypt SK")
	skFloat := make([]complex128, params.N())

	for i, s := range skInvNTT.Coeffs[0] {
		if s >= Q>>1 {
			skFloat[i] = -complex(float64(Q-s), 0)
		} else {
			skFloat[i] = complex(float64(s), 0)
		}
	}

	ringQ.InvMFormLvl(0, skInvNTT, skInvNTT)

	ctSk := make([]*ckks.Ciphertext, params.N()/params.Slots())
	ptSk := ckks.NewPlaintext(params, params.MaxLevel(), float64(params.Q()[params.MaxLevel()]))
	for i := range ctSk {
		encoder.Encode(ptSk, skFloat[i*params.Slots():(i+1)*params.Slots()], params.LogSlots())
		ctSk[i] = encryptor.EncryptNew(ptSk)
	}

	// Compute			   skLeft 				  skRight
	//    	   __________  	   _       __________ 	    _
	//   	  |		     |	  | |	  |		     |	   | |
	//   	  |		     |    | |	  |		     |     | |
	//   n    |	 ALeft   |  x | |  +  |  ARight  |  x  | | 	= 	A x sk
	//   	  |		     |	  | |	  |		     |     | |
	//   	  |__________|	  |_|	  |__________|     |_|
	//
	//			  N/n          1          N/n 		    1

	// Constructs matrix

	fmt.Println("Encode A")
	AVectors := make([][]complex128, params.Slots())
	for i := range AVectors {
		tmp := make([]complex128, params.N())
		for j := 0; j < params.N(); j++ {
			tmp[j] = complex(float64(lweReal[i].a.Coeffs[0][j]), float64(lweImag[i].a.Coeffs[0][j]))
		}

		AVectors[i] = tmp
	}

	// Diagonalize
	ptMatDiag := make([]*ckks.PtDiagMatrix, params.N()/params.Slots())
	AMatDiag := make(map[int][]complex128)
	for w := 0; w < params.N()/params.Slots(); w++ {
		for i := 0; i < params.Slots(); i++ {
			tmp := make([]complex128, params.Slots())
			for j := 0; j < params.Slots(); j++ {
				tmp[j] = AVectors[j][(j+i)%params.Slots()+w*params.Slots()]
			}

			AMatDiag[i] = tmp
		}

		ptMatDiag[w] = encoder.EncodeDiagMatrixBSGSAtLvl(params.MaxLevel(), AMatDiag, 1.0, 16.0, params.LogSlots())
	}

	fmt.Println("GenRepackKeys")
	rotKeyRepack := kgen.GenRotationKeysForRotations(params.RotationsForDiagMatrixMult(ptMatDiag[0]), false, sk)

	evalRepack := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: rotKeyRepack})

	fmt.Println("Start Repacking")
	startTotal := time.Now()
	ctAs := evalRepack.LinearTransform(ctSk[0], ptMatDiag[0])[0]
	for i := 1; i < params.N()/params.Slots(); i++ {

		startFrac := time.Now()
		evalRepack.Add(ctAs, evalRepack.LinearTransform(ctSk[i], ptMatDiag[i])[0], ctAs)
		fmt.Printf("Matrix[%d] : %s\n", i, time.Since(startFrac))
	}
	fmt.Printf("Done : %s\n", time.Since(startTotal))

	eval.Rescale(ctAs, 1.0, ctAs)
	eval.Add(ctAs, ptLWE, ctAs)

	ctAs.Scale = params.Scale()

	/*
		v := encoder.DecodePublic(decryptor.DecryptNew(ctAs), params.LogSlots(), 0)
		v2 := encoder.DecodePublic(ptLWE, params.LogSlots(), 0)

		for i := range v{
			fmt.Println(i, v[i], lweReal[i].b[0], lweImag[i].b[0], v2[i])
		}
	*/
}

func DecryptLWE(ringQ *ring.Ring, lwe RNSLWESample, scale float64, skInvNTT *ring.Poly) float64 {

	tmp := ringQ.NewPolyLvl(0)
	ringQ.MulCoeffsMontgomeryLvl(0, lwe.a, skInvNTT, tmp)
	qi := ringQ.Modulus[0]
	tmp0 := tmp.Coeffs[0]
	tmp1 := lwe.b[0]
	for j := 0; j < ringQ.N; j++ {
		tmp1 = ring.CRed(tmp1+tmp0[j], qi)
	}

	if tmp1 >= ringQ.Modulus[0]>>1 {
		tmp1 = ringQ.Modulus[0] - tmp1
		return -float64(tmp1) / scale
	} else {
		return float64(tmp1) / scale
	}
}

type RNSLWESample struct {
	b []uint64
	a *ring.Poly
}

func ExtractLWESamplesBitReversed(ct *ckks.Ciphertext, params ckks.Parameters) (LWEReal, LWEImag []RNSLWESample) {

	ringQ := params.RingQ()

	ringQ.InvNTTLvl(ct.Level(), ct.Value[0], ct.Value[0])
	ringQ.InvNTTLvl(ct.Level(), ct.Value[1], ct.Value[1])

	LWEReal = make([]RNSLWESample, params.Slots())
	LWEImag = make([]RNSLWESample, params.Slots())

	// Copies coefficients multiplied by X^{N-1} in reverse order :
	// a_{0} -a_{N-1} -a2_{N-2} ... -a_{1}
	acc := ringQ.NewPolyLvl(ct.Level())
	for i, qi := range ringQ.Modulus[:ct.Level()+1] {
		tmp0 := acc.Coeffs[i]
		tmp1 := ct.Value[1].Coeffs[i]
		tmp0[0] = tmp1[0]
		for j := 1; j < ringQ.N; j++ {
			tmp0[j] = qi - tmp1[ringQ.N-j]
		}
	}

	pol := ct.Value[0]

	gap := params.N() / (2 * params.Slots()) // Gap between plaintext coefficient if sparse packed

	// Real values
	for i, idx := 0, 0; i < params.Slots(); i, idx = i+1, idx+gap {

		iRev := utils.BitReverse64(uint64(i), uint64(params.LogSlots()))

		LWEReal[iRev].b = make([]uint64, len(pol.Coeffs))
		for j := range pol.Coeffs {
			LWEReal[iRev].b[j] = pol.Coeffs[j][idx]
		}

		LWEReal[iRev].a = acc.CopyNew()

		// Multiplies the accumulator by X^{N/(2*slots)}
		MulBySmallMonomial(ringQ, acc, gap)
	}

	// Imaginary values
	for i, idx := 0, 0; i < params.Slots(); i, idx = i+1, idx+gap {

		iRev := utils.BitReverse64(uint64(i), uint64(params.LogSlots()))

		LWEImag[iRev].b = make([]uint64, len(pol.Coeffs))
		for j := range pol.Coeffs {
			LWEImag[iRev].b[j] = pol.Coeffs[j][idx+(params.N()>>1)]
		}

		LWEImag[iRev].a = acc.CopyNew()

		// Multiplies the accumulator by X^{N/(2*slots)}
		MulBySmallMonomial(ringQ, acc, gap)
	}
	return
}

func MulBySmallMonomial(ringQ *ring.Ring, pol *ring.Poly, n int) {
	for i, qi := range ringQ.Modulus[:pol.Level()+1] {
		pol.Coeffs[i] = append(pol.Coeffs[i][ringQ.N-n:], pol.Coeffs[i][:ringQ.N-n]...)
		tmp := pol.Coeffs[i]
		for j := 0; j < n; j++ {
			tmp[j] = qi - tmp[j]
		}
	}
}
