package ckks

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// RotateHoisted takes an input Ciphertext and a list of rotations and returns a map of Ciphertext, where each element of the map is the input Ciphertext
// rotation by one element of the list. It is much faster than sequential calls to Rotate.
func (eval *evaluator) RotateHoisted(ctIn *Ciphertext, rotations []int) (cOut map[int]*Ciphertext) {

	level := ctIn.Level()

	eval.DecomposeNTT(level, ctIn.Value[1], eval.PoolDecompQ, eval.PoolDecompP)

	cOut = make(map[int]*Ciphertext)
	for _, i := range rotations {

		if i == 0 {
			cOut[i] = ctIn.CopyNew()
		} else {
			cOut[i] = NewCiphertext(eval.params, 1, level, ctIn.Scale)
			eval.permuteNTTHoisted(level, ctIn.Value[0], ctIn.Value[1], eval.PoolDecompQ, eval.PoolDecompP, i, cOut[i].Value[0], cOut[i].Value[1])
		}
	}

	return
}

// LinearTransform evaluates a linear transform on the ciphertext. The linearTransform can either be an (ordered) list of
// PtDiagMatrix or a single PtDiagMatrix. In either case a list of ciphertext is return (the second case returnign a list of
// containing a single ciphertext. A PtDiagMatrix is a diagonalized plaintext matrix contructed with an Encoder using
// the method encoder.EncodeDiagMatrixAtLvl(*).
func (eval *evaluator) LinearTransform(ctIn *Ciphertext, linearTransform interface{}) (ctOut []*Ciphertext) {

	switch element := linearTransform.(type) {
	case []*PtDiagMatrix:
		ctOut = make([]*Ciphertext, len(element))

		var maxLevel int
		for _, matrix := range element {
			maxLevel = utils.MaxInt(maxLevel, matrix.Level)
		}

		minLevel := utils.MinInt(maxLevel, ctIn.Level())

		eval.DecomposeNTT(minLevel, ctIn.Value[1], eval.PoolDecompQ, eval.PoolDecompP)

		for i, matrix := range element {
			ctOut[i] = NewCiphertext(eval.params, 1, minLevel, ctIn.Scale)

			if matrix.naive {
				eval.MultiplyByDiagMatrix(ctIn, matrix, eval.PoolDecompQ, eval.PoolDecompP, ctOut[i])
			} else {
				eval.MultiplyByDiagMatrixBSGS(ctIn, matrix, eval.PoolDecompQ, eval.PoolDecompP, ctOut[i])
			}
		}

	case *PtDiagMatrix:

		minLevel := utils.MinInt(element.Level, ctIn.Level())
		eval.DecomposeNTT(minLevel, ctIn.Value[1], eval.PoolDecompQ, eval.PoolDecompP)

		ctOut = []*Ciphertext{NewCiphertext(eval.params, 1, minLevel, ctIn.Scale)}

		if element.naive {
			eval.MultiplyByDiagMatrix(ctIn, element, eval.PoolDecompQ, eval.PoolDecompP, ctOut[0])
		} else {
			eval.MultiplyByDiagMatrixBSGS(ctIn, element, eval.PoolDecompQ, eval.PoolDecompP, ctOut[0])
		}
	}

	return
}

// InnerSumLog applies an optimized inner sum on the ciphetext (log2(n) + HW(n) rotations with double hoisting).
// The operation assumes that `ctIn` encrypts SlotCount/`batchSize` sub-vectors of size `batchSize` which it adds together (in parallel) by groups of `n`.
// It outputs in ctOut a ciphertext for which the "leftmost" sub-vector of each group is equal to the sum of the group.
// This method is faster than InnerSum when the number of rotations is large and uses log2(n) + HW(n) insteadn of 'n' keys.
func (eval *evaluator) InnerSumLog(ctIn *Ciphertext, batchSize, n int, ctOut *Ciphertext) {

	ringQ := eval.ringQ
	ringP := eval.ringP

	levelQ := ctIn.Level()

	//QiOverF := eval.params.QiOverflowMargin(levelQ)
	//PiOverF := eval.params.PiOverflowMargin()

	if n == 1 {
		if ctIn != ctOut {
			ring.CopyValuesLvl(levelQ, ctIn.Value[0], ctOut.Value[0])
			ring.CopyValuesLvl(levelQ, ctIn.Value[1], ctOut.Value[1])
		}
	} else {

		// Memory pool for ctIn = ctIn + rot(ctIn, 2^i) in Q
		tmpc0 := eval.poolQMul[0] // unused memory pool from evaluator
		tmpc1 := eval.poolQMul[1] // unused memory pool from evaluator
		tmpc2 := eval.poolQMul[2]

		// Accumulator outer loop for ctOut = ctOut + rot(ctIn, k) in QP
		ct0OutQ := eval.PoolQ[4]
		ct1OutQ := eval.PoolQ[5]
		ct0OutP := eval.PoolP[4]
		ct1OutP := eval.PoolP[5]

		// Memory pool for rot(ctIn, k)
		pool2Q := eval.PoolQ[2] // ctOut(c0', c1') from evaluator keyswitch memory pool
		pool3Q := eval.PoolQ[3] // ctOut(c0', c1') from evaluator keyswitch memory pool
		pool2P := eval.PoolP[2] // ctOut(c0', c1') from evaluator keyswitch memory pool
		pool3P := eval.PoolP[3] // ctOut(c0', c1') from evaluator keyswitch memory pool

		// Used by the key-switch
		// eval.poolQ[0]
		// eval.poolQ[1]
		// eval.poolP[0]
		// eval.poolP[1]

		state := false
		copy := true
		// Binary reading of the input n
		for i, j := 0, n; j > 0; i, j = i+1, j>>1 {

			// Starts by decomposing the input ciphertext
			if i == 0 {
				// If first iteration, then copies directly from the input ciphertext that hasn't been rotated
				eval.DecomposeNTT(levelQ, ctIn.Value[1], eval.PoolDecompQ, eval.PoolDecompP)
			} else {
				// Else copies from the rotated input ciphertext
				tmpc1.IsNTT = true
				eval.DecomposeNTT(levelQ, tmpc1, eval.PoolDecompQ, eval.PoolDecompP)
			}

			// If the binary reading scans a 1
			if j&1 == 1 {

				k := n - (n & ((2 << i) - 1))
				k *= batchSize

				// If the rotation is not zero
				if k != 0 {

					// Rotate((tmpc0, tmpc1), k)
					eval.permuteNTTHoistedNoModDown(levelQ, eval.PoolDecompQ, eval.PoolDecompP, k, pool2Q, pool3Q, pool2P, pool3P)

					// ctOut += Rotate((tmpc0, tmpc1), k)
					if copy {
						ring.CopyValuesLvl(levelQ, pool2Q, ct0OutQ)
						ring.CopyValuesLvl(levelQ, pool3Q, ct1OutQ)
						ring.CopyValues(pool2P, ct0OutP)
						ring.CopyValues(pool3P, ct1OutP)
						copy = false
					} else {
						ringQ.AddLvl(levelQ, ct0OutQ, pool2Q, ct0OutQ)
						ringQ.AddLvl(levelQ, ct1OutQ, pool3Q, ct1OutQ)
						ringP.Add(ct0OutP, pool2P, ct0OutP)
						ringP.Add(ct1OutP, pool3P, ct1OutP)
					}

					if i == 0 {
						ring.PermuteNTTWithIndexLvl(levelQ, ctIn.Value[0], eval.permuteNTTIndex[eval.params.GaloisElementForColumnRotationBy(k)], tmpc2)
					} else {
						ring.PermuteNTTWithIndexLvl(levelQ, tmpc0, eval.permuteNTTIndex[eval.params.GaloisElementForColumnRotationBy(k)], tmpc2)
					}

					ringQ.MulScalarBigintLvl(levelQ, tmpc2, ringP.ModulusBigint, tmpc2)
					ringQ.AddLvl(levelQ, ct0OutQ, tmpc2, ct0OutQ)

				} else {

					state = true

					// if n is not a power of two
					if n&(n-1) != 0 {
						eval.Baseconverter.ModDownSplitNTTPQ(levelQ, ct0OutQ, ct0OutP, ct0OutQ) // Division by P
						eval.Baseconverter.ModDownSplitNTTPQ(levelQ, ct1OutQ, ct1OutP, ct1OutQ) // Division by P

						// ctOut += (tmpc0, tmpc1)
						ringQ.AddLvl(levelQ, ct0OutQ, tmpc0, ctOut.Value[0])
						ringQ.AddLvl(levelQ, ct1OutQ, tmpc1, ctOut.Value[1])

					} else {
						ring.CopyValuesLvl(levelQ, tmpc0, ctOut.Value[0])
						ring.CopyValuesLvl(levelQ, tmpc1, ctOut.Value[1])
						ctOut.Value[0].IsNTT = true
						ctOut.Value[1].IsNTT = true
					}
				}
			}

			if !state {
				if i == 0 {
					eval.permuteNTTHoisted(levelQ, ctIn.Value[0], ctIn.Value[1], eval.PoolDecompQ, eval.PoolDecompP, (1<<i)*batchSize, tmpc0, tmpc1)

					ringQ.AddLvl(levelQ, tmpc0, ctIn.Value[0], tmpc0)
					ringQ.AddLvl(levelQ, tmpc1, ctIn.Value[1], tmpc1)
				} else {
					// (tmpc0, tmpc1) = Rotate((tmpc0, tmpc1), 2^i)
					eval.permuteNTTHoisted(levelQ, tmpc0, tmpc1, eval.PoolDecompQ, eval.PoolDecompP, (1<<i)*batchSize, pool2Q, pool3Q)
					ringQ.AddLvl(levelQ, tmpc0, pool2Q, tmpc0)
					ringQ.AddLvl(levelQ, tmpc1, pool3Q, tmpc1)
				}

			}
		}
	}
}

// InnerSum applies an naive inner sum on the ciphetext (n rotations with single hoisting).
// The operation assumes that `ctIn` encrypts SlotCount/`batchSize` sub-vectors of size `batchSize` which it adds together (in parallel) by groups of `n`.
// It outputs in ctOut a ciphertext for which the "leftmost" sub-vector of each group is equal to the sum of the group.
// This method is faster than InnerSumLog when the number of rotations is small but uses 'n' keys instead of log(n) + HW(n).
func (eval *evaluator) InnerSum(ctIn *Ciphertext, batchSize, n int, ctOut *Ciphertext) {

	ringQ := eval.ringQ
	ringP := eval.ringP

	levelQ := ctIn.Level()

	QiOverF := eval.params.QiOverflowMargin(levelQ) >> 1
	PiOverF := eval.params.PiOverflowMargin() >> 1

	// If sum with only the first element, then returns the input
	if n == 1 {

		// If the input-output points differ, copies on the output
		if ctIn != ctOut {
			ring.CopyValuesLvl(levelQ, ctIn.Value[0], ctOut.Value[0])
			ring.CopyValuesLvl(levelQ, ctIn.Value[1], ctOut.Value[1])
		}
		// If sum on at least two elements
	} else {

		// List of n-2 rotations
		rotations := []int{}
		for i := 1; i < n; i++ {
			rotations = append(rotations, i*batchSize)
		}

		// Memory pool
		tmpQ0 := eval.poolQMul[0] // unused memory pool from evaluator
		tmpQ1 := eval.poolQMul[1] // unused memory pool from evaluator

		pool2P := eval.PoolP[1] // ctOut(c0', c1') from evaluator keyswitch memory pool
		pool3P := eval.PoolP[2] // ctOut(c0', c1') from evaluator keyswitch memory pool

		// Basis decomposition
		eval.DecomposeNTT(levelQ, ctIn.Value[1], eval.PoolDecompQ, eval.PoolDecompP)

		// Pre-rotates all [1, ..., n-1] rotations
		// Hoisted rotation without division by P
		vecRotQ, vecRotP := eval.rotateHoistedNoModDown(ctIn, rotations, eval.PoolDecompQ, eval.PoolDecompP)

		// P*c0 -> tmpQ0
		ringQ.MulScalarBigintLvl(levelQ, ctIn.Value[0], ringP.ModulusBigint, tmpQ0)

		// Adds phi_k(P*c0) on each of the vecRotQ
		// Does not need to add on the vecRotP because mod P === 0
		for _, i := range rotations {
			if i != 0 {

				galEl := eval.params.GaloisElementForColumnRotationBy(i)

				_, generated := eval.rtks.Keys[galEl]
				if !generated {
					panic("switching key not available")
				}

				index := eval.permuteNTTIndex[galEl]

				ring.PermuteNTTWithIndexLvl(levelQ, tmpQ0, index, tmpQ1)  // phi(P*c0)
				ringQ.AddLvl(levelQ, vecRotQ[i][0], tmpQ1, vecRotQ[i][0]) // phi(d0_Q) += phi(P*c0)
			}
		}

		var reduce int
		// Sums elements [2, ..., n-1]
		for i := 1; i < n; i++ {

			j := i * batchSize

			if i == 1 {
				ring.CopyValuesLvl(levelQ, vecRotQ[j][0], tmpQ0)
				ring.CopyValuesLvl(levelQ, vecRotQ[j][1], tmpQ1)
				ring.CopyValues(vecRotP[j][0], pool2P)
				ring.CopyValues(vecRotP[j][1], pool3P)
			} else {
				ringQ.AddNoModLvl(levelQ, tmpQ0, vecRotQ[j][0], tmpQ0)
				ringQ.AddNoModLvl(levelQ, tmpQ1, vecRotQ[j][1], tmpQ1)
				ringP.AddNoMod(pool2P, vecRotP[j][0], pool2P)
				ringP.AddNoMod(pool3P, vecRotP[j][1], pool3P)
			}

			if reduce%QiOverF == QiOverF-1 {
				ringQ.ReduceLvl(levelQ, tmpQ0, tmpQ0)
				ringQ.ReduceLvl(levelQ, tmpQ1, tmpQ1)
			}

			if reduce%PiOverF == PiOverF-1 {
				ringP.Reduce(pool2P, pool2P)
				ringP.Reduce(pool3P, pool3P)
			}

			reduce++
		}

		if reduce%QiOverF != 0 {
			ringQ.ReduceLvl(levelQ, tmpQ0, tmpQ0)
			ringQ.ReduceLvl(levelQ, tmpQ1, tmpQ1)
		}

		if reduce%PiOverF != 0 {
			ringP.Reduce(pool2P, pool2P)
			ringP.Reduce(pool3P, pool3P)
		}

		// Division by P of sum(elements [2, ..., n-1] )
		eval.Baseconverter.ModDownSplitNTTPQ(levelQ, tmpQ0, pool2P, tmpQ0) // sum_{i=1, n-1}(phi(d0))/P
		eval.Baseconverter.ModDownSplitNTTPQ(levelQ, tmpQ1, pool3P, tmpQ1) // sum_{i=1, n-1}(phi(d1))/P

		// Adds element[1] (which did not require rotation)
		ringQ.AddLvl(levelQ, ctIn.Value[0], tmpQ0, ctOut.Value[0]) // sum_{i=1, n-1}(phi(d0))/P + ct0
		ringQ.AddLvl(levelQ, ctIn.Value[1], tmpQ1, ctOut.Value[1]) // sum_{i=1, n-1}(phi(d1))/P + ct1
	}
}

// ReplicateLog applies an optimized replication on the ciphetext (log2(n) + HW(n) rotations with double hoisting).
// It acts as the inverse of a inner sum (summing elements from left to right).
// The replication is parameterized by the size of the sub-vectors to replicate "batchSize" and
// the number of time "n" they need to be replicated.
// To ensure correctness, a gap of zero values of size batchSize * (n-1) must exist between
// two consecutive sub-vectors to replicate.
// This method is faster than Replicate when the number of rotations is large and uses log2(n) + HW(n) instead of 'n'.
func (eval *evaluator) ReplicateLog(ctIn *Ciphertext, batchSize, n int, ctOut *Ciphertext) {
	eval.InnerSumLog(ctIn, -batchSize, n, ctOut)
}

// Replicate applies naive replication on the ciphetext (n rotations with single hoisting).
// It acts as the inverse of a inner sum (summing elements from left to right).
// The replication is parameterized by the size of the sub-vectors to replicate "batchSize" and
// the number of time "n" they need to be replicated.
// To ensure correctness, a gap of zero values of size batchSize * (n-1) must exist between
// two consecutive sub-vectors to replicate.
// This method is faster than ReplicateLog when the number of rotations is small but uses 'n' keys instead of log2(n) + HW(n).
func (eval *evaluator) Replicate(ctIn *Ciphertext, batchSize, n int, ctOut *Ciphertext) {
	eval.InnerSum(ctIn, -batchSize, n, ctOut)
}

// MultiplyByDiagMatrix multiplies the ciphertext "ctIn" by the plaintext matrix "matrix" and returns the result on the ciphertext
// "ctOut". Memory pools for the decomposed ciphertext PoolDecompQ, PoolDecompP must be provided, those are list of poly of ringQ and ringP
// respectively, each of size params.Beta().
// The naive approach is used (single hoisting and no baby-step giant-step), which is faster than MultiplyByDiagMatrixBSGS
// for matrix of only a few non-zero diagonals but uses more keys.
func (eval *evaluator) MultiplyByDiagMatrix(ctIn *Ciphertext, matrix *PtDiagMatrix, PoolDecompQ, PoolDecompP []*ring.Poly, ctOut *Ciphertext) {

	ringQ := eval.ringQ
	ringP := eval.ringP

	levelQ := utils.MinInt(ctOut.Level(), utils.MinInt(ctIn.Level(), matrix.Level))
	levelP := eval.params.PCount() - 1

	QiOverF := eval.params.QiOverflowMargin(levelQ)
	PiOverF := eval.params.PiOverflowMargin()

	ksResP0 := eval.PoolP[0]  // Key-Switch ctOut[0] mod P
	ksResP1 := eval.PoolP[1]  // Key-Switch ctOut[1] mod P
	tmpP0 := eval.PoolP[2]    // Automorphism not-inplace pool res[0] mod P
	tmpP1 := eval.poolQMul[0] // Automorphism not-inplace pool res[1] mod P
	accP0 := eval.PoolP[3]    // Accumulator ctOut[0] mod P
	accP1 := eval.PoolP[4]    // Accumulator ctOut[1] mod P

	ct0TimesP := eval.PoolQ[0] // ct0 * P mod Q
	ksResQ0 := eval.PoolQ[1]   // Key-Switch ctOut[0] mod Q
	ksResQ1 := eval.PoolQ[2]   // Key-Switch ctOut[0] mod Q
	tmpQ0 := eval.PoolQ[3]     // Automorphism not-inplace pool ctOut[0] mod Q
	tmpQ1 := eval.PoolQ[4]     // Automorphism not-inplace pool ctOut[1] mod Q

	ringQ.MulScalarBigintLvl(levelQ, ctIn.Value[0], ringP.ModulusBigint, ct0TimesP) // P*c0

	var state bool
	var cnt int
	for k := range matrix.Vec {

		k &= int((ringQ.N >> 1) - 1)

		if k == 0 {
			state = true
		} else {

			galEl := eval.params.GaloisElementForColumnRotationBy(k)

			rtk, generated := eval.rtks.Keys[galEl]
			if !generated {
				panic("switching key not available")
			}

			index := eval.permuteNTTIndex[galEl]

			eval.KeyswitchHoistedNoModDown(levelQ, PoolDecompQ, PoolDecompP, rtk, ksResQ0, ksResQ1, ksResP0, ksResP1)

			ringQ.AddLvl(levelQ, ksResQ0, ct0TimesP, ksResQ0) // phi(d0_Q) += phi(P*c0)

			ring.PermuteNTTWithIndexLvl(levelQ, ksResQ0, index, tmpQ0) // phi(P*c0 + d0_Q)
			ring.PermuteNTTWithIndexLvl(levelQ, ksResQ1, index, tmpQ1) // phi(       d1_Q)

			ring.PermuteNTTWithIndexLvl(levelP, ksResP0, index, tmpP0) // phi(P*c0 + d0_P)
			ring.PermuteNTTWithIndexLvl(levelP, ksResP1, index, tmpP1) // phi(       d1_P)

			plaintextQ := matrix.Vec[k][0]
			plaintextP := matrix.Vec[k][1]

			if cnt == 0 {
				// keyswitch(c1_Q) = (d0_QP, d1_QP)
				ringQ.MulCoeffsMontgomeryLvl(levelQ, plaintextQ, tmpQ0, ctOut.Value[0]) // phi(P*c0 + d0_Q) * plaintext
				ringQ.MulCoeffsMontgomeryLvl(levelQ, plaintextQ, tmpQ1, ctOut.Value[1]) // phi(d1_Q) * plaintext
				ringP.MulCoeffsMontgomery(plaintextP, tmpP0, accP0)                     // phi(d0_P) * plaintext
				ringP.MulCoeffsMontgomery(plaintextP, tmpP1, accP1)                     // phi(d1_P) * plaintext
			} else {
				// keyswitch(c1_Q) = (d0_QP, d1_QP)
				ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, plaintextQ, tmpQ0, ctOut.Value[0]) // phi(P*c0 + d0_Q) * plaintext
				ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, plaintextQ, tmpQ1, ctOut.Value[1]) // phi(d1_Q) * plaintext
				ringP.MulCoeffsMontgomeryAndAdd(plaintextP, tmpP0, accP0)                     // phi(d0_P) * plaintext
				ringP.MulCoeffsMontgomeryAndAdd(plaintextP, tmpP1, accP1)                     // phi(d1_P) * plaintext
			}

			if cnt%QiOverF == QiOverF-1 {
				ringQ.ReduceLvl(levelQ, ctOut.Value[0], ctOut.Value[0])
				ringQ.ReduceLvl(levelQ, ctOut.Value[1], ctOut.Value[1])
			}

			if cnt%PiOverF == PiOverF-1 {
				ringP.Reduce(accP0, accP0)
				ringP.Reduce(accP1, accP1)
			}

			cnt++
		}
	}

	if cnt%QiOverF == 0 {
		ringQ.ReduceLvl(levelQ, ctOut.Value[0], ctOut.Value[0])
		ringQ.ReduceLvl(levelQ, ctOut.Value[1], ctOut.Value[1])
	}

	if cnt%PiOverF == 0 {
		ringP.Reduce(accP0, accP0)
		ringP.Reduce(accP1, accP1)
	}

	eval.Baseconverter.ModDownSplitNTTPQ(levelQ, ctOut.Value[0], accP0, ctOut.Value[0]) // sum(phi(c0 * P + d0_QP))/P
	eval.Baseconverter.ModDownSplitNTTPQ(levelQ, ctOut.Value[1], accP1, ctOut.Value[1]) // sum(phi(d1_QP))/P

	if state { // Rotation by zero
		ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[0][0], ctIn.Value[0], ctOut.Value[0]) // ctOut += c0_Q * plaintext
		ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[0][0], ctIn.Value[1], ctOut.Value[1]) // ctOut += c1_Q * plaintext
	}

	ctOut.Scale = matrix.Scale * ctIn.Scale
}

// MultiplyByDiagMatrixBSGS multiplies the ciphertext "ctIn" by the plaintext matrix "matrix" and returns the result on the ciphertext
// "ctOut". Memory pools for the decomposed ciphertext PoolDecompQ, PoolDecompP must be provided, those are list of poly of ringQ and ringP
// respectively, each of size params.Beta().
// The BSGS approach is used (double hoisting with baby-step giant-step), which is faster than MultiplyByDiagMatrix
// for matrix with more than a few non-zero diagonals and uses much less keys.
func (eval *evaluator) MultiplyByDiagMatrixBSGS(ctIn *Ciphertext, matrix *PtDiagMatrix, PoolDecompQ, PoolDecompP []*ring.Poly, ctOut *Ciphertext) {

	// N1*N2 = N
	N1 := matrix.N1

	ringQ := eval.ringQ
	ringP := eval.ringP

	levelQ := utils.MinInt(ctOut.Level(), utils.MinInt(ctIn.Level(), matrix.Level))
	levelP := eval.params.PCount() - 1

	QiOverF := eval.params.QiOverflowMargin(levelQ)
	PiOverF := eval.params.PiOverflowMargin()

	// Computes the rotations indexes of the non-zero rows of the diagonalized DFT matrix for the baby-step giang-step algorithm

	index, rotations := bsgsIndex(matrix.Vec, 1<<matrix.LogSlots, matrix.N1)

	// Pre-rotates ciphertext for the baby-step giant-step algorithm, does not divide by P yet
	vecRotQ, vecRotP := eval.rotateHoistedNoModDown(ctIn, rotations, eval.PoolDecompQ, eval.PoolDecompP)

	// Accumulator inner loop
	tmpQ0 := eval.poolQMul[0] // unused memory pool from evaluator
	tmpQ1 := eval.poolQMul[1] // unused memory pool from evaluator

	// Accumulator outer loop
	tmpQ2 := eval.poolQMul[2] // unused memory pool from evaluator
	tmpQ3 := eval.PoolQ[4]
	tmpP2 := eval.PoolP[3]
	tmpP3 := eval.PoolP[4]

	// Keyswitch accumulator
	pool2Q := eval.PoolQ[1] // ctOut(c0', c1') from evaluator keyswitch memory pool
	pool3Q := eval.PoolQ[2] // ctOut(c0', c1') from evaluator keyswitch memory pool
	pool2P := eval.PoolP[1] // ctOut(c0', c1') from evaluator keyswitch memory pool
	pool3P := eval.PoolP[2] // ctOut(c0', c1') from evaluator keyswitch memory pool

	// Do not use (used by switchKeysInPlaceNoModDown)
	// eval.PoolP[0]
	// eval.PoolQ[0]
	// eval.PoolQ[2]

	N1Rot := 0
	N2Rot := 0

	ringQ.MulScalarBigintLvl(levelQ, ctIn.Value[0], ringP.ModulusBigint, tmpQ0) // P*c0

	for _, i := range rotations {
		if i != 0 {

			galEl := eval.params.GaloisElementForColumnRotationBy(i)

			_, generated := eval.rtks.Keys[galEl]
			if !generated {
				panic("switching key not available")
			}

			index := eval.permuteNTTIndex[galEl]

			ring.PermuteNTTWithIndexLvl(levelQ, tmpQ0, index, tmpQ1)  // phi(P*c0)
			ringQ.AddLvl(levelQ, vecRotQ[i][0], tmpQ1, vecRotQ[i][0]) // phi(d0_Q) += phi(P*c0)
		}
	}

	// OUTER LOOP
	var cnt0 int
	for j := range index {

		if j != 0 {

			// INNER LOOP
			var state bool
			var cnt1 int
			for _, i := range index[j] {

				if i == 0 {
					state = true
				} else {

					N1Rot++

					plaintextQ := matrix.Vec[N1*j+i][0]
					plaintextP := matrix.Vec[N1*j+i][1]

					if cnt1 == 0 {
						ringQ.MulCoeffsMontgomeryConstantLvl(levelQ, plaintextQ, vecRotQ[i][0], tmpQ0) // phi(P*c0 + d0_Q) * plaintext
						ringQ.MulCoeffsMontgomeryConstantLvl(levelQ, plaintextQ, vecRotQ[i][1], tmpQ1) // phi(d1_Q) * plaintext
						ringP.MulCoeffsMontgomeryConstant(plaintextP, vecRotP[i][0], pool2P)           // phi(d0_P) * plaintext
						ringP.MulCoeffsMontgomeryConstant(plaintextP, vecRotP[i][1], pool3P)           // phi(d1_P) * plaintext
					} else {
						ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, plaintextQ, vecRotQ[i][0], tmpQ0) // phi(d0_Q) * plaintext
						ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, plaintextQ, vecRotQ[i][1], tmpQ1) // phi(d1_Q) * plaintext
						ringP.MulCoeffsMontgomeryAndAddNoMod(plaintextP, vecRotP[i][0], pool2P)                   // phi(d0_P) * plaintext
						ringP.MulCoeffsMontgomeryAndAddNoMod(plaintextP, vecRotP[i][1], pool3P)                   // phi(d1_P) * plaintext
					}

					if cnt1%(QiOverF>>1) == (QiOverF>>1)-1 {
						ringQ.ReduceLvl(levelQ, tmpQ0, tmpQ0)
						ringQ.ReduceLvl(levelQ, tmpQ1, tmpQ1)
					}

					if cnt1%(PiOverF>>1) == (PiOverF>>1)-1 {
						ringP.Reduce(pool2P, pool2P)
						ringP.Reduce(pool3P, pool3P)
					}

					cnt1++
				}
			}

			if cnt1%(QiOverF>>1) != 0 {
				ringQ.ReduceLvl(levelQ, tmpQ0, tmpQ0)
				ringQ.ReduceLvl(levelQ, tmpQ1, tmpQ1)
			}

			if cnt1%(PiOverF>>1) != 0 {
				ringP.Reduce(pool2P, pool2P)
				ringP.Reduce(pool3P, pool3P)
			}

			// Hoisting of the ModDown of sum(sum(phi(d0 + P*c0) * plaintext)) and sum(sum(phi(d1) * plaintext))
			eval.Baseconverter.ModDownSplitNTTPQ(levelQ, tmpQ0, pool2P, tmpQ0) // sum(phi(d0) * plaintext)/P
			eval.Baseconverter.ModDownSplitNTTPQ(levelQ, tmpQ1, pool3P, tmpQ1) // sum(phi(d1) * plaintext)/P

			// If i == 0
			if state {

				// If no loop before, then we copy the values on the accumulator instead of adding them
				if len(index[j]) == 1 {
					ringQ.MulCoeffsMontgomeryLvl(levelQ, matrix.Vec[N1*j][0], ctIn.Value[0], tmpQ0) // c0 * plaintext + sum(phi(d0) * plaintext)/P + phi(c0) * plaintext mod Q
					ringQ.MulCoeffsMontgomeryLvl(levelQ, matrix.Vec[N1*j][0], ctIn.Value[1], tmpQ1) // c1 * plaintext + sum(phi(d1) * plaintext)/P + phi(c1) * plaintext mod Q
				} else {
					ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[N1*j][0], ctIn.Value[0], tmpQ0) // c0 * plaintext + sum(phi(d0) * plaintext)/P + phi(c0) * plaintext mod Q
					ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[N1*j][0], ctIn.Value[1], tmpQ1) // c1 * plaintext + sum(phi(d1) * plaintext)/P + phi(c1) * plaintext mod Q
				}

				N1Rot++
			}

			galEl := eval.params.GaloisElementForColumnRotationBy(N1 * j)

			rtk, generated := eval.rtks.Keys[galEl]
			if !generated {
				panic("switching key not available")
			}

			index := eval.permuteNTTIndex[galEl]

			tmpQ1.IsNTT = true
			eval.SwitchKeysInPlaceNoModDown(levelQ, tmpQ1, rtk, pool2Q, pool2P, pool3Q, pool3P) // Switchkey(phi(tmpRes_1)) = (d0, d1) in base QP

			// Outer loop rotations
			ring.PermuteNTTWithIndexLvl(levelQ, tmpQ0, index, tmpQ1)    // phi(tmpRes_0)
			ringQ.AddLvl(levelQ, ctOut.Value[0], tmpQ1, ctOut.Value[0]) // ctOut += phi(tmpRes)

			N2Rot++

			if cnt0 == 0 {
				ring.PermuteNTTWithIndexLvl(levelQ, pool2Q, index, tmpQ2) // sum(phi(d0_Q))
				ring.PermuteNTTWithIndexLvl(levelQ, pool3Q, index, tmpQ3) // sum(phi(d1_Q))
				ring.PermuteNTTWithIndexLvl(levelP, pool2P, index, tmpP2) // sum(phi(d0_P))
				ring.PermuteNTTWithIndexLvl(levelP, pool3P, index, tmpP3) // sum(phi(d1_P))
			} else {
				ring.PermuteNTTWithIndexAndAddNoModLvl(levelQ, pool2Q, index, tmpQ2) // sum(phi(d0_Q))
				ring.PermuteNTTWithIndexAndAddNoModLvl(levelQ, pool3Q, index, tmpQ3) // sum(phi(d1_Q))
				ring.PermuteNTTWithIndexAndAddNoModLvl(levelP, pool2P, index, tmpP2) // sum(phi(d0_P))
				ring.PermuteNTTWithIndexAndAddNoModLvl(levelP, pool3P, index, tmpP3) // sum(phi(d1_P))
			}

			if cnt0%QiOverF == QiOverF-1 {
				ringQ.ReduceLvl(levelQ, tmpQ2, tmpQ2)
				ringQ.ReduceLvl(levelQ, tmpQ3, tmpQ3)
			}

			if cnt0%PiOverF == PiOverF-1 {
				ringP.Reduce(tmpP2, tmpP2)
				ringP.Reduce(tmpP3, tmpP3)
			}

			cnt0++
		}
	}

	if cnt0%QiOverF != 0 {
		ringQ.ReduceLvl(levelQ, tmpQ2, tmpQ2)
		ringQ.ReduceLvl(levelQ, tmpQ3, tmpQ3)
	}

	if cnt0%PiOverF != 0 {
		ringP.Reduce(tmpP2, tmpP2)
		ringP.Reduce(tmpP3, tmpP3)
	}

	// if j == 0 (N2 rotation by zero)
	var state bool
	var cnt1 int
	for _, i := range index[0] {

		if i == 0 {
			state = true
		} else {

			plaintextQ := matrix.Vec[i][0]
			plaintextP := matrix.Vec[i][1]
			N1Rot++
			// keyswitch(c1_Q) = (d0_QP, d1_QP)
			ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, plaintextQ, vecRotQ[i][0], tmpQ2) // phi(P*c0 + d0_Q) * plaintext
			ringQ.MulCoeffsMontgomeryConstantAndAddNoModLvl(levelQ, plaintextQ, vecRotQ[i][1], tmpQ3) // phi(d1_Q) * plaintext
			ringP.MulCoeffsMontgomeryAndAddNoMod(plaintextP, vecRotP[i][0], tmpP2)                    // phi(d0_P) * plaintext
			ringP.MulCoeffsMontgomeryAndAddNoMod(plaintextP, vecRotP[i][1], tmpP3)                    // phi(d1_P) * plaintext

			if cnt1%(QiOverF>>1) == (QiOverF>>1)-1 {
				ringQ.ReduceLvl(levelQ, tmpQ2, tmpQ2)
				ringQ.ReduceLvl(levelQ, tmpQ3, tmpQ3)
			}

			if cnt1%(PiOverF>>1) == (PiOverF>>1)-1 {
				ringP.Reduce(tmpP2, tmpP2)
				ringP.Reduce(tmpP3, tmpP3)
			}

			cnt1++
		}
	}

	if cnt1%(QiOverF>>1) != 0 {
		ringQ.ReduceLvl(levelQ, tmpQ2, tmpQ2)
		ringQ.ReduceLvl(levelQ, tmpQ3, tmpQ3)
	}

	if cnt1%(PiOverF>>1) != 0 {
		ringP.Reduce(tmpP2, tmpP2)
		ringP.Reduce(tmpP3, tmpP3)
	}

	eval.Baseconverter.ModDownSplitNTTPQ(levelQ, tmpQ2, tmpP2, tmpQ2) // sum(phi(c0 * P + d0_QP))/P
	eval.Baseconverter.ModDownSplitNTTPQ(levelQ, tmpQ3, tmpP3, tmpQ3) // sum(phi(d1_QP))/P

	ringQ.AddLvl(levelQ, ctOut.Value[0], tmpQ2, ctOut.Value[0]) // ctOut += sum(phi(c0 * P + d0_QP))/P
	ringQ.AddLvl(levelQ, ctOut.Value[1], tmpQ3, ctOut.Value[1]) // ctOut += sum(phi(d1_QP))/P

	if state { // Rotation by zero
		N1Rot++
		ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[0][0], ctIn.Value[0], ctOut.Value[0]) // ctOut += c0_Q * plaintext
		ringQ.MulCoeffsMontgomeryAndAddLvl(levelQ, matrix.Vec[0][0], ctIn.Value[1], ctOut.Value[1]) // ctOut += c1_Q * plaintext
	}

	ctOut.Scale = matrix.Scale * ctIn.Scale

	vecRotQ, vecRotP = nil, nil
}
