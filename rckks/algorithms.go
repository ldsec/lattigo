package rckks

import (
	"math/bits"
)

// PowerOf2 computes op^(2^logPow2), consuming logPow2 levels, and returns the result on opOut. Providing an evaluation
// key is necessary when logPow2 > 1.
func (eval *evaluator) PowerOf2(op *Ciphertext, logPow2 uint64, evakey *EvaluationKey, opOut *Ciphertext) {

	if logPow2 == 0 {

		if op != opOut {

			opOut.Copy(op.El())
		}

	} else {

		eval.MulRelin(op.El(), op.El(), evakey, opOut)

		eval.Rescale(opOut, eval.scale, opOut)

		for i := uint64(1); i < logPow2; i++ {

			eval.MulRelin(opOut.El(), opOut.El(), evakey, opOut)

			eval.Rescale(opOut, eval.scale, opOut)
		}
	}
}

// PowerNew computes op^degree, consuming log(degree) levels, and returns the result on a new element. Providing an evaluation
// key is necessary when degree > 2.
func (eval *evaluator) PowerNew(op *Ciphertext, degree uint64, evakey *EvaluationKey) (opOut *Ciphertext) {
	opOut = NewCiphertext(eval.params, 1, op.Level(), op.Scale())
	eval.Power(op, degree, evakey, opOut)
	return
}

// Power computes op^degree, consuming log(degree) levels, and returns the result on opOut. Providing an evaluation
// key is necessary when degree > 2.
func (eval *evaluator) Power(op *Ciphertext, degree uint64, evakey *EvaluationKey, opOut *Ciphertext) {

	tmpct0 := op.CopyNew()

	var logDegree, po2Degree uint64

	logDegree = uint64(bits.Len64(degree)) - 1
	po2Degree = 1 << logDegree

	eval.PowerOf2(tmpct0.Ciphertext(), logDegree, evakey, opOut)

	degree -= po2Degree

	for degree > 0 {

		logDegree = uint64(bits.Len64(degree)) - 1
		po2Degree = 1 << logDegree

		tmp := NewCiphertext(eval.params, 1, tmpct0.Level(), tmpct0.Scale())

		eval.PowerOf2(tmpct0.Ciphertext(), logDegree, evakey, tmp)

		eval.MulRelin(opOut.El(), tmp.El(), evakey, opOut)

		eval.Rescale(opOut, eval.scale, opOut)

		degree -= po2Degree
	}
}

// InverseNew computes 1/op and returns the result on a new element, iterating for n steps and consuming n levels. The algorithm requires the encrypted values to be in the range
// [-1.5 - 1.5i, 1.5 + 1.5i] or the result will be wrong. Each iteration increases the precision.
func (eval *evaluator) InverseNew(op *Ciphertext, steps uint64, evakey *EvaluationKey) (opOut *Ciphertext) {

	cbar := eval.NegNew(op)

	eval.AddConst(cbar, 1, cbar)

	tmp := eval.AddConstNew(cbar, 1)
	opOut = tmp.CopyNew().Ciphertext()

	for i := uint64(1); i < steps; i++ {

		eval.MulRelin(cbar.El(), cbar.El(), evakey, cbar.Ciphertext())

		eval.Rescale(cbar, eval.scale, cbar)

		tmp = eval.AddConstNew(cbar, 1)

		eval.MulRelin(tmp.El(), opOut.El(), evakey, tmp.Ciphertext())

		eval.Rescale(tmp, eval.scale, tmp)

		opOut = tmp.CopyNew().Ciphertext()
	}

	return opOut
}
