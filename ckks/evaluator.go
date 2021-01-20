package ckks

import (
	"errors"
	"math"
	"unsafe"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// Evaluator is an interface implementing the methodes to conduct homomorphic operations between ciphertext and/or plaintexts.
type Evaluator interface {
	Add(op0, op1 Operand, ctOut *Ciphertext)
	AddNoMod(op0, op1 Operand, ctOut *Ciphertext)
	AddNew(op0, op1 Operand) (ctOut *Ciphertext)
	AddNoModNew(op0, op1 Operand) (ctOut *Ciphertext)
	Sub(op0, op1 Operand, ctOut *Ciphertext)
	SubNoMod(op0, op1 Operand, ctOut *Ciphertext)
	SubNew(op0, op1 Operand) (ctOut *Ciphertext)
	SubNoModNew(op0, op1 Operand) (ctOut *Ciphertext)
	Neg(ct0 *Ciphertext, ctOut *Ciphertext)
	NegNew(ct0 *Ciphertext) (ctOut *Ciphertext)
	AddConstNew(ct0 *Ciphertext, constant interface{}) (ctOut *Ciphertext)
	AddConst(ct0 *Ciphertext, constant interface{}, ctOut *Ciphertext)
	MultByConstAndAdd(ct0 *Ciphertext, constant interface{}, ctOut *Ciphertext)
	MultByConstNew(ct0 *Ciphertext, constant interface{}) (ctOut *Ciphertext)
	MultByConst(ct0 *Ciphertext, constant interface{}, ctOut *Ciphertext)
	MultByGaussianInteger(ct0 *Ciphertext, cReal, cImag int64, ctOut *Ciphertext)
	MultByGaussianIntegerAndAdd(ct0 *Ciphertext, cReal, cImag int64, ctOut *Ciphertext)
	MultByiNew(ct0 *Ciphertext) (ctOut *Ciphertext)
	MultByi(ct0 *Ciphertext, ct1 *Ciphertext)
	DivByiNew(ct0 *Ciphertext) (ctOut *Ciphertext)
	DivByi(ct0 *Ciphertext, ct1 *Ciphertext)
	ScaleUpNew(ct0 *Ciphertext, scale float64) (ctOut *Ciphertext)
	ScaleUp(ct0 *Ciphertext, scale float64, ctOut *Ciphertext)
	SetScale(ct *Ciphertext, scale float64)
	MulByPow2New(ct0 *Ciphertext, pow2 uint64) (ctOut *Ciphertext)
	MulByPow2(ct0 *Element, pow2 uint64, ctOut *Element)
	ReduceNew(ct0 *Ciphertext) (ctOut *Ciphertext)
	Reduce(ct0 *Ciphertext, ctOut *Ciphertext) error
	DropLevelNew(ct0 *Ciphertext, levels uint64) (ctOut *Ciphertext)
	DropLevel(ct0 *Ciphertext, levels uint64)
	Rescale(ct0 *Ciphertext, threshold float64, c1 *Ciphertext) (err error)
	RescaleNew(ct0 *Ciphertext, threshold float64) (ctOut *Ciphertext, err error)
	RescaleMany(ct0 *Ciphertext, nbRescales uint64, c1 *Ciphertext) (err error)
	MulRelinNew(op0, op1 Operand, evakey *EvaluationKey) (ctOut *Ciphertext)
	MulRelin(op0, op1 Operand, evakey *EvaluationKey, ctOut *Ciphertext)
	RelinearizeNew(ct0 *Ciphertext, evakey *EvaluationKey) (ctOut *Ciphertext)
	Relinearize(ct0 *Ciphertext, evakey *EvaluationKey, ctOut *Ciphertext)
	SwitchKeysNew(ct0 *Ciphertext, switchingKey *SwitchingKey) (ctOut *Ciphertext)
	SwitchKeys(ct0 *Ciphertext, switchingKey *SwitchingKey, ctOut *Ciphertext)
	RotateNew(ct0 *Ciphertext, k uint64, evakey *RotationKeys) (ctOut *Ciphertext)
	Rotate(ct0 *Ciphertext, k uint64, evakey *RotationKeys, ctOut *Ciphertext)
	RotateHoisted(ctIn *Ciphertext, rotations []uint64, rotkeys *RotationKeys) (cOut map[uint64]*Ciphertext)
	ConjugateNew(ct0 *Ciphertext, evakey *RotationKeys) (ctOut *Ciphertext)
	Conjugate(ct0 *Ciphertext, evakey *RotationKeys, ctOut *Ciphertext)
	PowerOf2(el0 *Ciphertext, logPow2 uint64, evakey *EvaluationKey, elOut *Ciphertext)
	PowerNew(op *Ciphertext, degree uint64, evakey *EvaluationKey) (opOut *Ciphertext)
	Power(ct0 *Ciphertext, degree uint64, evakey *EvaluationKey, res *Ciphertext)
	InverseNew(ct0 *Ciphertext, steps uint64, evakey *EvaluationKey) (res *Ciphertext)
	EvaluatePoly(ct *Ciphertext, coeffs *Poly, evakey *EvaluationKey) (res *Ciphertext, err error)
	EvaluateCheby(ct *Ciphertext, cheby *ChebyshevInterpolation, evakey *EvaluationKey) (res *Ciphertext, err error)
}

// evaluator is a struct that holds the necessary elements to execute the homomorphic operations between Ciphertexts and/or Plaintexts.
// It also holds a small memory pool used to store intermediate computations.
type evaluator struct {
	params *Parameters
	scale  float64

	ringQ    *ring.Ring
	ringP    *ring.Ring
	poolQMul [3]*ring.Poly // Memory pool in order : for MForm(c0), MForm(c1), c2

	poolQ [4]*ring.Poly // Memory pool in order : Decomp(c2), for NTT^-1(c2), res(c0', c1')
	poolP [3]*ring.Poly // Memory pool in order : Decomp(c2), res(c0', c1')

	ctxpool *Ciphertext // Memory pool for ciphertext that need to be scaled up (to be removed eventually)

	baseconverter *ring.FastBasisExtender
	decomposer    *ring.Decomposer
}

// NewEvaluator creates a new Evaluator, that can be used to do homomorphic
// operations on the Ciphertexts and/or Plaintexts. It stores a small pool of polynomials
// and Ciphertexts that will be used for intermediate values.
func NewEvaluator(params *Parameters) Evaluator {

	var q, p *ring.Ring
	var err error
	if q, err = ring.NewRing(params.N(), params.qi); err != nil {
		panic(err)
	}

	if params.PiCount() != 0 {
		if p, err = ring.NewRing(params.N(), params.pi); err != nil {
			panic(err)
		}
	}

	var baseconverter *ring.FastBasisExtender
	var decomposer *ring.Decomposer
	var poolP [3]*ring.Poly
	if params.PiCount() != 0 {
		baseconverter = ring.NewFastBasisExtender(q, p)
		decomposer = ring.NewDecomposer(q.Modulus, p.Modulus)
		poolP = [3]*ring.Poly{p.NewPoly(), p.NewPoly(), p.NewPoly()}
	}

	return &evaluator{
		params:        params.Copy(),
		scale:         params.scale,
		ringQ:         q,
		ringP:         p,
		poolQMul:      [3]*ring.Poly{q.NewPoly(), q.NewPoly(), q.NewPoly()},
		poolQ:         [4]*ring.Poly{q.NewPoly(), q.NewPoly(), q.NewPoly(), q.NewPoly()},
		poolP:         poolP,
		ctxpool:       NewCiphertext(params, 1, params.MaxLevel(), params.scale),
		baseconverter: baseconverter,
		decomposer:    decomposer,
	}
}

func (eval *evaluator) getElemAndCheckBinary(op0, op1, opOut Operand, opOutMinDegree uint64) (el0, el1, elOut *Element) {
	if op0 == nil || op1 == nil || opOut == nil {
		panic("operands cannot be nil")
	}

	if op0.Degree()+op1.Degree() == 0 {
		panic("operands cannot be both plaintext")
	}

	if opOut.Degree() < opOutMinDegree {
		panic("receiver operand degree is too small")
	}
	el0, el1, elOut = op0.El(), op1.El(), opOut.El()
	return
}

func (eval *evaluator) newCiphertextBinary(op0, op1 Operand) (ctOut *Ciphertext) {

	maxDegree := utils.MaxUint64(op0.Degree(), op1.Degree())
	maxScale := utils.MaxFloat64(op0.Scale(), op1.Scale())
	minLevel := utils.MinUint64(op0.Level(), op1.Level())

	return NewCiphertext(eval.params, maxDegree, minLevel, maxScale)
}

// Add adds op0 to op1 and returns the result in ctOut.
func (eval *evaluator) Add(op0, op1 Operand, ctOut *Ciphertext) {
	el0, el1, elOut := eval.getElemAndCheckBinary(op0, op1, ctOut, utils.MaxUint64(op0.Degree(), op1.Degree()))
	eval.evaluateInPlace(el0, el1, elOut, eval.ringQ.AddLvl)
}

// AddNoMod adds op0 to op1 and returns the result in ctOut, without modular reduction.
func (eval *evaluator) AddNoMod(op0, op1 Operand, ctOut *Ciphertext) {
	el0, el1, elOut := eval.getElemAndCheckBinary(op0, op1, ctOut, utils.MaxUint64(op0.Degree(), op1.Degree()))
	eval.evaluateInPlace(el0, el1, elOut, eval.ringQ.AddNoModLvl)
}

// AddNew adds op0 to op1 and returns the result in a newly created element.
func (eval *evaluator) AddNew(op0, op1 Operand) (ctOut *Ciphertext) {
	ctOut = eval.newCiphertextBinary(op0, op1)
	eval.Add(op0, op1, ctOut)
	return
}

// AddNoModNew adds op0 to op1 without modular reduction, and returns the result in a newly created element.
func (eval *evaluator) AddNoModNew(op0, op1 Operand) (ctOut *Ciphertext) {
	ctOut = eval.newCiphertextBinary(op0, op1)
	eval.AddNoMod(op0, op1, ctOut)
	return
}

// Sub subtracts op1 from op0 and returns the result in ctOut.
func (eval *evaluator) Sub(op0, op1 Operand, ctOut *Ciphertext) {

	el0, el1, elOut := eval.getElemAndCheckBinary(op0, op1, ctOut, utils.MaxUint64(op0.Degree(), op1.Degree()))

	eval.evaluateInPlace(el0, el1, elOut, eval.ringQ.SubLvl)

	level := utils.MinUint64(utils.MinUint64(el0.Level(), el1.Level()), elOut.Level())

	if el0.Degree() < el1.Degree() {
		for i := el0.Degree() + 1; i < el1.Degree()+1; i++ {
			eval.ringQ.NegLvl(level, elOut.Value()[i], elOut.Value()[i])
		}
	}

}

// SubNoMod subtracts op1 from op0 and returns the result in ctOut, without modular reduction.
func (eval *evaluator) SubNoMod(op0, op1 Operand, ctOut *Ciphertext) {

	el0, el1, elOut := eval.getElemAndCheckBinary(op0, op1, ctOut, utils.MaxUint64(op0.Degree(), op1.Degree()))

	eval.evaluateInPlace(el0, el1, elOut, eval.ringQ.SubNoModLvl)

	level := utils.MinUint64(utils.MinUint64(el0.Level(), el1.Level()), elOut.Level())

	if el0.Degree() < el1.Degree() {
		for i := el0.Degree() + 1; i < el1.Degree()+1; i++ {
			eval.ringQ.NegLvl(level, elOut.Value()[i], elOut.Value()[i])
		}
	}

}

// SubNew subtracts op1 from op0 and returns the result in a newly created element.
func (eval *evaluator) SubNew(op0, op1 Operand) (ctOut *Ciphertext) {
	ctOut = eval.newCiphertextBinary(op0, op1)
	eval.Sub(op0, op1, ctOut)
	return
}

// SubNoModNew subtracts op1 from op0 without modular reduction, and returns the result in a newly created element.
func (eval *evaluator) SubNoModNew(op0, op1 Operand) (ctOut *Ciphertext) {
	ctOut = eval.newCiphertextBinary(op0, op1)
	eval.SubNoMod(op0, op1, ctOut)
	return
}

func (eval *evaluator) evaluateInPlace(c0, c1, ctOut *Element, evaluate func(uint64, *ring.Poly, *ring.Poly, *ring.Poly)) {

	var tmp0, tmp1 *Element

	level := utils.MinUint64(utils.MinUint64(c0.Level(), c1.Level()), ctOut.Level())

	maxDegree := utils.MaxUint64(c0.Degree(), c1.Degree())
	minDegree := utils.MinUint64(c0.Degree(), c1.Degree())

	// Else resizes the receiver element
	ctOut.Resize(eval.params, maxDegree)
	eval.DropLevel(ctOut.Ciphertext(), ctOut.Level()-utils.MinUint64(c0.Level(), c1.Level()))

	// Checks whether or not the receiver element is the same as one of the input elements
	// and acts accordingly to avoid unnecessary element creation or element overwriting,
	// and scales properly the element before the evaluation.
	if ctOut == c0 {

		if c0.Scale() > c1.Scale() {

			tmp1 = eval.ctxpool.El()

			if uint64(c0.Scale()/c1.Scale()) != 0 {
				eval.MultByConst(c1.Ciphertext(), uint64(c0.Scale()/c1.Scale()), tmp1.Ciphertext())
			}

		} else if c1.Scale() > c0.Scale() {

			if uint64(c1.Scale()/c0.Scale()) != 0 {
				eval.MultByConst(c0.Ciphertext(), uint64(c1.Scale()/c0.Scale()), c0.Ciphertext())
			}

			c0.SetScale(c1.Scale())

			tmp1 = c1

		} else {

			tmp1 = c1
		}

		tmp0 = c0

	} else if ctOut == c1 {

		if c1.Scale() > c0.Scale() {

			tmp0 = eval.ctxpool.El()
			if uint64(c1.Scale()/c0.Scale()) != 0 {
				eval.MultByConst(c0.Ciphertext(), uint64(c1.Scale()/c0.Scale()), tmp0.Ciphertext())
			}

		} else if c0.Scale() > c1.Scale() {

			if uint64(c0.Scale()/c1.Scale()) != 0 {
				eval.MultByConst(c1.Ciphertext(), uint64(c0.Scale()/c1.Scale()), ctOut.Ciphertext())
			}

			ctOut.SetScale(c0.Scale())

			tmp0 = c0

		} else {

			tmp0 = c0
		}

		tmp1 = c1

	} else {

		if c1.Scale() > c0.Scale() {

			tmp0 = eval.ctxpool.El()

			if uint64(c1.Scale()/c0.Scale()) != 0 {
				eval.MultByConst(c0.Ciphertext(), uint64(c1.Scale()/c0.Scale()), tmp0.Ciphertext())
			}

			tmp1 = c1

		} else if c0.Scale() > c1.Scale() {

			tmp1 = eval.ctxpool.El()

			if uint64(c0.Scale()/c1.Scale()) != 0 {
				eval.MultByConst(c1.Ciphertext(), uint64(c0.Scale()/c1.Scale()), tmp1.Ciphertext())
			}

			tmp0 = c0

		} else {
			tmp0 = c0
			tmp1 = c1
		}
	}

	for i := uint64(0); i < minDegree+1; i++ {
		evaluate(level, tmp0.Value()[i], tmp1.Value()[i], ctOut.Value()[i])
	}

	ctOut.SetScale(utils.MaxFloat64(c0.Scale(), c1.Scale()))

	// If the inputs degrees differ, it copies the remaining degree on the receiver.
	// Also checks that the receiver is not one of the inputs to avoid unnecessary work.

	if c0.Degree() > c1.Degree() && tmp0 != ctOut {
		for i := minDegree + 1; i < maxDegree+1; i++ {
			eval.ringQ.CopyLvl(level, tmp0.Value()[i], ctOut.Value()[i])
		}
	} else if c1.Degree() > c0.Degree() && tmp1 != ctOut {
		for i := minDegree + 1; i < maxDegree+1; i++ {
			eval.ringQ.CopyLvl(level, tmp1.Value()[i], ctOut.Value()[i])
		}
	}
}

// Neg negates the value of ct0 and returns the result in ctOut.
func (eval *evaluator) Neg(ct0 *Ciphertext, ctOut *Ciphertext) {

	level := utils.MinUint64(ct0.Level(), ctOut.Level())

	if ct0.Degree() != ctOut.Degree() {
		panic("cannot Negate: invalid receiver Ciphertext does not match input Ciphertext degree")
	}

	for i := range ct0.value {
		eval.ringQ.NegLvl(level, ct0.value[i], ctOut.Value()[i])
	}
}

// NegNew negates ct0 and returns the result in a newly created element.
func (eval *evaluator) NegNew(ct0 *Ciphertext) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.Scale())
	eval.Neg(ct0, ctOut)
	return
}

// AddConstNew adds the input constant (which can be a uint64, int64, float64 or complex128) to ct0 and returns the result in a new element.
func (eval *evaluator) AddConstNew(ct0 *Ciphertext, constant interface{}) (ctOut *Ciphertext) {
	ctOut = ct0.CopyNew().Ciphertext()
	eval.AddConst(ct0, constant, ctOut)
	return ctOut
}

func (eval *evaluator) getConstAndScale(level uint64, constant interface{}) (cReal, cImag, scale float64) {

	// Converts to float64 and determines if a scaling is required (which is the case if either real or imag have a rational part)
	scale = 1
	switch constant := constant.(type) {
	case complex128:
		cReal = real(constant)
		cImag = imag(constant)

		if cReal != 0 {
			valueInt := int64(cReal)
			valueFloat := cReal - float64(valueInt)

			if valueFloat != 0 {
				scale = float64(eval.ringQ.Modulus[level])
			}
		}

		if cImag != 0 {
			valueInt := int64(cImag)
			valueFloat := cImag - float64(valueInt)

			if valueFloat != 0 {
				scale = float64(eval.ringQ.Modulus[level])
			}
		}

	case float64:
		cReal = constant
		cImag = float64(0)

		if cReal != 0 {
			valueInt := int64(cReal)
			valueFloat := cReal - float64(valueInt)

			if valueFloat != 0 {
				scale = float64(eval.ringQ.Modulus[level])
			}
		}

	case uint64:
		cReal = float64(constant)
		cImag = float64(0)

	case int64:
		cReal = float64(constant)
		cImag = float64(0)

	case int:
		cReal = float64(constant)
		cImag = float64(0)
	}

	return
}

// AddConst adds the input constant (which can be a uint64, int64, float64 or complex128) to ct0 and returns the result in ctOut.
func (eval *evaluator) AddConst(ct0 *Ciphertext, constant interface{}, ctOut *Ciphertext) {

	var level = utils.MinUint64(ct0.Level(), ctOut.Level())
	var scaledConst, scaledConstReal, scaledConstImag uint64

	ringQ := eval.ringQ

	if ct0 != ctOut {
		if ct0.Degree() > ctOut.Degree() {
			panic("receiver operand degree is too small")
		} else if ctOut.Degree() > ct0.Degree() {
			tmp := ctOut.value
			tmp = tmp[:ct0.Degree()+1]
		}
	}

	for i := range ct0.Value()[1:] {
		ringQ.CopyLvl(level, ct0.Value()[i+1], ctOut.Value()[i+1])
	}

	cReal, cImag, _ := eval.getConstAndScale(level, constant)

	// Component wise addition of the following vector to the ciphertext:
	// [a + b*psi_qi^2, ....., a + b*psi_qi^2, a - b*psi_qi^2, ...., a - b*psi_qi^2] mod Qi
	// [{                  N/2                }{                N/2               }]
	// Which is equivalent outside of the NTT domain to adding a to the first coefficient of ct0 and b to the N/2-th coefficient of ct0.
	var qi uint64
	for i := uint64(0); i < level+1; i++ {
		scaledConstReal = 0
		scaledConstImag = 0
		scaledConst = 0

		qi = ringQ.Modulus[i]

		if cReal != 0 {
			scaledConstReal = scaleUpExact(cReal, ctOut.Scale(), qi)
			scaledConst = scaledConstReal
		}

		if cImag != 0 {
			scaledConstImag = ring.FastBRed(scaleUpExact(cImag, ctOut.Scale(), qi), ringQ.NttPsi[i][1], qi)
			scaledConst = ring.CRed(scaledConst+scaledConstImag, qi)
		}

		p0tmp := ct0.value[0].Coeffs[i]
		p1tmp := ctOut.Value()[0].Coeffs[i]

		for j := uint64(0); j < ringQ.N>>1; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))

			z[0] = ring.CRed(x[0]+scaledConst, qi)
			z[1] = ring.CRed(x[1]+scaledConst, qi)
			z[2] = ring.CRed(x[2]+scaledConst, qi)
			z[3] = ring.CRed(x[3]+scaledConst, qi)
			z[4] = ring.CRed(x[4]+scaledConst, qi)
			z[5] = ring.CRed(x[5]+scaledConst, qi)
			z[6] = ring.CRed(x[6]+scaledConst, qi)
			z[7] = ring.CRed(x[7]+scaledConst, qi)
		}

		if cImag != 0 {
			scaledConst = ring.CRed(scaledConstReal+(qi-scaledConstImag), qi)
		}

		for j := ringQ.N >> 1; j < ringQ.N; j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))

			z[0] = ring.CRed(x[0]+scaledConst, qi)
			z[1] = ring.CRed(x[1]+scaledConst, qi)
			z[2] = ring.CRed(x[2]+scaledConst, qi)
			z[3] = ring.CRed(x[3]+scaledConst, qi)
			z[4] = ring.CRed(x[4]+scaledConst, qi)
			z[5] = ring.CRed(x[5]+scaledConst, qi)
			z[6] = ring.CRed(x[6]+scaledConst, qi)
			z[7] = ring.CRed(x[7]+scaledConst, qi)
		}
	}

	for i := range ctOut.Value() {
		ctOut.Value()[i].Coeffs = ctOut.Value()[i].Coeffs[:level+1]
	}
}

// MultByConstAndAdd multiplies ct0 by the input constant, and adds it to the receiver element (it does not modify the input
// element), e.g., ctOut(x) = ctOut(x) + ct0(x) * (a+bi). This functions removes the need of storing the intermediate value c(x) * (a+bi).
// This function will modify the level and the scale of the receiver element depending on the level and the scale of the input
// element and the type of the constant. The level of the receiver element will be set to min(input.level, receiver.level).
// The scale of the receiver element will be set to the scale that the input element would have after the multiplication by the constant.
func (eval *evaluator) MultByConstAndAdd(ct0 *Ciphertext, constant interface{}, ctOut *Ciphertext) {

	var level = utils.MinUint64(ct0.Level(), ctOut.Level())

	// Forces a drop of ctOut level to ct0 level
	if ctOut.Level() > level {
		eval.DropLevel(ctOut, ctOut.Level()-level)
	}

	cReal, cImag, scale := eval.getConstAndScale(level, constant)

	var scaledConst, scaledConstReal, scaledConstImag uint64

	ringQ := eval.ringQ

	// If a scaling would be required to multiply by the constant,
	// it equalizes scales such that the scales match in the end.
	if scale != 1 {

		// If ctOut scaling is smaller than ct0's scale + the default scaling,
		// then brings ctOut scale to ct0's scale.
		if ctOut.Scale() < ct0.Scale()*scale {

			if uint64((scale*ct0.Scale())/ctOut.Scale()) != 0 {

				eval.MultByConst(ctOut, uint64((scale*ct0.Scale())/ctOut.Scale()), ctOut)

			}

			ctOut.SetScale(scale * ct0.Scale())

			// If ctOut.Scale() > ((a+bi)*scale)*ct0(x), then it sets the scale to
			// bring c(x)*scale to the level of ctOut(x) scale
		} else if ctOut.Scale() > ct0.Scale()*scale {
			scale = ctOut.Scale() / ct0.Scale()
		}

		// If no scaling is required, then it sets the appropriate scale such that
		// ct0(x)*scale matches ctOut(x) scale without modifying ct0(x) scale.
	} else {

		if ctOut.Scale() > ct0.Scale() {

			scale = ctOut.Scale() / ct0.Scale()

		} else if ct0.Scale() > ctOut.Scale() {

			if uint64(ct0.Scale()/ctOut.Scale()) != 0 {
				eval.MultByConst(ctOut, ct0.Scale()/ctOut.Scale(), ctOut)
			}

			ctOut.SetScale(ct0.Scale())
		}
	}

	// Component-wise multiplication of the following vector to the ciphertext:
	// [a + b*psi_qi^2, ....., a + b*psi_qi^2, a - b*psi_qi^2, ...., a - b*psi_qi^2] mod Qi
	// [{                  N/2                }{                N/2               }]
	// Which is equivalent outside of the NTT domain to adding a to the first coefficient of ct0 and b to the N/2-th coefficient of ct0.
	var scaledConstBred ring.FastBRedOperand
	for i := uint64(0); i < level+1; i++ {

		qi := ringQ.Modulus[i]

		scaledConstReal = 0
		scaledConstImag = 0
		scaledConst = 0

		if cReal != 0 {
			scaledConstReal = scaleUpExact(cReal, scale, qi)
			scaledConst = scaledConstReal
		}

		if cImag != 0 {
			scaledConstImag = scaleUpExact(cImag, scale, qi)
			scaledConstImag = ring.FastBRed(scaledConstImag, ringQ.NttPsi[i][1], qi)
			scaledConst = ring.CRed(scaledConst+scaledConstImag, qi)
		}

		scaledConstBred = ring.NewFastBRedOperand(scaledConst, qi)

		for u := range ct0.Value() {
			p0tmp := ct0.Value()[u].Coeffs[i]
			p1tmp := ctOut.Value()[u].Coeffs[i]

			for j := uint64(0); j < ringQ.N>>1; j = j + 8 {

				x := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))
				z := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))

				z[0] = ring.CRed(z[0]+ring.FastBRed(x[0], scaledConstBred, qi), qi)
				z[1] = ring.CRed(z[1]+ring.FastBRed(x[1], scaledConstBred, qi), qi)
				z[2] = ring.CRed(z[2]+ring.FastBRed(x[2], scaledConstBred, qi), qi)
				z[3] = ring.CRed(z[3]+ring.FastBRed(x[3], scaledConstBred, qi), qi)
				z[4] = ring.CRed(z[4]+ring.FastBRed(x[4], scaledConstBred, qi), qi)
				z[5] = ring.CRed(z[5]+ring.FastBRed(x[5], scaledConstBred, qi), qi)
				z[6] = ring.CRed(z[6]+ring.FastBRed(x[6], scaledConstBred, qi), qi)
				z[7] = ring.CRed(z[7]+ring.FastBRed(x[7], scaledConstBred, qi), qi)
			}
		}

		if cImag != 0 {
			scaledConst = ring.CRed(scaledConstReal+(qi-scaledConstImag), qi)
			scaledConstBred = ring.NewFastBRedOperand(scaledConst, qi)
		}

		for u := range ct0.Value() {
			p0tmp := ct0.Value()[u].Coeffs[i]
			p1tmp := ctOut.Value()[u].Coeffs[i]
			for j := ringQ.N >> 1; j < ringQ.N; j = j + 8 {

				x := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))
				z := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))

				z[0] = ring.CRed(z[0]+ring.FastBRed(x[0], scaledConstBred, qi), qi)
				z[1] = ring.CRed(z[1]+ring.FastBRed(x[1], scaledConstBred, qi), qi)
				z[2] = ring.CRed(z[2]+ring.FastBRed(x[2], scaledConstBred, qi), qi)
				z[3] = ring.CRed(z[3]+ring.FastBRed(x[3], scaledConstBred, qi), qi)
				z[4] = ring.CRed(z[4]+ring.FastBRed(x[4], scaledConstBred, qi), qi)
				z[5] = ring.CRed(z[5]+ring.FastBRed(x[5], scaledConstBred, qi), qi)
				z[6] = ring.CRed(z[6]+ring.FastBRed(x[6], scaledConstBred, qi), qi)
				z[7] = ring.CRed(z[7]+ring.FastBRed(x[7], scaledConstBred, qi), qi)
			}
		}
	}
}

// MultByConstNew multiplies ct0 by the input constant and returns the result in a newly created element.
// The scale of the output element will depend on the scale of the input element and the constant (if the constant
// needs to be scaled (its rational part is not zero)). The constant can be a uint64, int64, float64 or complex128.
func (eval *evaluator) MultByConstNew(ct0 *Ciphertext, constant interface{}) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.Scale())
	eval.MultByConst(ct0, constant, ctOut)
	return
}

// MultByConst multiplies ct0 by the input constant and returns the result in ctOut.
// The scale of the output element will depend on the scale of the input element and the constant (if the constant
// needs to be scaled (its rational part is not zero)). The constant can be a uint64, int64, float64 or complex128.
func (eval *evaluator) MultByConst(ct0 *Ciphertext, constant interface{}, ctOut *Ciphertext) {

	var level = utils.MinUint64(ct0.Level(), ctOut.Level())

	cReal, cImag, scale := eval.getConstAndScale(level, constant)

	// Component wise multiplication of the following vector with the ciphertext:
	// [a + b*psi_qi^2, ....., a + b*psi_qi^2, a - b*psi_qi^2, ...., a - b*psi_qi^2] mod Qi
	// [{                  N/2                }{                N/2               }]
	// Which is equivalent outside of the NTT domain to adding a to the first coefficient of ct0 and b to the N/2-th coefficient of ct0.
	ringQ := eval.ringQ
	var scaledConst, scaledConstReal, scaledConstImag uint64
	var scaledConstBred ring.FastBRedOperand
	for i := uint64(0); i < level+1; i++ {

		qi := ringQ.Modulus[i]

		scaledConstReal = 0
		scaledConstImag = 0
		scaledConst = 0

		if cReal != 0 {
			scaledConstReal = scaleUpExact(cReal, scale, qi)
			scaledConst = scaledConstReal
		}

		if cImag != 0 {
			scaledConstImag = scaleUpExact(cImag, scale, qi)
			scaledConstImag = ring.FastBRed(scaledConstImag, ringQ.NttPsi[i][1], qi)
			scaledConst = ring.CRed(scaledConst+scaledConstImag, qi)
		}

		scaledConstBred = ring.NewFastBRedOperand(scaledConst, qi)

		for u := range ct0.Value() {
			p0tmp := ct0.Value()[u].Coeffs[i]
			p1tmp := ctOut.Value()[u].Coeffs[i]

			for j := uint64(0); j < ringQ.N>>1; j = j + 8 {

				x := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))
				z := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))

				z[0] = ring.FastBRed(x[0], scaledConstBred, qi)
				z[1] = ring.FastBRed(x[1], scaledConstBred, qi)
				z[2] = ring.FastBRed(x[2], scaledConstBred, qi)
				z[3] = ring.FastBRed(x[3], scaledConstBred, qi)
				z[4] = ring.FastBRed(x[4], scaledConstBred, qi)
				z[5] = ring.FastBRed(x[5], scaledConstBred, qi)
				z[6] = ring.FastBRed(x[6], scaledConstBred, qi)
				z[7] = ring.FastBRed(x[7], scaledConstBred, qi)
			}
		}

		if cImag != 0 {
			scaledConst = ring.CRed(scaledConstReal+(qi-scaledConstImag), qi)
			scaledConstBred = ring.NewFastBRedOperand(scaledConst, qi)
		}

		for u := range ct0.Value() {
			p0tmp := ct0.Value()[u].Coeffs[i]
			p1tmp := ctOut.Value()[u].Coeffs[i]
			for j := ringQ.N >> 1; j < ringQ.N; j = j + 8 {

				x := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))
				z := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))

				z[0] = ring.FastBRed(x[0], scaledConstBred, qi)
				z[1] = ring.FastBRed(x[1], scaledConstBred, qi)
				z[2] = ring.FastBRed(x[2], scaledConstBred, qi)
				z[3] = ring.FastBRed(x[3], scaledConstBred, qi)
				z[4] = ring.FastBRed(x[4], scaledConstBred, qi)
				z[5] = ring.FastBRed(x[5], scaledConstBred, qi)
				z[6] = ring.FastBRed(x[6], scaledConstBred, qi)
				z[7] = ring.FastBRed(x[7], scaledConstBred, qi)
			}
		}
	}

	ctOut.SetScale(ct0.Scale() * scale)
}

func (eval *evaluator) MultByGaussianInteger(ct0 *Ciphertext, cReal, cImag int64, ctOut *Ciphertext) {

	ringQ := eval.ringQ

	level := utils.MinUint64(ct0.Level(), ctOut.Level())
	var scaledConst, scaledConstReal, scaledConstImag uint64
	var scaledConstBred ring.FastBRedOperand
	for i := uint64(0); i < level+1; i++ {

		qi := ringQ.Modulus[i]

		scaledConstReal = 0
		scaledConstImag = 0
		scaledConst = 0

		if cReal != 0 {
			if cReal < 0 {
				scaledConstReal = uint64(int64(qi) + cReal%int64(qi))
			} else {
				scaledConstReal = uint64(cReal)
			}
			scaledConst = scaledConstReal
		}

		if cImag != 0 {
			if cImag < 0 {
				scaledConstImag = uint64(int64(qi) + cImag%int64(qi))
			} else {
				scaledConstImag = uint64(cImag)
			}
			scaledConstImag = ring.FastBRed(scaledConstImag, ringQ.NttPsi[i][1], qi)
			scaledConst = ring.CRed(scaledConst+scaledConstImag, qi)
		}

		scaledConstBred = ring.NewFastBRedOperand(ring.BRedAdd(scaledConst, qi, ringQ.BredParams[i]), qi)

		for u := range ct0.Value() {
			p0tmp := ct0.Value()[u].Coeffs[i]
			p1tmp := ctOut.Value()[u].Coeffs[i]

			for j := uint64(0); j < ringQ.N>>1; j = j + 8 {

				x := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))
				z := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))

				z[0] = ring.FastBRed(x[0], scaledConstBred, qi)
				z[1] = ring.FastBRed(x[1], scaledConstBred, qi)
				z[2] = ring.FastBRed(x[2], scaledConstBred, qi)
				z[3] = ring.FastBRed(x[3], scaledConstBred, qi)
				z[4] = ring.FastBRed(x[4], scaledConstBred, qi)
				z[5] = ring.FastBRed(x[5], scaledConstBred, qi)
				z[6] = ring.FastBRed(x[6], scaledConstBred, qi)
				z[7] = ring.FastBRed(x[7], scaledConstBred, qi)
			}
		}

		if cImag != 0 {
			scaledConst = ring.BRedAdd(scaledConstReal+(qi-scaledConstImag), qi, ringQ.BredParams[i])
			scaledConstBred = ring.NewFastBRedOperand(scaledConst, qi)
		}

		for u := range ct0.Value() {
			p0tmp := ct0.Value()[u].Coeffs[i]
			p1tmp := ctOut.Value()[u].Coeffs[i]

			for j := ringQ.N >> 1; j < ringQ.N; j = j + 8 {

				x := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))
				z := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))

				z[0] = ring.FastBRed(x[0], scaledConstBred, qi)
				z[1] = ring.FastBRed(x[1], scaledConstBred, qi)
				z[2] = ring.FastBRed(x[2], scaledConstBred, qi)
				z[3] = ring.FastBRed(x[3], scaledConstBred, qi)
				z[4] = ring.FastBRed(x[4], scaledConstBred, qi)
				z[5] = ring.FastBRed(x[5], scaledConstBred, qi)
				z[6] = ring.FastBRed(x[6], scaledConstBred, qi)
				z[7] = ring.FastBRed(x[7], scaledConstBred, qi)
			}
		}
	}
}

func (eval *evaluator) MultByGaussianIntegerAndAdd(ct0 *Ciphertext, cReal, cImag int64, ctOut *Ciphertext) {

	ringQ := eval.ringQ

	level := utils.MinUint64(ct0.Level(), ctOut.Level())
	var scaledConst, scaledConstReal, scaledConstImag uint64
	var scaledConstBred ring.FastBRedOperand
	for i := uint64(0); i < level+1; i++ {

		qi := ringQ.Modulus[i]

		scaledConstReal = 0
		scaledConstImag = 0
		scaledConst = 0

		if cReal != 0 {
			if cReal < 0 {
				scaledConstReal = uint64(int64(qi) + cReal%int64(qi))
			} else {
				scaledConstReal = uint64(cReal)
			}
			scaledConst = scaledConstReal
		}

		if cImag != 0 {
			if cImag < 0 {
				scaledConstImag = uint64(int64(qi) + cImag%int64(qi))
			} else {
				scaledConstImag = uint64(cImag)
			}
			scaledConstImag = ring.FastBRed(scaledConstImag, ringQ.NttPsi[i][1], qi)
			scaledConst += scaledConstImag
		}

		scaledConstBred = ring.NewFastBRedOperand(ring.BRedAdd(scaledConst, qi, ringQ.BredParams[i]), qi)

		for u := range ct0.Value() {
			p0tmp := ct0.Value()[u].Coeffs[i]
			p1tmp := ctOut.Value()[u].Coeffs[i]

			for j := uint64(0); j < ringQ.N>>1; j = j + 8 {

				x := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))
				z := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))

				z[0] = ring.CRed(z[0]+ring.FastBRed(x[0], scaledConstBred, qi), qi)
				z[1] = ring.CRed(z[1]+ring.FastBRed(x[1], scaledConstBred, qi), qi)
				z[2] = ring.CRed(z[2]+ring.FastBRed(x[2], scaledConstBred, qi), qi)
				z[3] = ring.CRed(z[3]+ring.FastBRed(x[3], scaledConstBred, qi), qi)
				z[4] = ring.CRed(z[4]+ring.FastBRed(x[4], scaledConstBred, qi), qi)
				z[5] = ring.CRed(z[5]+ring.FastBRed(x[5], scaledConstBred, qi), qi)
				z[6] = ring.CRed(z[6]+ring.FastBRed(x[6], scaledConstBred, qi), qi)
				z[7] = ring.CRed(z[7]+ring.FastBRed(x[7], scaledConstBred, qi), qi)
			}
		}

		if cImag != 0 {
			scaledConst = ring.BRedAdd(scaledConstReal+(qi-scaledConstImag), qi, ringQ.BredParams[i])
			scaledConstBred = ring.NewFastBRedOperand(scaledConst, qi)
		}

		for u := range ct0.Value() {
			p0tmp := ct0.Value()[u].Coeffs[i]
			p1tmp := ctOut.Value()[u].Coeffs[i]

			for j := ringQ.N >> 1; j < ringQ.N; j = j + 8 {

				x := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))
				z := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))

				z[0] = ring.CRed(z[0]+ring.FastBRed(x[0], scaledConstBred, qi), qi)
				z[1] = ring.CRed(z[1]+ring.FastBRed(x[1], scaledConstBred, qi), qi)
				z[2] = ring.CRed(z[2]+ring.FastBRed(x[2], scaledConstBred, qi), qi)
				z[3] = ring.CRed(z[3]+ring.FastBRed(x[3], scaledConstBred, qi), qi)
				z[4] = ring.CRed(z[4]+ring.FastBRed(x[4], scaledConstBred, qi), qi)
				z[5] = ring.CRed(z[5]+ring.FastBRed(x[5], scaledConstBred, qi), qi)
				z[6] = ring.CRed(z[6]+ring.FastBRed(x[6], scaledConstBred, qi), qi)
				z[7] = ring.CRed(z[7]+ring.FastBRed(x[7], scaledConstBred, qi), qi)
			}
		}
	}
}

// MultByiNew multiplies ct0 by the imaginary number i, and returns the result in a newly created element.
// It does not change the scale.
func (eval *evaluator) MultByiNew(ct0 *Ciphertext) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, ct0.Level(), ct0.Scale())
	eval.MultByi(ct0, ctOut)
	return ctOut
}

// MultByi multiplies ct0 by the imaginary number i, and returns the result in ctOut.
// It does not change the scale.
func (eval *evaluator) MultByi(ct0 *Ciphertext, ctOut *Ciphertext) {

	var level = utils.MinUint64(ct0.Level(), ctOut.Level())

	ringQ := eval.ringQ

	var imag ring.FastBRedOperand

	// Equivalent to a product by the monomial x^(n/2) outside of the NTT domain
	for i := uint64(0); i < level+1; i++ {

		qi := ringQ.Modulus[i]

		imag = ringQ.NttPsi[i][1] // Psi^2

		for u := range ctOut.value {
			p0tmp := ct0.value[u].Coeffs[i]
			p1tmp := ctOut.value[u].Coeffs[i]

			for j := uint64(0); j < ringQ.N>>1; j = j + 8 {

				x := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))
				z := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))

				z[0] = ring.FastBRed(x[0], imag, qi)
				z[1] = ring.FastBRed(x[1], imag, qi)
				z[2] = ring.FastBRed(x[2], imag, qi)
				z[3] = ring.FastBRed(x[3], imag, qi)
				z[4] = ring.FastBRed(x[4], imag, qi)
				z[5] = ring.FastBRed(x[5], imag, qi)
				z[6] = ring.FastBRed(x[6], imag, qi)
				z[7] = ring.FastBRed(x[7], imag, qi)
			}
		}

		imag = ring.NewFastBRedOperand(qi-imag.Operand, qi)

		for u := range ctOut.value {
			p0tmp := ct0.value[u].Coeffs[i]
			p1tmp := ctOut.value[u].Coeffs[i]
			for j := ringQ.N >> 1; j < ringQ.N; j = j + 8 {

				x := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))
				z := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))

				z[0] = ring.FastBRed(x[0], imag, qi)
				z[1] = ring.FastBRed(x[1], imag, qi)
				z[2] = ring.FastBRed(x[2], imag, qi)
				z[3] = ring.FastBRed(x[3], imag, qi)
				z[4] = ring.FastBRed(x[4], imag, qi)
				z[5] = ring.FastBRed(x[5], imag, qi)
				z[6] = ring.FastBRed(x[6], imag, qi)
				z[7] = ring.FastBRed(x[7], imag, qi)

			}
		}
	}
}

// DivByiNew multiplies ct0 by the imaginary number 1/i = -i, and returns the result in a newly created element.
// It does not change the scale.
func (eval *evaluator) DivByiNew(ct0 *Ciphertext) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, ct0.Level(), ct0.Scale())
	eval.DivByi(ct0, ctOut)
	return
}

// DivByi multiplies ct0 by the imaginary number 1/i = -i, and returns the result in ctOut.
// It does not change the scale.
func (eval *evaluator) DivByi(ct0 *Ciphertext, ctOut *Ciphertext) {

	var level = utils.MinUint64(ct0.Level(), ctOut.Level())

	ringQ := eval.ringQ

	var imag ring.FastBRedOperand

	// Equivalent to a product by the monomial x^(3*n/2) outside of the NTT domain
	for i := uint64(0); i < level+1; i++ {

		qi := ringQ.Modulus[i]

		imag = ring.NewFastBRedOperand(qi-ringQ.NttPsi[i][1].Operand, qi) // -Psi^2

		for u := range ctOut.value {
			p0tmp := ct0.value[u].Coeffs[i]
			p1tmp := ctOut.value[u].Coeffs[i]
			for j := uint64(0); j < ringQ.N>>1; j = j + 8 {

				x := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))
				z := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))

				z[0] = ring.FastBRed(x[0], imag, qi)
				z[1] = ring.FastBRed(x[1], imag, qi)
				z[2] = ring.FastBRed(x[2], imag, qi)
				z[3] = ring.FastBRed(x[3], imag, qi)
				z[4] = ring.FastBRed(x[4], imag, qi)
				z[5] = ring.FastBRed(x[5], imag, qi)
				z[6] = ring.FastBRed(x[6], imag, qi)
				z[7] = ring.FastBRed(x[7], imag, qi)
			}
		}

		imag = ringQ.NttPsi[i][1] // Psi^2

		for u := range ctOut.value {
			p0tmp := ct0.value[u].Coeffs[i]
			p1tmp := ctOut.value[u].Coeffs[i]
			for j := ringQ.N >> 1; j < ringQ.N; j = j + 8 {

				x := (*[8]uint64)(unsafe.Pointer(&p0tmp[j]))
				z := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))

				z[0] = ring.FastBRed(x[0], imag, qi)
				z[1] = ring.FastBRed(x[1], imag, qi)
				z[2] = ring.FastBRed(x[2], imag, qi)
				z[3] = ring.FastBRed(x[3], imag, qi)
				z[4] = ring.FastBRed(x[4], imag, qi)
				z[5] = ring.FastBRed(x[5], imag, qi)
				z[6] = ring.FastBRed(x[6], imag, qi)
				z[7] = ring.FastBRed(x[7], imag, qi)
			}
		}
	}
}

// ScaleUpNew multiplies ct0 by 2^scale and sets its scale to its previous scale
// plus 2^n. It returns the result in a newly created element.
func (eval *evaluator) ScaleUpNew(ct0 *Ciphertext, scale float64) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.Scale())
	eval.ScaleUp(ct0, scale, ctOut)
	return
}

// ScaleUp multiplies ct0 by 2^scale and sets its scale to its previous scale
// plus 2^n. It returns the result in ctOut.
func (eval *evaluator) ScaleUp(ct0 *Ciphertext, scale float64, ctOut *Ciphertext) {
	eval.MultByConst(ct0, uint64(scale), ctOut)
	ctOut.SetScale(ct0.Scale() * scale)
}

// SetScale sets the scale of the ciphertext to the input scale (consumes a level)
func (eval *evaluator) SetScale(ct *Ciphertext, scale float64) {

	var tmp = eval.params.scale

	eval.scale = scale

	eval.MultByConst(ct, scale/ct.Scale(), ct)

	if err := eval.Rescale(ct, scale, ct); err != nil {
		panic(err)
	}

	ct.SetScale(scale)

	eval.scale = tmp
}

// MulByPow2New multiplies ct0 by 2^pow2 and returns the result in a newly created element.
func (eval *evaluator) MulByPow2New(ct0 *Ciphertext, pow2 uint64) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.Scale())
	eval.MulByPow2(ct0.El(), pow2, ctOut.El())
	return
}

// MulByPow2 multiplies ct0 by 2^pow2 and returns the result in ctOut.
func (eval *evaluator) MulByPow2(ct0 *Element, pow2 uint64, ctOut *Element) {
	var level = utils.MinUint64(ct0.Level(), ctOut.Level())
	for i := range ctOut.Value() {
		eval.ringQ.MulByPow2Lvl(level, ct0.value[i], pow2, ctOut.Value()[i])
	}
}

// ReduceNew applies a modular reduction to ct0 and returns the result in a newly created element.
// To be used in conjunction with functions that do not apply modular reduction.
func (eval *evaluator) ReduceNew(ct0 *Ciphertext) (ctOut *Ciphertext) {

	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.Scale())

	_ = eval.Reduce(ct0, ctOut)

	return ctOut
}

// Reduce applies a modular reduction to ct0 and returns the result in ctOut.
// To be used in conjunction with functions that do not apply modular reduction.
func (eval *evaluator) Reduce(ct0 *Ciphertext, ctOut *Ciphertext) error {

	if ct0.Degree() != ctOut.Degree() {
		return errors.New("cannot Reduce: degrees of receiver Ciphertext and input Ciphertext do not match")
	}

	for i := range ct0.value {
		eval.ringQ.ReduceLvl(utils.MinUint64(ct0.Level(), ctOut.Level()), ct0.value[i], ctOut.value[i])
	}

	return nil
}

// DropLevelNew reduces the level of ct0 by levels and returns the result in a newly created element.
// No rescaling is applied during this procedure.
func (eval *evaluator) DropLevelNew(ct0 *Ciphertext, levels uint64) (ctOut *Ciphertext) {
	ctOut = ct0.CopyNew().Ciphertext()
	eval.DropLevel(ctOut, levels)
	return
}

// DropLevel reduces the level of ct0 by levels and returns the result in ct0.
// No rescaling is applied during this procedure.
func (eval *evaluator) DropLevel(ct0 *Ciphertext, levels uint64) {
	level := ct0.Level()
	for i := range ct0.value {
		ct0.value[i].Coeffs = ct0.value[i].Coeffs[:level+1-levels]
	}
}

// RescaleNew divides ct0 by the last modulus in the moduli chain, and repeats this
// procedure (consuming one level each time) until the scale reaches the original scale or before it goes below it, and returns the result
// in a newly created element. Since all the moduli in the moduli chain are generated to be close to the
// original scale, this procedure is equivalent to dividing the input element by the scale and adding
// some error.
// Returns an error if "threshold <= 0", ct.Scale() = 0, ct.Level() = 0, ct.IsNTT() != true
func (eval *evaluator) RescaleNew(ct0 *Ciphertext, threshold float64) (ctOut *Ciphertext, err error) {

	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.Scale())

	return ctOut, eval.Rescale(ct0, threshold, ctOut)
}

// Rescale divides ct0 by the last modulus in the moduli chain, and repeats this
// procedure (consuming one level each time) until the scale reaches the original scale or before it goes below it, and returns the result
// in ctOut. Since all the moduli in the moduli chain are generated to be close to the
// original scale, this procedure is equivalent to dividing the input element by the scale and adding
// some error.
// Returns an error if "threshold <= 0", ct.Scale() = 0, ct.Level() = 0, ct.IsNTT() != true or if ct.Leve() != ctOut.Level()
func (eval *evaluator) Rescale(ct0 *Ciphertext, threshold float64, ctOut *Ciphertext) (err error) {

	ringQ := eval.ringQ

	if threshold <= 0 {
		return errors.New("cannot Rescale: threshold is 0")
	}

	if ct0.Scale() == 0 {
		return errors.New("cannot Rescale: ciphertext scale is 0")
	}

	if ct0.Level() == 0 {
		return errors.New("cannot Rescale: input Ciphertext already at level 0")
	}

	if ct0.Level() != ctOut.Level() {
		panic("cannot Rescale: degrees of receiver Ciphertext and input Ciphertext do not match")
	}

	if ct0.Scale() >= (threshold*float64(ringQ.Modulus[ctOut.Level()]))/2 {

		if !ct0.IsNTT() {
			panic("cannot Rescale: input Ciphertext not in NTT")
		}

		ctOut.Copy(ct0.El())

		for ctOut.Scale() >= (threshold*float64(ringQ.Modulus[ctOut.Level()]))/2 && ctOut.Level() != 0 {

			ctOut.DivScale(float64(ringQ.Modulus[ctOut.Level()]))

			for i := range ctOut.Value() {
				eval.ringQ.DivRoundByLastModulusNTT(ctOut.Value()[i])
			}

		}

	} else {
		ctOut.Copy(ct0.El())
	}

	return nil
}

// RescaleMany applies Rescale several times in a row on the input Ciphertext.
func (eval *evaluator) RescaleMany(ct0 *Ciphertext, nbRescales uint64, ctOut *Ciphertext) (err error) {

	if ct0.Level() < nbRescales {
		return errors.New("cannot RescaleMany: input Ciphertext level too low")
	}

	if ct0.Level() != ctOut.Level() {
		panic("cannot RescaleMany: degrees of receiver Ciphertext and input Ciphertext do not match")
	}

	if !ct0.IsNTT() {
		panic("cannot RescaleMany: input Ciphertext not in NTT")
	}

	ctOut.Copy(ct0.El())

	for i := uint64(0); i < nbRescales; i++ {
		ctOut.DivScale(float64(eval.ringQ.Modulus[ctOut.Level()-i]))
	}

	for i := range ctOut.Value() {
		eval.ringQ.DivRoundByLastModulusManyNTT(ctOut.Value()[i], nbRescales)
	}

	return nil
}

// MulRelinNew multiplies ct0 by ct1 and returns the result in a newly created element. The new scale is
// the multiplication between the scales of the input elements (addition when the scale is represented in log2). An evaluation
// key can be provided to apply a relinearization step to reduce the degree of the output element. This evaluation key is only
// required when the two input elements are Ciphertexts. If no evaluation key is provided and the input elements are two Ciphertexts,
// the resulting Ciphertext will be of degree two. This function only accepts Plaintexts (degree zero) and/or Ciphertexts of degree one.
func (eval *evaluator) MulRelinNew(op0, op1 Operand, evakey *EvaluationKey) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, utils.MinUint64(op0.Level(), op1.Level()), op0.Scale()+op1.Scale())
	eval.MulRelin(op0, op1, evakey, ctOut)

	return ctOut
}

// MulRelin multiplies ct0 by ct1 and returns the result in ctOut. The new scale is
// the multiplication between the scales of the input elements (addition when the scale is represented in log2). An evaluation
// key can be provided to apply a relinearization step to reduce the degree of the output element. This evaluation key is only
// required when the two input elements are Ciphertexts. If no evaluation key is provided and the input elements are two Ciphertexts,
// the resulting Ciphertext will be of degree two. This function only accepts Plaintexts (degree zero) and/or Ciphertexts of degree one.
func (eval *evaluator) MulRelin(op0, op1 Operand, evakey *EvaluationKey, ctOut *Ciphertext) {

	el0, el1, elOut := eval.getElemAndCheckBinary(op0, op1, ctOut, utils.MaxUint64(op0.Degree(), op1.Degree()))

	level := utils.MinUint64(utils.MinUint64(el0.Level(), el1.Level()), elOut.Level())

	if ctOut.Level() > level {
		eval.DropLevel(elOut.Ciphertext(), elOut.Level()-level)
	}

	if el0.Degree() > 1 || el1.Degree() > 1 {
		panic("cannot MulRelin: input elements must be of degree 0 or 1")
	}

	if !el0.IsNTT() {
		panic("cannot MulRelin: op0 must be in NTT")
	}

	if !el1.IsNTT() {
		panic("cannot MulRelin: op1 must be in NTT")
	}

	elOut.SetScale(el0.Scale() * el1.Scale())

	ringQ := eval.ringQ

	var c00, c01, c0, c1, c2 *ring.Poly

	// Case Ciphertext (x) Ciphertext
	if el0.Degree()+el1.Degree() == 2 {

		c00 = eval.poolQMul[0]
		c01 = eval.poolQMul[1]

		c0 = elOut.value[0]
		c1 = elOut.value[1]

		if evakey == nil {
			elOut.Resize(eval.params, 2)
			c2 = elOut.value[2]
		} else {
			c2 = eval.poolQMul[2]
		}

		// Avoid overwritting if the second input is the output
		var tmp0, tmp1 *Element
		if el1 == elOut {
			tmp0, tmp1 = el1, el0
		} else {
			tmp0, tmp1 = el0, el1
		}

		ringQ.MFormLvl(level, tmp0.value[0], c00)
		ringQ.MFormLvl(level, tmp0.value[1], c01)

		if el0 == el1 { // squaring case
			ringQ.MulCoeffsMontgomeryLvl(level, c00, tmp1.value[0], c0) // c0 = c[0]*c[0]
			ringQ.MulCoeffsMontgomeryLvl(level, c01, tmp1.value[1], c2) // c2 = c[1]*c[1]
			ringQ.MulCoeffsMontgomeryLvl(level, c00, tmp1.value[1], c1) // c1 = 2*c[0]*c[1]
			ringQ.AddLvl(level, c1, c1, c1)

		} else { // regular case
			ringQ.MulCoeffsMontgomeryLvl(level, c00, tmp1.value[0], c0) // c0 = c0[0]*c0[0]
			ringQ.MulCoeffsMontgomeryLvl(level, c01, tmp1.value[1], c2) // c2 = c0[1]*c1[1]
			ringQ.MulCoeffsMontgomeryLvl(level, c00, tmp1.value[1], c1)
			ringQ.MulCoeffsMontgomeryAndAddLvl(level, c01, tmp1.value[0], c1) // c1 = c0[0]*c1[1] + c0[1]*c1[0]
		}

		// Relinearize if a key was provided
		if evakey != nil {
			eval.switchKeysInPlace(level, c2, evakey.evakey, eval.poolQ[1], eval.poolQ[2])
			ringQ.AddLvl(level, c0, eval.poolQ[1], elOut.value[0])
			ringQ.AddLvl(level, c1, eval.poolQ[2], elOut.value[1])
		}

		// Case Plaintext (x) Ciphertext or Ciphertext (x) Plaintext
	} else {

		var tmp0, tmp1 *Element

		if el0.Degree() == 1 {
			tmp0, tmp1 = el1, el0
		} else {
			tmp0, tmp1 = el0, el1
		}

		c00 := eval.poolQMul[0]

		ringQ.MFormLvl(level, tmp0.value[0], c00)
		ringQ.MulCoeffsMontgomeryLvl(level, c00, tmp1.value[0], elOut.value[0])
		ringQ.MulCoeffsMontgomeryLvl(level, c00, tmp1.value[1], elOut.value[1])
	}
}

// RelinearizeNew applies the relinearization procedure on ct0 and returns the result in a newly
// created Ciphertext. The input Ciphertext must be of degree two.
func (eval *evaluator) RelinearizeNew(ct0 *Ciphertext, evakey *EvaluationKey) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1, ct0.Level(), ct0.Scale())
	eval.Relinearize(ct0, evakey, ctOut)
	return
}

// Relinearize applies the relinearization procedure on ct0 and returns the result in ctOut. The input Ciphertext must be of degree two.
func (eval *evaluator) Relinearize(ct0 *Ciphertext, evakey *EvaluationKey, ctOut *Ciphertext) {
	if ct0.Degree() != 2 {
		panic("cannot Relinearize: input Ciphertext is not of degree 2")
	}

	if ctOut != ct0 {
		ctOut.SetScale(ct0.Scale())
	}

	level := utils.MinUint64(ct0.Level(), ctOut.Level())
	ringQ := eval.ringQ

	eval.switchKeysInPlace(level, ct0.value[2], evakey.evakey, eval.poolQ[1], eval.poolQ[2])

	ringQ.AddLvl(level, ct0.value[0], eval.poolQ[1], ctOut.value[0])
	ringQ.AddLvl(level, ct0.value[1], eval.poolQ[2], ctOut.value[1])

	ctOut.Resize(eval.params, 1)
}

// SwitchKeysNew re-encrypts ct0 under a different key and returns the result in a newly created element.
// It requires a SwitchingKey, which is computed from the key under which the Ciphertext is currently encrypted,
// and the key under which the Ciphertext will be re-encrypted.
func (eval *evaluator) SwitchKeysNew(ct0 *Ciphertext, switchingKey *SwitchingKey) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.Scale())
	eval.SwitchKeys(ct0, switchingKey, ctOut)
	return
}

// SwitchKeys re-encrypts ct0 under a different key and returns the result in ctOut.
// It requires a SwitchingKey, which is computed from the key under which the Ciphertext is currently encrypted,
// and the key under which the Ciphertext will be re-encrypted.
func (eval *evaluator) SwitchKeys(ct0 *Ciphertext, switchingKey *SwitchingKey, ctOut *Ciphertext) {

	if ct0.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot SwitchKeys: input and output Ciphertext must be of degree 1")
	}

	level := utils.MinUint64(ct0.Level(), ctOut.Level())
	ringQ := eval.ringQ

	eval.switchKeysInPlace(level, ct0.value[1], switchingKey, eval.poolQ[1], eval.poolQ[2])

	ringQ.AddLvl(level, ct0.value[0], eval.poolQ[1], ctOut.value[0])
	ringQ.CopyLvl(level, eval.poolQ[2], ctOut.value[1])
}

// RotateNew rotates the columns of ct0 by k positions to the left, and returns the result in a newly created element.
// If the provided element is a Ciphertext, a key-switching operation is necessary and a rotation key for the specific rotation needs to be provided.
func (eval *evaluator) RotateNew(ct0 *Ciphertext, k uint64, evakey *RotationKeys) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.Scale())
	eval.Rotate(ct0, k, evakey, ctOut)
	return
}

// Rotate rotates the columns of ct0 by k positions to the left and returns the result in ctOut.
// If the provided element is a Ciphertext, a key-switching operation is necessary and a rotation key for the specific rotation needs to be provided.
func (eval *evaluator) Rotate(ct0 *Ciphertext, k uint64, evakey *RotationKeys, ctOut *Ciphertext) {

	if ct0.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot Rotate: input and output Ciphertext must be of degree 1")
	}

	k &= ((eval.ringQ.N >> 1) - 1)

	if k == 0 {
		ctOut.Copy(ct0.El())

	} else {

		ctOut.SetScale(ct0.Scale())

		// It checks in the RotationKeys if the corresponding rotation has been generated
		if evakey.evakeyRotColLeft[k] != nil {

			eval.permuteNTT(ct0, evakey.permuteNTTLeftIndex[k], evakey.evakeyRotColLeft[k], ctOut)

		} else {

			// If not, it checks if the left and right pow2 rotations have been generated
			hasPow2Rotations := true
			for i := uint64(1); i < eval.ringQ.N>>1; i <<= 1 {
				if evakey.evakeyRotColLeft[i] == nil || evakey.evakeyRotColRight[i] == nil {
					hasPow2Rotations = false
					break
				}
			}

			// If yes, it computes the least amount of rotation between left and right required to apply the demanded rotation
			if hasPow2Rotations {

				if utils.HammingWeight64(k) <= utils.HammingWeight64((eval.ringQ.N>>1)-k) {
					eval.rotateLPow2(ct0, k, evakey, ctOut)
				} else {
					eval.rotateRPow2(ct0, (eval.ringQ.N>>1)-k, evakey, ctOut)
				}

				// Otherwise, it returns an error indicating that the keys have not been generated
			} else {
				panic("cannot Rotate: specific rotation and pow2 rotations have not been generated")
			}
		}
	}
}

func (eval *evaluator) rotateLPow2(ct0 *Ciphertext, k uint64, evakey *RotationKeys, ctOut *Ciphertext) {
	eval.rotatePow2(ct0, k, evakey.permuteNTTLeftIndex, evakey.evakeyRotColLeft, ctOut)
}

func (eval *evaluator) rotateRPow2(ct0 *Ciphertext, k uint64, evakey *RotationKeys, ctOut *Ciphertext) {
	eval.rotatePow2(ct0, k, evakey.permuteNTTRightIndex, evakey.evakeyRotColRight, ctOut)
}

func (eval *evaluator) rotatePow2(ct0 *Ciphertext, k uint64, permuteNTTIndex map[uint64][]uint64, evakeyRotCol map[uint64]*SwitchingKey, ctOut *Ciphertext) {

	var evakeyIndex uint64

	evakeyIndex = 1

	level := utils.MinUint64(ct0.Level(), ctOut.Level())

	eval.ringQ.CopyLvl(level, ct0.value[0], ctOut.value[0])
	eval.ringQ.CopyLvl(level, ct0.value[1], ctOut.value[1])

	for k > 0 {

		if k&1 == 1 {

			eval.permuteNTT(ctOut, permuteNTTIndex[evakeyIndex], evakeyRotCol[evakeyIndex], ctOut)
		}

		evakeyIndex <<= 1
		k >>= 1
	}
}

// ConjugateNew conjugates ct0 (which is equivalent to a row rotation) and returns the result in a newly
// created element. If the provided element is a Ciphertext, a key-switching operation is necessary and a rotation key
// for the row rotation needs to be provided.
func (eval *evaluator) ConjugateNew(ct0 *Ciphertext, evakey *RotationKeys) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, ct0.Degree(), ct0.Level(), ct0.Scale())
	eval.Conjugate(ct0, evakey, ctOut)
	return
}

// Conjugate conjugates ct0 (which is equivalent to a row rotation) and returns the result in ctOut.
// If the provided element is a Ciphertext, a key-switching operation is necessary and a rotation key for the row rotation needs to be provided.
func (eval *evaluator) Conjugate(ct0 *Ciphertext, evakey *RotationKeys, ctOut *Ciphertext) {

	if ct0.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot Conjugate: input and output Ciphertext must be of degree 1")
	}

	if evakey.evakeyConjugate == nil {
		panic("cannot Conjugate: rows rotation key not generated")
	}

	ctOut.SetScale(ct0.Scale())

	eval.permuteNTT(ct0, evakey.permuteNTTConjugateIndex, evakey.evakeyConjugate, ctOut)
}

func (eval *evaluator) permuteNTT(ct0 *Ciphertext, index []uint64, rotKeys *SwitchingKey, ctOut *Ciphertext) {

	level := utils.MinUint64(ct0.Level(), ctOut.Level())

	pool2Q := eval.poolQ[1]
	pool3Q := eval.poolQ[2]

	eval.switchKeysInPlace(level, ct0.Value()[1], rotKeys, pool2Q, pool3Q)

	eval.ringQ.AddLvl(level, pool2Q, ct0.value[0], pool2Q)

	ring.PermuteNTTWithIndexLvl(level, pool2Q, index, ctOut.value[0])
	ring.PermuteNTTWithIndexLvl(level, pool3Q, index, ctOut.value[1])
}

func (eval *evaluator) switchKeysInPlaceNoModDown(level uint64, cx *ring.Poly, evakey *SwitchingKey, pool2Q, pool2P, pool3Q, pool3P *ring.Poly) {
	var reduce uint64

	ringQ := eval.ringQ
	ringP := eval.ringP

	// Pointers allocation
	c2QiQ := eval.poolQ[0]
	c2QiP := eval.poolP[0]

	c2 := eval.poolQ[3]

	evakey0Q := new(ring.Poly)
	evakey1Q := new(ring.Poly)
	evakey0P := new(ring.Poly)
	evakey1P := new(ring.Poly)

	// We switch the element on which the switching key operation will be conducted out of the NTT domain

	ringQ.InvNTTLvl(level, cx, c2)

	reduce = 0

	alpha := eval.params.Alpha()
	beta := uint64(math.Ceil(float64(level+1) / float64(alpha)))

	// Key switching with CRT decomposition for the Qi
	for i := uint64(0); i < beta; i++ {

		eval.decomposeAndSplitNTT(level, i, cx, c2, c2QiQ, c2QiP)

		evakey0Q.Coeffs = evakey.evakey[i][0].Coeffs[:level+1]
		evakey1Q.Coeffs = evakey.evakey[i][1].Coeffs[:level+1]
		evakey0P.Coeffs = evakey.evakey[i][0].Coeffs[len(ringQ.Modulus):]
		evakey1P.Coeffs = evakey.evakey[i][1].Coeffs[len(ringQ.Modulus):]

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

		//
		if reduce&3 == 3 {
			ringQ.ReduceConstantLvl(level, pool2Q, pool2Q)
			ringQ.ReduceConstantLvl(level, pool3Q, pool3Q)
			ringP.ReduceConstant(pool2P, pool2P)
			ringP.ReduceConstant(pool3P, pool3P)
		}

		reduce++
	}

	ringQ.ReduceLvl(level, pool2Q, pool2Q)
	ringQ.ReduceLvl(level, pool3Q, pool3Q)
	ringP.Reduce(pool2P, pool2P)
	ringP.Reduce(pool3P, pool3P)
}

// switchKeysInPlace applies the general key-switching procedure of the form [c0 + cx*evakey[0], c1 + cx*evakey[1]]
func (eval *evaluator) switchKeysInPlace(level uint64, cx *ring.Poly, evakey *SwitchingKey, p0, p1 *ring.Poly) {

	eval.switchKeysInPlaceNoModDown(level, cx, evakey, p0, eval.poolP[1], p1, eval.poolP[2])

	eval.baseconverter.ModDownSplitNTTPQ(level, p0, eval.poolP[1], p0)
	eval.baseconverter.ModDownSplitNTTPQ(level, p1, eval.poolP[2], p1)
}

// decomposeAndSplitNTT decomposes the input polynomial into the target CRT basis.
func (eval *evaluator) decomposeAndSplitNTT(level, beta uint64, c2NTT, c2InvNTT, c2QiQ, c2QiP *ring.Poly) {

	ringQ := eval.ringQ
	ringP := eval.ringP

	eval.decomposer.DecomposeAndSplit(level, beta, c2InvNTT, c2QiQ, c2QiP)

	p0idxst := beta * eval.params.Alpha()
	p0idxed := p0idxst + eval.decomposer.Xalpha()[beta]

	// c2_qi = cx mod qi mod qi
	for x := uint64(0); x < level+1; x++ {

		qi := ringQ.Modulus[x]
		nttPsi := ringQ.NttPsi[x]

		if p0idxst <= x && x < p0idxed {
			p0tmp := c2NTT.Coeffs[x]
			p1tmp := c2QiQ.Coeffs[x]
			for j := uint64(0); j < ringQ.N; j++ {
				p1tmp[j] = p0tmp[j]
			}
		} else {
			ring.NTTLazy(c2QiQ.Coeffs[x], c2QiQ.Coeffs[x], ringQ.N, nttPsi, qi)
		}
	}
	// c2QiP = c2 mod qi mod pj
	ringP.NTTLazy(c2QiP, c2QiP)
}

// RotateHoisted takes an input Ciphertext and a list of rotations and returns a map of Ciphertext, where each element of the map is the input Ciphertext
// rotation by one element of the list. It is much faster than sequential calls to Rotate.
func (eval *evaluator) RotateHoisted(ct0 *Ciphertext, rotations []uint64, rotkeys *RotationKeys) (cOut map[uint64]*Ciphertext) {

	// Pre-computation for rotations using hoisting
	ringQ := eval.ringQ
	ringP := eval.ringP

	c2NTT := ct0.value[1]
	c2InvNTT := ringQ.NewPoly()
	ringQ.InvNTTLvl(ct0.Level(), c2NTT, c2InvNTT)

	alpha := eval.params.Alpha()
	beta := uint64(math.Ceil(float64(ct0.Level()+1) / float64(alpha)))

	c2QiQDecomp := make([]*ring.Poly, beta)
	c2QiPDecomp := make([]*ring.Poly, beta)

	for i := uint64(0); i < beta; i++ {
		c2QiQDecomp[i] = ringQ.NewPoly()
		c2QiPDecomp[i] = ringP.NewPoly()
		eval.decomposeAndSplitNTT(ct0.Level(), i, c2NTT, c2InvNTT, c2QiQDecomp[i], c2QiPDecomp[i])
	}

	cOut = make(map[uint64]*Ciphertext)

	for _, i := range rotations {

		i &= ((ringQ.N >> 1) - 1)

		if i == 0 {
			cOut[i] = ct0.CopyNew().Ciphertext()
		} else {
			cOut[i] = NewCiphertext(eval.params, 1, ct0.Level(), ct0.Scale())
			eval.permuteNTTHoisted(ct0, c2QiQDecomp, c2QiPDecomp, i, rotkeys, cOut[i])
		}
	}

	return
}

func (eval *evaluator) permuteNTTHoisted(ct0 *Ciphertext, c2QiQDecomp, c2QiPDecomp []*ring.Poly, k uint64, rotKeys *RotationKeys, ctOut *Ciphertext) {

	if ct0.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot switchKeyHoisted: input and output Ciphertext must be of degree 1")
	}

	k &= 2*eval.ringQ.N - 1

	if rotKeys.permuteNTTLeftIndex[k] == nil {
		panic("cannot switchKeyHoisted: specific rotation has not been generated")
	}

	ctOut.SetScale(ct0.Scale())

	pool2Q := eval.poolQ[0]
	pool3Q := eval.poolQ[1]

	pool2P := eval.poolP[0]
	pool3P := eval.poolP[1]

	level := ctOut.Level()

	eval.keyswitchHoisted(level, c2QiQDecomp, c2QiPDecomp, rotKeys.evakeyRotColLeft[k], pool2Q, pool3Q, pool2P, pool3P)

	eval.ringQ.AddLvl(level, pool2Q, ct0.value[0], pool2Q)

	ring.PermuteNTTWithIndexLvl(level, pool2Q, rotKeys.permuteNTTLeftIndex[k], ctOut.value[0])
	ring.PermuteNTTWithIndexLvl(level, pool3Q, rotKeys.permuteNTTLeftIndex[k], ctOut.value[1])
}

func (eval *evaluator) keyswitchHoisted(level uint64, c2QiQDecomp, c2QiPDecomp []*ring.Poly, evakey *SwitchingKey, pool2Q, pool3Q, pool2P, pool3P *ring.Poly) {

	eval.keyswitchHoistedNoModDown(level, c2QiQDecomp, c2QiPDecomp, evakey, pool2Q, pool3Q, pool2P, pool3P)

	// Computes pool2Q = pool2Q/pool2P and pool3Q = pool3Q/pool3P
	eval.baseconverter.ModDownSplitNTTPQ(level, pool2Q, pool2P, pool2Q)
	eval.baseconverter.ModDownSplitNTTPQ(level, pool3Q, pool3P, pool3Q)
}

func (eval *evaluator) keyswitchHoistedNoModDown(level uint64, c2QiQDecomp, c2QiPDecomp []*ring.Poly, evakey *SwitchingKey, pool2Q, pool3Q, pool2P, pool3P *ring.Poly) {

	ringQ := eval.ringQ
	ringP := eval.ringP

	alpha := eval.params.Alpha()
	beta := uint64(math.Ceil(float64(level+1) / float64(alpha)))

	evakey0Q := new(ring.Poly)
	evakey1Q := new(ring.Poly)
	evakey0P := new(ring.Poly)
	evakey1P := new(ring.Poly)

	// Key switching with CRT decomposition for the Qi
	var reduce uint64
	for i := uint64(0); i < beta; i++ {

		evakey0Q.Coeffs = evakey.evakey[i][0].Coeffs[:level+1]
		evakey1Q.Coeffs = evakey.evakey[i][1].Coeffs[:level+1]
		evakey0P.Coeffs = evakey.evakey[i][0].Coeffs[len(ringQ.Modulus):]
		evakey1P.Coeffs = evakey.evakey[i][1].Coeffs[len(ringQ.Modulus):]

		if i == 0 {
			ringQ.MulCoeffsMontgomeryLvl(level, evakey0Q, c2QiQDecomp[i], pool2Q)
			ringQ.MulCoeffsMontgomeryLvl(level, evakey1Q, c2QiQDecomp[i], pool3Q)
			ringP.MulCoeffsMontgomery(evakey0P, c2QiPDecomp[i], pool2P)
			ringP.MulCoeffsMontgomery(evakey1P, c2QiPDecomp[i], pool3P)
		} else {
			ringQ.MulCoeffsMontgomeryAndAddNoModLvl(level, evakey0Q, c2QiQDecomp[i], pool2Q)
			ringQ.MulCoeffsMontgomeryAndAddNoModLvl(level, evakey1Q, c2QiQDecomp[i], pool3Q)
			ringP.MulCoeffsMontgomeryAndAddNoMod(evakey0P, c2QiPDecomp[i], pool2P)
			ringP.MulCoeffsMontgomeryAndAddNoMod(evakey1P, c2QiPDecomp[i], pool3P)
		}

		if reduce&7 == 1 {
			ringQ.ReduceLvl(level, pool2Q, pool2Q)
			ringQ.ReduceLvl(level, pool3Q, pool3Q)
			ringP.Reduce(pool2P, pool2P)
			ringP.Reduce(pool3P, pool3P)
		}

		reduce++
	}

	if (reduce-1)&7 != 1 {
		ringQ.ReduceLvl(level, pool2Q, pool2Q)
		ringQ.ReduceLvl(level, pool3Q, pool3Q)
		ringP.Reduce(pool2P, pool2P)
		ringP.Reduce(pool3P, pool3P)
	}
}
