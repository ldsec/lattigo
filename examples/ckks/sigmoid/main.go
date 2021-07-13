package main

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"
)

func chebyshevinterpolation() {

	var err error

	// This example packs random 8192 float64 values in the range [-8, 8]
	// and approximates the function 1/(exp(-x) + 1) over the range [-8, 8].
	// The result is then parsed and compared to the expected result.

	// Scheme params are taken directly from the proposed defaults
	params, err := ckks.NewParametersFromLiteral(ckks.PN14QP438)
	if err != nil {
		panic(err)
	}

	encoder := ckks.NewEncoder(params)

	// Keys
	kgen := ckks.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPair()

	// Relinearization key
	rlk := kgen.GenRelinearizationKey(sk, 2)

	// Encryptor
	encryptor := ckks.NewEncryptor(params, pk)

	// Decryptor
	decryptor := ckks.NewDecryptor(params, sk)

	// Evaluator
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})

	// Values to encrypt
	values := make([]complex128, params.Slots())
	for i := range values {
		values[i] = complex(utils.RandFloat64(-8, 8), 0)
	}

	fmt.Printf("CKKS parameters: logN = %d, logQ = %d, levels = %d, scale= %f, sigma = %f \n",
		params.LogN(), params.LogQP(), params.MaxLevel()+1, params.Scale(), params.Sigma())

	fmt.Println()
	fmt.Printf("Values     : %6f %6f %6f %6f...\n",
		round(values[0]), round(values[1]), round(values[2]), round(values[3]))
	fmt.Println()

	// Plaintext creation and encoding process
	plaintext := encoder.EncodeNew(values, params.LogSlots())

	// Encryption process
	var ciphertext *ckks.Ciphertext
	ciphertext = encryptor.EncryptNew(plaintext)

	fmt.Println("Evaluation of the function 1/(exp(-x)+1) in the range [-8, 8] (degree of approximation: 32)")

	// Evaluation process
	// We approximate f(x) in the range [-8, 8] with a Chebyshev interpolant of 33 coefficients (degree 32).
	chebyapproximation := ckks.Approximate(f, -8, 8, 33)

	a := chebyapproximation.A()
	b := chebyapproximation.B()

	// Change of variable
	evaluator.MultByConst(ciphertext, 2/(b-a), ciphertext)
	evaluator.AddConst(ciphertext, (-a-b)/(b-a), ciphertext)
	if err := evaluator.Rescale(ciphertext, params.Scale(), ciphertext); err != nil {
		panic(err)
	}

	// We evaluate the interpolated Chebyshev interpolant on the ciphertext
	if ciphertext, err = evaluator.EvaluateCheby(ciphertext, chebyapproximation, ciphertext.Scale); err != nil {
		panic(err)
	}

	fmt.Println("Done... Consumed levels:", params.MaxLevel()-ciphertext.Level())

	// Computation of the reference values
	for i := range values {
		values[i] = f(values[i])
	}

	// Print results and comparison
	printDebug(params, ciphertext, values, decryptor, encoder)

}

func f(x complex128) complex128 {
	return 1 / (cmplx.Exp(-x) + 1)
}

func round(x complex128) complex128 {
	var factor float64 = 100000000
	a := math.Round(real(x)*factor) / factor
	b := math.Round(imag(x)*factor) / factor
	return complex(a, b)
}

func printDebug(params ckks.Parameters, ciphertext *ckks.Ciphertext, valuesWant []complex128, decryptor ckks.Decryptor, encoder ckks.Encoder) (valuesTest []complex128) {

	valuesTest = encoder.Decode(decryptor.DecryptNew(ciphertext), params.LogSlots())

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))
	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale))
	fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant: %6.10f %6.10f %6.10f %6.10f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])
	fmt.Println()

	precStats := ckks.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, params.LogSlots(), 0)

	fmt.Println(precStats.String())

	return
}

func main() {
	chebyshevinterpolation()
}
