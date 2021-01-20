package main

import (
	"fmt"
	"math"

	"github.com/ldsec/lattigo/v2/ckks"
	"github.com/ldsec/lattigo/v2/utils"
)

func chebyshevinterpolation() {

	// This example packs random 8192 float64 values in the range [-8, 8]
	// and approximates the function 1/(exp(-x) + 1) over the range [-8, 8].
	// The result is then parsed and compared to the expected result.

	// Scheme params
	LogN := uint64(14)
	LogSlots := uint64(13)

	LogModuli := ckks.LogModuli{
		LogQi: []uint64{55, 40, 40, 40, 40},
		LogPi: []uint64{45},
	}

	Scale := float64(1 << 40)

	params, err := ckks.NewParametersFromLogModuli(LogN, &LogModuli)
	if err != nil {
		panic(err)
	}
	params.SetScale(Scale)
	params.SetLogSlots(LogSlots)

	encoder := ckks.NewEncoder(params)

	// Keys
	kgen := ckks.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPair()

	// Relinearization key
	rlk := kgen.GenRelinKey(sk)

	// Encryptor
	encryptor := ckks.NewEncryptorFromPk(params, pk)

	// Decryptor
	decryptor := ckks.NewDecryptor(params, sk)

	// Evaluator
	evaluator := ckks.NewEvaluator(params)

	// Values to encrypt
	values0 := make([]complex128, params.Slots())
	values1 := make([]complex128, params.Slots())
	for i := range values0 {
		values0[i] = complex(math.Floor(utils.RandFloat64(0, 3)), 0)
		values1[i] = complex(math.Floor(utils.RandFloat64(0, 3)), 0)
	}

	fmt.Printf("CKKS parameters: logN = %d, logQ = %d, levels = %d, scale= %f, sigma = %f \n",
		params.LogN(), params.LogQP(), params.MaxLevel()+1, params.Scale(), params.Sigma())

	fmt.Println()
	fmt.Printf("Values0     : %6f %6f %6f %6f...\n",
		round(values0[0]), round(values0[1]), round(values0[2]), round(values0[3]))
	fmt.Printf("Values1     : %6f %6f %6f %6f...\n",
		round(values1[0]), round(values1[1]), round(values1[2]), round(values1[3]))
	fmt.Println()

	// Plaintext creation and encoding process
	plaintext0 := encoder.EncodeNew(values0, params.LogSlots())
	plaintext1 := encoder.EncodeNew(values1, params.LogSlots())

	// Encryption process
	var ciphertext0, ciphertext1 *ckks.Ciphertext
	ciphertext0 = encryptor.EncryptNew(plaintext0)
	ciphertext1 = encryptor.EncryptNew(plaintext1)

	fmt.Println("Evaluation of the function 1-(a-b)^2 in the range [-1, 1]")

	evaluator.Sub(ciphertext0, ciphertext1, ciphertext0)

	x2 := evaluator.MulRelinNew(ciphertext0, ciphertext0, rlk)
	evaluator.Rescale(x2, params.Scale(), x2)
	x4 := evaluator.MulRelinNew(x2, x2, rlk)
	evaluator.Rescale(x4, params.Scale(), x4)

	evaluator.MultByConst(x2, 5, x2)
	evaluator.Sub(x4, x2, x4)

	x4.SetScale(x4.Scale() * 4)

	evaluator.AddConst(x4, 1, x4)

	fmt.Println(params.Qi()[0], params.Qi()[1], params.Qi()[2])

	fmt.Println("Done... Consumed levels:", params.MaxLevel()-x4.Level())

	// Computation of the reference values
	for i := range values0 {
		x := (values0[i] - values1[i])
		x2 := x * x
		x4 := x2 * x2
		values0[i] = 0.25 * (x4 - 5*x2 + 4)
	}

	// Print results and comparison
	printDebug(params, x4, values0, decryptor, encoder)

}

func round(x complex128) complex128 {
	var factor float64 = 100000000
	a := math.Round(real(x)*factor) / factor
	b := math.Round(imag(x)*factor) / factor
	return complex(a, b)
}

func printDebug(params *ckks.Parameters, ciphertext *ckks.Ciphertext, valuesWant []complex128, decryptor ckks.Decryptor, encoder ckks.Encoder) (valuesTest []complex128) {

	valuesTest = encoder.DecodePublic(decryptor.DecryptNew(ciphertext), params.LogSlots(), 0)

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))
	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale()))
	fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant: %6.10f %6.10f %6.10f %6.10f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])
	fmt.Println()

	precStats := ckks.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, 0)

	fmt.Println(precStats.String())

	return
}

func main() {
	chebyshevinterpolation()
}
