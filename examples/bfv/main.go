package main

import (
	"fmt"
	"math"
	"math/bits"

	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"

	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
)

func obliviousRiding() {

	// This example simulates a situation where an anonymous rider
	// wants to find the closest available rider within a given area.
	// The application is inspired by the paper https://oride.epfl.ch/
	//
	// 		A. Pham, I. Dacosta, G. Endignoux, J. Troncoso-Pastoriza,
	//		K. Huguenin, and J.-P. Hubaux. ORide: A Privacy-Preserving
	//		yet Accountable Ride-Hailing Service. In Proceedings of the
	//		26th USENIX Security Symposium, Vancouver, BC, Canada, August 2017.
	//
	// Each area is represented as a rectangular grid where each driver
	// anyonymously signs in (i.e. the server only knows the driver is located
	// in the area).
	//
	// First, the rider generates an ephemeral key pair (riderSk, riderPk), which she
	// uses to encrypt her coordinates. She then sends the tuple (riderPk, enc(coordinates))
	// to the server handling the area she is in.
	//
	// Once the public key and the encrypted rider coordinates of the rider
	// have been received by the server, the rider's public key is transferred
	// to all the drivers within the area, with a randomized different index
	// for each of them, that indicates in which coefficient each driver must
	// encode her coordinates.
	//
	// Each driver encodes her coordinates in the designated coefficient and
	// uses the received public key to encrypt her encoded coordinates.
	// She then sends back the encrypted coordinates to the server.
	//
	// Once the encrypted coordinates of the drivers have been received, the server
	// homomorphically computes the squared distance: (x0 - x1)^2 + (y0 - y1)^2 between
	// the rider and each of the drivers, and sends back the encrypted result to the rider.
	//
	// The rider decrypts the result and chooses the closest driver.

	// Number of drivers in the area
	nbDrivers := 2048 //max is N

	// BFV parameters (128 bit security) with plaintext modulus 65929217
	paramDef := bfv.PN13QP218
	paramDef.T = 0x3ee0001

	params, err := bfv.NewParametersFromLiteral(paramDef)
	if err != nil {
		panic(err)
	}

	encoder := bfv.NewEncoder(params)

	// Rider's keygen
	kgen := bfv.NewKeyGenerator(params)

	riderSk, riderPk := kgen.GenKeyPair()

	decryptor := bfv.NewDecryptor(params, riderSk)

	encryptorRiderPk := bfv.NewEncryptor(params, riderPk)

	encryptorRiderSk := bfv.NewEncryptor(params, riderSk)

	evaluator := bfv.NewEvaluator(params, rlwe.EvaluationKey{})

	fmt.Println("============================================")
	fmt.Println("Homomorphic computations on batched integers")
	fmt.Println("============================================")
	fmt.Println()
	fmt.Printf("Parameters : N=%d, T=%d, Q = %d bits, sigma = %f \n",
		1<<params.LogN(), params.T(), params.LogQP(), params.Sigma())
	fmt.Println()

	maxvalue := uint64(math.Sqrt(float64(params.T()))) // max values = floor(sqrt(plaintext modulus))
	mask := uint64(1<<bits.Len64(maxvalue) - 1)        // binary mask upper-bound for the uniform sampling

	fmt.Printf("Generating %d driversData and 1 Rider randomly positioned on a grid of %d x %d units \n",
		nbDrivers, maxvalue, maxvalue)
	fmt.Println()

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	// Rider coordinates [x, y, x, y, ....., x, y]
	riderPosX, riderPosY := ring.RandUniform(prng, maxvalue, mask), ring.RandUniform(prng, maxvalue, mask)

	Rider := make([]uint64, 1<<params.LogN())
	for i := 0; i < nbDrivers; i++ {
		Rider[(i << 1)] = riderPosX
		Rider[(i<<1)+1] = riderPosY
	}

	riderPlaintext := bfv.NewPlaintext(params)
	encoder.EncodeUint(Rider, riderPlaintext)

	// driversData coordinates [0, 0, ..., x, y, ..., 0, 0]
	driversData := make([][]uint64, nbDrivers)

	driversPlaintexts := make([]*bfv.Plaintext, nbDrivers)
	for i := 0; i < nbDrivers; i++ {
		driversData[i] = make([]uint64, 1<<params.LogN())
		driversData[i][(i << 1)] = ring.RandUniform(prng, maxvalue, mask)
		driversData[i][(i<<1)+1] = ring.RandUniform(prng, maxvalue, mask)
		driversPlaintexts[i] = bfv.NewPlaintext(params)
		encoder.EncodeUint(driversData[i], driversPlaintexts[i])
	}

	fmt.Printf("Encrypting %d driversData (x, y) and 1 Rider (%d, %d) \n",
		nbDrivers, riderPosX, riderPosY)
	fmt.Println()

	RiderCiphertext := encryptorRiderSk.EncryptNew(riderPlaintext)

	DriversCiphertexts := make([]*bfv.Ciphertext, nbDrivers)
	for i := 0; i < nbDrivers; i++ {
		DriversCiphertexts[i] = encryptorRiderPk.EncryptNew(driversPlaintexts[i])
	}

	fmt.Println("Computing encrypted distance = ((CtD1 + CtD2 + CtD3 + CtD4...) - CtR)^2 ...")
	fmt.Println()

	evaluator.Neg(RiderCiphertext, RiderCiphertext)
	for i := 0; i < nbDrivers; i++ {
		evaluator.Add(RiderCiphertext, DriversCiphertexts[i], RiderCiphertext)
	}

	result := encoder.DecodeUintNew(decryptor.DecryptNew(evaluator.MulNew(RiderCiphertext, RiderCiphertext)))

	minIndex, minPosX, minPosY, minDist := 0, params.T(), params.T(), params.T()

	errors := 0

	for i := 0; i < nbDrivers; i++ {

		driverPosX, driverPosY := driversData[i][i<<1], driversData[i][(i<<1)+1]

		computedDist := result[i<<1] + result[(i<<1)+1]
		expectedDist := distance(driverPosX, driverPosY, riderPosX, riderPosY)

		if computedDist == expectedDist {
			if computedDist < minDist {
				minIndex = i
				minPosX, minPosY = driverPosX, driverPosY
				minDist = computedDist
			}
		} else {
			errors++
		}

		if i < 4 || i > nbDrivers-5 {
			fmt.Printf("Distance with Driver %d : %8d = (%4d - %4d)^2 + (%4d - %4d)^2 --> correct: %t\n",
				i, computedDist, driverPosX, riderPosX, driverPosY, riderPosY, computedDist == expectedDist)
		}

		if i == nbDrivers>>1 {
			fmt.Println("...")
		}
	}

	fmt.Printf("\nFinished with %.2f%% errors\n\n", 100*float64(errors)/float64(nbDrivers))
	fmt.Printf("Closest Driver to Rider is n°%d (%d, %d) with a distance of %d units\n",
		minIndex, minPosX, minPosY, int(math.Sqrt(float64(minDist))))
}

func distance(a, b, c, d uint64) uint64 {
	if a > c {
		a, c = c, a
	}
	if b > d {
		b, d = d, b
	}
	x, y := a-c, b-d
	return x*x + y*y
}

func main() {
	obliviousRiding()
}
