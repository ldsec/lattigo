package rckks

import (
	"fmt"
	"math"
	"sort"
)

// PrecisionStats is a struct storing statistic about the precision of a CKKS plaintext
type PrecisionStats struct {
	MaxDelta        float64
	MinDelta        float64
	MaxPrecision    float64
	MinPrecision    float64
	MeanDelta       float64
	MeanPrecision   float64
	MedianDelta     float64
	MedianPrecision float64

	Dist []struct {
		Prec  float64
		Count int
	}

	cdfResol int
}

func (prec PrecisionStats) String() string {
	return fmt.Sprintf("\nMinimum precision : (%.2f) bits \n", prec.MinPrecision) +
		fmt.Sprintf("Maximum precision : (%.2f) bits \n", prec.MaxPrecision) +
		fmt.Sprintf("Mean    precision : (%.2f) bits \n", prec.MeanPrecision) +
		fmt.Sprintf("Median  precision : (%.2f) bits \n", prec.MedianPrecision)
}

// GetPrecisionStats generates a PrecisionStats struct from the reference values and the decrypted values
func GetPrecisionStats(params *Parameters, encoder Encoder, decryptor Decryptor, valuesWant []float64, element interface{}) (prec PrecisionStats) {

	var valuesTest []float64

	slots := params.Slots()

	switch element.(type) {
	case *Ciphertext:
		valuesTest = encoder.Decode(decryptor.DecryptNew(element.(*Ciphertext)), slots)
	case *Plaintext:
		valuesTest = encoder.Decode(element.(*Plaintext), slots)
	case []float64:
		valuesTest = element.([]float64)
	}

	var delta float64

	diff := make([]float64, params.Slots())

	prec.MaxDelta = 0
	prec.MinDelta = 1

	prec.MeanDelta = 0

	prec.cdfResol = 500

	prec.Dist = make([]struct {
		Prec  float64
		Count int
	}, prec.cdfResol)

	precision := make([]float64, len(valuesWant))

	for i := range valuesWant {

		delta = math.Abs(valuesTest[i] - valuesWant[i])
		precision[i] = math.Log2(1 / delta)

		diff[i] += delta

		prec.MeanDelta += diff[i]

		if delta > prec.MaxDelta {
			prec.MaxDelta = delta
		}

		if delta < prec.MinDelta {
			prec.MinDelta = delta
		}
	}

	prec.calcCDF(precision, prec.Dist)

	prec.MinPrecision = deltaToPrecision(prec.MaxDelta)
	prec.MaxPrecision = deltaToPrecision(prec.MinDelta)
	prec.MeanDelta /= float64(params.Slots())
	prec.MeanPrecision = deltaToPrecision(prec.MeanDelta)
	prec.MedianDelta = calcmedian(diff)
	prec.MedianPrecision = deltaToPrecision(prec.MedianDelta)

	return prec
}

func deltaToPrecision(c float64) float64 {
	return math.Log2(1 / c)
}

func (prec *PrecisionStats) calcCDF(precs []float64, res []struct {
	Prec  float64
	Count int
}) {
	sortedPrecs := make([]float64, len(precs))
	copy(sortedPrecs, precs)
	sort.Float64s(sortedPrecs)
	minPrec := sortedPrecs[0]                  //math.Log2(1/real(prec.MaxDelta))
	maxPrec := sortedPrecs[len(sortedPrecs)-1] //math.Log2(1/real(prec.MinDelta))
	for i := 0; i < prec.cdfResol; i++ {
		curPrec := minPrec + float64(i)*(maxPrec-minPrec)/float64(prec.cdfResol)
		for countSmaller, p := range sortedPrecs {
			if p >= curPrec {
				res[i].Prec = curPrec
				res[i].Count = countSmaller
				break
			}
		}
	}
}

func calcmedian(values []float64) (median float64) {

	sort.Float64s(values)

	index := len(values) / 2

	if len(values)&1 == 1 {
		return values[index]
	}

	if index+1 == len(values) {
		return values[index]
	}

	return (values[index] + values[index+1]) / 2
}
