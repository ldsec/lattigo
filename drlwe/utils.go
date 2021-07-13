package drlwe

func extendBasisSmallNormAndCenter(Q uint64, modulusP []uint64, coeffsQ []uint64, coeffsP [][]uint64) {
	QHalf := Q >> 1
	var sign uint64
	for j, c := range coeffsQ {

		sign = 1
		if c > QHalf {
			c = Q - c
			sign = 0
		}

		for i, pi := range modulusP {
			coeffsP[i][j] = (c * sign) | (pi-c)*(sign^1)
		}
	}
}
