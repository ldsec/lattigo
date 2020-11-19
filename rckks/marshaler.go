package rckks

import (
	"encoding/binary"
	"errors"
	"math"

	"github.com/ldsec/lattigo/v2/ring"
)

// GetDataLen returns the length in bytes of the target Ciphertext.
func (ciphertext *Ciphertext) GetDataLen(WithMetaData bool) (dataLen uint64) {
	if WithMetaData {
		dataLen += 11
	}

	for _, el := range ciphertext.value {
		dataLen += el.GetDataLen(WithMetaData)
	}

	return dataLen
}

// MarshalBinary encodes a Ciphertext on a byte slice. The total size
// in byte is 4 + 8* N * numberModuliQ * (degree + 1).
func (ciphertext *Ciphertext) MarshalBinary() (data []byte, err error) {

	data = make([]byte, ciphertext.GetDataLen(true))

	data[0] = uint8(ciphertext.Degree() + 1)

	binary.LittleEndian.PutUint64(data[1:9], math.Float64bits(ciphertext.Scale()))

	if ciphertext.isNTT {
		data[10] = 1
	}

	var pointer, inc uint64

	pointer = 11

	for _, el := range ciphertext.value {

		if inc, err = el.WriteTo(data[pointer:]); err != nil {
			return nil, err
		}

		pointer += inc
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled Ciphertext on the target Ciphertext.
// The target Ciphertext must be of the appropriate format and size, it can be created with the
// method NewCiphertext(uint64).
func (ciphertext *Ciphertext) UnmarshalBinary(data []byte) (err error) {
	if len(data) < 11 {
		return errors.New("too small bytearray")
	}

	ciphertext.Element = new(Element)

	ciphertext.value = make([]*ring.Poly, uint8(data[0]))

	ciphertext.scale = math.Float64frombits(binary.LittleEndian.Uint64(data[1:9]))

	if uint8(data[10]) == 1 {
		ciphertext.isNTT = true
	}

	var pointer, inc uint64
	pointer = 11

	for i := range ciphertext.value {

		ciphertext.value[i] = new(ring.Poly)

		if inc, err = ciphertext.value[i].DecodePolyNew(data[pointer:]); err != nil {
			return err
		}

		pointer += inc
	}

	if pointer != uint64(len(data)) {
		return errors.New("remaining unparsed data")
	}

	return nil
}

// GetDataLen returns the length in bytes of the target SecretKey.
func (sk *SecretKey) GetDataLen(WithMetaData bool) (dataLen uint64) {
	return sk.sk.GetDataLen(WithMetaData)
}

// MarshalBinary encodes a SecretKey in a byte slice.
func (sk *SecretKey) MarshalBinary() (data []byte, err error) {

	data = make([]byte, sk.GetDataLen(true))

	if _, err = sk.sk.WriteTo(data); err != nil {
		return nil, err
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled SecretKey on the target secret-key.
func (sk *SecretKey) UnmarshalBinary(data []byte) (err error) {

	sk.sk = new(ring.Poly)

	if _, err = sk.sk.DecodePolyNew(data); err != nil {
		return err
	}

	return nil
}

// GetDataLen returns the length in bytes of the target PublicKey.
func (pk *PublicKey) GetDataLen(WithMetaData bool) (dataLen uint64) {

	for _, el := range pk.pk {
		dataLen += el.GetDataLen(WithMetaData)
	}

	return
}

// MarshalBinary encodes a PublicKey in a byte slice.
func (pk *PublicKey) MarshalBinary() (data []byte, err error) {

	dataLen := pk.GetDataLen(true)

	data = make([]byte, dataLen)

	var pointer, inc uint64

	if inc, err = pk.pk[0].WriteTo(data[pointer:]); err != nil {
		return nil, err
	}

	if _, err = pk.pk[1].WriteTo(data[pointer+inc:]); err != nil {
		return nil, err
	}

	return data, err

}

// UnmarshalBinary decodes a previously marshaled PublicKey in the target PublicKey.
func (pk *PublicKey) UnmarshalBinary(data []byte) (err error) {

	var pointer, inc uint64

	pk.pk[0] = new(ring.Poly)
	pk.pk[1] = new(ring.Poly)

	if inc, err = pk.pk[0].DecodePolyNew(data[pointer:]); err != nil {
		return err
	}

	if _, err = pk.pk[1].DecodePolyNew(data[pointer+inc:]); err != nil {
		return err
	}

	return nil
}

// GetDataLen returns the length in bytes of the target EvaluationKey.
func (evaluationkey *EvaluationKey) GetDataLen(WithMetaData bool) (dataLen uint64) {
	return evaluationkey.evakey.GetDataLen(WithMetaData)
}

// MarshalBinary encodes an evaluation key in a byte slice.
func (evaluationkey *EvaluationKey) MarshalBinary() (data []byte, err error) {

	var pointer uint64

	dataLen := evaluationkey.evakey.GetDataLen(true)

	data = make([]byte, dataLen)

	if _, err = evaluationkey.evakey.encode(pointer, data); err != nil {
		return nil, err
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled evaluation-key in the target evaluation-key.
func (evaluationkey *EvaluationKey) UnmarshalBinary(data []byte) (err error) {
	evaluationkey.evakey = new(SwitchingKey)
	if _, err = evaluationkey.evakey.decode(data); err != nil {
		return err
	}
	return nil
}

// GetDataLen returns the length in bytes of the target SwitchingKey.
func (switchkey *SwitchingKey) GetDataLen(WithMetaData bool) (dataLen uint64) {

	if WithMetaData {
		dataLen++
	}

	for j := uint64(0); j < uint64(len(switchkey.evakey)); j++ {
		dataLen += switchkey.evakey[j][0].GetDataLen(WithMetaData)
		dataLen += switchkey.evakey[j][1].GetDataLen(WithMetaData)
	}

	return
}

// MarshalBinary encodes an SwitchingKey in a byte slice.
func (switchkey *SwitchingKey) MarshalBinary() (data []byte, err error) {

	data = make([]byte, switchkey.GetDataLen(true))

	if _, err = switchkey.encode(0, data); err != nil {
		return nil, err
	}

	return data, nil
}

// UnmarshalBinary decode a previously marshaled SwitchingKey in the target SwitchingKey.
func (switchkey *SwitchingKey) UnmarshalBinary(data []byte) (err error) {

	if _, err = switchkey.decode(data); err != nil {
		return err
	}

	return nil
}

func (switchkey *SwitchingKey) encode(pointer uint64, data []byte) (uint64, error) {

	var err error

	var inc uint64

	data[pointer] = uint8(len(switchkey.evakey))

	pointer++

	for j := uint64(0); j < uint64(len(switchkey.evakey)); j++ {

		if inc, err = switchkey.evakey[j][0].WriteTo(data[pointer:]); err != nil {
			return pointer, err
		}

		pointer += inc

		if inc, err = switchkey.evakey[j][1].WriteTo(data[pointer:]); err != nil {
			return pointer, err
		}

		pointer += inc
	}

	return pointer, nil
}

func (switchkey *SwitchingKey) decode(data []byte) (pointer uint64, err error) {

	decomposition := uint64(data[0])

	pointer = uint64(1)

	switchkey.evakey = make([][2]*ring.Poly, decomposition)

	var inc uint64

	for j := uint64(0); j < decomposition; j++ {

		switchkey.evakey[j][0] = new(ring.Poly)
		if inc, err = switchkey.evakey[j][0].DecodePolyNew(data[pointer:]); err != nil {
			return pointer, err
		}
		pointer += inc

		switchkey.evakey[j][1] = new(ring.Poly)
		if inc, err = switchkey.evakey[j][1].DecodePolyNew(data[pointer:]); err != nil {
			return pointer, err
		}
		pointer += inc

	}

	return pointer, nil
}

// GetDataLen returns the length in bytes of the target RotationKeys.
func (rotationkey *RotationKeys) GetDataLen(WithMetaData bool) (dataLen uint64) {
	for i := range rotationkey.evakeyRotColLeft {
		if WithMetaData {
			dataLen += 4
		}
		dataLen += rotationkey.evakeyRotColLeft[i].GetDataLen(WithMetaData)
	}

	for i := range rotationkey.evakeyRotColRight {
		if WithMetaData {
			dataLen += 4
		}
		dataLen += rotationkey.evakeyRotColRight[i].GetDataLen(WithMetaData)
	}

	return
}

// MarshalBinary encodes a RotationKeys structure in a byte slice.
func (rotationkey *RotationKeys) MarshalBinary() (data []byte, err error) {

	data = make([]byte, rotationkey.GetDataLen(true))

	mappingColL := []uint64{}
	mappingColR := []uint64{}

	for i := range rotationkey.evakeyRotColLeft {
		mappingColL = append(mappingColL, i)
	}

	for i := range rotationkey.evakeyRotColRight {
		mappingColR = append(mappingColR, i)
	}

	var pointer uint64

	for _, i := range mappingColL {

		binary.BigEndian.PutUint32(data[pointer:pointer+4], uint32(i))
		data[pointer] = uint8(RotationLeft)
		pointer += 4

		pointer, _ = rotationkey.evakeyRotColLeft[i].encode(pointer, data)
	}

	for _, i := range mappingColR {

		binary.BigEndian.PutUint32(data[pointer:pointer+4], uint32(i))
		data[pointer] = uint8(RotationRight)
		pointer += 4

		pointer, _ = rotationkey.evakeyRotColRight[i].encode(pointer, data)
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled RotationKeys in the target RotationKeys.
func (rotationkey *RotationKeys) UnmarshalBinary(data []byte) (err error) {

	var rotationType int
	var rotationNumber uint64

	pointer := uint64(0)
	var inc uint64

	dataLen := len(data)

	for dataLen > 0 {

		rotationType = int(data[pointer])
		rotationNumber = (uint64(data[pointer+1]) << 16) | (uint64(data[pointer+2]) << 8) | (uint64(data[pointer+3]))

		pointer += 4

		if rotationType == RotationLeft {

			if rotationkey.evakeyRotColLeft == nil {
				rotationkey.evakeyRotColLeft = make(map[uint64]*SwitchingKey)
			}

			if rotationkey.permuteNTTLeftIndex == nil {
				rotationkey.permuteNTTLeftIndex = make(map[uint64][]uint64)
			}

			rotationkey.evakeyRotColLeft[rotationNumber] = new(SwitchingKey)
			if inc, err = rotationkey.evakeyRotColLeft[rotationNumber].decode(data[pointer:]); err != nil {
				return err
			}

			N := uint64(len(rotationkey.evakeyRotColLeft[rotationNumber].evakey[0][0].Coeffs[0]))

			rotationkey.permuteNTTLeftIndex[rotationNumber] = ring.PermuteNTTIndex(GaloisGen, rotationNumber, N, 4*N)

		} else if rotationType == RotationRight {

			if rotationkey.evakeyRotColRight == nil {
				rotationkey.evakeyRotColRight = make(map[uint64]*SwitchingKey)
			}

			rotationkey.evakeyRotColRight[rotationNumber] = new(SwitchingKey)
			if inc, err = rotationkey.evakeyRotColRight[rotationNumber].decode(data[pointer:]); err != nil {
				return err
			}

			if rotationkey.permuteNTTRightIndex == nil {
				rotationkey.permuteNTTRightIndex = make(map[uint64][]uint64)
			}

			N := uint64(len(rotationkey.evakeyRotColRight[rotationNumber].evakey[0][0].Coeffs[0]))

			rotationkey.permuteNTTRightIndex[rotationNumber] = ring.PermuteNTTIndex(GaloisGen, (4*N)-rotationNumber, N, 4*N)

		} else {

			return err
		}

		pointer += inc

		dataLen -= int(4 + inc)
	}

	return nil
}
