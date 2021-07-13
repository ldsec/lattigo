package ckks

import (
	"encoding/binary"
	"errors"
	"math"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
)

// GetDataLen returns the length in bytes of the target Ciphertext.
func (ciphertext *Ciphertext) GetDataLen(WithMetaData bool) (dataLen int) {
	// MetaData is :
	// 1 byte : Degree
	// 9 byte : Scale
	// 1 byte : isNTT
	if WithMetaData {
		dataLen += 10
	}

	for _, el := range ciphertext.Value {
		dataLen += el.GetDataLen(WithMetaData)
	}

	return dataLen
}

// MarshalBinary encodes a Ciphertext on a byte slice. The total size
// in byte is 4 + 8* N * numberModuliQ * (degree + 1).
func (ciphertext *Ciphertext) MarshalBinary() (data []byte, err error) {

	data = make([]byte, ciphertext.GetDataLen(true))

	data[0] = uint8(ciphertext.Degree() + 1)

	binary.LittleEndian.PutUint64(data[1:9], math.Float64bits(ciphertext.Scale))

	var pointer, inc int

	pointer = 10

	for _, el := range ciphertext.Value {

		if inc, err = el.WriteTo(data[pointer:]); err != nil {
			return nil, err
		}

		pointer += inc
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled Ciphertext on the target Ciphertext.
func (ciphertext *Ciphertext) UnmarshalBinary(data []byte) (err error) {
	if len(data) < 10 { // cf. ciphertext.GetDataLen()
		return errors.New("too small bytearray")
	}

	ciphertext.Ciphertext = new(rlwe.Ciphertext)

	ciphertext.Value = make([]*ring.Poly, uint8(data[0]))

	ciphertext.Scale = math.Float64frombits(binary.LittleEndian.Uint64(data[1:9]))

	var pointer, inc int
	pointer = 10

	for i := range ciphertext.Value {

		ciphertext.Value[i] = new(ring.Poly)

		if inc, err = ciphertext.Value[i].DecodePolyNew(data[pointer:]); err != nil {
			return err
		}

		pointer += inc
	}

	if pointer != len(data) {
		return errors.New("remaining unparsed data")
	}

	return nil
}
