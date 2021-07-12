package rlwe

import (
	"encoding/binary"

	"github.com/ldsec/lattigo/v2/ring"
)

// SecretKey is a type for generic RLWE secret keys.
type SecretKey struct {
	Value [2]*ring.Poly
}

// PublicKey is a type for generic RLWE public keys.
type PublicKey struct {
	Value [2][2]*ring.Poly
}

// SwitchingKey is a type for generic RLWE public switching keys.
type SwitchingKey struct {
	Value [][2][2]*ring.Poly
}

// RelinearizationKey is a type for generic RLWE public relinearization keys. It stores a slice with a
// switching key per relinearizable degree. The switching key at index i is used to relinearize a degree
// i+2 ciphertexts back to a degree i + 1 one.
type RelinearizationKey struct {
	Keys []*SwitchingKey
}

// RotationKeySet is a type for storing generic RLWE public rotation keys. It stores a map indexed by the
// galois element defining the automorphism.
type RotationKeySet struct {
	Keys map[uint64]*SwitchingKey
}

// EvaluationKey is a type for storing generic RLWE public evaluation keys. An evaluation key is a union
// of a relinearization key and a set of rotation keys.
type EvaluationKey struct {
	Rlk  *RelinearizationKey
	Rtks *RotationKeySet
}

// NewSecretKey generates a new SecretKey with zero values.
func NewSecretKey(params Parameters) *SecretKey {
	sk := new(SecretKey)
	sk.Value = [2]*ring.Poly{ring.NewPoly(params.N(), params.QCount()), ring.NewPoly(params.N(), params.PCount())}
	return sk
}

// NewPublicKey returns a new PublicKey with zero values.
func NewPublicKey(params Parameters) (pk *PublicKey) {
	return &PublicKey{
		Value: [2][2]*ring.Poly{
			{ring.NewPoly(params.N(), params.QCount()), ring.NewPoly(params.N(), params.PCount())},
			{ring.NewPoly(params.N(), params.QCount()), ring.NewPoly(params.N(), params.PCount())}},
	}
}

// Equals checks two PublicKey struct for equality.
func (pk *PublicKey) Equals(other *PublicKey) bool {
	if pk == other {
		return true
	}
	nilVal := [2][2]*ring.Poly{}
	return pk.Value != nilVal && other.Value != nilVal &&
		pk.Value[0][0].Equals(other.Value[0][0]) &&
		pk.Value[0][1].Equals(other.Value[0][1]) &&
		pk.Value[1][0].Equals(other.Value[1][0]) &&
		pk.Value[1][1].Equals(other.Value[1][1])
}

// NewRotationKeySet returns a new RotationKeySet with pre-allocated switching keys for each distinct galoisElement value.
func NewRotationKeySet(params Parameters, galoisElement []uint64) (rotKey *RotationKeySet) {
	rotKey = new(RotationKeySet)
	rotKey.Keys = make(map[uint64]*SwitchingKey, len(galoisElement))
	for _, galEl := range galoisElement {
		rotKey.Keys[galEl] = NewSwitchingKey(params)
	}
	return
}

// GetRotationKey return the rotation key for the given galois element or nil if such key is not in the set. The
// second argument is true  iff the first one is non-nil.
func (rtks *RotationKeySet) GetRotationKey(galoisEl uint64) (*SwitchingKey, bool) {
	rotKey, inSet := rtks.Keys[galoisEl]
	return rotKey, inSet
}

// NewSwitchingKey returns a new public switching key with pre-allocated zero-value
func NewSwitchingKey(params Parameters) *SwitchingKey {
	swk := new(SwitchingKey)
	swk.Value = make([][2][2]*ring.Poly, params.Beta())
	for i := 0; i < params.Beta(); i++ {
		swk.Value[i][0][0] = ring.NewPoly(params.N(), params.QCount())
		swk.Value[i][0][1] = ring.NewPoly(params.N(), params.PCount())
		swk.Value[i][1][0] = ring.NewPoly(params.N(), params.QCount())
		swk.Value[i][1][1] = ring.NewPoly(params.N(), params.PCount())
	}

	return swk
}

// NewRelinKey creates a new EvaluationKey with zero values.
func NewRelinKey(params Parameters, maxRelinDegree int) (evakey *RelinearizationKey) {

	evakey = new(RelinearizationKey)

	evakey.Keys = make([]*SwitchingKey, maxRelinDegree)

	for d := 0; d < maxRelinDegree; d++ {
		evakey.Keys[d] = NewSwitchingKey(params)
	}

	return
}

// GetDataLen returns the length in bytes of the target SecretKey.
func (sk *SecretKey) GetDataLen(WithMetadata bool) (dataLen int) {
	return sk.Value[0].GetDataLen(WithMetadata) + sk.Value[1].GetDataLen(WithMetadata)
}

// MarshalBinary encodes a secret key in a byte slice.
func (sk *SecretKey) MarshalBinary() (data []byte, err error) {

	data = make([]byte, sk.GetDataLen(true))

	var inc int
	if inc, err = sk.Value[0].WriteTo(data); err != nil {
		return nil, err
	}

	if _, err = sk.Value[1].WriteTo(data[inc:]); err != nil {
		return nil, err
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled SecretKey in the target SecretKey.
func (sk *SecretKey) UnmarshalBinary(data []byte) (err error) {

	sk.Value[0] = new(ring.Poly)
	sk.Value[1] = new(ring.Poly)

	var inc int
	if inc, err = sk.Value[0].DecodePolyNew(data); err != nil {
		return err
	}

	if _, err = sk.Value[1].DecodePolyNew(data[inc:]); err != nil {
		return err
	}

	return nil
}

// CopyNew creates a deep copy of the receiver secret key and returns it.
func (sk *SecretKey) CopyNew() *SecretKey {
	if sk == nil || sk.Value[0] == nil {
		return nil
	}
	return &SecretKey{[2]*ring.Poly{sk.Value[0].CopyNew(), sk.Value[1].CopyNew()}}
}

// GetDataLen returns the length in bytes of the target PublicKey.
func (pk *PublicKey) GetDataLen(WithMetadata bool) (dataLen int) {

	for _, el := range pk.Value {
		dataLen += el[0].GetDataLen(WithMetadata)
		dataLen += el[1].GetDataLen(WithMetadata)
	}

	return
}

// MarshalBinary encodes a PublicKey in a byte slice.
func (pk *PublicKey) MarshalBinary() (data []byte, err error) {

	dataLen := pk.GetDataLen(true)

	data = make([]byte, dataLen)

	var pointer, inc int

	if inc, err = pk.Value[0][0].WriteTo(data[pointer:]); err != nil {
		return nil, err
	}
	pointer += inc

	if inc, err = pk.Value[0][1].WriteTo(data[pointer:]); err != nil {
		return nil, err
	}
	pointer += inc

	if inc, err = pk.Value[1][0].WriteTo(data[pointer:]); err != nil {
		return nil, err
	}
	pointer += inc

	if _, err = pk.Value[1][1].WriteTo(data[pointer:]); err != nil {
		return nil, err
	}

	return data, err

}

// UnmarshalBinary decodes a previously marshaled PublicKey in the target PublicKey.
func (pk *PublicKey) UnmarshalBinary(data []byte) (err error) {

	var pointer, inc int

	pk.Value[0][0] = new(ring.Poly)
	pk.Value[0][1] = new(ring.Poly)
	pk.Value[1][0] = new(ring.Poly)
	pk.Value[1][1] = new(ring.Poly)

	if inc, err = pk.Value[0][0].DecodePolyNew(data[pointer:]); err != nil {
		return err
	}
	pointer += inc

	if inc, err = pk.Value[0][1].DecodePolyNew(data[pointer:]); err != nil {
		return err
	}
	pointer += inc

	if inc, err = pk.Value[1][0].DecodePolyNew(data[pointer:]); err != nil {
		return err
	}
	pointer += inc

	if _, err = pk.Value[1][1].DecodePolyNew(data[pointer:]); err != nil {
		return err
	}

	return nil
}

// CopyNew creates a deep copy of the receiver PublicKey and returns it.
func (pk *PublicKey) CopyNew() *PublicKey {
	if pk == nil || pk.Value[0][0] == nil || pk.Value[0][1] == nil || pk.Value[1][0] == nil || pk.Value[1][1] == nil {
		return nil
	}
	return &PublicKey{[2][2]*ring.Poly{{pk.Value[0][0].CopyNew(), pk.Value[0][1].CopyNew()}, {pk.Value[1][0].CopyNew(), pk.Value[1][1].CopyNew()}}}
}

// Equals checks two RelinearizationKeys for equality.
func (rlk *RelinearizationKey) Equals(other *RelinearizationKey) bool {
	if rlk == other {
		return true
	}
	if (rlk == nil) != (other == nil) {
		return false
	}
	if len(rlk.Keys) != len(other.Keys) {
		return false
	}
	for i := range rlk.Keys {
		if !rlk.Keys[i].Equals(other.Keys[i]) {
			return false
		}
	}
	return true
}

// GetDataLen returns the length in bytes of the target EvaluationKey.
func (rlk *RelinearizationKey) GetDataLen(WithMetadata bool) (dataLen int) {

	if WithMetadata {
		dataLen++
	}

	for _, evakey := range rlk.Keys {
		dataLen += (*SwitchingKey)(evakey).GetDataLen(WithMetadata)
	}

	return
}

// MarshalBinary encodes an EvaluationKey key in a byte slice.
func (rlk *RelinearizationKey) MarshalBinary() (data []byte, err error) {

	var pointer int

	dataLen := rlk.GetDataLen(true)

	data = make([]byte, dataLen)

	data[0] = uint8(len(rlk.Keys))

	pointer++

	for _, evakey := range rlk.Keys {

		if pointer, err = (*SwitchingKey)(evakey).encode(pointer, data); err != nil {
			return nil, err
		}
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled EvaluationKey in the target EvaluationKey.
func (rlk *RelinearizationKey) UnmarshalBinary(data []byte) (err error) {

	deg := int(data[0])

	rlk.Keys = make([]*SwitchingKey, deg)

	pointer := 1
	var inc int
	for i := 0; i < deg; i++ {
		rlk.Keys[i] = new(SwitchingKey)
		if inc, err = rlk.Keys[i].decode(data[pointer:]); err != nil {
			return err
		}
		pointer += inc
	}

	return nil
}

// CopyNew creates a deep copy of the receiver RelinearizationKey and returns it.
func (rlk *RelinearizationKey) CopyNew() *RelinearizationKey {
	if rlk == nil || len(rlk.Keys) == 0 {
		return nil
	}
	rlkb := &RelinearizationKey{Keys: make([]*SwitchingKey, len(rlk.Keys))}
	for i, swk := range rlk.Keys {
		rlkb.Keys[i] = swk.CopyNew()
	}
	return rlkb
}

// Equals checks two SwitchingKeys for equality.
func (swk *SwitchingKey) Equals(other *SwitchingKey) bool {
	if swk == other {
		return true
	}
	if (swk == nil) != (other == nil) {
		return false
	}
	if len(swk.Value) != len(other.Value) {
		return false
	}
	for i := range swk.Value {
		if !(swk.Value[i][0][0].Equals(other.Value[i][0][0]) &&
			swk.Value[i][0][1].Equals(other.Value[i][0][1]) &&
			swk.Value[i][1][1].Equals(other.Value[i][1][1]) &&
			swk.Value[i][1][1].Equals(other.Value[i][1][1])) {
			return false
		}
	}
	return true
}

// GetDataLen returns the length in bytes of the target SwitchingKey.
func (swk *SwitchingKey) GetDataLen(WithMetadata bool) (dataLen int) {

	if WithMetadata {
		dataLen++
	}

	for j := uint64(0); j < uint64(len(swk.Value)); j++ {
		dataLen += swk.Value[j][0][0].GetDataLen(WithMetadata)
		dataLen += swk.Value[j][0][1].GetDataLen(WithMetadata)
		dataLen += swk.Value[j][1][0].GetDataLen(WithMetadata)
		dataLen += swk.Value[j][1][1].GetDataLen(WithMetadata)
	}

	return
}

// MarshalBinary encodes an SwitchingKey in a byte slice.
func (swk *SwitchingKey) MarshalBinary() (data []byte, err error) {

	data = make([]byte, swk.GetDataLen(true))

	if _, err = swk.encode(0, data); err != nil {
		return nil, err
	}

	return data, nil
}

// UnmarshalBinary decode a previously marshaled SwitchingKey in the target SwitchingKey.
func (swk *SwitchingKey) UnmarshalBinary(data []byte) (err error) {

	if _, err = swk.decode(data); err != nil {
		return err
	}

	return nil
}

// CopyNew creates a deep copy of the receiver SwitchingKey and returns it.
func (swk *SwitchingKey) CopyNew() *SwitchingKey {
	if swk == nil || len(swk.Value) == 0 {
		return nil
	}
	swkb := &SwitchingKey{Value: make([][2][2]*ring.Poly, len(swk.Value))}
	for i, el := range swk.Value {
		swkb.Value[i] = [2][2]*ring.Poly{{el[0][0].CopyNew(), el[0][1].CopyNew()}, {el[1][0].CopyNew(), el[1][1].CopyNew()}}
	}
	return swkb
}

func (swk *SwitchingKey) encode(pointer int, data []byte) (int, error) {

	var err error
	var inc int

	data[pointer] = uint8(len(swk.Value))

	pointer++

	for j := 0; j < len(swk.Value); j++ {

		if inc, err = swk.Value[j][0][0].WriteTo(data[pointer : pointer+swk.Value[j][0][0].GetDataLen(true)]); err != nil {
			return pointer, err
		}

		pointer += inc

		if inc, err = swk.Value[j][0][1].WriteTo(data[pointer : pointer+swk.Value[j][0][1].GetDataLen(true)]); err != nil {
			return pointer, err
		}

		pointer += inc

		if inc, err = swk.Value[j][1][0].WriteTo(data[pointer : pointer+swk.Value[j][1][0].GetDataLen(true)]); err != nil {
			return pointer, err
		}

		pointer += inc

		if inc, err = swk.Value[j][1][1].WriteTo(data[pointer : pointer+swk.Value[j][1][1].GetDataLen(true)]); err != nil {
			return pointer, err
		}

		pointer += inc
	}

	return pointer, nil
}

func (swk *SwitchingKey) decode(data []byte) (pointer int, err error) {

	decomposition := int(data[0])

	pointer = 1

	swk.Value = make([][2][2]*ring.Poly, decomposition)

	var inc int

	for j := 0; j < decomposition; j++ {

		swk.Value[j][0][0] = new(ring.Poly)
		if inc, err = swk.Value[j][0][0].DecodePolyNew(data[pointer:]); err != nil {
			return pointer, err
		}
		pointer += inc

		swk.Value[j][0][1] = new(ring.Poly)
		if inc, err = swk.Value[j][0][1].DecodePolyNew(data[pointer:]); err != nil {
			return pointer, err
		}
		pointer += inc

		swk.Value[j][1][0] = new(ring.Poly)
		if inc, err = swk.Value[j][1][0].DecodePolyNew(data[pointer:]); err != nil {
			return pointer, err
		}
		pointer += inc

		swk.Value[j][1][1] = new(ring.Poly)
		if inc, err = swk.Value[j][1][1].DecodePolyNew(data[pointer:]); err != nil {
			return pointer, err
		}
		pointer += inc

	}

	return pointer, nil
}

// Equals checks to RotationKeySets for equality.
func (rtks *RotationKeySet) Equals(other *RotationKeySet) bool {
	if rtks == other {
		return true
	}
	if (rtks == nil) || (other == nil) {
		return false
	}
	if len(rtks.Keys) != len(other.Keys) {
		return false
	}
	return rtks.Includes(other)
}

// Includes checks whether the receiver RotationKeySet includes the given other RotationKeySet.
func (rtks *RotationKeySet) Includes(other *RotationKeySet) bool {
	if (rtks == nil) || (other == nil) {
		return false
	}
	for galEl := range other.Keys {
		if _, inSet := rtks.Keys[galEl]; !inSet {
			return false
		}
	}
	return true
}

// GetDataLen returns the length in bytes of the target RotationKeys.
func (rtks *RotationKeySet) GetDataLen(WithMetaData bool) (dataLen int) {
	for _, k := range rtks.Keys {
		if WithMetaData {
			dataLen += 4
		}
		dataLen += k.GetDataLen(WithMetaData)
	}
	return
}

// MarshalBinary encodes a RotationKeys struct in a byte slice.
func (rtks *RotationKeySet) MarshalBinary() (data []byte, err error) {

	data = make([]byte, rtks.GetDataLen(true))

	pointer := int(0)

	for galEL, key := range rtks.Keys {

		binary.BigEndian.PutUint32(data[pointer:pointer+4], uint32(galEL))
		pointer += 4

		if pointer, err = key.encode(pointer, data); err != nil {
			return nil, err
		}
	}

	return data, nil
}

// UnmarshalBinary decodes a previously marshaled RotationKeys in the target RotationKeys.
func (rtks *RotationKeySet) UnmarshalBinary(data []byte) (err error) {

	rtks.Keys = make(map[uint64]*SwitchingKey)

	for len(data) > 0 {

		galEl := uint64(binary.BigEndian.Uint32(data))
		data = data[4:]

		swk := new(SwitchingKey)
		var inc int
		if inc, err = swk.decode(data); err != nil {
			return err
		}
		data = data[inc:]
		rtks.Keys[galEl] = swk

	}

	return nil
}
