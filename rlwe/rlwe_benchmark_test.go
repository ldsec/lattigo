package rlwe

import (
	"encoding/json"
	"runtime"
	"testing"
)

func BenchmarkRLWE(b *testing.B) {
	defaultParams := TestParams
	if testing.Short() {
		defaultParams = TestParams[:2]
	}
	if *flagParamString != "" {
		var jsonParams ParametersLiteral
		json.Unmarshal([]byte(*flagParamString), &jsonParams)
		defaultParams = []ParametersLiteral{jsonParams} // the custom test suite reads the parameters from the -params flag
	}

	for _, defaultParam := range defaultParams {
		params, err := NewParametersFromLiteral(defaultParam)
		if err != nil {
			panic(err)
		}

		kgen := NewKeyGenerator(params)
		keySwitcher := NewKeySwitcher(params)

		for _, testSet := range []func(kgen KeyGenerator, keySwitcher *KeySwitcher, b *testing.B){
			benchHoistedKeySwitch,
		} {
			testSet(kgen, keySwitcher, b)
			runtime.GC()
		}
	}
}

func benchHoistedKeySwitch(kgen KeyGenerator, keySwitcher *KeySwitcher, b *testing.B) {

	params := kgen.(*keyGenerator).params
	skIn := kgen.GenSecretKey()
	skOut := kgen.GenSecretKey()
	plaintext := NewPlaintext(params, params.MaxLevel())
	plaintext.Value.IsNTT = true
	encryptor := NewEncryptor(params, skIn)
	ciphertext := NewCiphertextNTT(params, 1, plaintext.Level())
	encryptor.Encrypt(plaintext, ciphertext)

	swk := kgen.GenSwitchingKey(skIn, skOut)

	b.Run(testString(params, "DecomposeNTT/"), func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			keySwitcher.DecomposeNTT(ciphertext.Level(), ciphertext.Value[1], keySwitcher.PoolDecompQ, keySwitcher.PoolDecompP)
		}
	})

	b.Run(testString(params, "KeySwitchHoisted/"), func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			keySwitcher.KeyswitchHoisted(ciphertext.Level(), keySwitcher.PoolDecompQ, keySwitcher.PoolDecompP, swk, ciphertext.Value[0], ciphertext.Value[1], keySwitcher.PoolP[1], keySwitcher.PoolP[2])
		}
	})
}
