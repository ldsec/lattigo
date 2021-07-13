
# Changelog
All notable changes to this project will be documented in this file. 

## Unreleased

- Added SECURITY.md
- ALL: when possible, public functions now use `int` instead of `uint64` as parameters and return values.
- ALL: `ring.Ring` are not instantiated once in the parameters and read only. They are then accessed by other structs, like the encryptor or evaluator.
- RING: removed `MulyPoly` and its related tests.
- RING: `ring.Ring` is now read only and thread safe.
- RING: RNS rescaling API is now inplace and can take a different poly as output.
- RING: added `ReadFromDistLvl` and `ReadAndAddFromDistLvl` to Gaussian sampler API.
- RING: added `IsNTT` and `IsMForm` flags in the `ring.Poly` type. For now, these flags are never checked or changed by the `ring` package.
- RLWE: added a new `rlwe` package as common implementation base for the lattigo RLWE schemes.
- RLWE: extracted the `rlwe.Parameters` type as common base for BFV and CKKS parameters.
- RLWE: extracted the `rlwe.KeyGenerator` type as common key-generator for BFV and CKKS.
- RLWE: extracted the `rlwe.Ciphertext` type as common base for BFV and CKKS ciphertexts.
- RLWE: extracted the `rlwe.Plaintext` type as common base for BFV and CKKS plaintext.
- RLWE: extracted the `rlwe.Encryptor`  type as common base for BFV and CKKS encryptors.
- RLWE: extracted the `rlwe.Decryptor`  type as common base for BFV and CKKS decryptors.
- RLWE: extracted the `rlwe.KeySwitcher` type as a common key-switching implementation for BFV and CKKS evaluators.
- RLWE: renamed the `Parameters.Copy()` method to `Parameters.CopyNew()` for consistency.
- RLWE: added `Parameter` struct now stores the relevant `ring.Ring` instances and has getter methods to access them.
- RLWE: added equality and inclusion check methods for the `rlwe.RotatationKeySet` type.
- RLWE: added tests for encryption, decryption, key-generation and key-switching.
- RLWE: moved keys related marshalling tests of `bfv` and `ckks` packages the `rlwe` package. 
- DRLWE: added a new `drlwe` package as a common implementation base for the lattigo multiparty RLWE schemes.
- DRLWE: added tests for the protocols.
- DRLWE: moved keys related marshalling tests of `dbfv` and `dckks` packages the `drlwe` package. 
- BFV/CKKS: the schemes are now using a common implementation for their keys.
- BFV/CKKS: the rotation-keys are now indexed by their corresponding Galois automorphism.
- BFV/CKKS: the `Evaluator` interface now has a single method for all column rotations and one method for the row-rotation/conjugate.
- BFV/CKKS: the relinearization and rotation keys are now passed to the `Evaluator` constructor methods (and no longer to the operations methods).
- BFV/CKKS: added the ParameterLiteral type for literally specifying scheme parameters in Go programs.
- BFV/CKKS: removed the now obsolete `Moduli` and `LogModuli` types and their associated `Parameters` constructors.
- BFV/CKKS: `Parameters` types are now passed by value in most situations.
- BFV/CKKS: added `encoding/json`-compatible JSON serialisers and deserialisers for the `Parameters` types.
- BFV/CKKS: removed the scheme-specific key types.
- BFV/CKKS: added a `-params=[params json]` flag for all test and bench suites for specifying parameters from the command line.
- DBFV/DCKKS: added a common interface and implementation for each multiparty protocol.
- DBFV/DCKKS: added standalone Encryption-To-Shares (`E2SProtocol`) and Shares-To-Encryption (`S2EProtocol`) protocols for encrypted vs secret-shared domain switching.
- DBFV/DCKKS: generalized the Refresh-and-permute protocol into generic `MaskedTransformProtocol` that accepts an arbitrary linear function.
- DCKKS: public-refresh now takes a target desired output scale, which allows to refresh the ciphertext to the default scale.
- BFV: the moduli of `ringQMul` are not generated based on `N` and`Q`.
- CKKS: added `Parameter` methods computing the required rotations for relevant `Evaluator` operations.
- CKKS: added methods for operating linear-transformation and improved several aspects listed below:
- CKKS: improved the tests for `CoeffsToSlots` and `SlotsToCoeffs`.

#### CKKS Bootstrapping
- The procedure now allows for a more granular parameterization.
- Added flag in bootstrapping parameters for bit-reversed inputs (with bit-reversed output) CoeffsToSlots and SlotsToCoeffs.
- Added optional Arcsine.
- The procedure now uses the new linear-transformation API.
- `CoeffsToSlots` and `SlotsToCoeffs` are now standalone public functions.

#### New CKKS Evaluator methods 
- `RotateHoisted`: evaluate several rotations on a single ciphertext.
- `LinearTransform`: evaluate one or more `PtDiagMatrix` on a ciphertext using `MultiplyByDiagMatrix` or `MultiplyByDiagMatrixBSGS` according to the encoding of `PtDiagMatrix`.
- `MultiplyByDiagMatrix`: multiplies a ciphertext with a `PtDiagMatrix` using n rotations with single hoisting.
- `MultiplyByDiagMatrixBSGS`: multiplies a ciphertext with a `PtDiagMatrix` using 2sqrt(n) rotations with double-hoisting.
- `InnerSumLog`: optimal log approach that works for any value (not only powers of two) and can be parameterized to inner sum batches of values (sub-vectors).
- `InnerSum`: naive approach that is faster for small values but needs more keys.
- `ReplicateLog`: optimal log approach that works for any value (not only powers of two) and can be parameterized to replicate batches of values (sub-vectors).
- `Replicate`: naive approach that is faster for small values but needs more keys.

#### New CKKS Encoder methods
- `PtDiagMatrix`: struct that represents a linear transformation.
- `EncodeDiagMatrixBSGSAtLvl`: encodes a `PtDiagMatrix` at a given level, with a given scale for the BSGS algorithm.
- `EncodeDiagMatrixAtLvl`: encodes a `PtDiagMatrix` at a given level, with a given scale for a naive evaluation.
- `DecodePublic`: adds Gaussian noise of variance floor(sigma * sqrt(2*pi)) before the decoding step (see SECURITY.md).
- `DecodeCoeffsPublic`: adds Gaussian noise of variance floor(sigma * sqrt(2*pi)) before the decoding step (see SECURITY.md).
- `GetErrSTDFreqDom` : get the error standard deviation in the frequency domain (slots).
- `GetErrSTDTimeDom`: get the error standard deviation in the time domain (coefficients).

#### CKKS Fixes
- `MultByi` now correctly sets the output ciphertext scale.
- `Relinearize` now correctly sets the output ciphertext level.
- matrix-vector multiplication now correctly manages ciphertext of higher level than the plaintext matrix.
- matrix-vector encoding now properly works for negative diagonal indexes.

#### Others
- PrecisionStats now includes the standard deviation of the error in the slots and coefficients domains.

## [2.1.1] - 2020-12-23

### Added
- BFV/CKKS: added a check for minimum polynomial degree when creating parameters.
- BFV: added the `bfv.Element.Level` method.
- RING: test for sparse ternary sampler.

### Changed
- BFV/CKKS: pk is now (-as + e, a) instead of (-(as + e), a).
- BFV: harmonized the EvaluationKey setter from `SetRelinKeys` to `Set`.
- CKKS: renamed `BootstrappParams` into `BootstrappingParameters`.
- CKKS: the `Evaluator.DropLevel`, `Parameters.SetLogSlots` and `Element.Copy` methods no longer return errors.
- RING: minimum poly degree modulus is 16 to ensure the NTT correctness.
- RING: isPrime has been replaced by big.ProbablyPrime, which is deterministic for integers < 2^64.

### Fixed
- ALL: reduced cyclomatic complexity of several functions.
- ALL: fixed all instances reported by staticcheck and gosec excluding G103 (audit the use of unsafe).
- ALL: test vectors are now generated using the crypto/rand instead of math/rand package.
- ALL: fixed some unhandled errors.
- BFV/CKKS:  improved the documentation: documented several hard-coded values and fixed typos.
- RING: fixed bias in sparse ternary sampling for some parameters.
- RING: tests for the modular reduction algorithms are now deterministic.

## [2.1.0] - 2020-12-11

### Added
- BFV: special-purpose plaintext types (`PlaintextRingT` or `PlaintextMul`) for optimized ct-pt operations. See bfv/encoder.go and bfv/plaintext.go.
- BFV: allocation-free `Encoder` methods.
- RING: `GenNTTPrimes` now takes the value `Nth` (for Nth primitive root) as input rather than `logN`.

### Changed
- BFV: the `Encoder.DecodeUint64` and `Encoder.DecodeInt64` methods now take the output slice as argument.
- CKKS: API of `Evaluator.RotateColumns` becomes `Evaluator.Rotate`.
- CKKS: the change of variable in `Evaluator.EvaluateCheby` isn't done automatically anymore and the user must do it before calling the function to ensure correctness.
- CKKS: when encoding, the number of slots must now be given in log2 basis. This is to prevent errors that would induced by zero values or non power of two values.
- CKKS: new encoder API : `EncodeAtLvlNew` and `EncodeNTTAtLvlNew`, which allow a user to encode a plaintext at a specific level.

### Removed
- CKKS: removed method `Evaluator.EvaluateChebySpecial`.
- BFV: removed `QiMul` field from `bfv.Parameters`. It is now automatically generated.

## [2.0.0] - 2020-10-07

### Performance
- Global 1.5x speed-up across all arithmetic (this does not include sampling).

### Added
- BFV/CKKS: Added fast encryption (directly in Q without the rescaling by P).
- CKKS: Added full-RNS scale-invariant bootstrapping (<https://eprint.iacr.org/2020/1203>).
- CKKS: Added parameterized tests for a range of experiments.
- CKKS: Added arbitrary precision encoding/decoding.
- CKKS: Added scale invariant polynomial evaluation.
- CKKS: Added encoding/decoding for coefficient packing.
- CKKS The user can now choose to encode a plaintext in or out of the NTT domain (the latter option leads to slightly faster encryptions).
- CKKS: Added secret-key gen with error distribution.
- DBFV: Added collective refresh with arbitrary permutation/linear transformation.
- DCKKS: Added collective refresh with arbitrary permutation/linear transformation.
- RING: Added arbitrary precision complex arithmetic, including cos and sin functions.
- RING: Added polynomial interpolation.
- RING: Added polynomial inversion.
- RING: Extracted interface type Scaler for polynomial coefficient scaling.
- RING: Added type RNSScaler as an efficient, cross-platform implementation of the Scaler interface.

### Changed 
- ALL: all tests now use "require".
- BFV/CKKS: Now parameters without P can be used, but the key-switching is disabled.
- BFV/CKKS: Now parameters do not use methods to access internal values.
- BFV/CKKS: New rotations keys optimized for hoisting rotations of the form (-phi^{-1}(s1)a + phi(s0) + e, a).
- BFV: The Decoder uses the RNSScaler implementation of the Scaler interface to perform the t/Q rescaling.
- CKKS: Simplified the code of the hybrid key-switching (does not affect user experience).
- CKKS: The encoding/decoding operations at level 0 are now 500% faster.
- CKKS: The encoder now accepts slices of complex values with length equal to or smaller than the specified number of slots.
- RING: Improved primes finding.
- RING: All Gaussian sampling now uses Ziggurat sampling.
- RING: Revamped polynomial samplers to make them more memory efficient, consistent user friendly, and to enable parallel sampling.
- RING: The SimpleScaler type now use slightly slower but cross-platform big.Int/Float.
- UTILS: Complete revamp of the PRNG (Blake2b XOF), to make it more user friendly and consistent.

### Removed
- BFV/CKKS: Parameters API generation GenFromLogModuli() and GenFromModuli() have been removed and replaced by Gen().
- CKKS: EvaluatePolyFast(.) and EvaluatePolyEco(.) are replaced by EvaluatePoly(.).
- CKKS: EvaluateChebyFast(.) and EvaluateChebyEco(.) are replaced by EvaluatePolyCheby(.).
- CKKS: EvaluateChebyEcoSpecial(.) and EvaluateChebyFastSpecial(.) are replaced by EvaluatePolyChebySpecial(.).
- RING: The Float128 type was removed due to cross-platform incompatibility.

### Fixes
- BFV: Fixed multiplication that was failing when #Qi != #QMul.
- BFV: Fixed a mempool corruption when encrypting from SK.
- CKKS: The function mulrelin now always returns a fully reduced polynomial.
- CKKS: The encoder now correctly checks that the number of slots is a power of two.
- RING: Prevented a rare case of uint64 overflow during prime sampling.
- RING: Prevented a rare case where two identical primes could be returned when sampling primes.

## [1.3.1] - 2020-02-26
### Added
- BFV/CKKS: Added API for encrypting using a CRP (common reference polynomial).
- BFV/CKKS: Added API for encrypting faster (encrypts zero directly in Q instead of QP and does not need to divide by P).
- BFV/CKKS: Parameters can now be created without the modulus P. Doing so disables all key-switching operations.
- CKKS: Added tests for hoisted rotations.
- RING: Added benchmarks for a NTT using purely Barrett reduction for comparison purposes.
### Changed 
- BFV/CKKS: Changed the switching keys from (-as1 + (s0-s1) + e, a) to (-as1 + s0 + e, a). This does not affect the user experience as it only changes the internal behavior, which is kept consistent. However, Rotation and KeySwitching keys generated with older releases will induce wrong results and will need to be re-generated.
### Fixes
- BFV: Fixed EncryptFromSK that was not correctly wiping the memory pool before using it, which lead to back encryptions.
- BFV: Fixed an index out of bound error that would happen during the multiplication if #QMul > #Qi.
- CKKS: Removed some redundant operations in the hoisted rotations.
- CKKS: MulRelin now always returns a fully reduced ciphertext.
- DCKKS: PCKS and CKS now correctly set the scale of the output ciphertext to the scale of the input ciphertext.
- RING: Fixed GenerateNTTPrimes that could return twice the same prime if the initial value was prime.
- RING: The function context.UniformPoly now samples based on the number of moduli of the context rather than based on the input polynomial.

## [1.3.0] - 2019-12-20
### Added
- All schemes: New switching-keys and key-switching algorithm based on the concept presented in https://eprint.iacr.org/2019/688.pdf.
- All schemes: New marshaling interface for all structures.
- BFV/CKKS: New Parameters structs and API that enable a better customization and fine tuning for specific applications.
- CKKS: New API for hoisted rotations, faster than sequential rotations.
- DBFV/DCKKS: Added collective refresh of a ciphertext (decentralized bootstrapping).
- RING: Added Ziggurat sampling, available from the context.
- RING: Enabled dense and sparse ternary polynomials sampling directly from the context.
- RING: New API enabling "level"-wise polynomial arithmetic.
- RING: New API for modulus switching with flooring and rounding.
- UTILS: The package utils now regroups all the utility methods which were previously duplicated among packages.
### Removed
- BFV/CKKS/DBFV/DCKKS: Removed their respective context. Ring context remains public.
- All schemes: Removed key-switching with bit decomposition. This option will however be re-introduced at a later stage since applications using small parameters can be impacted by this change.
- BFV/CKKS/RING: Removed redundant/irrelevant tests and benchmarks.
- BFV: Removed context QP as it is now not used in the multiplication.
- BFV: Removed int encoder, now only batch encoding is supported.
- CKKS: Modulus switching is moved to the Ring package.
- RING: Removed the algorithms that needed Float128 during the BFV multiplication.
- RING: Removed most wrapping methods for bigInt, which are now replaced by the native math/big package.
- RING: Removed ternary sampler, which is now part of the context.
### Changed
- All schemes: Encryptor, Decryptor, Encoder, Evaluator, KeyGenerator are now interface types.
- All schemes: Improved Godoc and error strings.
- ALl schemes: Greatly reduced the number of methods that could return an error.
- All schemes: New tests and benchmarks with fully supported regex.
- All schemes: Coefficient-wise arithmetic using double slices is now substantially faster.
- BFV/CKKS: Changed the name of the underlying ring contexts. Q now represents the ciphertext modulus (with QMul being the extended ciphertext modulus for BFV) and QP represents modulus of the keys (where P is the special primes used during the new key-switching).
- BFV/CKKS/DBFV/DCKKS: The structures are now created using the parameters instead of the context.
- BFV: Quantization during multiplication does not use Float128 any more, resulting in a substantial speed improvement.
- BFV: BatchEncoder has been renamed Encoder.
- CKKS: The scale is now stored as a float64 instead of a power of 2.
- CKKS: Rounding is applied instead of flooring when a real value is converted to an integer value. This change affects the rescaling and the encoding.
- CKKS: Use of one context for all levels, instead of requiring one ring context per level.
- CKKS: New baby-step giant-step algorithm for evaluating polynomials in standard and Chebyshev basis.
- CKKS: Reduced the number of NTT needed during the encryption.
- CKKS: API for MultConst is now MultByConst.
- BFV/CKKS: New API for the rotation-keys generation.
- DBFV/DCKKS: Complete revamp of the API and interfaces enabling a much easier integration into larger systems.
- DBFV/DCKKS: Improved PCKS and CKS using the concept of the new key-switching technique which enables to reduces the added noise.
- DCKKS: All protocols work for ciphertexts at any levels.
- RING: Faster MulScalarBigint (now similar to MulScalar).
- UTILS: PRNG must be keyed to be forward secure.
### Fixes
- All packages: Corrected typos, godoc and golint.
- CKKS: ciphertext rotation now correctly sets the scale of the output ciphertext.
- DBFV/DCKKS: Correctness is now ensured when the same protocol instance is used to generate multiples shares.

## [1.2.0] - 2019-12-01
Internal version, merged with 1.3.0.

## [1.1.0] - 2019-10-01
### Added
- CHANGELOG.md file.
- BFV: New methods on bfvcontext to access information/stored variables.
- BFV: When creating a new batch encoder, an error will now be returned if the plaintext modulus does not allow NTT.
- CKKS: New methods on ckkscontext to access information/stored variables.
- CKKS: Added marshalling and tests for marshalling.
- CKKS: Added default parameters and parameters marshalling.
- CKKS (API change): encryption can now also be done with the secret-key.
- CKKS (API change): new separate struct for the encoder, that will store a small memory pool for temporary elements.
- BFV/CKKS: Operand interface.
- RING: New method MultByVector.
- RING (API change): New Ternary Sampler, which enables to specify the key distribution {-1, 0, 1} -> [(1-p)/2, p, (1-p)/2]; it is faster than the previous implementation.
- GoDoc for BFV, CKKS, Ring, DBFV, and DCKKS.
- README for BFV and CKKS.
- Code cleaning for all packages.
- Minor optimizations for all packages.

### Changed
- Updated README.md.
- BFV (API change): bfvcontext now only accepts as input a struct of the type Parameters, similar to the one used for DefaultParams.
- BFV (API change): Removed bfvcontext from ciphertexts and plaintexts.
- BFV (API change): Encryption can now also be done with the secret-key.
- BFV: The value logQ in the bfvcontext now stores the bit size of the product of the ciphertext's moduli.
- BFV: The information printed by the tests now better conveys the parameters.
- BFV: Updated and optimized default parameters.
- BFV: Updated example with secure parameters.
- CKKS (API change): ckkscontext now only accepts as input a struct of the type Parameters, similar to the one used for DefaultParams.
- CKKS (API change): ckkscontext now requires as input the moduli chain (in bit size), instead of a generic logQ and levels. This allows more fine-grained control on the rescaling and levels.
- CKKS (API change): removed ckkscontext from ciphertexts and plaintexts.
- CKKS: Updated the value logQ in the ckkscontext, which now stores the bit size of the product of the ciphertext's moduli.
- CKKS: Updated the information printed by the tests to more reflect the parameters.
- BFV/CKKS: Greatly simplified the code related to the rotations.
- BFV/CKKS: Reduced the number of instances where an error could be returned and updated information of the returning errors.
- RING (API change): Changed the API of the validation method of the context to better reflect its purpose.
- BFV/CKKS/RING: The copy function will now copy the input poly on the target poly (previously the target was copied on the input).
- Updates on all packages and tests to comply with the API changes in BFV and CKKS.

### Removed
- The evaluator of both BFV and CKKS cannot operate on two plaintexts anymore. They now always return an element of type Ciphertext.
- The contexts of BFV and CKKS will not store their checksum anymore, nor will the evaluator check for context consistency of the input and output elements.

### Fixed
- Fixed overflow occurring in the basis extension when small and large moduli are used together.

## [1.0.0] - 2019-08-17
### Added
- First public release.
