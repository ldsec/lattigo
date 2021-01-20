# Security of Approximate-Numbers Homomorphic Encryption

Classified as an _approximate decryption_ scheme, the CKKS scheme is secure as long as the decryption is only revealed to entities with knowledge of the secret-key. However, some applications may need to share the decrypted value and doing so without taking the appropriate measures may lead to an efficient key-recovery attack as reported by  [Li and Micciancio](https://eprint.iacr.org/2020/1533). 

# Mitigation of the key-recovery attack
Lattigo implements tools to mitigate _Li and Micciancio_'s attack. In particular, the decoding step of the CKKS (and its real number variant R-CKKS) allows the user add a key-independent error **_e_** of standard deviation **_σ_** to the decrypted plaintext before decoding.
Estimating **_σ_** must be done carefully and we provide the following iterative process to do so :
 1. Given a circuit **_C_** taking as inputs vectors **_ω_** of size **_n_** following a distribution **_χ_** and a security parameter **_λ_**, select the appropriate parameters allowing the homomorphic evaluation of **_d(ω)_**, denominated by **_H(d(ω))_**, which includes the encoding, encryption, evaluation, decryption and decoding.
 2. Sample inputs vectors **_ω_** from the distribution **_χ_** and compute the standard deviation **_σ_** in the time domain (coefficient domain) of **_e = d(ω) - H(d(ω))_**. This can be done using the encoder method **GetErrSTDTimeDom(_d(ω)_, _H(d(ω))_, _Δ_)**, where **_Δ_** is the scale of the plaintext after the decryption.
 3. Use the encoder method **DecodePublic** with the parameter **_σ_** to decode plaintexts that will be published. **DecodePublic** will add an error **_e_** of standard deviation **_σ_** bounded by **B = _σ • (2π)^0.5_**.