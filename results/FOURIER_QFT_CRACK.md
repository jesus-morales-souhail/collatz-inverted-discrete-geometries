# DFT, inverse DFT, and simulated QFT

## Steps

1. Build orbit signals under \(f\) and \(g\) (parity or \(\log(1+X_t)\)).  
2. Classical DFT; recover the signal with the inverse DFT.  
3. Encode a unit amplitude vector and apply the unitary QFT matrix of size \(2^n\).  
4. Compare the two spectral measures (total variation, \(L^2\)).  
5. Optional: weights \(w(k;\theta)=\mathrm{softmax}(\theta)\) maximising \(\mathbb{E}_w[|p_g-p_f|]\).

## Typical numbers (default script)

- DFT reconstruction error on the order of \(10^{-15}\).  
- QFT round-trip fidelity equal to 1 within numerical error.  
- Clear separation of classical parity power spectra between \(f\) and \(g\).  
- A separability score in \([0,1]\) summarised as `crack_probability` in the JSON (spectral only).

## Note

The QFT is implemented as dense linear algebra on a classical machine. No quantum device is involved.

Figures: `figures/fourier/`.
