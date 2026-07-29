# Fourier / QFT crack detector

## Pipeline
1. Classical DFT and inverse DFT on Collatz orbit signals.
2. Encode signal as quantum amplitudes |ψ⟩.
3. Apply unitary QFT (simulated) and read Born probabilities P(k).
4. Compare f vs g spectra → crack_probability score.
5. Variational w(k;θ)=softmax(θ) maximizing E_w[p_g−p_f].

## Numbers (default run)
- DFT roundtrip error ~ 9.84e-16
- QFT fidelity ~ 1.000000
- QFT TV(f,g) = 0.1958
- crack_probability = 0.6038
- variational J* = 0.0579

## Boundary
Not a hardware quantum claim. Not a trading claim.
Not a proof of universal inverted divergence.

Figures: `/home/ashpokemon/Proyectos/collatz-inverted-discrete-geometries/figures/fourier`