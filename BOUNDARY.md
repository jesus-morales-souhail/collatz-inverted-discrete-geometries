# Limits of the claims

## In scope

- Definitions of $f$ (normal Collatz) and $g$ (parity-swapped map) on integers.
- Discrete spaces that imitate Euclidean, spherical and hyperbolic geometry.
- Empirical comparisons of orbits under $f$ and $g$ (growth, path lengths, spectral separation).
- The Markov property for the deterministic iteration $X_{n+1}=g(X_n)$.
- A theorem for the family $x_0=2k+\frac{1}{2}$ under the floor-parity rule $3x+1$ when the integer part is even.
- Pedagogy: binary plots, elementary cellular automata, train diagrams.
- The discrete Fourier/Ohm product on scalar fields defined on the cell line (constant $\kappa$).

## Out of scope

- A proof of the classical Collatz conjecture on all of $\mathbb{Z}^+$.
- A claim that every inverted integer orbit diverges (counterexamples: finite cycles).
- A canonical notion of parity on a smooth Riemannian manifold without discretisation.
- Deriving arithmetic divergence from Euclid’s postulates alone.
- Identifying $f$ or $g$ with linear heat or current transport.
- Any financial application or investment claim.
- Quantum computational advantage; the QFT used here is a classical unitary matrix of size $2^n$.

## Short table

| Statement | Status |
|-----------|--------|
| $f$ reaches 4-2-1 for large tested ranges | Empirical; global proof open |
| $g$ explores more than $f$ on our monitors | Empirical, reproducible here |
| Even half-integers $\to\infty$ under the stated rule | Theorem for that family |
| All orbits of $g$ go to $\infty$ | False |
| Same statements on arbitrary continuous manifolds | Not defined without extra structure |
