# Conclusion

On the discrete models used in this repository (integer line, modular ring, hyperbolic tree levels, and their product monitors), the inverted map

\[
g(n)=\begin{cases}3n+1 & n\text{ even}\\ n/2 & n\text{ odd}\end{cases}
\]

tends to **explore more** of the state space than the normal Collatz map \(f\). The difference shows up in several independent checks: visit rates to 1 on \(\mathbb{Z}/N\mathbb{Z}\), growth and path monitors on non-compact components, and spectral separation after DFT/QFT of orbit signals.

This does **not** mean that every inverted orbit diverges. Finite cycles exist. It also does not settle the classical Collatz conjecture for \(f\).

Related snapshots: `two_maps_euclidean_spherical.json`, `discrete_geometry_ca.json`, `fourier_qft_crack.json`, `summary_reproduced.json`.
