# Conclusion

On the discrete models used in this repository (integer line, modular ring, hyperbolic tree levels, and their product monitors), the inverted map

$$
g(n)=\begin{cases}
3n+1 & \text{if }n\text{ is even}\\
n/2 & \text{if }n\text{ is odd}
\end{cases}
$$

tends to **explore more** of the state space than the normal Collatz map $f$. The difference shows up in several independent checks: visit rates to $1$ on $\mathbb{Z}/N\mathbb{Z}$, growth and path monitors on non-compact components, and spectral separation after DFT/QFT of orbit signals.

This does **not** mean that every inverted orbit diverges. Finite cycles exist. It also does not settle the classical Collatz conjecture for $f$.

Bitstreams extracted from $f$ and $g$ on lattice planes (parity, LSBs, $\mathbb{Z}/N\mathbb{Z}$, growth, product $\mathbb{Z}^2$/$\mathbb{Z}^3$, log-floor) mostly **fail** a self-contained statistical battery at $\alpha=0.01$. That is a negative / mixed result, consistent with the usual fate of chaos-map PRNG candidates. It is not a cryptographic evaluation and not a full NIST/TestU01/PractRand campaign.

Related snapshots: `two_maps_euclidean_spherical.json`, `discrete_geometry_ca.json`, `fourier_qft_crack.json`, `constitutive_planes.json`, `prng_statistical_quality.json`, `summary_reproduced.json`.
