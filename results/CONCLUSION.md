# Conclusion

## Main result: coverage residual

The central product of this repository is **budgeted coverage**, not a proof of the Collatz conjecture.

1. On the cycle $4\to 2\to 1$ under $f$, each node has a unique successor: **choice entropy $0$**.
2. The inverse tree from $1$ (generators $n\mapsto 2n$ and, when legal, $n\mapsto(n-1)/3$) **grows** with depth; mean branching temperature $T_{\mathrm{eff}}=\log(N_{d+1}/N_d)$ is **positive** in the recorded run.
3. The residual
   $$
   R(D)=1-\frac{\#\{\text{covered states in the window}\}}{\#\{\text{window}\}}
   $$
   is the observable. At depth $D=28$ on $\{1,\ldots,8000\}$, $R\approx 0.91$ (partial coverage under budget). **$R=0$ at finite $D$ is not claimed.**
4. Under the same forward budget, $g$ covers a **slightly larger** fraction of that integer window than $f$.

Full write-up and flags: [`COVERAGE_RESIDUAL.md`](COVERAGE_RESIDUAL.md).  
Reproduce: `python scripts/coverage_residual.py`.

That is the strongest structural result here: a measurable coverage product (nucleus + generators + budget + residual) that works without solving the classical conjecture.

---

## Supporting results

On the discrete models in this repository (integer line, modular ring, hyperbolic tree levels, product monitors), $g$ tends to explore more than $f$ in several independent checks (visits to $1$ on $\mathbb{Z}/N\mathbb{Z}$, growth proxies, spectral separation). That is empirical. Finite cycles of $g$ exist; universal divergence under $g$ is false. The classical conjecture for $f$ remains open.

Bitstreams from $f$ and $g$ on lattice planes mostly **fail** a self-contained statistical battery. Broader integer coverage is not statistical randomness.

Even half-integers under the stated floor rule diverge (theorem for that family only).

Table: [`CONVERGENCE_VS_EXPLORATION.md`](CONVERGENCE_VS_EXPLORATION.md).

Related snapshots: `coverage_residual.json`, `convergence_vs_exploration.json`, `two_maps_euclidean_spherical.json`, `discrete_geometry_ca.json`, `fourier_qft_crack.json`, `constitutive_planes.json`, `prng_statistical_quality.json`, `summary_reproduced.json`.
