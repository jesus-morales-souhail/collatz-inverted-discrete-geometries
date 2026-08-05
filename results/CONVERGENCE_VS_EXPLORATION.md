# Convergence vs exploration

State variable: one integer $X_t=n\ge 0$.

$$
f(n)=\begin{cases}
n/2 & \text{if }n\text{ is even}\\
3n+1 & \text{if }n\text{ is odd}
\end{cases}
\qquad
g(n)=\begin{cases}
3n+1 & \text{if }n\text{ is even}\\
n/2 & \text{if }n\text{ is odd}
\end{cases}
$$

## What the numbers say

| Check | $f$ | $g$ |
|-------|-----|-----|
| Reach $\{1,2,4\}$ (large sample on $\mathbb{Z}^+$) | 1 | not the classical attractor |
| End near attractor proxy (fixed-step batch) | 0.715 | 0.385 |
| Visit 1 or 4-2-1 (same batch) | 0.74 | 0.385 |
| Visit class 1 on $\mathbb{Z}/N\mathbb{Z}$ | 0.9844 | 0.4375 |
| Mean $\log(1+x_E)$ (product monitor) | 1.821 | 2.218 |
| Mean hyperbolic level (product) | 0.3375 | 0.45 |
| $x_{10}>x_0$ under $g$ | — | 0.4996 |
| $x_{100}>x_0$ under $g$ | — | 0.0155 |
| Hit a cycle within 200 steps under $g$ (sample) | — | 1 |
| Mean PRNG pass rate ($\alpha=0.01$, 9 tests) | 0.04444 | 0.01111 |
| Spectral separability score | — | crack_probability = 0.6038 (DFT TV = 0.7115) |

Log-orbit batch (40 odd seeds, 512 steps), $Y=\log(1+|X|)$:

| | mean $\Delta Y$ | median $\Delta Y$ | mean OLS $b$ in $\Delta Y=a+bY$ |
|--|-----------------|-------------------|----------------------------------|
| $f$ | -0.004724 | -0.004938 | -0.6086 |
| $g$ | -0.002986 | -0.002622 | -1.279 |

## Reading (no extra story)

- On $\mathbb{Z}^+$, $f$ is the map people study for the 4-2-1 attractor. In the samples here it reaches that cycle for essentially all tested seeds.
- $g$ visits 1 less often on the ring, and the product monitors show larger log size and slightly higher tree level under the same step count.
- That is **exploration relative to $f$ on these monitors**, not a theorem that every $g$-orbit goes to infinity. Cycles of $g$ are known; in one sample every seed hit a cycle inside 200 steps. Short-horizon growth under $g$ drops as the horizon grows (see figure 02).
- On $\mathbb{Z}/N\mathbb{Z}$ there is no infinity. Orbits are finite. The question is coverage and visits to the class of 1, not escape.
- Bitstreams from both maps fail almost all tests in the self-contained PRNG battery. Broader integer coverage is not the same as statistical randomness.
- Under $g$, even always maps to odd ($3n+1$ is odd). Parity is not an iid coin.

## Sources

`summary_reproduced.json`, `discrete_geometry_ca.json`, `two_maps_euclidean_spherical.json`,
`fourier_qft_crack.json`, `prng_statistical_quality.json`, plus the log batch above.

Figures: `figures/explore/`.

Reproduce:

```bash
python scripts/convergence_vs_exploration.py
```
