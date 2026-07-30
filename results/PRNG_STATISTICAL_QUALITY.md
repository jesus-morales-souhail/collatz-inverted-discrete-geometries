# PRNG statistical quality on lattice planes

Self-contained battery inspired by NIST SP 800-22 (monobit, block frequency,
runs, serial, poker, autocorrelation). **Not** a full NIST, TestU01 or PractRand
campaign. **No** cryptographic suitability is claimed.

Seeds: 48 odd integers. Steps per seed: 2048. Significance: $\alpha=0.01$.

## Pass rates (fraction of 9 tests with $p\ge 0.01$)

| Plane | pass rate $f$ | pass rate $g$ | bits $f$ | bits $g$ |
|-------|---------------|---------------|----------|----------|
| `E_parity_Z` | 0.00 | 0.00 | 98352 | 98352 |
| `E_lsb1` | 0.00 | 0.00 | 98352 | 98352 |
| `E_lsb4` | 0.33 | 0.00 | 393408 | 393408 |
| `S_mod256_low8` | 0.00 | 0.00 | 786816 | 786816 |
| `S_mod64_low6` | 0.00 | 0.00 | 590112 | 590112 |
| `H_growth` | 0.00 | 0.00 | 98304 | 98304 |
| `Z2_parity_x_mod` | 0.00 | 0.00 | 196704 | 196704 |
| `Z3_parity_mod_growth` | 0.00 | 0.11 | 295056 | 295056 |
| `log2_floor_lsb4` | 0.11 | 0.00 | 393408 | 393408 |
| `delta_parity` | 0.00 | 0.00 | 98304 | 98304 |

Mean pass rate: $f$ = 0.04, $g$ = 0.01.

## Log-orbit drift (Lagarias-type contact with diffusion)

Set $Y_t=\log(1+|X_t|)$. Empirical mean step $\mathbb{E}[\Delta Y]$ and OLS
slope $b$ in $\Delta Y=a+bY$ (a crude mean-reversion probe; OU would have $b<0$).

- seed 0 path under $f$: mean ΔY = 0.000109, b = -1.474
- seed 0 path under $g$: mean ΔY = -0.0006769, b = -1.994

Under fair independent parity the crude one-step log factors of $f$ and $g$ are
equal after swapping branch labels. Empirical differences come from parity
dependence (under $g$, even always maps to odd) and multi-divisions by 2.
See Lagarias surveys on the $3x+1$ problem.

## Lattice planes (dimensions added progressively)

- `E_parity_Z`: 1D lattice Z, parity coordinate
- `E_lsb1`: 1D, LSB of X_t
- `E_lsb4`: 1D, 4 LSBs of X_t (more dimensions from same integer)
- `S_mod256_low8`: ring lattice Z/256Z, 8-bit word
- `S_mod64_low6`: ring lattice Z/64Z
- `H_growth`: growth direction bit (tree-level proxy)
- `Z2_parity_x_mod`: product lattice Z^2 bits
- `Z3_parity_mod_growth`: product lattice Z^3 bits
- `log2_floor_lsb4`: scale coordinate floor(log2(1+X)) LSBs
- `delta_parity`: parity of absolute step size

## Named frameworks (pointers only)

- **Skew product / RDS** (Arnold): formal home of the word “product” $(x,y)\mapsto(Tx,S_x(y))$.
- **PDMP** (Davis 1984): deterministic flow between random jumps — not fitted here.
- **MSM**: discrete states from continuous data (e.g. molecular trajectories) — not built here.
- **Log diffusion**: the honest contact of OU-type language with Collatz-type maps.

Figures: `figures/prng/`.
