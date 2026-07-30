# Collatz inverted vs normal on discrete geometries

Jesús Morales Souhail  
[ORCID 0009-0000-7637-1818](https://orcid.org/0009-0000-7637-1818) · [github.com/jesus-morales-souhail](https://github.com/jesus-morales-souhail)  
July 2026 · independent notes · not peer reviewed

I keep two arithmetic maps and a few discrete spaces that imitate Euclidean, spherical and hyperbolic geometry. Parity is a property of integers. It is not assumed to sit on a smooth manifold without an explicit discretisation.

| | |
|:--|:--|
| [`START_HERE.md`](START_HERE.md) | reading order |
| [`BOUNDARY.md`](BOUNDARY.md) | limits of the claims |
| [`docs/CLAIMS.md`](docs/CLAIMS.md) | claim table |
| [`results/CONCLUSION.md`](results/CONCLUSION.md) | short conclusion |
| [`results/CONVERGENCE_VS_EXPLORATION.md`](results/CONVERGENCE_VS_EXPLORATION.md) | attractor vs coverage on the monitors |

---

## Maps

**Normal $f$** (classical Collatz; global convergence to $4\to 2\to 1$ remains a conjecture on $\mathbb{Z}^+$):

$$
f(n)=\begin{cases}
n/2 & \text{if }n\text{ is even}\\
3n+1 & \text{if }n\text{ is odd}
\end{cases}
$$

**Inverted $g$** (branches swapped):

$$
g(n)=\begin{cases}
3n+1 & \text{if }n\text{ is even}\\
n/2 & \text{if }n\text{ is odd}
\end{cases}
$$

$g$ is **not** claimed to send every orbit to $+\infty$. Cycles exist (for example $0\leftrightarrow 1$ and $6\to 19\to 9\to 4\to 13\to 6$).

---

## Main observations

### Exploration under $g$ versus $f$

On $\mathbb{Z}$, on $\mathbb{Z}/N\mathbb{Z}$, on a discrete hyperbolic tree (value + level), and on a product of path monitors, $g$ tends to explore more of the state space than $f$: less concentration on $1$ on the ring, larger growth proxies on non-compact components when the same number of steps is used. That is an empirical statement about the models implemented here, reproducible from the scripts.

Table of checks: [`results/CONVERGENCE_VS_EXPLORATION.md`](results/CONVERGENCE_VS_EXPLORATION.md).

```bash
python scripts/convergence_vs_exploration.py
```

### Markov structure

$$
X_{n+1}=g(X_n)
$$

with probability one. The future depends only on the present state. Successive parities are not independent fair coins: under $g$, an even integer is always followed by an odd one because $3n+1$ is odd.

### Even half-integers

If $x_0=2k+\frac{1}{2}$ for integer $k\ge 0$ and the rule “floor even $\to 3x+1$” is used, the form is preserved and the sequence is strictly increasing, so $x_n\to+\infty$. This is a theorem for that family only.

### Spherical model

On $\mathbb{Z}/N\mathbb{Z}$ every orbit stays in a finite set. There is no analogue of divergence to infinity. Attraction toward the class of $1$ still makes sense and is stronger under $f$ than under $g$ in the runs recorded here.

---

## Discrete models

| Model | Space | Notes |
|-------|--------|--------|
| Euclidean | $\mathbb{Z}^+$ | classical arithmetic; $n\to\infty$ possible |
| Spherical | $\mathbb{Z}/N\mathbb{Z}$ | compact; no infinity |
| Hyperbolic | value + tree level | capacity grows like $2^{\mathrm{level}}$ |
| Product | $(x_E,x_S,x_H^{\mathrm{val}},x_H^{\mathrm{lvl}})$ | $S$ remains in $\{0,\ldots,N-1\}$ |

---

## Fourier / QFT (simulated)

Orbit signals (parity or $\log(1+X_t)$) are transformed with the classical DFT and recovered with the inverse DFT. A unitary QFT matrix of size $2^n$ is applied to a normalised amplitude vector (classical linear algebra). Separation between the spectral measures of $f$ and $g$ is summarised by a score in $[0,1]$. A weight $w(k;\theta)=\mathrm{softmax}(\theta)$ can be fitted to emphasise bins where the two measures differ.

This is not a claim about quantum hardware.

```bash
python scripts/fourier_qft_crack.py
```

---

## Constitutive product on cell fields

The linear constitutive forms of heat conduction and electrical conduction,

$$
q=-k\nabla T,\qquad J=\sigma E=-\sigma\nabla V,
$$

share the same discrete skeleton on a 1D grid: $\mathrm{flux}=-\kappa\,(d*u)$. With constant $\kappa$ and a circular embedding the DFT factors this as

$$
F\{\mathrm{flux}\}_k=-\kappa\,H_k\,F\{u\}_k,\qquad
H_k=\frac{e^{2\pi i k/N}-1}{\Delta x}.
$$

Orbits under $f$ and $g$ define scalar fields on the cell line (Euclidean log, modular value, level-weighted log, and the cell index). The same flux law applies to each field. The maps themselves are not claimed to be Ohm or Fourier laws. If $\kappa$ depends on parity, the product becomes Hadamard in space and is no longer a single factor $H_k$ in frequency.

```bash
python scripts/constitutive_planes.py
```

Figures: `figures/constitutive/`. Short note: [`results/CONSTITUTIVE_PLANES.md`](results/CONSTITUTIVE_PLANES.md).

---

## PRNG statistical quality on lattice planes

Bitstreams are taken from orbits under $f$ and $g$ on a sequence of lattice planes: parity on $\mathbb{Z}$, LSBs of $X_t$, words on the ring $\mathbb{Z}/N\mathbb{Z}$, growth bits, product coordinates in $\mathbb{Z}^{2}$ and $\mathbb{Z}^{3}$, and LSBs of $\lfloor\log_2(1+|X|)\rfloor$.

A self-contained battery inspired by NIST SP 800-22 (monobit, block frequency, runs, serial, poker, autocorrelation) is run at $\alpha=0.01$. This is **not** a full NIST, TestU01 or PractRand campaign, and **no** cryptographic suitability is claimed. Chaos-map generators often fail serious batteries; the expected scientific outcome is mixed or negative.

On the sample used in this repository the mean pass rate is low for both maps (see [`results/PRNG_STATISTICAL_QUALITY.md`](results/PRNG_STATISTICAL_QUALITY.md)).

The same script records a log-orbit diagnostic $Y_t=\log(1+|X_t|)$: empirical drift $\mathbb{E}[\Delta Y]$ and OLS slope $b$ in $\Delta Y=a+bY$. That is the honest contact with OU-type language (mean reversion as a regression coefficient), not an external coupling of unrelated processes. Skew products, PDMPs and MSMs are named only as frameworks; they are not fitted here.

```bash
python scripts/prng_statistical_quality.py
```

Figures: `figures/prng/`.

---

## Setup

```bash
git clone https://github.com/jesus-morales-souhail/collatz-inverted-discrete-geometries.git
cd collatz-inverted-discrete-geometries
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python tests/test_core.py
python scripts/run_all_checks.py
```

Figure generators (optional):

```bash
python scripts/markov_dependent.py
python scripts/discrete_geometry_ca.py
python scripts/two_maps_euclidean_spherical.py
python scripts/fourier_qft_crack.py
python scripts/constitutive_planes.py
python scripts/prng_statistical_quality.py
python scripts/convergence_vs_exploration.py
```

---

## Layout

```
src/        maps, geometries, Fourier, constitutive, PRNG lattice bits
scripts/    reproducible numerics and figures
tests/      small tests of the maps and theorems
results/    JSON and short notes from runs
figures/    plots and a short train animation
docs/       claim table
```

---

## Contact

jmskjym@gmail.com
