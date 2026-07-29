# Collatz inverted vs normal on discrete geometries

Jesús Morales Souhail  
[ORCID 0009-0000-7637-1818](https://orcid.org/0009-0000-7637-1818) · [github.com/jesus-morales-souhail](https://github.com/jesus-morales-souhail)  
July 2026 · independent work · not peer reviewed

Two arithmetic maps, measured on **discrete** models that imitate Euclidean, spherical, and hyperbolic geometry—and on their product—without pretending that parity is well-defined on a smooth manifold without discretization.

| | |
|:--|:--|
| [`START_HERE.md`](START_HERE.md) | how to read this folder |
| [`BOUNDARY.md`](BOUNDARY.md) | what this repo is not allowed to claim |
| [`docs/CLAIMS.md`](docs/CLAIMS.md) | short list of claims with status |
| [`results/CONCLUSION_CANONICA.md`](results/CONCLUSION_CANONICA.md) | one-sentence lab conclusion |

---

## The two equations

**Normal \(f\)** (classical Collatz — attractor 4→2→1 is a *conjecture* on \(\mathbb{Z}^+\)):

\[
f(n)=\begin{cases} n/2 & n\text{ even}\\ 3n+1 & n\text{ odd}\end{cases}
\]

**Inverted \(g\)** (parity-swapped — *more exploration*; **not** always \(\to\infty\)):

\[
g(n)=\begin{cases} 3n+1 & n\text{ even}\\ n/2 & n\text{ odd}\end{cases}
\]

---

## What we measured (and stand behind)

### 1. Geometric exploration crack (empirical, multi-model)

In discrete Euclidean \(\mathbb{Z}\), spherical \(\mathbb{Z}/N\mathbb{Z}\), hyperbolic (value + tree level), and a product monitor \((E,S,H)\):

- \(f\) pushes toward the finite attractor pattern (4-2-1 / visit 1).
- \(g\) explores more: longer paths in non-compact monitors, less attraction to 1 on the ring, higher work \(\int F\,ds\) when growth dominates.

**Canonical one-liner:**  
*The inverted map explores more than the normal map across the discrete geometry models we defined—not only in crypto-style samples and not only in one plane.*

### 2. Markov dependence (proved by construction + checked)

\[
X_{n+1}=g(X_n)\quad\text{with probability }1
\]

\[
P(X_{n+1}\mid X_n,\ldots,X_0)=P(X_{n+1}\mid X_n)
\]

If steps were fully independent of the present, there would be no dynamical rule and no trajectories. Successive parities under \(g\) are **not** i.i.d. coins (e.g. even \(\mapsto\) always odd via \(3n+1\)).

### 3. Theorem: even half-integers under a floor-parity inverted rule

For \(x_0=2k+\frac12\) (\(k\ge 0\)) and “floor even \(\to 3x+1\)”: the form is preserved and \(x_{n+1}>x_n\), hence \(x_n\to+\infty\).

This is a **family theorem**, not a claim about all integer orbits of \(g\).

### 4. Spherical discrete model kills “infinity”

On \(\mathbb{Z}/N\mathbb{Z}\) every orbit stays in a finite set. The analogue of “\(\to 1\)” still makes sense; the analogue of “\(\to\infty\)” does not.

---

## Discrete geometries (no continuous parity cheat)

| Model | Space | Infinity? | Role of \(f\) / \(g\) |
|-------|--------|-----------|----------------------|
| Euclidean | \(\mathbb{Z}^+\) | \(n\to\infty\) possible | Classical arithmetic |
| Spherical | \(\mathbb{Z}/N\mathbb{Z}\) | **No** | Compact attractor statistics |
| Hyperbolic | (value, tree level) | value + capacity \(2^{\mathrm{level}}\) | Growth aligns with expansion under \(g\) |
| 4D product | \((x_E,x_S,x_H^{\mathrm{val}},x_H^{\mathrm{lvl}})\) | \(S\) always bounded | Mode \(f\to 1\) vs mode \(g\to\) escape on non-compact coords |

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

Optional figure scripts (matplotlib):

```bash
python scripts/markov_dependent.py
python scripts/discrete_geometry_ca.py
python scripts/two_maps_euclidean_spherical.py
python scripts/fourier_qft_crack.py   # DFT → iDFT → QFT → crack probability
```

### Fourier / QFT crack channel

Classical **DFT** and **inverse DFT** (numerical roundtrip ~1e-15), then a **simulated unitary QFT** on amplitudes encoded from orbit signals. Separation of spectral measures of \(f\) vs \(g\) yields a `crack_probability` score; a variational weight \(w(k;\theta)=\mathrm{softmax}(\theta)\) maximizes \(\mathbb{E}_w[|p_g-p_f|]\).

See `results/FOURIER_QFT_CRACK.md` and `BOUNDARY.md` (no hardware quantum claim).

---

## Layout

```
src/           maps, Markov, discrete geometries, even-½ theorem
scripts/       reproducible checks and figure generators
tests/         minimal claim tests
results/       JSON/MD snapshots of runs
figures/       binary, train, Markov, discrete_geo, riemann monitors
docs/          claims table
BOUNDARY.md    claim fence
```

---

## What is *not* in this repo

Trading bots, broker code, or profitability claims. Those were lab-only and are out of scope for a mathematics-first public package.

---

## Contact

jmskjym@gmail.com
