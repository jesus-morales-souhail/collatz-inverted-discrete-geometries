# Two jumps on the coverage product

Not a proof of the Collatz conjecture. Raises the level of the residual product already in the repo.

The number $R\approx 0.91$ at depth $28$ on $\{1,\ldots,8000\}$ is residual under a **finite budget**.
It is not an ocean to empty in one symbolic step. Two standard moves still help:

1. Treat coverage as a function of depth (discrete asymptotic / sum over layers).
2. Move the same process to another space where “infinity” is not the question.

## Path 1 — residual as a function of depth

Grow the inverse tree from $1$ (generators $n\mapsto 2n$ and, when legal,
$n\mapsto(n-1)/3$). At each depth $d$ count how much of the window is covered and set

$$
R(d)=1-\frac{N_{\mathrm{cov}}(d)}{N_{\mathrm{win}}}.
$$

Window $\{1,\ldots,8000\}$. Max depth $36$.

| depth $d$ | $R(d)$ |
|-----------|--------|
| $0$ | $0.9999$ |
| $8$ | $0.9979$ |
| $16$ | $0.9848$ |
| $24$ | $0.9440$ |
| $28$ | $0.9085$ |
| $32$ | $0.8684$ |
| $36$ | $0.8219$ |

Mean $T_{\mathrm{eff}}=\log(N_{d+1}/N_d)$ over layers: **0.2223**.
Empirical slope of $\log N_d$ on mid depths: **0.2345** (fit only, not a theorem).

What this does: shows how residual falls when the budget grows. What it does not do: claim
$\lim_{d\to\infty} R(d)=0$ on $\mathbb{Z}^+$. That would be Collatz-hard in residual language.

Figures: `figures/coverage/03_R_of_D_path1.png`.

## Path 2 — change of space

Same forward budget (seeds $1..300$, $100$ steps) on an integer window, versus all residues
on $\mathbb{Z}/128\mathbb{Z}$.

| space | coverage $f$ | coverage $g$ | residual $f$ | residual $g$ |
|-------|--------------|--------------|--------------|--------------|
| $\mathbb{Z}$ window $1..8000$ | 0.0759 | 0.0789 | 0.9241 | 0.9211 |
| $\mathbb{Z}/128\mathbb{Z}$ (all residues) | 1.0000 | 1.0000 | 0.0000 | 0.0000 |

On the ring the state space is finite. Coverage can saturate; there is no $n\to\infty$.
The question changed: visit classes, not escape to infinity. That is the useful content of a
“change of set” here — not a stereographic trick that erases the residual by definition.

Figure: `figures/coverage/04_Z_vs_mod_path2.png`.

## What is still open

- Shape of $R(d)$ for larger $d$ and larger windows (cost grows with the inverse tree).
- Other spaces (e.g. $2$-adic extensions) as further domain changes; not implemented here.
- Any claim that either path settles the classical conjecture.

## Reproduce

```bash
python scripts/coverage_two_jumps.py
```

JSON: `results/coverage_two_jumps.json`.
