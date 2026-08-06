# Coverage residual (main result of this repository)

Not a proof of the Collatz conjecture. Among the constructions in this repo, this is the primary structural result.

Product: nucleus + generators + budget $\to$ coverage and residual $R$.

## Nucleus (choice entropy)

On the cycle $4\to 2\to 1\to 4$ under $f$, each node has a unique successor. Branching factor $1$. Choice entropy $0$ bits.

## Inverse tree from $1$

Generators: $n\mapsto 2n$, and $n\mapsto(n-1)/3$ when that value is a positive odd integer.

| quantity | value (this run) |
|----------|------------------|
| max depth $D$ | $28$ |
| window | $\{1,\ldots,8000\}$ |
| covered in window | $732$ |
| coverage fraction | $0.0915$ |
| residual $R$ | $0.9085$ |
| mean $T_{\mathrm{eff}}=\log(N_{d+1}/N_d)$ | $0.219$ |

First layer sizes $N_d$: $1,1,1,1,1,2,2,4,4,6,6,8,\ldots$

$R=0$ at finite $D$ is not claimed.

## Forward coverage (same budget)

Seeds $1,\ldots,300$. Steps $100$. Window $\{1,\ldots,8000\}$. Modular ring $N=128$ with all residues as seeds.

| map | coverage on $\mathbb{Z}$ window | residual $R$ | coverage on $\mathbb{Z}/N\mathbb{Z}$ |
|-----|----------------------------------|--------------|--------------------------------------|
| $f$ | $0.0759$ | $0.9241$ | $1.0000$ |
| $g$ | $0.0789$ | $0.9211$ | $1.0000$ |

On the ring with a full set of seeds both maps saturate. On the integer window, $g$ covers slightly more than $f$ under this budget.

## Core checks (this run)

| check | result |
|-------|--------|
| cycle choice entropy $0$ | yes |
| inverse tree grows | yes |
| mean $T_{\mathrm{eff}}>0$ | yes |
| $g$ covers more than $f$ on the $\mathbb{Z}$ window | yes |
| inverse $R<1/2$ at $D=28$, $X=8000$ | no (not required) |
| $g$ covers more than $f$ on the full modular base | no (both $1$) |

Score on the full flag list: $4/6$. The four core product checks above hold.

## Reproduce

```bash
python scripts/coverage_residual.py
```

Figures: `figures/coverage/`. JSON: `results/coverage_residual.json`.
