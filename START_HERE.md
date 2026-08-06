# Reading order

1. [`BOUNDARY.md`](BOUNDARY.md) — what is and is not claimed.  
2. [`results/COVERAGE_RESIDUAL.md`](results/COVERAGE_RESIDUAL.md) — **main result** (budgeted residual $R(D)$).  
3. [`results/COVERAGE_TWO_JUMPS.md`](results/COVERAGE_TWO_JUMPS.md) — $R(d)$ vs depth; $\mathbb{Z}$ vs $\mathbb{Z}/N\mathbb{Z}$.  
4. [`results/CONCLUSION.md`](results/CONCLUSION.md).  
5. [`docs/CLAIMS.md`](docs/CLAIMS.md).  
6. `python tests/test_core.py`, `python scripts/coverage_residual.py`, `python scripts/coverage_two_jumps.py`.  
7. Other notes under `results/` and plots under `figures/`.

Main claim of the repo: coverage under a fixed budget (nucleus, inverse tree, residual $R$, forward $f$ vs $g$) works as a measurable product. That is not a proof of the classical Collatz conjecture and not a claim that every inverted orbit diverges.
