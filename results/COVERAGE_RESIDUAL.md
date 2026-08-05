# Coverage residual

Not a Collatz proof. Product: nucleus + generators + budget → coverage and residual \(R\).

## Nucleus (choice entropy)

Cycle \(4\to 2\to 1\to 4\) under \(f\): branching factor 1, choice entropy **0 bits**.

## Inverse tree from 1

- depth max = 28
- window \({1,\ldots,8000}\)
- nodes in window covered: **732**
- coverage fraction: **0.0915**
- residual \(R\): **0.9085**
- mean \(T_{\mathrm{eff}}\): **0.219**

\(N_d\) (first depths): [1, 1, 1, 1, 1, 2, 2, 4, 4, 6, 6, 8]

## Forward coverage (budgeted)

Seeds \(1\ldots300\), steps = 100, \(X_{\max}=8000\).

| | coverage on \(\mathbb{Z}\) window | residual \(R\) | coverage on \(\mathbb{Z}/N\mathbb{Z}\) |
|--|--------------------------------------|----------------|-----------------------------------------------|
| \(f\) | 0.0759 | 0.9241 | 1.0000 |
| \(g\) | 0.0789 | 0.9211 | 1.0000 |

Modular: all residues as seeds, \(N=128\).

## Quality flags (coverage product)

{
  "cycle_choice_entropy_zero": true,
  "inverse_tree_grows": true,
  "inverse_residual_below_half": false,
  "g_covers_more_mod_than_f": false,
  "g_covers_more_Z_than_f": true,
  "mean_T_eff_positive": true
}

Score: **4/6**

## Reading

- Low choice entropy at the cycle = concentrated nucleus.
- Inverse tree expands microstates (\(N_d\), \(T_{\mathrm{eff}}\)).
- Residual \(R\) is the observable; \(R=0\) would be full coverage (Collatz-hard).
- Forward \(g\) vs \(f\) compares coverage under the same budget on \(\mathbb{Z}\) and on the ring.

Figures: `figures/coverage/`.

```bash
python scripts/coverage_residual.py
```
