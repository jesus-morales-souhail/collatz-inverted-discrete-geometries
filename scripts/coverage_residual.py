#!/usr/bin/env python3
"""
Coverage residual diagnostics (local product). Not a Collatz proof.

  python scripts/coverage_residual.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.coverage_residual import compare_coverage_bundle, inverse_tree_levels

RES = ROOT / "results"
OUT = ROOT / "figures" / "coverage"
RES.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)


def main() -> int:
    rep = compare_coverage_bundle(
        X_max=8000,
        inv_depth=28,
        forward_steps=100,
        n_seeds=300,
        mod_N=128,
    )
    inv = rep["inverse_from_1"]
    ff, fg = rep["forward_Z"]["f"], rep["forward_Z"]["g"]
    mf, mg = rep["forward_mod_all_residues"]["f"], rep["forward_mod_all_residues"]["g"]

    (RES / "coverage_residual.json").write_text(json.dumps(rep, indent=2))

    # figures
    Nd = inv["N_d"]
    d = np.arange(len(Nd))
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))

    axes[0].semilogy(d, np.maximum(Nd, 1), "o-", color="#2f9e44")
    axes[0].set_xlabel("depth d from 1")
    axes[0].set_ylabel(r"$N_d$ (new nodes)")
    axes[0].set_title("Inverse tree growth")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(d, inv["covered_cumulative_in_1_X"], "o-", color="#1c7ed6")
    axes[1].axhline(inv["X_max"], color="gray", ls=":", lw=0.8)
    axes[1].set_xlabel("depth d")
    axes[1].set_ylabel(f"covered in 1..{inv['X_max']}")
    axes[1].set_title(f"Coverage cumulative  R={inv['residual_R']:.3f}")
    axes[1].grid(True, alpha=0.3)

    labels = ["Z cov f", "Z cov g", "mod cov f", "mod cov g"]
    vals = [
        ff["coverage_frac"],
        fg["coverage_frac"],
        mf["coverage_frac_mod"],
        mg["coverage_frac_mod"],
    ]
    axes[2].bar(range(4), vals, color=["#1c7ed6", "#e03131", "#4dabf7", "#ff6b6b"])
    axes[2].set_xticks(range(4))
    axes[2].set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
    axes[2].set_ylim(0, 1.05)
    axes[2].set_title("Forward coverage under budget")
    axes[2].grid(True, axis="y", alpha=0.3)

    fig.suptitle("Coverage product (not a Collatz proof)", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "01_coverage_residual.png", dpi=140, facecolor="white")
    plt.close(fig)

    # T_eff
    te = inv["T_eff_log_branch"]
    fig, ax = plt.subplots(figsize=(7, 3.2))
    ax.plot(range(len(te)), te, "s-", color="#e8590c")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("depth")
    ax.set_ylabel(r"$T_{\mathrm{eff}}=\log(N_{d+1}/N_d)$")
    ax.set_title("Branching temperature of inverse tree")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "02_T_eff_branching.png", dpi=140, facecolor="white")
    plt.close(fig)

    md = f"""# Coverage residual

Not a Collatz proof. Product: nucleus + generators + budget → coverage and residual \(R\).

## Nucleus (choice entropy)

Cycle \(4\\to 2\\to 1\\to 4\) under \(f\): branching factor 1, choice entropy **0 bits**.

## Inverse tree from 1

- depth max = {inv['max_depth']}
- window \({{1,\\ldots,{inv['X_max']}}}\)
- nodes in window covered: **{inv['covered_final']}**
- coverage fraction: **{inv['coverage_frac']:.4f}**
- residual \(R\): **{inv['residual_R']:.4f}**
- mean \(T_{{\\mathrm{{eff}}}}\): **{inv['mean_T_eff']:.4g}**

\(N_d\) (first depths): {inv['N_d'][:12]}

## Forward coverage (budgeted)

Seeds \(1\\ldots{ff['n_seeds']}\), steps = {ff['steps']}, \(X_{{\\max}}={ff.get('X_max')}\).

| | coverage on \(\\mathbb{{Z}}\) window | residual \(R\) | coverage on \(\\mathbb{{Z}}/N\\mathbb{{Z}}\) |
|--|--------------------------------------|----------------|-----------------------------------------------|
| \(f\) | {ff['coverage_frac']:.4f} | {ff['residual_R']:.4f} | {mf['coverage_frac_mod']:.4f} |
| \(g\) | {fg['coverage_frac']:.4f} | {fg['residual_R']:.4f} | {mg['coverage_frac_mod']:.4f} |

Modular: all residues as seeds, \(N={mf['mod_N']}\).

## Quality flags (coverage product)

{json.dumps(rep['quality_flags'], indent=2)}

Score: **{rep['quality_score']}**

## Reading

- Low choice entropy at the cycle = concentrated nucleus.
- Inverse tree expands microstates (\(N_d\), \(T_{{\\mathrm{{eff}}}}\)).
- Residual \(R\) is the observable; \(R=0\) would be full coverage (Collatz-hard).
- Forward \(g\) vs \(f\) compares coverage under the same budget on \(\\mathbb{{Z}}\) and on the ring.

Figures: `figures/coverage/`.

```bash
python scripts/coverage_residual.py
```
"""
    (RES / "COVERAGE_RESIDUAL.md").write_text(md)

    print("quality", rep["quality_score"], rep["quality_flags"])
    print("inverse R", inv["residual_R"], "mean T_eff", inv["mean_T_eff"])
    print("Z cov f/g", ff["coverage_frac"], fg["coverage_frac"])
    print("mod cov f/g", mf["coverage_frac_mod"], mg["coverage_frac_mod"])
    print("wrote", RES / "coverage_residual.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
