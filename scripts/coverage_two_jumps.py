#!/usr/bin/env python3
"""
Two ways to raise the coverage product (not a Collatz proof):

  1) Residual as a function of depth: R(d), N_d, T_eff(d)
  2) Same budget on Z vs Z/NZ (change of space)

  python scripts/coverage_two_jumps.py
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

from src.coverage_residual import forward_coverage, residual_curve

RES = ROOT / "results"
OUT = ROOT / "figures" / "coverage"
RES.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)


def write_note(out: dict, X_max: int, max_depth: int, mod_N: int) -> None:
    p1 = out["path1_residual_curve"]
    p2z = out["path2_change_of_space"]["Z_window"]
    p2m = out["path2_change_of_space"]["Z_mod_N"]
    sel = p1["R_at_selected_D"]
    rows = "\n".join("| ${}$ | ${:.4f}$ |".format(k, v) for k, v in sel.items())

    # Keep math out of f-strings where possible so $$ and braces stay intact.
    body = """# Two jumps on the coverage product

Not a proof of the Collatz conjecture. Raises the level of the residual product already in the repo.

The number $R\\approx 0.91$ at depth $28$ on $\\{1,\\ldots,8000\\}$ is residual under a **finite budget**.
It is not an ocean to empty in one symbolic step. Two standard moves still help:

1. Treat coverage as a function of depth (discrete asymptotic / sum over layers).
2. Move the same process to another space where “infinity” is not the question.

## Path 1 — residual as a function of depth

Grow the inverse tree from $1$ (generators $n\\mapsto 2n$ and, when legal,
$n\\mapsto(n-1)/3$). At each depth $d$ count how much of the window is covered and set

$$
R(d)=1-\\frac{N_{\\mathrm{cov}}(d)}{N_{\\mathrm{win}}}.
$$

Window $\\{1,\\ldots,XMAX\\}$. Max depth $DMAX$.

| depth $d$ | $R(d)$ |
|-----------|--------|
ROWS

Mean $T_{\\mathrm{eff}}=\\log(N_{d+1}/N_d)$ over layers: **TEFF**.
Empirical slope of $\\log N_d$ on mid depths: **SLOPE** (fit only, not a theorem).

What this does: shows how residual falls when the budget grows. What it does not do: claim
$\\lim_{d\\to\\infty} R(d)=0$ on $\\mathbb{Z}^+$. That would be Collatz-hard in residual language.

Figures: `figures/coverage/03_R_of_D_path1.png`.

## Path 2 — change of space

Same forward budget (seeds $1..300$, $100$ steps) on an integer window, versus all residues
on $\\mathbb{Z}/MODN\\mathbb{Z}$.

| space | coverage $f$ | coverage $g$ | residual $f$ | residual $g$ |
|-------|--------------|--------------|--------------|--------------|
| $\\mathbb{Z}$ window $1..XMAX$ | COV_ZF | COV_ZG | R_ZF | R_ZG |
| $\\mathbb{Z}/MODN\\mathbb{Z}$ (all residues) | COV_MF | COV_MG | R_MF | R_MG |

On the ring the state space is finite. Coverage can saturate; there is no $n\\to\\infty$.
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
"""
    body = (
        body.replace("XMAX", str(X_max))
        .replace("DMAX", str(max_depth))
        .replace("MODN", str(mod_N))
        .replace("ROWS", rows)
        .replace("TEFF", f"{p1['mean_T_eff']:.4g}")
        .replace("SLOPE", f"{p1['slope_log_N_d_mid']:.4g}")
        .replace("COV_ZF", f"{p2z['coverage_f']:.4f}")
        .replace("COV_ZG", f"{p2z['coverage_g']:.4f}")
        .replace("R_ZF", f"{p2z['R_f']:.4f}")
        .replace("R_ZG", f"{p2z['R_g']:.4f}")
        .replace("COV_MF", f"{p2m['coverage_f']:.4f}")
        .replace("COV_MG", f"{p2m['coverage_g']:.4f}")
        .replace("R_MF", f"{p2m['R_f']:.4f}")
        .replace("R_MG", f"{p2m['R_g']:.4f}")
    )
    (RES / "COVERAGE_TWO_JUMPS.md").write_text(body)


def main() -> int:
    X_max = 8000
    max_depth = 36
    mod_N = 128
    seeds = list(range(1, 301))
    steps = 100

    curve = residual_curve(max_depth=max_depth, X_max=X_max, root=1)
    d = np.array(curve["depth"], dtype=float)
    R = np.array(curve["R_of_d"], dtype=float)
    Nd = np.array(curve["N_d"], dtype=float)

    fwd_f = forward_coverage(seeds, steps, inverted=False, X_max=X_max, mod_N=mod_N)
    fwd_g = forward_coverage(seeds, steps, inverted=True, X_max=X_max, mod_N=mod_N)
    mod_f = forward_coverage(list(range(mod_N)), steps, inverted=False, mod_N=mod_N)
    mod_g = forward_coverage(list(range(mod_N)), steps, inverted=True, mod_N=mod_N)

    mid = (d >= 8) & (d <= max_depth) & (Nd > 0)
    if mid.sum() >= 4:
        slope_logN = float(np.polyfit(d[mid], np.log(Nd[mid]), 1)[0])
    else:
        slope_logN = float("nan")

    out = {
        "path1_residual_curve": {
            "X_max": X_max,
            "max_depth": max_depth,
            "R_final": curve["R_final"],
            "mean_T_eff": curve["mean_T_eff"],
            "slope_log_N_d_mid": slope_logN,
            "R_at_selected_D": {
                str(k): curve["R_of_d"][k]
                for k in (0, 8, 16, 24, 28, 32, max_depth)
                if k <= max_depth
            },
            "N_d_head": curve["N_d"][:16],
        },
        "path2_change_of_space": {
            "Z_window": {
                "X_max": X_max,
                "seeds": len(seeds),
                "steps": steps,
                "coverage_f": fwd_f["coverage_frac"],
                "coverage_g": fwd_g["coverage_frac"],
                "R_f": fwd_f["residual_R"],
                "R_g": fwd_g["residual_R"],
            },
            "Z_mod_N": {
                "N": mod_N,
                "seeds": "all residues",
                "steps": steps,
                "coverage_f": mod_f["coverage_frac_mod"],
                "coverage_g": mod_g["coverage_frac_mod"],
                "R_f": mod_f["residual_R_mod"],
                "R_g": mod_g["residual_R_mod"],
            },
        },
        "note": (
            "R(D) at finite D is not lim R. Modular coverage is a different question "
            "(finite space). Neither path proves the Collatz conjecture."
        ),
    }
    (RES / "coverage_two_jumps.json").write_text(json.dumps(out, indent=2))

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.6))
    axes[0].plot(d, R, "o-", color="#1c7ed6", ms=3)
    axes[0].set_xlabel("depth d")
    axes[0].set_ylabel(r"$R(d)$ on $\{1,\ldots,X\}$")
    axes[0].set_title(rf"Path 1: residual vs budget ($X={X_max}$)")
    axes[0].set_ylim(-0.02, 1.05)
    axes[0].grid(True, alpha=0.3)

    axes[1].semilogy(d, np.maximum(Nd, 1), "s-", color="#2f9e44", ms=3)
    axes[1].set_xlabel("depth d")
    axes[1].set_ylabel(r"$N_d$ (new nodes)")
    axes[1].set_title("Inverse tree layer size")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(range(len(curve["T_eff"])), curve["T_eff"], "^-", color="#e8590c", ms=3)
    axes[2].axhline(0, color="k", lw=0.6)
    axes[2].set_xlabel("depth")
    axes[2].set_ylabel(r"$T_{\mathrm{eff}}=\log(N_{d+1}/N_d)$")
    axes[2].set_title("Branching rate")
    axes[2].grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "03_R_of_D_path1.png", dpi=140, facecolor="white")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    labels = [
        r"$\mathbb{Z}$ window $f$",
        r"$\mathbb{Z}$ window $g$",
        r"$\mathbb{Z}/N\mathbb{Z}$ $f$",
        r"$\mathbb{Z}/N\mathbb{Z}$ $g$",
    ]
    vals = [
        fwd_f["coverage_frac"],
        fwd_g["coverage_frac"],
        mod_f["coverage_frac_mod"],
        mod_g["coverage_frac_mod"],
    ]
    ax.bar(range(4), vals, color=["#1c7ed6", "#e03131", "#4dabf7", "#ff6b6b"])
    ax.set_xticks(range(4))
    ax.set_xticklabels(labels, rotation=15, ha="right", fontsize=8)
    ax.set_ylim(0, 1.08)
    ax.set_ylabel("coverage fraction")
    ax.set_title(r"Path 2: same budget, different space")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "04_Z_vs_mod_path2.png", dpi=140, facecolor="white")
    plt.close(fig)

    write_note(out, X_max, max_depth, mod_N)

    print("R_final", curve["R_final"], "mean_T_eff", curve["mean_T_eff"])
    print("R selected", out["path1_residual_curve"]["R_at_selected_D"])
    print("mod cov f/g", mod_f["coverage_frac_mod"], mod_g["coverage_frac_mod"])
    print("wrote", RES / "coverage_two_jumps.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
