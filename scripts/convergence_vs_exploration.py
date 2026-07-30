#!/usr/bin/env python3
"""
Convergence vs exploration for f and g on the discrete planes used here.

  python scripts/convergence_vs_exploration.py

Pulls numbers from existing result JSON when present, recomputes a few
log-orbit drifts, writes results/CONVERGENCE_VS_EXPLORATION.md and a figure.
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

from src.prng_lattice import log_drift_stats, orbit_array

RES = ROOT / "results"
OUT = ROOT / "figures" / "explore"
OUT.mkdir(parents=True, exist_ok=True)


def jload(name: str) -> dict | None:
    p = RES / name
    if not p.exists():
        return None
    return json.loads(p.read_text())


def main() -> int:
    summary = jload("summary_reproduced.json") or {}
    disc = jload("discrete_geometry_ca.json") or {}
    two = jload("two_maps_euclidean_spherical.json") or {}
    prng = jload("prng_statistical_quality.json") or {}
    four = jload("fourier_qft_crack.json") or {}
    const = jload("constitutive_planes.json") or {}

    # spherical visit-1
    sph_f = summary.get("spherical_frac_visit_1", {}).get("normal_f")
    sph_g = summary.get("spherical_frac_visit_1", {}).get("inverted_g")
    if sph_f is None and "spherical" in disc:
        sph_f = disc["spherical"]["normal_f"]["frac_visit_1"]
        sph_g = disc["spherical"]["inverted_g"]["frac_visit_1"]

    # product monitors
    prod = summary.get("product_4d", {})
    logE_f = prod.get("mean_log1p_E_normal")
    logE_g = prod.get("mean_log1p_E_inverted")
    lvl_f = prod.get("mean_H_level_normal")
    lvl_g = prod.get("mean_H_level_inverted")

    # euclidean attractor proxies (fixed-step sample in discrete_geometry)
    euc = disc.get("euclidean", {})
    euc_f_attr = euc.get("normal_f", {}).get("frac_end_near_attractor")
    euc_g_attr = euc.get("inverted_g", {}).get("frac_end_near_attractor")
    euc_f_vis = euc.get("normal_f", {}).get("frac_visit_1_or_421")
    euc_g_vis = euc.get("inverted_g", {}).get("frac_visit_1_or_421")

    # inverted short-horizon growth and long-horizon cycles
    inv_h = (two.get("discrete") or {}).get("inverted_horizons") or {}
    inv_cyc = (two.get("discrete") or {}).get("inverted_cycles") or {}
    h10 = (inv_h.get("h10") or {}).get("frac_x_h_gt_x0")
    h100 = (inv_h.get("h100") or {}).get("frac_x_h_gt_x0")
    frac_cycle = inv_cyc.get("frac_hit_cycle_within_200_steps")

    # f reach 421
    reach_f = (two.get("discrete") or {}).get("original", {}).get("frac_reach_421")
    if reach_f is None:
        reach_f = summary.get("original_reach_421_frac_0_to_n")

    # spectral
    crack = (four.get("report") or {}).get("crack_probability")
    dft_tv = (four.get("report") or {}).get("dft_tv")

    # prng means
    mean_prng_f = prng.get("mean_pass_rate_f")
    mean_prng_g = prng.get("mean_pass_rate_g")
    if mean_prng_f is None and prng.get("pass_rate_table"):
        mean_prng_f = float(np.mean([r["pass_rate_f"] for r in prng["pass_rate_table"]]))
        mean_prng_g = float(np.mean([r["pass_rate_g"] for r in prng["pass_rate_table"]]))

    # constitutive flux energies (n0=10)
    cf = (const.get("compare_f_g") or {}).get("flux_energy_const_f") or {}
    cg = (const.get("compare_f_g") or {}).get("flux_energy_const_g") or {}

    # log drift over a batch of odd seeds
    seeds = list(range(3, 3 + 40 * 2, 2))
    steps = 512
    dY_f, dY_g, b_f, b_g = [], [], [], []
    for s in seeds:
        sf = log_drift_stats(orbit_array(s, steps, False))
        sg = log_drift_stats(orbit_array(s, steps, True))
        dY_f.append(sf["mean_dY"])
        dY_g.append(sg["mean_dY"])
        b_f.append(sf["reversion_b"])
        b_g.append(sg["reversion_b"])

    log_batch = {
        "n_seeds": len(seeds),
        "steps": steps,
        "mean_dY_f": float(np.mean(dY_f)),
        "mean_dY_g": float(np.mean(dY_g)),
        "median_dY_f": float(np.median(dY_f)),
        "median_dY_g": float(np.median(dY_g)),
        "mean_reversion_b_f": float(np.mean(b_f)),
        "mean_reversion_b_g": float(np.mean(b_g)),
    }

    out = {
        "state_variable": "X_t = n in Z_>=0 (one integer)",
        "maps": {"f": "even->n/2, odd->3n+1", "g": "even->3n+1, odd->n/2"},
        "euclidean": {
            "frac_reach_421_f": reach_f,
            "frac_end_near_attractor_f": euc_f_attr,
            "frac_end_near_attractor_g": euc_g_attr,
            "frac_visit_1_or_421_f": euc_f_vis,
            "frac_visit_1_or_421_g": euc_g_vis,
            "inverted_frac_x10_gt_x0": h10 if h10 is not None else summary.get("inverted_frac_x10_gt_x0"),
            "inverted_frac_x100_gt_x0": h100,
            "inverted_frac_hit_cycle_within_200": frac_cycle,
        },
        "spherical_Z_mod_N": {
            "frac_visit_1_f": sph_f,
            "frac_visit_1_g": sph_g,
            "note": "finite space: no n->infinity",
        },
        "hyperbolic_product": {
            "mean_log1p_E_f": logE_f,
            "mean_log1p_E_g": logE_g,
            "mean_H_level_f": lvl_f,
            "mean_H_level_g": lvl_g,
        },
        "spectral": {"dft_tv": dft_tv, "crack_probability": crack},
        "prng_mean_pass_rate": {"f": mean_prng_f, "g": mean_prng_g, "alpha": 0.01},
        "constitutive_flux_energy_n0_10": {"f": cf, "g": cg},
        "log_orbit_batch": log_batch,
        "reading": {
            "convergence": "f concentrates on the 4-2-1 cycle (global proof open)",
            "exploration": "g covers more on the monitors here; not universal divergence",
            "cycles_g": "finite cycles exist; many seeds hit a cycle within 200 steps in the sample",
            "compact": "on Z/NZ every orbit is eventually periodic",
            "prng": "bitstreams from both maps mostly fail the self-contained battery",
        },
    }
    (RES / "convergence_vs_exploration.json").write_text(json.dumps(out, indent=2))

    # figure: three bars
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8))

    # 1) spherical visit 1
    ax = axes[0]
    if sph_f is not None:
        ax.bar([0, 1], [sph_f, sph_g], color=["#1c7ed6", "#e03131"])
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["f", "g"])
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("fraction")
        ax.set_title(r"$\mathbb{Z}/N\mathbb{Z}$: visit class 1")
        ax.grid(True, axis="y", alpha=0.3)

    # 2) product log E
    ax = axes[1]
    if logE_f is not None:
        ax.bar([0, 1], [logE_f, logE_g], color=["#1c7ed6", "#e03131"])
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["f", "g"])
        ax.set_title(r"product: mean $\log(1+x_E)$")
        ax.grid(True, axis="y", alpha=0.3)

    # 3) mean dY
    ax = axes[2]
    ax.bar([0, 1], [log_batch["mean_dY_f"], log_batch["mean_dY_g"]], color=["#1c7ed6", "#e03131"])
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["f", "g"])
    ax.set_title(r"batch mean $\Delta Y$, $Y=\log(1+|X|)$")
    ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle(r"Convergence ($f$) vs broader coverage ($g$) — same step budgets", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "01_convergence_vs_exploration.png", dpi=140, facecolor="white")
    plt.close(fig)

    # second figure: inverted horizons if present
    if inv_h:
        keys = sorted(inv_h.keys(), key=lambda k: int(k[1:]) if k.startswith("h") else 0)
        xs = [int(k[1:]) for k in keys]
        ys = [inv_h[k]["frac_x_h_gt_x0"] for k in keys]
        fig, ax = plt.subplots(figsize=(7, 3.5))
        ax.plot(xs, ys, "o-", color="#e03131")
        ax.set_xlabel("horizon (steps)")
        ax.set_ylabel(r"fraction with $x_h > x_0$ under $g$")
        ax.set_title("Inverted map: short growth is not permanent escape")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(OUT / "02_g_horizons_growth_fraction.png", dpi=140, facecolor="white")
        plt.close(fig)

    # markdown — short, plain
    def fmt(x, nd=4):
        if x is None:
            return "—"
        if isinstance(x, float):
            return f"{x:.{nd}g}"
        return str(x)

    md = f"""# Convergence vs exploration

State variable: one integer $X_t = n \\ge 0$.

$$
f:\\ \\text{{even}}\\to n/2,\\ \\text{{odd}}\\to 3n+1
\\qquad
g:\\ \\text{{even}}\\to 3n+1,\\ \\text{{odd}}\\to n/2
$$

## What the numbers say

| Check | $f$ | $g$ |
|-------|-----|-----|
| Reach $\\{{1,2,4\\}}$ (large sample on $\\mathbb{{Z}}^+$) | {fmt(reach_f, 4)} | not the classical attractor |
| End near attractor proxy (fixed-step batch) | {fmt(euc_f_attr)} | {fmt(euc_g_attr)} |
| Visit 1 or 4-2-1 (same batch) | {fmt(euc_f_vis)} | {fmt(euc_g_vis)} |
| Visit class 1 on $\\mathbb{{Z}}/N\\mathbb{{Z}}$ | {fmt(sph_f)} | {fmt(sph_g)} |
| Mean $\\log(1+x_E)$ (product monitor) | {fmt(logE_f)} | {fmt(logE_g)} |
| Mean hyperbolic level (product) | {fmt(lvl_f)} | {fmt(lvl_g)} |
| $x_{{10}}>x_0$ under $g$ | — | {fmt(h10 if h10 is not None else summary.get("inverted_frac_x10_gt_x0"))} |
| $x_{{100}}>x_0$ under $g$ | — | {fmt(h100)} |
| Hit a cycle within 200 steps under $g$ (sample) | — | {fmt(frac_cycle)} |
| Mean PRNG pass rate ($\\alpha=0.01$, 9 tests) | {fmt(mean_prng_f)} | {fmt(mean_prng_g)} |
| Spectral separability score | — | crack_probability = {fmt(crack)} (DFT TV = {fmt(dft_tv)}) |

Log-orbit batch ({log_batch['n_seeds']} odd seeds, {log_batch['steps']} steps), $Y=\\log(1+|X|)$:

| | mean $\\Delta Y$ | median $\\Delta Y$ | mean OLS $b$ in $\\Delta Y=a+bY$ |
|--|-----------------|-------------------|----------------------------------|
| $f$ | {log_batch['mean_dY_f']:.4g} | {log_batch['median_dY_f']:.4g} | {log_batch['mean_reversion_b_f']:.4g} |
| $g$ | {log_batch['mean_dY_g']:.4g} | {log_batch['median_dY_g']:.4g} | {log_batch['mean_reversion_b_g']:.4g} |

## Reading (no extra story)

- On $\\mathbb{{Z}}^+$, $f$ is the map people study for the 4-2-1 attractor. In the samples here it reaches that cycle for essentially all tested seeds.
- $g$ visits 1 less often on the ring, and the product monitors show larger log size and slightly higher tree level under the same step count.
- That is **exploration relative to $f$ on these monitors**, not a theorem that every $g$-orbit goes to infinity. Cycles of $g$ are known; in one sample every seed hit a cycle inside 200 steps. Short-horizon growth under $g$ drops as the horizon grows (see figure 02).
- On $\\mathbb{{Z}}/N\\mathbb{{Z}}$ there is no infinity. Orbits are finite. The question is coverage and visits to the class of 1, not escape.
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
"""
    (RES / "CONVERGENCE_VS_EXPLORATION.md").write_text(md)
    print("wrote", RES / "CONVERGENCE_VS_EXPLORATION.md")
    print("sph visit1 f/g", sph_f, sph_g)
    print("log batch", log_batch)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
