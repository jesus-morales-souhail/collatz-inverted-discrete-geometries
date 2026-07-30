#!/usr/bin/env python3
"""
PRNG-style statistical quality of bitstreams from f and g on lattice planes.

  python scripts/prng_statistical_quality.py

Self-contained subset inspired by NIST SP 800-22. Not full NIST / TestU01 /
PractRand. No cryptographic claim. Mixed or negative results expected.
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

from src.prng_lattice import (
    EXTRACTORS,
    compare_f_g_prng,
    log_drift_stats,
    orbit_array,
    run_battery,
    stream_from_seeds,
)

OUT = ROOT / "figures" / "prng"
RES = ROOT / "results"
OUT.mkdir(parents=True, exist_ok=True)


def main() -> int:
    print("PRNG / statistical quality on lattice planes (f vs g)")

    seeds = list(range(3, 3 + 48 * 2, 2))  # 48 odd seeds
    steps = 2048
    report = compare_f_g_prng(seeds=seeds, steps=steps)

    # ── figures ──────────────────────────────────────────────────────────
    planes = [row["plane"] for row in report["pass_rate_table"]]
    pf = [row["pass_rate_f"] for row in report["pass_rate_table"]]
    pg = [row["pass_rate_g"] for row in report["pass_rate_table"]]

    fig, ax = plt.subplots(figsize=(11, 5))
    x = np.arange(len(planes))
    w = 0.38
    ax.bar(x - w / 2, pf, w, label="f", color="#1c7ed6")
    ax.bar(x + w / 2, pg, w, label="g", color="#e03131")
    ax.set_xticks(x)
    ax.set_xticklabels(planes, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("fraction of tests with p ≥ 0.01")
    ax.set_ylim(0, 1.05)
    ax.axhline(1.0, color="gray", ls=":", lw=0.8)
    ax.set_title("Pass rate of self-contained battery (α=0.01) by lattice plane")
    ax.legend()
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "01_pass_rates_by_plane.png", dpi=140, facecolor="white")
    plt.close(fig)

    # fail counts heatmap-like
    fig, ax = plt.subplots(figsize=(11, 4.5))
    fails_f = [row["n_fail_f"] for row in report["pass_rate_table"]]
    fails_g = [row["n_fail_g"] for row in report["pass_rate_table"]]
    ax.plot(x, fails_f, "o-", label="f fails", color="#1c7ed6")
    ax.plot(x, fails_g, "s-", label="g fails", color="#e03131")
    ax.set_xticks(x)
    ax.set_xticklabels(planes, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("number of failed tests (of 9)")
    ax.set_title("Failures by plane (higher = worse as a PRNG candidate)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "02_failures_by_plane.png", dpi=140, facecolor="white")
    plt.close(fig)

    # log drift comparison across seeds
    drifts_f, drifts_g, revs_f, revs_g = [], [], [], []
    for s in seeds[:32]:
        of = orbit_array(s, steps, False)
        og = orbit_array(s, steps, True)
        sf, sg = log_drift_stats(of), log_drift_stats(og)
        drifts_f.append(sf["mean_dY"])
        drifts_g.append(sg["mean_dY"])
        revs_f.append(sf["reversion_b"])
        revs_g.append(sg["reversion_b"])

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].hist(drifts_f, bins=15, alpha=0.7, label="f", color="#1c7ed6")
    axes[0].hist(drifts_g, bins=15, alpha=0.7, label="g", color="#e03131")
    axes[0].axvline(0, color="k", lw=0.8)
    axes[0].set_xlabel(r"mean $\Delta Y$, $Y=\log(1+|X|)$")
    axes[0].set_title("Empirical log-orbit drift")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].hist(revs_f, bins=15, alpha=0.7, label="f", color="#1c7ed6")
    axes[1].hist(revs_g, bins=15, alpha=0.7, label="g", color="#e03131")
    axes[1].axvline(0, color="k", lw=0.8)
    axes[1].set_xlabel(r"OLS slope $b$ in $\Delta Y = a + b Y$")
    axes[1].set_title("Mean-reversion probe (OU would have $b<0$)")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "03_log_drift_reversion.png", dpi=140, facecolor="white")
    plt.close(fig)

    # progressive dimensions: show pass rate as we add planes in order
    order = [
        "E_parity_Z",
        "E_lsb1",
        "E_lsb4",
        "S_mod64_low6",
        "S_mod256_low8",
        "H_growth",
        "Z2_parity_x_mod",
        "Z3_parity_mod_growth",
        "log2_floor_lsb4",
        "delta_parity",
    ]
    rate_f = {r["plane"]: r["pass_rate_f"] for r in report["pass_rate_table"]}
    rate_g = {r["plane"]: r["pass_rate_g"] for r in report["pass_rate_table"]}
    fig, ax = plt.subplots(figsize=(9, 4))
    xf = np.arange(len(order))
    ax.step(xf, [rate_f[p] for p in order], where="mid", label="f", color="#1c7ed6")
    ax.step(xf, [rate_g[p] for p in order], where="mid", label="g", color="#e03131")
    ax.set_xticks(xf)
    ax.set_xticklabels(order, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("pass rate")
    ax.set_title("Adding lattice dimensions / planes (left → right)")
    ax.set_ylim(-0.05, 1.05)
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "04_progressive_planes.png", dpi=140, facecolor="white")
    plt.close(fig)

    # poster
    mean_pass_f = float(np.mean(pf))
    mean_pass_g = float(np.mean(pg))
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 8)
    ax.text(5, 7.3, "Bitstreams from Collatz maps on lattice planes", ha="center", fontsize=13, fontweight="bold")
    ax.text(
        5, 5.5,
        "Battery (α=0.01): monobit, block freq, runs, serial 2/3,\n"
        "poker m=4, autocorr lags 1/8/32\n"
        f"Mean pass rate  f={mean_pass_f:.2f}   g={mean_pass_g:.2f}\n\n"
        "Planes: parity on Z → LSBs → Z/NZ → growth → product Z2, Z3 → log-floor\n"
        "Not NIST/TestU01/PractRand full suite. No crypto claim.",
        ha="center", va="center", fontsize=10, family="monospace",
        bbox=dict(boxstyle="round", facecolor="#f8f9fa", edgecolor="#495057", lw=1.5),
    )
    ax.text(
        5, 2.3,
        r"Log-orbit: $Y_t=\log(1+|X_t|)$. Empirical drift and OLS slope $b$ in $\Delta Y=a+bY$."
        "\n"
        "Contact with OU-type modelling is this diagnostic, not an external coupling.\n"
        "Skew product / PDMP / MSM: named frameworks only; not fitted here.",
        ha="center", fontsize=9,
        bbox=dict(boxstyle="round", facecolor="#fff9db", edgecolor="#868e96"),
    )
    ax.text(
        5, 0.7,
        f"seeds={report['n_seeds']} odd, steps={steps}  ·  figures/prng/",
        ha="center", fontsize=9, family="monospace",
    )
    fig.tight_layout()
    fig.savefig(OUT / "00_POSTER_prng.png", dpi=150, facecolor="white")
    plt.close(fig)

    # compact JSON (drop heavy per-test detail for some planes if needed — keep full)
    # strip nested full test lists to a compact form for readability
    compact = {
        "alpha": report["alpha"],
        "n_seeds": report["n_seeds"],
        "steps": report["steps"],
        "battery": report["battery"],
        "note": report["note"],
        "pass_rate_table": report["pass_rate_table"],
        "lattice_planes": report["lattice_planes"],
        "framework_pointers": report["framework_pointers"],
        "theoretical_parity": report["theoretical_parity"],
        "log_drift_seed0_f": report["f"]["log_drift_seed0"],
        "log_drift_seed0_g": report["g"]["log_drift_seed0"],
        "mean_pass_rate_f": mean_pass_f,
        "mean_pass_rate_g": mean_pass_g,
        "per_plane_detail": {},
    }
    for name in planes:
        compact["per_plane_detail"][name] = {
            "f": {
                "n_bits": report["f"]["planes"][name]["n_bits"],
                "pass_rate": report["f"]["planes"][name]["pass_rate"],
                "tests": [
                    {
                        "name": t["name"],
                        "p_value": t["p_value"],
                        "pass": t["pass_at_0.01"],
                    }
                    for t in report["f"]["planes"][name]["tests"]
                ],
            },
            "g": {
                "n_bits": report["g"]["planes"][name]["n_bits"],
                "pass_rate": report["g"]["planes"][name]["pass_rate"],
                "tests": [
                    {
                        "name": t["name"],
                        "p_value": t["p_value"],
                        "pass": t["pass_at_0.01"],
                    }
                    for t in report["g"]["planes"][name]["tests"]
                ],
            },
        }

    (RES / "prng_statistical_quality.json").write_text(json.dumps(compact, indent=2))

    # markdown note
    lines = [
        "# PRNG statistical quality on lattice planes",
        "",
        "Self-contained battery inspired by NIST SP 800-22 (monobit, block frequency,",
        "runs, serial, poker, autocorrelation). **Not** a full NIST, TestU01 or PractRand",
        "campaign. **No** cryptographic suitability is claimed.",
        "",
        f"Seeds: {report['n_seeds']} odd integers. Steps per seed: {steps}. Significance: $\\alpha=0.01$.",
        "",
        "## Pass rates (fraction of 9 tests with $p\\ge 0.01$)",
        "",
        "| Plane | pass rate $f$ | pass rate $g$ | bits $f$ | bits $g$ |",
        "|-------|---------------|---------------|----------|----------|",
    ]
    for row in report["pass_rate_table"]:
        lines.append(
            f"| `{row['plane']}` | {row['pass_rate_f']:.2f} | {row['pass_rate_g']:.2f} | "
            f"{row['n_bits_f']} | {row['n_bits_g']} |"
        )
    lines += [
        "",
        f"Mean pass rate: $f$ = {mean_pass_f:.2f}, $g$ = {mean_pass_g:.2f}.",
        "",
        "## Log-orbit drift (Lagarias-type contact with diffusion)",
        "",
        r"Set $Y_t=\log(1+|X_t|)$. Empirical mean step $\mathbb{E}[\Delta Y]$ and OLS",
        r"slope $b$ in $\Delta Y=a+bY$ (a crude mean-reversion probe; OU would have $b<0$).",
        "",
        f"- seed 0 path under $f$: mean ΔY = {report['f']['log_drift_seed0']['mean_dY']:.4g}, "
        f"b = {report['f']['log_drift_seed0']['reversion_b']:.4g}",
        f"- seed 0 path under $g$: mean ΔY = {report['g']['log_drift_seed0']['mean_dY']:.4g}, "
        f"b = {report['g']['log_drift_seed0']['reversion_b']:.4g}",
        "",
        "Under fair independent parity the crude one-step log factors of $f$ and $g$ are",
        "equal after swapping branch labels. Empirical differences come from parity",
        "dependence (under $g$, even always maps to odd) and multi-divisions by 2.",
        "See Lagarias surveys on the $3x+1$ problem.",
        "",
        "## Lattice planes (dimensions added progressively)",
        "",
    ]
    for k, v in report["lattice_planes"].items():
        lines.append(f"- `{k}`: {v}")
    lines += [
        "",
        "## Named frameworks (pointers only)",
        "",
        "- **Skew product / RDS** (Arnold): formal home of the word “product” $(x,y)\\mapsto(Tx,S_x(y))$.",
        "- **PDMP** (Davis 1984): deterministic flow between random jumps — not fitted here.",
        "- **MSM**: discrete states from continuous data (e.g. molecular trajectories) — not built here.",
        "- **Log diffusion**: the honest contact of OU-type language with Collatz-type maps.",
        "",
        "Figures: `figures/prng/`.",
        "",
    ]
    (RES / "PRNG_STATISTICAL_QUALITY.md").write_text("\n".join(lines))

    print(f"mean pass rate f={mean_pass_f:.3f} g={mean_pass_g:.3f}")
    print("wrote", RES / "prng_statistical_quality.json")
    print("figures →", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
