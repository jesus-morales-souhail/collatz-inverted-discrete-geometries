#!/usr/bin/env python3
"""
Constitutive flux (Fourier/Ohm form) on Collatz cell fields E / S / H / base.

  python scripts/constitutive_planes.py
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

from src.constitutive_discrete import (
    branch_conductivity,
    build_planes,
    collatz_orbit,
    compare_f_g_constitutive,
    constitutive_flux,
    dft_product_identity,
    discrete_gradient,
    potential_log,
)

OUT = ROOT / "figures" / "constitutive"
RES = ROOT / "results"
OUT.mkdir(parents=True, exist_ok=True)


def main() -> int:
    print("constitutive flux on cell fields (E/S/H/base)")

    x = np.linspace(0, 2 * np.pi, 128, endpoint=False)
    u = np.sin(3 * x) + 0.3 * np.cos(5 * x)
    idn = dft_product_identity(u, kappa=2.5, dx=x[1] - x[0])
    print(f"DFT product residual (sin, circular model): {idn['relative_err_circular_model']:.3e}")

    rep = compare_f_g_constitutive(n0=10, steps=64, kappa=1.0, N_mod=64)
    print("flux energy const κ  f:", rep["flux_energy_const_f"])
    print("flux energy const κ  g:", rep["flux_energy_const_g"])
    print("flux energy branch κ f:", rep["flux_energy_branch_f"])
    print("flux energy branch κ g:", rep["flux_energy_branch_g"])

    of = collatz_orbit(10, 64, inverted=False)
    og = collatz_orbit(10, 64, inverted=True)
    pf, pg = build_planes(of), build_planes(og)

    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    planes_f = [pf.cell, pf.euclidean, pf.spherical, pf.hyperbolic]
    planes_g = [pg.cell, pg.euclidean, pg.spherical, pg.hyperbolic]
    titles = ["cell t", "E: log(1+|X|)", "S: X mod N", "H: level-weighted log"]
    for ax, tf, tg, title in zip(axes.ravel(), planes_f, planes_g, titles):
        ax.plot(tf, label="f", color="#1c7ed6")
        ax.plot(tg, label="g", color="#e03131")
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7)
    fig.suptitle(r"Fields $\varphi$ on the cell line from one orbit ($n_0=10$)")
    fig.tight_layout()
    fig.savefig(OUT / "01_four_plane_potentials.png", dpi=140, facecolor="white")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, planes, title in [
        (axes[0], pf, r"flux $=-\kappa\nabla\varphi$  ($f$, $\kappa=1$)"),
        (axes[1], pg, r"flux $=-\kappa\nabla\varphi$  ($g$, $\kappa=1$)"),
    ]:
        for name, u, col in [
            ("E", planes.euclidean, "#1c7ed6"),
            ("S", planes.spherical, "#ae3ec9"),
            ("H", planes.hyperbolic, "#e03131"),
        ]:
            fl = constitutive_flux(u, 1.0)
            ax.plot(fl, label=name, color=col, lw=1.2)
        ax.set_title(title)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.set_xlabel("cell t")
    fig.tight_layout()
    fig.savefig(OUT / "02_flux_const_kappa.png", dpi=140, facecolor="white")
    plt.close(fig)

    u = pf.euclidean
    info = dft_product_identity(u, kappa=1.0)
    k = np.arange(len(u))
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.semilogy(k, np.abs(info["Q"]) + 1e-15, label=r"$|\mathcal{F}\{\mathrm{flux}\}|$", color="#1c7ed6")
    ax.semilogy(
        k,
        np.abs(info["Q_pred"]) + 1e-15,
        "--",
        label=r"$|-\kappa H_k \mathcal{F}\{\varphi\}|$",
        color="#e03131",
    )
    ax.set_title(
        rf"Circular DFT product  relative residual $= {info['relative_err_circular_model']:.2e}$"
    )
    ax.set_xlabel("frequency bin k")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "03_dft_product_identity.png", dpi=140, facecolor="white")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5.2))
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7)
    ax.text(5, 6.4, r"Fourier / Ohm form on discrete cell fields", ha="center", fontsize=13, fontweight="bold")
    ax.text(
        5, 4.5,
        r"$q=-k\nabla T$" + "     " + r"$J=\sigma E=-\sigma\nabla V$" + "\n"
        r"$\mathrm{flux}=-\kappa\,(d*\varphi)$" + "\n"
        r"$\mathcal{F}\{\mathrm{flux}\}_k=-\kappa\,H_k\,\mathcal{F}\{\varphi\}_k$"
        + r"  ($\kappa$ constant, circular)" + "\n\n"
        r"Same product on $\varphi^E$, $\varphi^S$ (mod $N$), $\varphi^H$, and cell index $t$.",
        ha="center", va="center", fontsize=11, family="monospace",
        bbox=dict(boxstyle="round", facecolor="#f8f9fa", edgecolor="#495057", lw=1.5),
    )
    ax.text(
        5, 2.0,
        r"Orbits under $f$ and $g$ define $\varphi$; they are not themselves Ohm/Fourier laws." + "\n"
        r"Parity-dependent $\kappa_t$ is a Hadamard product in space (no single $H_k$ factor).",
        ha="center", fontsize=9,
        bbox=dict(boxstyle="round", facecolor="#fff9db", edgecolor="#868e96"),
    )
    ax.text(
        5, 0.6,
        rf"$n_0=10$: flux energy on E  $f={rep['flux_energy_const_f']['E_euclidean']:.3g}$"
        rf"  $g={rep['flux_energy_const_g']['E_euclidean']:.3g}$",
        ha="center", fontsize=9, family="monospace",
    )
    fig.tight_layout()
    fig.savefig(OUT / "00_POSTER_constitutive.png", dpi=150, facecolor="white")
    plt.close(fig)

    out = {
        "sin_test_dft_identity_err": idn["relative_err_circular_model"],
        "compare_f_g": rep,
        "note": (
            "With constant kappa the flux is -kappa * grad(phi). Under a circular DFT "
            "this factors as -kappa * H_k * F{phi}_k. Branch-dependent kappa is "
            "spatial Hadamard multiplication. Maps f and g only supply the fields."
        ),
    }
    (RES / "constitutive_planes.json").write_text(json.dumps(out, indent=2))
    (RES / "CONSTITUTIVE_PLANES.md").write_text(
        "\n".join(
            [
                "# Constitutive product on discrete planes",
                "",
                r"Fourier: $q=-k\nabla T$. Ohm: $J=\sigma E=-\sigma\nabla V$.",
                "",
                r"On a 1D grid the gradient is convolution with $d$. With constant $\kappa$",
                r"and a circular embedding the DFT gives",
                "",
                "$$",
                r"F\{\mathrm{flux}\}_k=-\kappa H_k F\{u\}_k.",
                "$$",
                "",
                r"The cell line carries fields $\varphi^E$, $\varphi^S$, $\varphi^H$ and the",
                r"base index. Each admits the same flux form. The maps $f$ and $g$ build",
                r"the potentials; they are not linear Ohm/Fourier laws.",
                "",
                "Figures: `figures/constitutive/`.",
                "",
            ]
        )
    )

    print("wrote", RES / "constitutive_planes.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
