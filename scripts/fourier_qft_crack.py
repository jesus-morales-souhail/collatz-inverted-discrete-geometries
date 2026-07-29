#!/usr/bin/env python3
"""
DFT → inverse DFT → QFT (simulated) on Collatz f vs g signals.
Derive crack probability + variational frequency weight.

  python scripts/fourier_qft_crack.py
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

from src.fourier_crack import (
    analyze_crack,
    apply_qft,
    dft,
    dft_roundtrip_error,
    encode_amplitudes,
    idft,
    log_signal,
    parity_signal,
    variational_separation,
)

OUT = ROOT / "figures" / "fourier"
RES = ROOT / "results"
OUT.mkdir(parents=True, exist_ok=True)
RES.mkdir(parents=True, exist_ok=True)


def main() -> int:
    seeds = list(range(1, 65))
    steps = 64
    n_qubits = 6  # 64 bins

    print("=" * 70)
    print(" DFT ↔ iDFT → QFT  |  spectral separation f vs g")
    print("=" * 70)

    report, ex = analyze_crack(seeds, steps=steps, n_qubits=n_qubits)
    var = variational_separation(ex["probs_f_qft"], ex["probs_g_qft"], n_iter=250)

    print(f"DFT reconstruction max err f/g: {report.reconstruction_err_f:.2e} / {report.reconstruction_err_g:.2e}")
    print(f"QFT roundtrip fidelity f/g:     {report.qft_fidelity_f:.6f} / {report.qft_fidelity_g:.6f}")
    print(f"DFT TV(f,g):  {report.dft_tv:.4f}   L2: {report.dft_l2:.4f}")
    print(f"QFT TV(f,g):  {report.qft_tv:.4f}   L2: {report.qft_l2:.4f}")
    print(f"QFT entropy f/g: {report.qft_entropy_f:.3f} / {report.qft_entropy_g:.3f}")
    print(f"Separability score: {report.crack_probability:.4f}")
    print(f"Variational J (weighted Δp): {var['J_final']:.4f}")

    # ── figures ──
    # 1) classical DFT power
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    freqs = ex["freqs_rfft"]
    axes[0].semilogy(freqs, ex["Pf_classical"] + 1e-15, label="f normal", color="#1c7ed6")
    axes[0].semilogy(freqs, ex["Pg_classical"] + 1e-15, label="g inverted", color="#e03131")
    axes[0].set_title("Classical power spectrum (mean over seeds)")
    axes[0].set_xlabel("frequency")
    axes[0].set_ylabel("|DFT|² mean")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    k = np.arange(len(ex["probs_f_qft"]))
    axes[1].plot(k, ex["probs_f_qft"], label="QFT Born f", color="#1c7ed6")
    axes[1].plot(k, ex["probs_g_qft"], label="QFT Born g", color="#e03131")
    axes[1].set_title(f"QFT probabilities (n={n_qubits} qubits simulated)")
    axes[1].set_xlabel("basis index k")
    axes[1].set_ylabel("P(k)")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    fig.suptitle(
        f"Spectral separation f vs g: crack_p={report.crack_probability:.3f}  |  QFT TV={report.qft_tv:.3f}",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(OUT / "01_dft_vs_qft_spectra.png", dpi=140, facecolor="white")
    plt.close(fig)

    # 2) roundtrip classical
    x = log_signal(27, 64, inverted=False)
    X = dft(x)
    x2 = idft(X).real
    fig, ax = plt.subplots(figsize=(9, 3.5))
    ax.plot(x, label="original log orbit f", lw=2)
    ax.plot(x2, "--", label="iDFT(DFT(x))", lw=1.5)
    ax.set_title(f"Classical invertibility  max|err|={dft_roundtrip_error(x):.2e}")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "02_dft_inverse_roundtrip.png", dpi=140, facecolor="white")
    plt.close(fig)

    # 3) single seed QFT probs
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
    for ax, inv, title, col in [
        (axes[0], False, "QFT |ψ_f⟩ from log-orbit n0=10", "#1c7ed6"),
        (axes[1], True, "QFT |ψ_g⟩ from log-orbit n0=10", "#e03131"),
    ]:
        s = log_signal(10, 64, inverted=inv)
        amp, nq = encode_amplitudes(s, n_qubits)
        _, pr = apply_qft(amp, nq)
        ax.bar(np.arange(len(pr)), pr, color=col, width=1.0)
        ax.set_title(title)
        ax.set_xlabel("k")
        ax.set_ylabel("P(k)")
    fig.tight_layout()
    fig.savefig(OUT / "03_qft_single_seed.png", dpi=140, facecolor="white")
    plt.close(fig)

    # 4) variational weights
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].plot(var["J_history"], color="#7048e8")
    axes[0].set_title("Variational ascent of J = E_w[p_g − p_f]")
    axes[0].set_xlabel("iteration")
    axes[0].set_ylabel("J")
    axes[0].grid(True, alpha=0.3)
    axes[1].plot(var["weights"], label="w(k)", color="#7048e8")
    axes[1].plot(var["delta"], label="Δp = p_g−p_f", color="#e03131", alpha=0.7)
    axes[1].set_title("Learned frequency weights (variable function)")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "04_variational_weight.png", dpi=140, facecolor="white")
    plt.close(fig)

    # 5) poster
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.text(5, 5.5, "DFT → iDFT → QFT → separability score", ha="center", fontsize=13, fontweight="bold")
    ax.text(
        5, 4.2,
        "1) Señal de órbita f o g (paridad / log)\n"
        "2) DFT clásica + inversa (recuperación exacta numérica)\n"
        "3) Codificar |ψ⟩ y aplicar QFT unitaria simulada → P(k)=|φ_k|²\n"
        "4) Distancia TV/L2 entre espectros f y g → crack_probability\n"
        "5) Función variable w(k;θ)=softmax(θ) maximiza E_w[p_g−p_f]",
        ha="center", va="center", fontsize=10, family="monospace",
        bbox=dict(boxstyle="round", facecolor="#e7f5ff", edgecolor="#1c7ed6", lw=2),
    )
    ax.text(
        5, 1.5,
        f"crack_p = {report.crack_probability:.3f}   |   QFT TV = {report.qft_tv:.3f}   |   "
        f"J* = {var['J_final']:.3f}\n"
        "QFT implemented as dense linear algebra on a classical machine.",
        ha="center", fontsize=9,
        bbox=dict(boxstyle="round", facecolor="#fff3bf", edgecolor="#f59f00"),
    )
    fig.tight_layout()
    fig.savefig(OUT / "00_POSTER_fourier_qft.png", dpi=150, facecolor="white")
    plt.close(fig)

    out = {
        "report": {
            "dft_tv": report.dft_tv,
            "dft_l2": report.dft_l2,
            "qft_tv": report.qft_tv,
            "qft_l2": report.qft_l2,
            "qft_entropy_f": report.qft_entropy_f,
            "qft_entropy_g": report.qft_entropy_g,
            "reconstruction_err_f": report.reconstruction_err_f,
            "reconstruction_err_g": report.reconstruction_err_g,
            "qft_fidelity_f": report.qft_fidelity_f,
            "qft_fidelity_g": report.qft_fidelity_g,
            "crack_probability": report.crack_probability,
            "notes": report.notes,
        },
        "variational": {
            "J_final": var["J_final"],
            "top_bins": var["top_bins"],
            "top_weights": var["top_weights"],
            "top_delta": var["top_delta"],
            "variable_function": var["variable_function"],
        },
        "params": {"seeds": f"1..{seeds[-1]}", "steps": steps, "n_qubits": n_qubits},
    }
    (RES / "fourier_qft_crack.json").write_text(json.dumps(out, indent=2))
    (RES / "FOURIER_QFT_CRACK.md").write_text(
        "\n".join(
            [
                "# Fourier / QFT crack detector",
                "",
                "## Pipeline",
                "1. Classical DFT and inverse DFT on Collatz orbit signals.",
                "2. Encode signal as quantum amplitudes |ψ⟩.",
                "3. Apply unitary QFT (simulated) and read Born probabilities P(k).",
                "4. Compare f vs g spectra → crack_probability score.",
                "5. Variational w(k;θ)=softmax(θ) maximizing E_w[p_g−p_f].",
                "",
                "## Numbers (default run)",
                f"- DFT roundtrip error ~ {report.reconstruction_err_f:.2e}",
                f"- QFT fidelity ~ {report.qft_fidelity_f:.6f}",
                f"- QFT TV(f,g) = {report.qft_tv:.4f}",
                f"- crack_probability = {report.crack_probability:.4f}",
                f"- variational J* = {var['J_final']:.4f}",
                "",
                "## Boundary",
                "Simulated QFT only; not a claim of quantum hardware advantage.",
                "Not a proof of universal inverted divergence.",
                "",
                f"Figures: `{OUT}`",
            ]
        )
    )
    print("wrote", RES / "fourier_qft_crack.json")
    print("figures in", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
