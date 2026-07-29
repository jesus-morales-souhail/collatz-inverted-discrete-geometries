from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
#!/usr/bin/env python3
"""
Modelos DISCRETOS que imitan cada geometría + Collatz f/g + CA elementales.

No se pega paridad en una variedad continua sin discretizar.
All constructions are discrete; parity is never assumed on a continuum without a grid.

Ecuaciones:
  f (normal):     par → n/2,      impar → 3n+1
  g (invertida):  par → 3n+1,     impar → n/2

Geometrías discretas:
  E  euclidiano  : Z+  (valor n → ∞ posible)
  S  esférico    : Z/NZ  (compacto, sin ∞)
  H  hiperbólico : (valor, nivel en árbol binario regular) — capacidad exp.
  M  modo        : 0=f, 1=g

Producto 4D:
  X = (x_E, x_S, x_H_val, x_H_level)  + modo
  →1:  aplicar f en todas
  →∞:  aplicar g en E y H (no compactas); S siempre mod N

Autómatas celulares elementales (Wolfram) + grieta de paridad Collatz.

  python experiments/collatz_geometry/discrete_geometry_ca.py
"""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).resolve().parent / "figures" / "discrete_geo"
RES = Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
RES.mkdir(parents=True, exist_ok=True)


# ═══════════════════════════════════════════════════════════════════════════
# Collatz f / g
# ═══════════════════════════════════════════════════════════════════════════

def f_normal(n: int) -> int:
    n = int(n)
    if n < 0:
        n = abs(n)
    return n // 2 if n % 2 == 0 else 3 * n + 1


def g_inverted(n: int) -> int:
    n = int(n)
    if n < 0:
        n = abs(n)
    return (3 * n + 1) if n % 2 == 0 else n // 2


def f_mod(n: int, N: int) -> int:
    return f_normal(n) % N


def g_mod(n: int, N: int) -> int:
    return g_inverted(n) % N


# ═══════════════════════════════════════════════════════════════════════════
# 1. Euclidiano: Z+
# ═══════════════════════════════════════════════════════════════════════════

def euclidean_orbit(n0: int, steps: int, inverted: bool) -> np.ndarray:
    fn = g_inverted if inverted else f_normal
    x = max(int(n0), 0)
    seq = [x]
    for _ in range(steps):
        x = fn(x)
        seq.append(x)
    return np.array(seq, dtype=np.int64)


def euclidean_stats(seeds: List[int], steps: int = 40) -> dict:
    out = {}
    for inv, name in [(False, "normal_f"), (True, "inverted_g")]:
        hit_1 = 0
        hit_cycle_421 = 0
        diverged = 0  # max > 1e6 or final > 100*start
        cycled_small = 0
        for s in seeds:
            if s == 0 and not inv:
                hit_1 += 1
                continue
            seq = euclidean_orbit(s, steps, inv)
            if 1 in seq or (4 in seq and 2 in seq):
                hit_cycle_421 += 1
            if seq[-1] == 1 or (len(set(seq[-6:])) <= 3 and 1 in seq[-6:]):
                hit_1 += 1
            if seq.max() > 10**6 or (s > 0 and seq[-1] > 100 * s and seq[-1] > seq[0]):
                diverged += 1
            # cycle detection in last part
            if len(set(seq[-10:])) < 8 and seq.max() < 10**5:
                cycled_small += 1
        n = len(seeds)
        out[name] = {
            "frac_visit_1_or_421": hit_cycle_421 / n,
            "frac_end_near_attractor": hit_1 / n,
            "frac_diverged_proxy": diverged / n,
            "frac_small_cycle_proxy": cycled_small / n,
        }
    return out


# ═══════════════════════════════════════════════════════════════════════════
# 2. Esférico: Z/NZ
# ═══════════════════════════════════════════════════════════════════════════

def spherical_orbit(n0: int, N: int, steps: int, inverted: bool) -> np.ndarray:
    fn = g_mod if inverted else f_mod
    x = int(n0) % N
    seq = [x]
    for _ in range(steps):
        x = fn(x, N)
        seq.append(x)
    return np.array(seq, dtype=np.int64)


def spherical_stats(N: int = 128, steps: int = 80) -> dict:
    out = {}
    for inv, name in [(False, "normal_f"), (True, "inverted_g")]:
        visit_1 = 0
        n_cycles = []
        unique_frac = []
        for s in range(N):
            seq = spherical_orbit(s, N, steps, inv)
            if 1 in seq:
                visit_1 += 1
            unique_frac.append(len(set(seq)) / N)
            # cycle length from rho algorithm on deterministic map
            # detect period at end
            period = 1
            for p in range(1, min(N + 1, 40)):
                if seq[-1] == seq[-1 - p]:
                    # check full period
                    if all(seq[-1 - i] == seq[-1 - i - p] for i in range(p)):
                        period = p
                        break
            n_cycles.append(period)
        out[name] = {
            "N": N,
            "frac_visit_1": visit_1 / N,
            "mean_unique_frac": float(np.mean(unique_frac)),
            "mean_period_proxy": float(np.mean(n_cycles)),
            "no_infinity": True,
            "note": "Espacio finito: no existe n→∞",
        }
    return out


# ═══════════════════════════════════════════════════════════════════════════
# 3. Hiperbólico: valor + nivel en árbol (crecimiento exp de capacidad)
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class HypState:
    value: int
    level: int  # profundidad en árbol binario regular (capacidad ~ 2^level)

    def capacity(self) -> int:
        return 2 ** max(self.level, 0)


def hyperbolic_step(state: HypState, inverted: bool) -> HypState:
    """
    Valor sigue f o g.
    Nivel: si el valor crece, sube de nivel (más espacio); si decrece, baja (mín 0).
    Imita expansión hiperbólica: más vecinos/capacidad al alejarse del origen.
    """
    fn = g_inverted if inverted else f_normal
    v0 = state.value
    v1 = fn(v0)
    if v1 > v0:
        lvl = state.level + 1
    elif v1 < v0:
        lvl = max(0, state.level - 1)
    else:
        lvl = state.level
    return HypState(value=v1, level=lvl)


def hyperbolic_orbit(n0: int, steps: int, inverted: bool) -> List[HypState]:
    s = HypState(value=max(int(n0), 0), level=0)
    seq = [s]
    for _ in range(steps):
        s = hyperbolic_step(s, inverted)
        seq.append(s)
    return seq


def hyperbolic_stats(seeds: List[int], steps: int = 30) -> dict:
    out = {}
    for inv, name in [(False, "normal_f"), (True, "inverted_g")]:
        final_levels = []
        final_vals = []
        max_caps = []
        for s0 in seeds:
            orb = hyperbolic_orbit(s0, steps, inv)
            final_levels.append(orb[-1].level)
            final_vals.append(orb[-1].value)
            max_caps.append(max(h.capacity() for h in orb))
        out[name] = {
            "mean_final_level": float(np.mean(final_levels)),
            "mean_final_value": float(np.mean(final_vals)),
            "mean_log1p_final_value": float(np.mean(np.log1p(final_vals))),
            "mean_max_capacity": float(np.mean(max_caps)),
            "mean_log2_max_capacity": float(np.mean(np.log2(np.maximum(max_caps, 1)))),
        }
    return out


# ═══════════════════════════════════════════════════════════════════════════
# 4. Producto 4D (E, S, H_val, H_level) + modo
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class State4D:
    x_E: int
    x_S: int
    x_H_val: int
    x_H_level: int
    mode: int  # 0 = normal f, 1 = inverted g


def step_4d(X: State4D, N: int, force_mode: int | None = None) -> State4D:
    """
    force_mode: None → usa X.mode; 0 fuerza f; 1 fuerza g.
    Esférico siempre mod N (acotado).
    E y H usan f o g según modo.
    """
    m = X.mode if force_mode is None else force_mode
    inv = m == 1
    fn = g_inverted if inv else f_normal
    fnm = g_mod if inv else f_mod

    e = fn(X.x_E)
    s = fnm(X.x_S, N)
    h = hyperbolic_step(HypState(X.x_H_val, X.x_H_level), inv)
    return State4D(e, s, h.value, h.level, m)


def orbit_4d(n0: int, N: int, steps: int, mode: int) -> List[State4D]:
    X = State4D(n0, n0 % N, n0, 0, mode)
    seq = [X]
    for _ in range(steps):
        X = step_4d(X, N, force_mode=mode)
        seq.append(X)
    return seq


def product_4d_stats(seeds: List[int], N: int = 64, steps: int = 35) -> dict:
    out = {}
    for mode, name in [(0, "bundle_to_1_normal"), (1, "bundle_to_inf_inverted")]:
        e_final, s_final, h_level, h_val = [], [], [], []
        for s0 in seeds:
            orb = orbit_4d(s0, N, steps, mode)
            last = orb[-1]
            e_final.append(last.x_E)
            s_final.append(last.x_S)
            h_level.append(last.x_H_level)
            h_val.append(last.x_H_val)
        out[name] = {
            "mean_log1p_E": float(np.mean(np.log1p(e_final))),
            "mean_S": float(np.mean(s_final)),  # always in 0..N-1
            "max_S": int(np.max(s_final)),
            "mean_H_level": float(np.mean(h_level)),
            "mean_log1p_H_val": float(np.mean(np.log1p(h_val))),
            "frac_E_eq_1": float(np.mean([e == 1 for e in e_final])),
            "frac_E_gt_1e4": float(np.mean([e > 1e4 for e in e_final])),
            "S_always_bounded": True,
        }
    return out


# ═══════════════════════════════════════════════════════════════════════════
# 5. Autómatas celulares elementales (Wolfram)
# ═══════════════════════════════════════════════════════════════════════════

def wolfram_rule_table(rule: int) -> Dict[Tuple[int, int, int], int]:
    """Byte de regla → mapa vecindad 3 → bit."""
    table = {}
    for i, neigh in enumerate(
        [(1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0), (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)]
    ):
        bit = (rule >> (7 - i)) & 1
        table[neigh] = bit
    return table


def ca_evolve(rule: int, size: int = 121, steps: int = 60, seed_center: bool = True) -> np.ndarray:
    """Evoluciona CA 1D; filas = tiempo, columnas = espacio. Condiciones periódicas."""
    table = wolfram_rule_table(rule)
    grid = np.zeros((steps, size), dtype=np.uint8)
    if seed_center:
        grid[0, size // 2] = 1
    else:
        rng = np.random.default_rng(0)
        grid[0] = rng.integers(0, 2, size)
    for t in range(steps - 1):
        for i in range(size):
            L = grid[t, (i - 1) % size]
            C = grid[t, i]
            R = grid[t, (i + 1) % size]
            grid[t + 1, i] = table[(int(L), int(C), int(R))]
    return grid


def collatz_parity_tape(n0: int, steps: int, inverted: bool) -> np.ndarray:
    """
    Grieta: cinta 1D de paridades a lo largo de la órbita Collatz
    (Markov temporal, no espacial — comparada con CA).
    """
    fn = g_inverted if inverted else f_normal
    x = int(n0)
    bits = []
    for _ in range(steps):
        bits.append(x % 2)
        x = fn(x)
    return np.array(bits, dtype=np.uint8)


def parity_ca_hybrid(n0: int, rule: int, width: int = 31, steps: int = 40) -> np.ndarray:
    """
    Grieta experimental:
    - fila t: ventana de paridades alrededor del estado Collatz invertido
      (bits de n, n>>1, n>>2, ... embebidos en la línea)
    - luego un paso de regla Wolfram sobre esa configuración
    Muestra interacción dependencia aritmética + dependencia espacial.
    """
    grid = np.zeros((steps, width), dtype=np.uint8)
    x = int(n0)
    for t in range(steps):
        # rellenar con bits de x (LSB al centro)
        bits = [(x >> k) & 1 for k in range(width)]
        # centrar LSB
        row = np.zeros(width, dtype=np.uint8)
        mid = width // 2
        for k, b in enumerate(bits[: mid + 1]):
            if mid - k >= 0:
                row[mid - k] = b
        grid[t] = row
        # un paso CA
        table = wolfram_rule_table(rule)
        new = np.zeros(width, dtype=np.uint8)
        for i in range(width):
            L, C, R = row[(i - 1) % width], row[i], row[(i + 1) % width]
            new[i] = table[(int(L), int(C), int(R))]
        # actualizar x con Collatz invertida (aritmética)
        x = g_inverted(x)
        # mezclar: siguiente fila se reconstruye de x, no de new
        # (la CA se visualiza como capa paralela)
        if t + 1 < steps:
            grid[t] = np.maximum(grid[t], new)  # mostrar huella CA
    return grid


# ═══════════════════════════════════════════════════════════════════════════
# Figuras
# ═══════════════════════════════════════════════════════════════════════════

def plot_euclidean(seeds_demo=(6, 10, 27)):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, inv, title in [
        (axes[0], False, "Euclidiano Z — NORMAL f → atractor"),
        (axes[1], True, "Euclidiano Z — INVERTIDA g → crece/cicla"),
    ]:
        for n0 in seeds_demo:
            seq = euclidean_orbit(n0, 40, inv)
            ax.semilogy(np.maximum(seq, 1), "o-", ms=2, lw=1, label=f"n0={n0}")
        ax.set_xlabel("t")
        ax.set_ylabel("n (log)")
        ax.set_title(title)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "01_euclidean_Z.png", dpi=140, facecolor="white")
    plt.close(fig)


def plot_spherical(N=64):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, inv, title in [
        (axes[0], False, f"Esférico Z/{N}Z — NORMAL f mod N"),
        (axes[1], True, f"Esférico Z/{N}Z — INVERTIDA g mod N"),
    ]:
        for n0 in [1, 3, 7, 15, 31]:
            seq = spherical_orbit(n0, N, 50, inv)
            ax.plot(seq, "o-", ms=2, lw=1, label=f"n0={n0}")
        ax.axhline(1, color="k", ls="--", lw=0.8, alpha=0.5)
        ax.set_ylim(-1, N)
        ax.set_xlabel("t")
        ax.set_ylabel(f"n mod {N}")
        ax.set_title(title + "\n(sin ∞: espacio finito)")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "02_spherical_modN.png", dpi=140, facecolor="white")
    plt.close(fig)


def plot_hyperbolic():
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, inv, title in [
        (axes[0], False, "Hiperbólico — NORMAL: valor + nivel"),
        (axes[1], True, "Hiperbólico — INVERTIDA: valor + nivel"),
    ]:
        for n0 in [6, 10, 12, 20]:
            orb = hyperbolic_orbit(n0, 25, inv)
            vals = [h.value for h in orb]
            lvls = [h.level for h in orb]
            ax.plot(vals, lvls, "o-", ms=3, lw=1, label=f"n0={n0}")
        ax.set_xlabel("valor n")
        ax.set_ylabel("nivel (capacidad 2^nivel)")
        ax.set_title(title)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
        if inv:
            ax.set_xscale("log")
    fig.tight_layout()
    fig.savefig(OUT / "03_hyperbolic_tree.png", dpi=140, facecolor="white")
    plt.close(fig)


def plot_product_4d(N=64, steps=30):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    seeds = [6, 10, 12, 15, 21]
    for mode, color, name in [(0, "#1c7ed6", "→1 normal f"), (1, "#e03131", "→∞ inverted g")]:
        for n0 in seeds:
            orb = orbit_4d(n0, N, steps, mode)
            E = [x.x_E for x in orb]
            S = [x.x_S for x in orb]
            Hl = [x.x_H_level for x in orb]
            Hv = [x.x_H_val for x in orb]
            axes[0, 0].plot(np.log1p(E), color=color, alpha=0.5, lw=1)
            axes[0, 1].plot(S, color=color, alpha=0.5, lw=1)
            axes[1, 0].plot(Hl, color=color, alpha=0.5, lw=1)
            axes[1, 1].plot(np.log1p(Hv), color=color, alpha=0.5, lw=1)
    axes[0, 0].set_title("x_E euclidiano log1p")
    axes[0, 1].set_title(f"x_S esférico mod {N} (siempre acotado)")
    axes[1, 0].set_title("x_H nivel hiperbólico")
    axes[1, 1].set_title("x_H valor log1p")
    for ax in axes.ravel():
        ax.grid(True, alpha=0.3)
        ax.set_xlabel("t")
    # leyenda manual
    axes[0, 0].plot([], [], color="#1c7ed6", label="modo f →1")
    axes[0, 0].plot([], [], color="#e03131", label="modo g →∞")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Producto 4D discreto (E, S, H_level, H_val) — dos modos", fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT / "04_product_4D.png", dpi=140, facecolor="white")
    plt.close(fig)


def plot_ca_rules():
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    for ax, rule, title in [
        (axes[0, 0], 30, "Rule 30 — caos (clase III)"),
        (axes[0, 1], 90, "Rule 90 — Sierpinski (clase III)"),
        (axes[1, 0], 110, "Rule 110 — estructuras (clase IV)"),
        (axes[1, 1], 184, "Rule 184 — tráfico"),
    ]:
        g = ca_evolve(rule, size=101, steps=50)
        ax.imshow(g, cmap="binary", interpolation="nearest", aspect="auto")
        ax.set_title(title)
        ax.set_xlabel("espacio")
        ax.set_ylabel("tiempo ↓")
    fig.suptitle("Autómatas celulares elementales (Wolfram) — dependencia local Markov", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "05_cellular_automata.png", dpi=140, facecolor="white")
    plt.close(fig)


def plot_collatz_parity_vs_ca():
    """Grieta: cinta de paridades Collatz vs CA Rule 90 desde misma semilla bit."""
    fig, axes = plt.subplots(2, 2, figsize=(11, 6))
    # parity tapes
    for ax, inv, title in [
        (axes[0, 0], False, "Paridades NORMAL f (semillas 1..40)"),
        (axes[0, 1], True, "Paridades INVERTIDA g (semillas 1..40)"),
    ]:
        M = np.array([collatz_parity_tape(s, 50, inv) for s in range(1, 41)])
        ax.imshow(M, cmap="binary", aspect="auto", interpolation="nearest")
        ax.set_xlabel("t")
        ax.set_ylabel("semilla")
        ax.set_title(title)
    # CA
    axes[1, 0].imshow(ca_evolve(90, 81, 40), cmap="binary", aspect="auto")
    axes[1, 0].set_title("Rule 90 (espacial, Markov local)")
    axes[1, 1].imshow(ca_evolve(30, 81, 40), cmap="binary", aspect="auto")
    axes[1, 1].set_title("Rule 30 (caos espacial)")
    for ax in axes[1]:
        ax.set_xlabel("espacio")
        ax.set_ylabel("tiempo")
    fig.suptitle(
        "Collatz: dependencia aritmética en el tiempo  |  CA: dependencia espacial en la línea",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(OUT / "06_parity_collatz_vs_CA.png", dpi=140, facecolor="white")
    plt.close(fig)


def plot_poster(eu, sp, hy, pr):
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.text(5, 9.5, "GEOMETRÍAS DISCRETAS + COLLATZ f/g + CA", ha="center", fontsize=14, fontweight="bold")
    ax.text(
        5, 8.6,
        "No se pega paridad en variedad continua sin discretizar.\n"
        "Modelos: Z (E) · Z/NZ (S) · árbol (H) · producto 4D · autómatas de Wolfram.",
        ha="center", fontsize=9,
        bbox=dict(boxstyle="round", facecolor="#e7f5ff", edgecolor="#1c7ed6"),
    )
    boxes = [
        (0.3, 5.2, "EUCLIDEO Z", f"f→421: {eu['normal_f']['frac_visit_1_or_421']:.2f}\n"
         f"g diverg~: {eu['inverted_g']['frac_diverged_proxy']:.2f}"),
        (2.7, 5.2, "ESFERICO Z/N", f"f visit 1: {sp['normal_f']['frac_visit_1']:.2f}\n"
         f"g visit 1: {sp['inverted_g']['frac_visit_1']:.2f}\n(sin ∞)"),
        (5.1, 5.2, "HIPERBOLICO", f"f nivel: {hy['normal_f']['mean_final_level']:.1f}\n"
         f"g nivel: {hy['inverted_g']['mean_final_level']:.1f}\n"
         f"g log val: {hy['inverted_g']['mean_log1p_final_value']:.1f}"),
        (7.5, 5.2, "PRODUCTO 4D", f"→1 log E: {pr['bundle_to_1_normal']['mean_log1p_E']:.2f}\n"
         f"→∞ log E: {pr['bundle_to_inf_inverted']['mean_log1p_E']:.2f}\n"
         f"S acotada siempre"),
    ]
    for x, y, title, body in boxes:
        ax.add_patch(plt.Rectangle((x, y), 2.2, 2.8, facecolor="#f8f9fa", edgecolor="#343a40", lw=1.5))
        ax.text(x + 1.1, y + 2.4, title, ha="center", fontsize=9, fontweight="bold")
        ax.text(x + 1.1, y + 1.2, body, ha="center", va="center", fontsize=8, family="monospace")
    ax.text(
        5, 4.3, "Markov: en todos los modelos X_{t+1} depende solo del presente (local en CA, aritmético en Collatz).",
        ha="center", fontsize=9,
    )
    ax.text(
        5, 2.8,
        "Grieta: en S no hay ∞; en E/H la invertida empuja escape; el producto 4D separa\n"
        "atracción a 1 (f) vs escape (g) manteniendo la componente esférica finita.",
        ha="center", fontsize=9,
        bbox=dict(boxstyle="round", facecolor="#fff3bf", edgecolor="#f59f00", lw=2),
    )
    ax.text(
        5, 1.2,
        "Enlace lab Morales-Souhail (github): celdas medibles, dependencia local,\n"
        "Parity is only used after discretisation.",
        ha="center", fontsize=8, style="italic", color="#495057",
    )
    fig.tight_layout()
    fig.savefig(OUT / "00_POSTER_discrete_geometry.png", dpi=150, facecolor="white")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

def main() -> int:
    print("=" * 72)
    print(" GEOMETRÍAS DISCRETAS + COLLATZ f/g + CA ELEMENTALES")
    print("=" * 72)

    seeds = list(range(1, 201))
    print("[1] Euclidiano…")
    eu = euclidean_stats(seeds, steps=45)
    print("   ", eu)

    print("[2] Esférico…")
    sp = spherical_stats(N=128, steps=100)
    print("   ", sp)

    print("[3] Hiperbólico…")
    hy = hyperbolic_stats(seeds[:80], steps=30)
    print("   ", hy)

    print("[4] Producto 4D…")
    pr = product_4d_stats(seeds[:60], N=64, steps=35)
    print("   ", pr)

    print("[5] Figuras…")
    plot_euclidean()
    plot_spherical(64)
    plot_hyperbolic()
    plot_product_4d(64, 30)
    plot_ca_rules()
    plot_collatz_parity_vs_ca()
    plot_poster(eu, sp, hy, pr)

    report = {
        "principle": "No continuous parity on Riemann manifold without discretization. Discrete models only.",
        "equations": {
            "f_normal": "even→n/2, odd→3n+1",
            "g_inverted": "even→3n+1, odd→n/2",
        },
        "euclidean": eu,
        "spherical": sp,
        "hyperbolic": hy,
        "product_4d": pr,
        "wolfram_rules_shown": [30, 90, 110, 184],
        "related": "https://github.com/jesus-morales-souhail",
        "conclusions_es": [
            "Euclidiano Z: f empuja a 4-2-1 (conjetura); g mezcla divergencia y ciclos (~no siempre ∞).",
            "Esférico Z/N: no existe ∞; f visita 1 más que g; g tiene más exploración en el anillo.",
            "Hiperbólico: g sube de nivel (capacidad exp) cuando el valor crece; f baja de nivel hacia el atractor.",
            "Producto 4D: modo f → componentes E/H acotadas hacia 1; modo g → E/H crecen; S siempre finita.",
            "CA elementales: Markov espacial local; misma lógica que Collatz (dependencia del presente).",
            "Grieta: paridad Collatz = cinta temporal; CA = cinta espacial; juntas muestran dependencia ≠ independencia.",
        ],
    }
    (RES / "discrete_geometry_ca.json").write_text(json.dumps(report, indent=2))
    (RES / "DISCRETE_GEOMETRY_CA.md").write_text(
        "\n".join(
            [
                "# Geometrías discretas + Collatz + CA",
                "",
                "## Principio",
                "No se define par/impar en una variedad continua sin discretizar.",
                "",
                "## Modelos",
                "- E: Z+",
                "- S: Z/NZ",
                "- H: (valor, nivel árbol)",
                "- 4D: (E, S, H_val, H_level)",
                "- CA: reglas 30, 90, 110, 184",
                "",
                "## Conclusiones",
                *[f"- {c}" for c in report["conclusions_es"]],
                "",
                f"Figuras: `{OUT}`",
            ]
        )
    )
    print("\nJSON:", RES / "discrete_geometry_ca.json")
    print("DONE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
