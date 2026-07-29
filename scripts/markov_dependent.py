from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
#!/usr/bin/env python3
"""
Cadena de Markov DEPENDIENTE para Collatz Invertida (y original).

  X_{n+1} = f(X_n)   con probabilidad 1  (transición degenerada)

Propiedad de Markov:
  P(X_{n+1} | X_n, X_{n-1}, …, X_0) = P(X_{n+1} | X_n)

Si fuera independiente del presente:
  P(X_{n+1} | X_n) = P(X_{n+1})
→ no habría regla f, ni trayectorias, ni “∞ vs acotado” desde un origen.

Uso:
  python experiments/collatz_geometry/markov_dependent.py
"""

from __future__ import annotations

import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Callable, Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "figures" / "markov"
RES = Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
RES.mkdir(parents=True, exist_ok=True)


# ═══════════════════════════════════════════════════════════════════════════
# 1. Definición formal f
# ═══════════════════════════════════════════════════════════════════════════

def last_digit_floor(x: float) -> int:
    return abs(int(math.floor(abs(x) + 1e-15))) % 10


def f_inverted(x: float) -> float:
    """
    f(x) = 3x+1  si último dígito de floor(|x|) es par
         = x/2   si es impar
    """
    d = last_digit_floor(x)
    if d % 2 == 0:
        return 3.0 * float(x) + 1.0
    return float(x) / 2.0


def f_inverted_int(n: int) -> int:
    """Versión entera clásica del lab: par→3n+1, impar→n//2."""
    n = int(n)
    return (3 * n + 1) if n % 2 == 0 else n // 2


def f_original_int(n: int) -> int:
    n = int(n)
    return n // 2 if n % 2 == 0 else 3 * n + 1


def transition_kernel(x: float, inverted: bool = True) -> Dict[float, float]:
    """
    Núcleo de Markov degenerado: un solo sucesor con masa 1.
    P(y | x) = 1_{y = f(x)}
    """
    y = f_inverted(x) if inverted else float(f_original_int(int(x)))
    return {y: 1.0}


# ═══════════════════════════════════════════════════════════════════════════
# 2. Propiedad de Markov — verificación operativa
# ═══════════════════════════════════════════════════════════════════════════

def orbit(x0: float, steps: int, f: Callable) -> List[float]:
    x = float(x0)
    seq = [x]
    for _ in range(steps):
        x = float(f(x))
        seq.append(x)
    return seq


def verify_markov_deterministic(
    seeds: List[int], steps: int = 30, inverted: bool = True
) -> dict:
    """
    Comprueba: dado el mismo X_n = s, el sucesor es siempre el mismo
    (independiente de cómo se llegó a s).
    Eso es Markov + determinismo.
    """
    f = f_inverted_int if inverted else f_original_int
    # mapa estado → conjunto de sucesores observados (con historial distinto)
    state_to_next: Dict[int, set] = defaultdict(set)
    state_to_predecessors: Dict[int, set] = defaultdict(set)

    for s0 in seeds:
        x = int(s0)
        for _ in range(steps):
            y = f(x)
            state_to_next[x].add(y)
            state_to_predecessors[y].add(x)
            x = y

    multi_next = {s: list(ns) for s, ns in state_to_next.items() if len(ns) > 1}
    # cada estado visitado tiene exactamente 1 sucesor
    n_states = len(state_to_next)
    n_unique_succ = sum(1 for ns in state_to_next.values() if len(ns) == 1)

    return {
        "inverted": inverted,
        "n_states_visited": n_states,
        "n_states_unique_successor": n_unique_succ,
        "markov_deterministic_ok": len(multi_next) == 0,
        "multi_successor_violations": multi_next,  # debe ser {}
        "note": "Si multi_successor vacío: P(X_{n+1}|X_n)=δ_{f(X_n)} y no depende del camino previo.",
    }


def verify_independence_fails(seeds: List[int], steps: int = 20) -> dict:
    """
    Si fuera independiente: P(X_{n+1}|X_n)=P(X_{n+1}).
    En un modelo independiente, la correlación corr(X_n, X_{n+1}) ≈ 0
    y la paridad bit0_n no predeciría bit0_{n+1} vía f.

    Aquí medimos:
    - correlación de paridades sucesivas bajo f (dependiente)
    - vs paridades i.i.d. Bernoulli(0.5) (independiente)
    """
    f = f_inverted_int
    par_pairs = []  # (bit0_n, bit0_{n+1}) under dynamics
    for s0 in seeds:
        x = int(s0)
        for _ in range(steps):
            y = f(x)
            par_pairs.append((x % 2, y % 2))
            x = y

    # P(paridad_{n+1} | paridad_n) empírica bajo f
    cnt = Counter(par_pairs)
    total = sum(cnt.values())
    # condicionales
    cond = {}
    for p0 in (0, 1):
        sub = {k: v for k, v in cnt.items() if k[0] == p0}
        s = sum(sub.values()) or 1
        cond[p0] = {p1: sub.get((p0, p1), 0) / s for p1 in (0, 1)}

    # bajo independencia de paridades 50-50, P(next|now)=0.5 siempre
    # distancia total variation media a 0.5
    tv = 0.0
    for p0 in (0, 1):
        for p1 in (0, 1):
            tv += abs(cond[p0][p1] - 0.5)
    tv /= 4.0

    # simulación independiente: paridades i.i.d.
    rng = np.random.default_rng(0)
    ind_pairs = list(zip(rng.integers(0, 2, total), rng.integers(0, 2, total)))
    cnt_i = Counter(ind_pairs)
    cond_i = {}
    for p0 in (0, 1):
        sub = {k: v for k, v in cnt_i.items() if k[0] == p0}
        s = sum(sub.values()) or 1
        cond_i[p0] = {p1: sub.get((p0, p1), 0) / s for p1 in (0, 1)}
    tv_i = 0.0
    for p0 in (0, 1):
        for p1 in (0, 1):
            tv_i += abs(cond_i[p0][p1] - 0.5)
    tv_i /= 4.0

    return {
        "n_transitions": total,
        "P_parity_next_given_now_under_f": cond,
        "mean_abs_dev_from_half_under_f": tv,
        "mean_abs_dev_from_half_iid": float(tv_i),
        "interpretation": (
            "Bajo f, la paridad siguiente NO es ~0.5 condicional al presente "
            "(TV>>0). Bajo i.i.d. sí ~0.5. Por tanto las paridades NO son monedas independientes."
        ),
    }


# ═══════════════════════════════════════════════════════════════════════════
# 3. Grafo de transiciones (estados finitos: dígitos / módulo)
# ═══════════════════════════════════════════════════════════════════════════

def digit_transition_graph(inverted: bool = True) -> dict:
    """
    Proyectamos a último dígito 0..9 (cadena en espacio finito).
    Ojo: f en enteros no es solo función del dígito (3n+1 cambia más),
    pero la DECISIÓN de rama sí lo es; el dígito siguiente depende del valor.
    Para la versión entera f_inverted_int, el siguiente dígito depende de n completo.
    Aquí mostramos: decisión de rama 100% función del dígito/paridad.
    """
    branch = {}
    for d in range(10):
        if inverted:
            branch[d] = "3n+1" if d % 2 == 0 else "n/2"
        else:
            branch[d] = "n/2" if d % 2 == 0 else "3n+1"
    return {
        "space": "last_digit_decides_branch",
        "inverted": inverted,
        "branch_by_digit": branch,
        "markov_on_full_state": "X_{n+1}=f(X_n) exacto en el valor completo",
    }


def integer_transition_sample(n_max: int = 50) -> dict:
    """Lista x → f(x) para enteros 0..n_max (invertida)."""
    edges = []
    for x in range(n_max + 1):
        y = f_inverted_int(x)
        edges.append({"from": x, "to": y, "prob": 1.0, "branch": "3n+1" if x % 2 == 0 else "n/2"})
    return {"edges": edges, "deterministic": True}


# ═══════════════════════════════════════════════════════════════════════════
# 4. Figuras
# ═══════════════════════════════════════════════════════════════════════════

def plot_markov_schema():
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # izquierda: dependiente
    ax = axes[0]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")
    ax.set_title("Markov DEPENDIENTE (Collatz Invertida)", fontsize=12, fontweight="bold")
    # nodos X0 X1 X2
    positions = [(1.5, 3), (5, 3), (8.5, 3)]
    labels = [r"$X_0$", r"$X_1=f(X_0)$", r"$X_2=f(X_1)$"]
    for (x, y), lab in zip(positions, labels):
        ax.add_patch(plt.Circle((x, y), 0.7, facecolor="#d0ebff", edgecolor="#1c7ed6", lw=2))
        ax.text(x, y, lab, ha="center", va="center", fontsize=11)
    ax.annotate("", xy=(4.2, 3), xytext=(2.3, 3),
                arrowprops=dict(arrowstyle="->", color="#1c7ed6", lw=2))
    ax.annotate("", xy=(7.7, 3), xytext=(5.8, 3),
                arrowprops=dict(arrowstyle="->", color="#1c7ed6", lw=2))
    ax.text(5, 1.2,
            r"$P(X_{n+1}\mid X_n,X_{n-1},\ldots)=P(X_{n+1}\mid X_n)=\mathbf{1}_{f(X_n)}$"
            "\nFuturo ← solo presente. Pasado remoto irrelevante dado $X_n$.",
            ha="center", fontsize=9,
            bbox=dict(boxstyle="round", facecolor="#e7f5ff", edgecolor="#1c7ed6"))

    # derecha: independiente (lo que NO es)
    ax = axes[1]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")
    ax.set_title("Si fuera INDEPENDIENTE (NO es el modelo)", fontsize=12, fontweight="bold", color="#c92a2a")
    for i, (x, y) in enumerate([(1.5, 3), (5, 3), (8.5, 3)]):
        ax.add_patch(plt.Circle((x, y), 0.7, facecolor="#ffe3e3", edgecolor="#c92a2a", lw=2))
        ax.text(x, y, rf"$Y_{i}$", ha="center", va="center", fontsize=12)
    # sin flechas entre ellos — ruido
    ax.text(5, 4.5, "tiradas sueltas", ha="center", color="#c92a2a", fontsize=10)
    ax.text(5, 1.2,
            r"$P(Y_{n+1}\mid Y_n)=P(Y_{n+1})$"
            "\nNo hay regla $f$, no hay trayectoria, no hay ∞ vs acotado."
            "\nEl modelo perdería su razón de ser.",
            ha="center", fontsize=9,
            bbox=dict(boxstyle="round", facecolor="#fff5f5", edgecolor="#c92a2a"))

    fig.tight_layout()
    p = OUT / "01_markov_dependiente_vs_independiente.png"
    fig.savefig(p, dpi=140, facecolor="white")
    plt.close(fig)
    print("wrote", p)


def plot_parity_dependence(cond: dict):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    # heatmap cond under f
    M = np.array([[cond[0][0], cond[0][1]], [cond[1][0], cond[1][1]]])
    im = axes[0].imshow(M, vmin=0, vmax=1, cmap="Blues")
    axes[0].set_xticks([0, 1])
    axes[0].set_yticks([0, 1])
    axes[0].set_xticklabels(["next PAR", "next IMPAR"])
    axes[0].set_yticklabels(["now PAR", "now IMPAR"])
    for i in range(2):
        for j in range(2):
            axes[0].text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", color="black")
    axes[0].set_title(r"P(paridad$_{n+1}$ | paridad$_n$) bajo $f$ invertida")
    fig.colorbar(im, ax=axes[0], fraction=0.046)

    # iid ~ 0.5
    M2 = np.full((2, 2), 0.5)
    im2 = axes[1].imshow(M2, vmin=0, vmax=1, cmap="Reds")
    axes[1].set_xticks([0, 1])
    axes[1].set_yticks([0, 1])
    axes[1].set_xticklabels(["next PAR", "next IMPAR"])
    axes[1].set_yticklabels(["now PAR", "now IMPAR"])
    for i in range(2):
        for j in range(2):
            axes[1].text(j, i, "0.50", ha="center", va="center")
    axes[1].set_title("Si monedas i.i.d. (independiente)")
    fig.colorbar(im2, ax=axes[1], fraction=0.046)

    fig.suptitle("Las paridades sucesivas NO son independientes — están ligadas por f", fontsize=11)
    fig.tight_layout()
    p = OUT / "02_parity_transition_not_iid.png"
    fig.savefig(p, dpi=140, facecolor="white")
    plt.close(fig)
    print("wrote", p)


def plot_small_chain_graph():
    """Grafo x→f(x) para 0..15."""
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(-1, 16)
    ax.set_ylim(-1, 12)
    ax.axis("off")
    ax.set_title(r"Cadena invertida entera: cada estado un único sucesor ($P=1$)", fontsize=12)

    # layout: estados en círculo o grid
    for x in range(16):
        y = f_inverted_int(x)
        # nodos
        ax.add_patch(plt.Circle((x, 8), 0.35, facecolor="#e7f5ff", edgecolor="#1c7ed6", lw=1.5))
        ax.text(x, 8, str(x), ha="center", va="center", fontsize=8, fontweight="bold")
        # flecha a y (colocar y abajo si y>15 mostrar número)
        ax.annotate(
            "",
            xy=(min(y, 15) if y <= 15 else 15, 3 if y <= 15 else 3),
            xytext=(x, 7.5),
            arrowprops=dict(arrowstyle="->", color="#495057", lw=1),
        )
        ax.text(x, 5.5, f"→{y}", ha="center", fontsize=7, color="#c92a2a" if x % 2 == 0 else "#2f9e44")

    ax.text(7.5, 1.2, "rojo etiqueta: 3n+1 (par)   verde: n/2 (impar)   ·   un solo destino por estado",
            ha="center", fontsize=9,
            bbox=dict(boxstyle="round", facecolor="#f8f9fa"))
    # ciclo 0-1 highlight
    ax.annotate("", xy=(1, 9.2), xytext=(0, 9.2),
                arrowprops=dict(arrowstyle="->", color="#e03131", lw=2))
    ax.annotate("", xy=(0, 9.5), xytext=(1, 9.5),
                arrowprops=dict(arrowstyle="->", color="#e03131", lw=2))
    ax.text(0.5, 10.3, "ciclo 0 ↔ 1", ha="center", color="#e03131", fontsize=10, fontweight="bold")

    fig.tight_layout()
    p = OUT / "03_transition_graph_0_15.png"
    fig.savefig(p, dpi=140, facecolor="white")
    plt.close(fig)
    print("wrote", p)


def plot_trajectory_markov():
    """Trayectoria: solo el estado actual basta para dibujar el siguiente."""
    seq = orbit(6, 12, f_inverted_int)
    fig, ax = plt.subplots(figsize=(10, 3.5))
    ax.plot(seq, "o-", color="#1c7ed6", ms=8, lw=2)
    for i, v in enumerate(seq):
        ax.text(i, v + 0.8, str(v), ha="center", fontsize=8)
    ax.set_xlabel("n (tiempo)")
    ax.set_ylabel(r"$X_n$")
    ax.set_title(r"Trayectoria $X_{n+1}=f(X_n)$ desde 6 — Markov dependiente del presente")
    ax.grid(True, alpha=0.3)
    ax.text(
        0.98, 0.95,
        "Dado $X_n=4$, el siguiente es\nsiempre 13, da igual el pasado.",
        transform=ax.transAxes, ha="right", va="top", fontsize=9,
        bbox=dict(boxstyle="round", facecolor="#fff3bf"),
    )
    fig.tight_layout()
    p = OUT / "04_trajectory_from_6.png"
    fig.savefig(p, dpi=140, facecolor="white")
    plt.close(fig)
    print("wrote", p)


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

def main() -> int:
    print("=" * 70)
    print(" CADENA DE MARKOV DEPENDIENTE — Collatz Invertida")
    print("=" * 70)

    seeds = list(range(0, 500))
    markov = verify_markov_deterministic(seeds, steps=40, inverted=True)
    markov_o = verify_markov_deterministic(seeds, steps=40, inverted=False)
    indep = verify_independence_fails(seeds, steps=25)
    graph = digit_transition_graph(True)
    edges = integer_transition_sample(30)

    print("\n[Markov determinista invertida]", markov["markov_deterministic_ok"],
          "estados", markov["n_states_visited"])
    print("[Markov determinista original ]", markov_o["markov_deterministic_ok"])
    print("[Paridades bajo f] P(next|now):", indep["P_parity_next_given_now_under_f"])
    print("[TV vs 0.5] f=", indep["mean_abs_dev_from_half_under_f"],
          "iid=", indep["mean_abs_dev_from_half_iid"])

    plot_markov_schema()
    plot_parity_dependence(indep["P_parity_next_given_now_under_f"])
    plot_small_chain_graph()
    plot_trajectory_markov()

    report = {
        "definition": {
            "state": "X_n",
            "kernel": "P(X_{n+1}=y | X_n=x) = 1 if y=f(x) else 0",
            "f_inverted": "even last digit / parity → 3x+1; odd → x/2",
            "markov_property": "P(X_{n+1}|X_n,...,X_0)=P(X_{n+1}|X_n)",
        },
        "why_not_independent": (
            "Independencia total ⇒ P(X_{n+1}|X_n)=P(X_{n+1}): no hay f, "
            "no hay trayectorias, el modelo pierde su razón de ser."
        ),
        "verification": {
            "inverted": markov,
            "original": markov_o,
            "parity_dependence": indep,
        },
        "branch_by_digit": graph,
        "sample_edges_0_30": edges,
        "conclusions_es": [
            "La Collatz Invertida es una cadena de Markov determinista: un solo sucesor con probabilidad 1.",
            "El futuro depende del presente; el pasado remoto no aporta una vez conocido X_n.",
            "Las paridades sucesivas NO son monedas i.i.d.; el ~50% global es promedio sobre semillas, no independencia paso a paso.",
            "Si los eventos fueran totalmente independientes, no existiría relación dinámica con el pasado ni con el presente vía f.",
        ],
    }

    out = RES / "markov_dependent.json"
    out.write_text(json.dumps(report, indent=2))
    md = RES / "MARKOV_DEPENDENT.md"
    md.write_text(
        "\n".join(
            [
                "# Cadena de Markov dependiente — Collatz Invertida",
                "",
                "## Definición",
                "",
                r"$$X_{n+1}=f(X_n),\quad P=1$$",
                "",
                r"$$f(x)=\begin{cases}3x+1 & \text{paridad/dígito par}\\ x/2 & \text{impar}\end{cases}$$",
                "",
                "## Propiedad de Markov",
                "",
                r"$$P(X_{n+1}\mid X_n,\ldots,X_0)=P(X_{n+1}\mid X_n)$$",
                "",
                "## Verificación",
                "",
                f"- Determinista OK (invertida): **{markov['markov_deterministic_ok']}**",
                f"- Estados visitados: {markov['n_states_visited']}",
                f"- Desviación de paridades vs 0.5 bajo f: **{indep['mean_abs_dev_from_half_under_f']:.3f}** (i.i.d. ~{indep['mean_abs_dev_from_half_iid']:.3f})",
                "",
                "## Conclusión",
                "",
                "Dependiente del **presente**. No independiente. El ~50% no es moneda paso a paso.",
                "",
                f"Figuras en `{OUT}`",
            ]
        )
    )
    print(f"\nJSON: {out}")
    print(f"MD:   {md}")
    print("Figuras:", list(OUT.glob("*.png")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
