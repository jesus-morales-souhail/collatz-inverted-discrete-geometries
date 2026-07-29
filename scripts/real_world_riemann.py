from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
#!/usr/bin/env python3
"""
Las DOS ecuaciones Collatz en TODOS los planos/espacios modelo
relevantes — ancladas al MUNDO REAL (precios de mercado).

Espacios (geometrías modelo de curvatura constante + métricas riemannianas):
  1) Euclidiano plano R^d          K = 0
  2) Esfera S^{d-1}                K = +1
  3) Hiperbólico H^2 (Poincaré)    K = -1
  4) Riemann pullback:
       - métrica euclídea g = I
       - métrica log-precio (escala invariante)
       - métrica conformal 1/r^2
       - métrica tipo Fisher en scores

Semillas REALES: último dígito / enteros escalados de ticks y barras
(crypto, bolsa, oro, petróleo, tierras raras).

  python experiments/collatz_geometry/real_world_riemann.py
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


from src.collatz_maps import f_normal as collatz_original_step, g_inverted as inverted_collatz_step
from src.geometry import project_to_sphere, spherical_geodesic_distance


# ═══════════════════════════════════════════════════════════════════════════
# Mapas
# ═══════════════════════════════════════════════════════════════════════════

def step(n: int, inverted: bool) -> int:
    n = int(n)
    return inverted_collatz_step(n) if inverted else collatz_original_step(n)


def orbit_from(n0: int, n_steps: int, inverted: bool) -> np.ndarray:
    x = int(abs(n0))
    out = [x]
    for _ in range(n_steps):
        x = step(x, inverted)
        out.append(abs(x))
    return np.asarray(out, dtype=np.int64)


# ═══════════════════════════════════════════════════════════════════════════
# Carga mundo real
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class RealSeries:
    name: str
    asset_class: str
    prices: np.ndarray
    seeds_digit: np.ndarray
    seeds_scaled: np.ndarray


def load_binance_prices(path: Path, max_rows: int = 15000) -> np.ndarray:
    import json as _json

    prices = []
    with path.open() as f:
        for i, line in enumerate(f):
            if max_rows and i >= max_rows:
                break
            line = line.strip()
            if not line:
                continue
            j = line.find('"p":"')
            if j < 0:
                prices.append(float(_json.loads(line)["p"]))
            else:
                k = line.find('"', j + 5)
                prices.append(float(line[j + 5 : k]))
    return np.asarray(prices, dtype=np.float64)


def load_bar_closes(path: Path) -> np.ndarray:
    import pandas as pd

    df = pd.read_csv(path)
    cols = {c.lower(): c for c in df.columns}
    cname = None
    for k, v in cols.items():
        if "close" in k:
            cname = v
            break
    if cname is None:
        raise ValueError(f"no close in {path}")
    return df[cname].astype(float).to_numpy()


def series_from_prices(name: str, asset_class: str, prices: np.ndarray, pip: int = 2) -> RealSeries:
    prices = prices[np.isfinite(prices) & (prices > 0)]
    digits = np.array([last_digit(float(p), pip) for p in prices], dtype=np.int64)
    # parte entera escalada (restringida a rango manejable)
    scaled = np.rint(prices * (10**pip)).astype(np.int64)
    scaled = np.abs(scaled) % 10_000_000  # cap
    return RealSeries(name, asset_class, prices, digits, scaled)


def load_universe(max_ticks: int = 12000) -> List[RealSeries]:
    data_raw = ROOT / "data" / "raw"
    data_bars = ROOT / "data" / "bars"
    out: List[RealSeries] = []

    binance = [
        ("BTCUSDT", "crypto", 2),
        ("ETHUSDT", "crypto", 2),
        ("SOLUSDT", "crypto", 3),
        ("PAXGUSDT", "commodity", 2),
        ("DOGEUSDT", "crypto", 5),
    ]
    for sym, cls, pip in binance:
        p = data_raw / f"{sym}_aggTrades_latest.jsonl"
        if p.exists():
            prices = load_binance_prices(p, max_ticks)
            out.append(series_from_prices(sym, cls, prices, pip))
            print(f"  loaded {sym}: {len(prices)} ticks")

    bars = [
        ("SPY", "equity", "SPY_1m_7d.csv", 2),
        ("QQQ", "equity", "QQQ_1m_7d.csv", 2),
        ("AAPL", "equity", "AAPL_1m_7d.csv", 2),
        ("NVDA", "equity", "NVDA_1m_7d.csv", 2),
        ("GLD", "commodity", "GLD_1m_7d.csv", 2),
        ("USO", "commodity", "USO_1m_7d.csv", 2),
        ("REMX", "rare_earth", "REMX_1m_7d.csv", 2),
        ("MP", "rare_earth", "MP_1m_7d.csv", 2),
        ("UUUU", "rare_earth", "UUUU_1d_2y.csv", 2),
    ]
    for name, cls, fname, pip in bars:
        p = data_bars / fname
        if p.exists():
            prices = load_bar_closes(p)
            out.append(series_from_prices(name, cls, prices, pip))
            print(f"  loaded {name}: {len(prices)} bars")
    return out


# ═══════════════════════════════════════════════════════════════════════════
# Embeddings geométricos
# ═══════════════════════════════════════════════════════════════════════════

def emb_euclidean_phase(seq: np.ndarray) -> np.ndarray:
    """R^2: (log1p x_t, log1p x_{t+1}) — plano log."""
    s = np.log1p(seq.astype(np.float64))
    if len(s) < 2:
        return np.zeros((1, 2))
    return np.column_stack([s[:-1], s[1:]])


def emb_euclidean_delay3(seq: np.ndarray) -> np.ndarray:
    s = np.log1p(seq.astype(np.float64))
    if len(s) < 3:
        return np.zeros((1, 3))
    return np.column_stack([s[:-2], s[1:-1], s[2:]])


def emb_spherical(seq: np.ndarray, window: int = 8) -> np.ndarray:
    s = np.log1p(seq.astype(np.float64))
    if len(s) < window:
        v = np.zeros(window)
        v[: len(s)] = s
        return project_to_sphere(v).reshape(1, -1)
    rows = [project_to_sphere(s[i : i + window]) for i in range(len(s) - window + 1)]
    return np.asarray(rows)


def to_poincare_disk(xy: np.ndarray, eps: float = 1e-9) -> np.ndarray:
    """
    Mapea puntos de R^2 al disco de Poincaré (modelo K=-1).
    stereo-like: v → v / (1 + ||v||)  garantiza ||p|| < 1.
    """
    r = np.linalg.norm(xy, axis=1, keepdims=True)
    return xy / (1.0 + r + eps)


def hyperbolic_distance_poincare(a: np.ndarray, b: np.ndarray, eps: float = 1e-12) -> float:
    """Distancia hiperbólica en el disco de Poincaré."""
    a = np.asarray(a, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    # clamp inside disk
    def clip(p):
        n = np.linalg.norm(p)
        if n >= 1.0:
            p = p * (1.0 - 1e-6) / (n + eps)
        return p

    a, b = clip(a), clip(b)
    num = np.linalg.norm(a - b) ** 2
    den = (1.0 - np.dot(a, a)) * (1.0 - np.dot(b, b))
    arg = 1.0 + 2.0 * num / max(den, eps)
    return float(np.arccosh(min(max(arg, 1.0), 1e8)))


# ─── Métricas riemannianas a lo largo de una curva en R^d ─────────────────

def riemann_length(
    pts: np.ndarray,
    metric: str = "euclidean",
) -> dict:
    """
    Longitud de curva poligonal bajo distintas métricas riemannianas g:

      euclidean:   ds^2 = dx·dx                 (K=0 flat)
      log_scale:   ds^2 = sum (dx_i / (x_i+ε))^2  (invariante de escala)
      conformal:   ds^2 = dx·dx / (r^2+ε)       (conformal a euclídeo)
      fisher:      ds^2 = sum dx_i^2 / (p_i(1-p_i)+ε) tras softmax
    """
    pts = np.asarray(pts, dtype=np.float64)
    if len(pts) < 2:
        return {"length": 0.0, "mean_speed": 0.0, "energy": 0.0}

    dseg = []
    for i in range(len(pts) - 1):
        x, y = pts[i], pts[i + 1]
        dx = y - x
        if metric == "euclidean":
            ds = float(np.linalg.norm(dx))
        elif metric == "log_scale":
            # pullback de log coords
            mid = 0.5 * (np.abs(x) + np.abs(y)) + 1e-9
            ds = float(np.linalg.norm(dx / mid))
        elif metric == "conformal":
            r2 = float(np.dot(x, x) + 1e-9)
            ds = float(np.linalg.norm(dx) / np.sqrt(r2))
        elif metric == "fisher":
            # interpreta coords como logits → probs
            def soft(v):
                z = v - np.max(v)
                e = np.exp(z)
                return e / (e.sum() + 1e-12)

            p = soft(x)
            q = soft(y)
            # approx Fisher-Rao on simplex via Hellinger
            ds = float(np.linalg.norm(np.sqrt(np.clip(p, 1e-12, 1)) - np.sqrt(np.clip(q, 1e-12, 1))))
        else:
            raise ValueError(metric)
        dseg.append(ds)

    dseg = np.asarray(dseg)
    return {
        "length": float(dseg.sum()),
        "mean_speed": float(dseg.mean()),
        "energy": float(np.sum(dseg**2)),  # acción discreta
        "max_speed": float(dseg.max()),
    }


def path_metrics_all_geometries(seq: np.ndarray) -> dict:
    """Calcula longitudes/distancias en todos los espacios para una órbita."""
    seq = np.asarray(seq, dtype=np.int64)
    eu = emb_euclidean_phase(seq)
    e3 = emb_euclidean_delay3(seq)
    sp = emb_spherical(seq, window=min(8, max(3, len(seq) // 2)))
    # hiperbólico: Poincaré sobre fase euclídea
    po = to_poincare_disk(eu)

    # esférico: suma geodésicas
    sph_len = 0.0
    if len(sp) >= 2:
        sph_len = float(
            sum(spherical_geodesic_distance(sp[i], sp[i + 1]) for i in range(len(sp) - 1))
        )
        sph_end = spherical_geodesic_distance(sp[0], sp[-1])
    else:
        sph_end = 0.0

    hyp_len = 0.0
    if len(po) >= 2:
        hyp_len = float(
            sum(hyperbolic_distance_poincare(po[i], po[i + 1]) for i in range(len(po) - 1))
        )
        hyp_end = hyperbolic_distance_poincare(po[0], po[-1])
    else:
        hyp_end = 0.0

    # euclídeo R^2
    if len(eu) >= 2:
        eu_end = float(np.linalg.norm(eu[-1] - eu[0]))
        eu_len = float(np.sum(np.linalg.norm(np.diff(eu, axis=0), axis=1)))
    else:
        eu_end = eu_len = 0.0

    out = {
        "euclidean_R2": {
            "curvature_model": 0.0,
            "path_length": eu_len,
            "end_distance": eu_end,
            **{f"riemann_{k}": v for k, v in riemann_length(eu, "euclidean").items()},
        },
        "euclidean_R3_delay": {
            "curvature_model": 0.0,
            **riemann_length(e3, "euclidean"),
        },
        "sphere_S": {
            "curvature_model": 1.0,
            "path_length_geodesic": sph_len,
            "end_geodesic": sph_end,
        },
        "hyperbolic_H2": {
            "curvature_model": -1.0,
            "path_length_geodesic": hyp_len,
            "end_geodesic": hyp_end,
        },
        "riemann_log_scale": riemann_length(eu, "log_scale"),
        "riemann_conformal": riemann_length(eu, "conformal"),
        "riemann_fisher": riemann_length(eu, "fisher"),
        # crecimiento discreto
        "discrete": {
            "x0": int(seq[0]),
            "x_final": int(seq[-1]),
            "grew": bool(seq[-1] > seq[0]),
            "log_growth": float(np.log1p(seq[-1]) - np.log1p(max(seq[0], 1))),
            "max_val": int(seq.max()),
        },
    }
    return out


# ═══════════════════════════════════════════════════════════════════════════
# Batería sobre mundo real
# ═══════════════════════════════════════════════════════════════════════════

def aggregate_metrics(list_of_dicts: List[dict], path: Tuple[str, ...]) -> dict:
    vals = []
    for d in list_of_dicts:
        cur = d
        ok = True
        for k in path:
            if not isinstance(cur, dict) or k not in cur:
                ok = False
                break
            cur = cur[k]
        if ok and isinstance(cur, (int, float)) and np.isfinite(cur):
            vals.append(float(cur))
    if not vals:
        return {"n": 0}
    a = np.asarray(vals)
    return {
        "n": len(a),
        "mean": float(a.mean()),
        "median": float(np.median(a)),
        "std": float(a.std()),
        "p90": float(np.percentile(a, 90)),
    }


def run_on_series(
    series: RealSeries,
    n_steps: int = 12,
    max_seeds: int = 800,
    seed_mode: str = "digit",
) -> dict:
    """
    Para cada precio real → semilla entera → órbita original e invertida
    → métricas en todos los espacios.
    """
    seeds = series.seeds_digit if seed_mode == "digit" else series.seeds_scaled
    # subsample uniforme
    idx = np.linspace(0, len(seeds) - 1, num=min(max_seeds, len(seeds)), dtype=int)
    seeds = seeds[idx]

    orig_m, inv_m = [], []
    grew_o = grew_i = 0
    for s in seeds:
        o = orbit_from(int(s), n_steps, inverted=False)
        i = orbit_from(int(s), n_steps, inverted=True)
        mo = path_metrics_all_geometries(o)
        mi = path_metrics_all_geometries(i)
        orig_m.append(mo)
        inv_m.append(mi)
        if mo["discrete"]["grew"]:
            grew_o += 1
        if mi["discrete"]["grew"]:
            grew_i += 1

    def pack(ms: List[dict]) -> dict:
        return {
            "euclidean_path": aggregate_metrics(ms, ("euclidean_R2", "path_length")),
            "euclidean_end": aggregate_metrics(ms, ("euclidean_R2", "end_distance")),
            "sphere_path": aggregate_metrics(ms, ("sphere_S", "path_length_geodesic")),
            "sphere_end": aggregate_metrics(ms, ("sphere_S", "end_geodesic")),
            "hyperbolic_path": aggregate_metrics(ms, ("hyperbolic_H2", "path_length_geodesic")),
            "hyperbolic_end": aggregate_metrics(ms, ("hyperbolic_H2", "end_geodesic")),
            "riemann_log_length": aggregate_metrics(ms, ("riemann_log_scale", "length")),
            "riemann_conformal_length": aggregate_metrics(ms, ("riemann_conformal", "length")),
            "riemann_fisher_length": aggregate_metrics(ms, ("riemann_fisher", "length")),
            "log_growth": aggregate_metrics(ms, ("discrete", "log_growth")),
        }

    n = len(seeds)
    return {
        "name": series.name,
        "asset_class": series.asset_class,
        "n_prices": int(len(series.prices)),
        "n_seeds": n,
        "n_steps": n_steps,
        "seed_mode": seed_mode,
        "frac_grew_original": grew_o / max(n, 1),
        "frac_grew_inverted": grew_i / max(n, 1),
        "original": pack(orig_m),
        "inverted": pack(inv_m),
        # ratios inverted/original de longitudes medias (grieta geométrica)
        "ratio_inv_over_orig": {
            "euclidean_path": _safe_ratio(pack(inv_m)["euclidean_path"], pack(orig_m)["euclidean_path"]),
            "sphere_path": _safe_ratio(pack(inv_m)["sphere_path"], pack(orig_m)["sphere_path"]),
            "hyperbolic_path": _safe_ratio(pack(inv_m)["hyperbolic_path"], pack(orig_m)["hyperbolic_path"]),
            "riemann_log": _safe_ratio(pack(inv_m)["riemann_log_length"], pack(orig_m)["riemann_log_length"]),
            "riemann_fisher": _safe_ratio(pack(inv_m)["riemann_fisher_length"], pack(orig_m)["riemann_fisher_length"]),
        },
    }


def _safe_ratio(a: dict, b: dict) -> Optional[float]:
    if a.get("n", 0) == 0 or b.get("n", 0) == 0:
        return None
    bm = b.get("mean", 0.0)
    if abs(bm) < 1e-15:
        return None
    return float(a["mean"] / bm)


def class_summary(results: List[dict]) -> dict:
    by: Dict[str, List[dict]] = {}
    for r in results:
        by.setdefault(r["asset_class"], []).append(r)
    out = {}
    for cls, rows in by.items():
        def mean_key(path_orig_inv):
            vals = []
            for r in rows:
                # path_orig_inv like ("ratio_inv_over_orig", "sphere_path")
                cur = r
                for k in path_orig_inv:
                    cur = cur.get(k, {}) if isinstance(cur, dict) else None
                    if cur is None:
                        break
                if isinstance(cur, (int, float)) and np.isfinite(cur):
                    vals.append(float(cur))
            return float(np.mean(vals)) if vals else None

        out[cls] = {
            "n_assets": len(rows),
            "mean_frac_grew_inverted": float(np.mean([r["frac_grew_inverted"] for r in rows])),
            "mean_frac_grew_original": float(np.mean([r["frac_grew_original"] for r in rows])),
            "mean_ratio_sphere_inv_orig": mean_key(("ratio_inv_over_orig", "sphere_path")),
            "mean_ratio_euclid_inv_orig": mean_key(("ratio_inv_over_orig", "euclidean_path")),
            "mean_ratio_hyper_inv_orig": mean_key(("ratio_inv_over_orig", "hyperbolic_path")),
            "mean_ratio_riemann_log": mean_key(("ratio_inv_over_orig", "riemann_log")),
            "assets": [r["name"] for r in rows],
        }
    return out


# ═══════════════════════════════════════════════════════════════════════════
# Figuras
# ═══════════════════════════════════════════════════════════════════════════

def make_figures(results: List[dict], out_dir: Path) -> List[str]:
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    paths = []

    # 1) frac grew by asset
    names = [r["name"] for r in results]
    fo = [r["frac_grew_original"] for r in results]
    fi = [r["frac_grew_inverted"] for r in results]
    x = np.arange(len(names))
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.bar(x - 0.2, fo, 0.4, label="Original crece", color="#1f77b4")
    ax.bar(x + 0.2, fi, 0.4, label="Invertida crece", color="#d62728")
    ax.axhline(0.5, color="k", ls="--", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("fracción x_final > x_0")
    ax.set_title("Mundo real: crecimiento de órbitas Collatz (semillas = dígitos de precio)")
    ax.legend()
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    p = fig_dir / "real_growth_by_asset.png"
    fig.savefig(p, dpi=130)
    plt.close(fig)
    paths.append(str(p))

    # 2) ratios inv/orig por geometría
    geos = ["euclidean_path", "sphere_path", "hyperbolic_path", "riemann_log", "riemann_fisher"]
    fig, ax = plt.subplots(figsize=(10, 4))
    for g in geos:
        vals = [r["ratio_inv_over_orig"].get(g) or np.nan for r in results]
        ax.plot(names, vals, "-o", ms=4, label=g)
    ax.axhline(1.0, color="k", ls="--", lw=0.8)
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("longitud_inv / longitud_orig")
    ax.set_title("Ratio de exploración Invertida/Original por geometría (mundo real)")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    p = fig_dir / "real_geometry_ratios.png"
    fig.savefig(p, dpi=130)
    plt.close(fig)
    paths.append(str(p))

    # 3) ejemplo órbita real BTC dígito en 3 geometrías
    # se regenera una órbita demo
    if results:
        demo_seed = 7
        for inv, lab, col in [(False, "orig", "#1f77b4"), (True, "inv", "#d62728")]:
            pass
        fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))
        for seed in [2, 4, 5, 7, 9]:
            so = orbit_from(seed, 20, False)
            si = orbit_from(seed, 20, True)
            eo, ei = emb_euclidean_phase(so), emb_euclidean_phase(si)
            axes[0].plot(eo[:, 0], eo[:, 1], lw=0.8, alpha=0.7)
            axes[0].plot(ei[:, 0], ei[:, 1], lw=0.8, alpha=0.7, ls="--")
        axes[0].set_title("Euclidiano R² log-fase")
        axes[0].grid(True, alpha=0.3)
        # sphere PCA-ish: first 2 coords
        for seed in [2, 4, 5, 7, 9]:
            si = emb_spherical(orbit_from(seed, 25, True), 6)
            so = emb_spherical(orbit_from(seed, 25, False), 6)
            axes[1].plot(so[:, 0], so[:, 1], lw=0.8, alpha=0.7)
            axes[1].plot(si[:, 0], si[:, 1], lw=0.8, alpha=0.7, ls="--")
        axes[1].set_title("Esfera (coords 0,1)")
        axes[1].grid(True, alpha=0.3)
        for seed in [2, 4, 5, 7, 9]:
            eo = to_poincare_disk(emb_euclidean_phase(orbit_from(seed, 20, False)))
            ei = to_poincare_disk(emb_euclidean_phase(orbit_from(seed, 20, True)))
            axes[2].plot(eo[:, 0], eo[:, 1], lw=0.8, alpha=0.7)
            axes[2].plot(ei[:, 0], ei[:, 1], lw=0.8, alpha=0.7, ls="--")
        circ = plt.Circle((0, 0), 1, fill=False, color="k", lw=0.8)
        axes[2].add_patch(circ)
        axes[2].set_aspect("equal")
        axes[2].set_title("Hiperbólico (disco Poincaré)")
        axes[2].grid(True, alpha=0.3)
        fig.suptitle("Sólido=original, punteado=invertida (semillas dígito 2–9)")
        fig.tight_layout()
        p = fig_dir / "three_geometries_demo.png"
        fig.savefig(p, dpi=130)
        plt.close(fig)
        paths.append(str(p))

    return paths


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", type=int, default=12)
    ap.add_argument("--max-seeds", type=int, default=600)
    ap.add_argument("--max-ticks", type=int, default=12000)
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "results")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    hw = type('H', (), {'backend':'numpy','device_name':'cpu','cpu_threads':1})()
    print("=" * 72)
    print(" COLLATZ × GEOMETRÍAS × MUNDO REAL")
    print(" Euclídeo (K=0) · Esfera (K=+1) · Hiperbólico (K=-1) · Riemann pullback")
    print("=" * 72)
    print(f"hardware: {hw.device_name} | {hw.backend}")

    print("\n[1] Cargando series reales…")
    universe = load_universe(args.max_ticks)
    if not universe:
        print("No hay datos. Abort.")
        return 1

    print(f"\n[2] Órbitas + métricas ({len(universe)} activos, steps={args.steps})…")
    results = []
    for s in universe:
        print(f"  → {s.name} ({s.asset_class})…")
        results.append(run_on_series(s, n_steps=args.steps, max_seeds=args.max_seeds))

    print("\n[3] Resumen por clase…")
    by_cls = class_summary(results)

    print("\n[4] Figuras…")
    figs = make_figures(results, args.out.parent)

    # conclusiones
    conclusions = []
    conclusions.append(
        "GEOMETRÍAS USADAS: R^d (K=0), S^{n} (K=+1), H^2 Poincaré (K=-1), "
        "y longitudes riemannianas (euclídea, log-scale, conformal, Fisher) "
        "sobre el pullback de la órbita en el plano log-fase."
    )
    conclusions.append(
        "MUNDO REAL: semillas = último dígito (y opcionalmente entero escalado) "
        "de precios reales crypto/bolsa/commodities/tierras raras."
    )
    # media global
    fi = float(np.mean([r["frac_grew_inverted"] for r in results]))
    fo = float(np.mean([r["frac_grew_original"] for r in results]))
    conclusions.append(
        f"CRECIMIENTO GLOBAL: original crece en {100*fo:.1f}% de semillas-reales; "
        f"invertida en {100*fi:.1f}% (comparar con ~50% abstracto a h=10)."
    )
    # ratios
    ratios_s = [r["ratio_inv_over_orig"].get("sphere_path") for r in results]
    ratios_s = [x for x in ratios_s if x is not None]
    if ratios_s:
        conclusions.append(
            f"ESFERA: longitud geodésica inv/orig media = {np.mean(ratios_s):.3f} "
            f"(>1 ⇒ la invertida recorre más la cúpula sobre datos reales)."
        )
    ratios_e = [r["ratio_inv_over_orig"].get("euclidean_path") for r in results]
    ratios_e = [x for x in ratios_e if x is not None]
    if ratios_e:
        conclusions.append(
            f"EUCLÍDEO: path inv/orig media = {np.mean(ratios_e):.3f}."
        )
    ratios_h = [r["ratio_inv_over_orig"].get("hyperbolic_path") for r in results]
    ratios_h = [x for x in ratios_h if x is not None]
    if ratios_h:
        conclusions.append(
            f"HIPERBÓLICO: path inv/orig media = {np.mean(ratios_h):.3f}."
        )
    conclusions.append(
        "RIEMANN: las longitudes con g log-scale / conformal / Fisher miden la misma "
        "curva con distinta noción de distancia (escala-invariante, conformal, simplex). "
        "No hace falta curvatura variable a lo largo del camino para comparar mapas: "
        "los espacios modelo K∈{-1,0,+1} ya separan contracción vs exploración."
    )
    conclusions.append(
        "APLICACIÓN REAL: si en un activo la invertida muestra ratio>>1 en esfera/H^2 "
        "y frac_grew≈0.5, hay 'grieta geométrica' de exploración; eso NO implica edge "
        "de trading hasta superar costes (como vimos en el lab paper)."
    )

    report = {
        "spaces": {
            "euclidean_R2_R3": {"curvature_K": 0, "metric": "dx·dx"},
            "sphere_S": {"curvature_K": 1, "metric": "geodesic arccos <u,v>"},
            "hyperbolic_H2": {"curvature_K": -1, "metric": "Poincaré disk"},
            "riemann_pullback": {
                "metrics": ["euclidean", "log_scale", "conformal", "fisher"],
                "note": "longitud ∫ sqrt(g_ij dx^i dx^j) a lo largo de la órbita embebida",
            },
        },
        "equations": {
            "original": "even→n/2, odd→3n+1",
            "inverted": "even→3n+1, odd→n/2",
        },
        "hardware": {"device": hw.device_name, "backend": hw.backend},
        "by_asset": results,
        "by_class": by_cls,
        "figures": figs,
        "conclusions_es": conclusions,
    }

    out_json = args.out / "real_world_riemann.json"
    out_json.write_text(json.dumps(report, indent=2))

    md_lines = [
        "# Collatz en geometrías (incl. Riemann) × mundo real",
        "",
        "## Espacios",
        "",
        "| Espacio | Curvatura K | Qué mide |",
        "|---------|-------------|----------|",
        "| Euclidiano R²/R³ | 0 | longitud / fin-inicio en fase log |",
        "| Esfera S^{n} | +1 | geodésica esférica |",
        "| Hiperbólico H² | −1 | geodésica Poincaré |",
        "| Riemann pullback | (según g) | log-scale, conformal, Fisher |",
        "",
        "## Por clase de activo (semillas = dígitos de precio real)",
        "",
    ]
    for cls, info in sorted(by_cls.items()):
        md_lines.append(f"### {cls}")
        md_lines.append(
            f"- crece original={info['mean_frac_grew_original']:.3f}, "
            f"invertida={info['mean_frac_grew_inverted']:.3f}"
        )
        md_lines.append(
            f"- ratio path inv/orig: euclid={info['mean_ratio_euclid_inv_orig']}, "
            f"sphere={info['mean_ratio_sphere_inv_orig']}, "
            f"hyper={info['mean_ratio_hyper_inv_orig']}, "
            f"riemann_log={info['mean_ratio_riemann_log']}"
        )
        md_lines.append(f"- assets: {', '.join(info['assets'])}")
        md_lines.append("")
    md_lines.append("## Conclusiones")
    md_lines.append("")
    for c in conclusions:
        md_lines.append(f"- {c}")
    md_path = args.out / "REAL_WORLD_RIEMANN_REPORT.md"
    md_path.write_text("\n".join(md_lines))

    print("\n" + "=" * 72)
    print(" CONCLUSIÓN")
    print("=" * 72)
    for c in conclusions:
        print(" •", c)
    print("\nPor clase:")
    for cls, info in sorted(by_cls.items()):
        print(
            f"  [{cls}] grew_inv={info['mean_frac_grew_inverted']:.3f} "
            f"ratio_sph={info['mean_ratio_sphere_inv_orig']} "
            f"ratio_E={info['mean_ratio_euclid_inv_orig']} "
            f"ratio_H={info['mean_ratio_hyper_inv_orig']}"
        )
    print(f"\nJSON: {out_json}")
    print(f"MD:   {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
