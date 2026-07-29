from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
#!/usr/bin/env python3
"""
Test matemático de las DOS ecuaciones Collatz
en espacio euclidiano plano y en geometría esférica.

Ecuación 1 — Original (atractor 4→2→1):
    par   → n/2
    impar → 3n+1

Ecuación 2 — Invertida (crecimiento / “→∞” a horizonte corto):
    par   → 3n+1
    impar → n/2

Uso:
  python experiments/collatz_geometry/two_maps_euclidean_spherical.py
  python experiments/collatz_geometry/two_maps_euclidean_spherical.py --n-max 50000
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

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
    return inverted_collatz_step(n) if inverted else collatz_original_step(n)


def orbit(n0: int, max_steps: int, inverted: bool) -> np.ndarray:
    x = int(n0)
    out = [x]
    for _ in range(max_steps):
        x = step(x, inverted)
        out.append(x)
        if not inverted and x in (1, 2, 4) and len(out) > 3:
            # atractor conocido; seguir un poco para el ciclo
            pass
    return np.asarray(out, dtype=np.int64)


def detect_cycle(seq: np.ndarray) -> Tuple[bool, int, int]:
    """Detecta ciclo simple; devuelve (found, start_idx, length)."""
    seen: Dict[int, int] = {}
    for i, v in enumerate(seq.tolist()):
        if v in seen:
            return True, seen[v], i - seen[v]
        seen[v] = i
    return False, -1, 0


def reaches_attractor_421(n0: int, max_steps: int = 5000) -> Tuple[bool, int, int]:
    """¿La Collatz original cae en {1,2,4}? Devuelve (ok, steps, max_val)."""
    x = int(n0)
    mx = x
    for t in range(max_steps):
        if x in (1, 2, 4):
            return True, t, mx
        x = collatz_original_step(x)
        if x > mx:
            mx = x
        if x < 0:  # no esperado
            return False, t, mx
    return x in (1, 2, 4), max_steps, mx


# ═══════════════════════════════════════════════════════════════════════════
# Estadística discreta (0 … n_max-1)
# ═══════════════════════════════════════════════════════════════════════════

def discrete_battery(n_max: int, horizons: List[int], max_steps_attr: int = 10000) -> dict:
    """Batería sobre enteros no negativos."""
    seeds = np.arange(n_max, dtype=np.int64)

    # --- Original: fracción que llega a 4-2-1 ---
    hit = 0
    max_excursions = []
    steps_to_hit = []
    for n in seeds:
        if n == 0:
            # 0 → 0 ciclo trivial en original (par → 0)
            hit += 1
            steps_to_hit.append(0)
            max_excursions.append(0)
            continue
        ok, st, mx = reaches_attractor_421(int(n), max_steps_attr)
        if ok:
            hit += 1
            steps_to_hit.append(st)
            max_excursions.append(mx)

    orig = {
        "n_max": n_max,
        "frac_reach_421": hit / n_max,
        "n_reach_421": hit,
        "mean_steps_to_421": float(np.mean(steps_to_hit)) if steps_to_hit else None,
        "median_steps_to_421": float(np.median(steps_to_hit)) if steps_to_hit else None,
        "mean_max_excursion": float(np.mean(max_excursions)) if max_excursions else None,
        "p99_max_excursion": float(np.percentile(max_excursions, 99)) if max_excursions else None,
    }

    # --- Invertida: crecimiento a horizontes fijos + ciclos ---
    inv_stats = {}
    for h in horizons:
        grow = 0
        log_ratio = []
        final_vals = []
        for n in seeds:
            x = int(n)
            x0 = max(x, 1)
            for _ in range(h):
                x = inverted_collatz_step(x)
            if x > n:
                grow += 1
            log_ratio.append(np.log1p(abs(x)) - np.log1p(x0))
            final_vals.append(x)
        inv_stats[f"h{h}"] = {
            "frac_x_h_gt_x0": grow / n_max,
            "mean_log1p_growth": float(np.mean(log_ratio)),
            "median_log1p_growth": float(np.median(log_ratio)),
            "frac_final_gt_1e6": float(np.mean(np.abs(final_vals) > 1e6)),
            "frac_final_gt_1e12": float(np.mean(np.abs(final_vals) > 1e12)),
        }

    # ciclos en invertida (órbitas cortas)
    cycle_count = 0
    for n in range(min(n_max, 5000)):
        seq = orbit(n, 200, inverted=True)
        found, _, _ = detect_cycle(seq)
        if found:
            cycle_count += 1
    inv_cycles = {
        "sample": min(n_max, 5000),
        "frac_hit_cycle_within_200_steps": cycle_count / min(n_max, 5000),
        "note": "Muchas semillas de la invertida caen en ciclos finitos a largo plazo; "
        "el ~50% es crecimiento a horizonte corto, no teorema de escape a ∞.",
    }

    return {"original": orig, "inverted_horizons": inv_stats, "inverted_cycles": inv_cycles}


# ═══════════════════════════════════════════════════════════════════════════
# Espacio euclidiano plano R^d
# ═══════════════════════════════════════════════════════════════════════════

def embed_orbit_euclidean(seq: np.ndarray, mode: str = "phase") -> np.ndarray:
    """
    Embebe la órbita en R^2 o R^3:
      phase: (x_t, x_{t+1})
      log_phase: (log1p|x_t|, log1p|x_{t+1}|)
      delay3: (x_t, x_{t+1}, x_{t+2}) normalizado por escala
    """
    s = seq.astype(np.float64)
    if mode == "phase":
        if len(s) < 2:
            return np.zeros((1, 2))
        return np.column_stack([s[:-1], s[1:]])
    if mode == "log_phase":
        if len(s) < 2:
            return np.zeros((1, 2))
        a = np.log1p(np.abs(s[:-1]))
        b = np.log1p(np.abs(s[1:]))
        return np.column_stack([a, b])
    if mode == "delay3":
        if len(s) < 3:
            return np.zeros((1, 3))
        return np.column_stack([s[:-2], s[1:-1], s[2:]])
    raise ValueError(mode)


def euclidean_orbit_metrics(seq: np.ndarray) -> dict:
    """Métricas euclídeas de la trayectoria embebida."""
    pts = embed_orbit_euclidean(seq, "log_phase")
    if len(pts) < 2:
        return {"path_length": 0.0, "mean_step": 0.0, "end_radius": 0.0, "expansion": 0.0}
    d = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    r0 = float(np.linalg.norm(pts[0]) + 1e-12)
    r1 = float(np.linalg.norm(pts[-1]) + 1e-12)
    # radio en el plano log-phase ~ tamaño de la órbita
    return {
        "path_length": float(d.sum()),
        "mean_step": float(d.mean()),
        "end_radius": r1,
        "start_radius": r0,
        "expansion": float(r1 / r0),
        "max_radius": float(np.max(np.linalg.norm(pts, axis=1))),
    }


def euclidean_battery(sample_seeds: List[int], steps: int = 40) -> dict:
    out = {"original": [], "inverted": []}
    for inv, key in [(False, "original"), (True, "inverted")]:
        expansions = []
        path_lens = []
        for n in sample_seeds:
            seq = orbit(n, steps, inverted=inv)
            m = euclidean_orbit_metrics(seq)
            expansions.append(m["expansion"])
            path_lens.append(m["path_length"])
        out[key] = {
            "n_seeds": len(sample_seeds),
            "steps": steps,
            "mean_expansion": float(np.mean(expansions)),
            "median_expansion": float(np.median(expansions)),
            "p90_expansion": float(np.percentile(expansions, 90)),
            "mean_path_length": float(np.mean(path_lens)),
            "frac_expansion_gt_1": float(np.mean(np.array(expansions) > 1.0)),
            "frac_expansion_gt_2": float(np.mean(np.array(expansions) > 2.0)),
        }
    return out


# ═══════════════════════════════════════════════════════════════════════════
# Geometría esférica S^{d-1}
# ═══════════════════════════════════════════════════════════════════════════

def embed_orbit_spherical(seq: np.ndarray, window: int = 8) -> np.ndarray:
    """
    Ventanas deslizantes de log1p|x| normalizadas a la esfera unitaria.
    Cada punto es un estado en S^{window-1}.
    """
    s = np.log1p(np.abs(seq.astype(np.float64)))
    if len(s) < window:
        v = np.zeros(window)
        v[: len(s)] = s
        return project_to_sphere(v).reshape(1, -1)
    rows = []
    for i in range(len(s) - window + 1):
        rows.append(project_to_sphere(s[i : i + window]))
    return np.asarray(rows)


def spherical_orbit_metrics(seq: np.ndarray, window: int = 8) -> dict:
    pts = embed_orbit_spherical(seq, window)
    if len(pts) < 2:
        return {
            "geodesic_path_length": 0.0,
            "mean_geodesic_step": 0.0,
            "end_to_start_geodesic": 0.0,
            "mean_radius_eucl_before_proj": 0.0,
        }
    geos = [
        spherical_geodesic_distance(pts[i], pts[i + 1]) for i in range(len(pts) - 1)
    ]
    # "concentración" hacia un polo: media del primer componente
    mean_pole = float(np.mean(pts[:, 0]))
    return {
        "geodesic_path_length": float(np.sum(geos)),
        "mean_geodesic_step": float(np.mean(geos)),
        "end_to_start_geodesic": spherical_geodesic_distance(pts[0], pts[-1]),
        "mean_first_coord": mean_pole,
        "n_sphere_points": len(pts),
    }


def spherical_battery(sample_seeds: List[int], steps: int = 40, window: int = 8) -> dict:
    out = {}
    for inv, key in [(False, "original"), (True, "inverted")]:
        g_paths = []
        g_end = []
        mean_steps = []
        for n in sample_seeds:
            seq = orbit(n, steps, inverted=inv)
            m = spherical_orbit_metrics(seq, window)
            g_paths.append(m["geodesic_path_length"])
            g_end.append(m["end_to_start_geodesic"])
            mean_steps.append(m["mean_geodesic_step"])
        out[key] = {
            "n_seeds": len(sample_seeds),
            "steps": steps,
            "sphere_dim": window,
            "mean_geodesic_path": float(np.mean(g_paths)),
            "mean_end_to_start_geodesic": float(np.mean(g_end)),
            "mean_geodesic_step": float(np.mean(mean_steps)),
            "median_end_to_start_geodesic": float(np.median(g_end)),
        }
    return out


def cloud_on_sphere(seeds: List[int], steps: int, inverted: bool, window: int = 6) -> np.ndarray:
    """Último estado esférico de cada semilla → nube en S^{w-1} (para PCA 2D)."""
    rows = []
    for n in seeds:
        seq = orbit(n, steps, inverted=inverted)
        pts = embed_orbit_spherical(seq, window)
        rows.append(pts[-1])
    return np.asarray(rows)


def pca_2d(x: np.ndarray) -> np.ndarray:
    x = x - x.mean(axis=0, keepdims=True)
    # SVD
    try:
        _, _, vt = np.linalg.svd(x, full_matrices=False)
        return x @ vt[:2].T
    except Exception:
        return x[:, :2]


# ═══════════════════════════════════════════════════════════════════════════
# Lyapunov-like / expansión local
# ═══════════════════════════════════════════════════════════════════════════

def local_expansion_stats(n_max: int = 5000, inverted: bool = True) -> dict:
    """
    Ratio |f(n)| / max(|n|,1) y log — proxy de expansión local del mapa.
    """
    ratios = []
    for n in range(n_max):
        y = step(n, inverted)
        ratios.append(abs(y) / max(abs(n), 1))
    r = np.asarray(ratios, dtype=float)
    # paridad
    even = r[0::2] if n_max > 1 else r
    odd = r[1::2] if n_max > 1 else r
    return {
        "inverted": inverted,
        "mean_ratio": float(r.mean()),
        "median_ratio": float(np.median(r)),
        "mean_log_ratio": float(np.mean(np.log(r + 1e-15))),
        "even_mean_ratio": float(even.mean()) if len(even) else None,
        "odd_mean_ratio": float(odd.mean()) if len(odd) else None,
        "theory_note": (
            "Original: par→1/2, impar→~3; Invertida: par→~3, impar→1/2. "
            "Media geométrica de factores de paridad aleatoria: "
            "orig ≈ sqrt(0.5*3)=sqrt(1.5)>1 en tramos impares pero el mapa cae al ciclo; "
            "inv ≈ misma media local pero con paridad intercambiada y dinámicas enteras distintas."
        ),
    }


# ═══════════════════════════════════════════════════════════════════════════
# Figuras
# ═══════════════════════════════════════════════════════════════════════════

def make_figures(out_dir: Path, sample: List[int]) -> List[str]:
    paths = []
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # 1) Órbitas ejemplo
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, inv, title in [
        (axes[0], False, "Original → atractor"),
        (axes[1], True, "Invertida → crecimiento/ciclos"),
    ]:
        for n in [6, 7, 27, 100, 999]:
            seq = orbit(n, 30, inverted=inv)
            ax.plot(np.log1p(np.abs(seq)), label=f"n={n}", lw=1.0)
        ax.set_xlabel("paso t")
        ax.set_ylabel("log1p |x_t|")
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7)
    fig.tight_layout()
    p = fig_dir / "orbits_log.png"
    fig.savefig(p, dpi=130)
    plt.close(fig)
    paths.append(str(p))

    # 2) Fase euclídea log
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, inv, title in [
        (axes[0], False, "Euclidiano: fase log (original)"),
        (axes[1], True, "Euclidiano: fase log (invertida)"),
    ]:
        for n in [7, 27, 100, 255, 1000]:
            seq = orbit(n, 50, inverted=inv)
            pts = embed_orbit_euclidean(seq, "log_phase")
            ax.plot(pts[:, 0], pts[:, 1], "-o", ms=2, lw=0.7, label=f"n={n}")
        ax.set_xlabel("log1p|x_t|")
        ax.set_ylabel("log1p|x_{t+1}|")
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7)
    fig.tight_layout()
    p = fig_dir / "euclidean_log_phase.png"
    fig.savefig(p, dpi=130)
    plt.close(fig)
    paths.append(str(p))

    # 3) Nubes esféricas PCA 2D
    seeds = sample[:400]
    cloud_o = cloud_on_sphere(seeds, steps=25, inverted=False, window=6)
    cloud_i = cloud_on_sphere(seeds, steps=25, inverted=True, window=6)
    po, pi = pca_2d(cloud_o), pca_2d(cloud_i)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].scatter(po[:, 0], po[:, 1], s=8, alpha=0.5, c="#1f77b4")
    axes[0].set_title("Esfera S^5 → PCA2 (original, t final)")
    axes[0].set_aspect("equal", adjustable="datalim")
    axes[0].grid(True, alpha=0.3)
    axes[1].scatter(pi[:, 0], pi[:, 1], s=8, alpha=0.5, c="#d62728")
    axes[1].set_title("Esfera S^5 → PCA2 (invertida, t final)")
    axes[1].set_aspect("equal", adjustable="datalim")
    axes[1].grid(True, alpha=0.3)
    fig.tight_layout()
    p = fig_dir / "spherical_pca_clouds.png"
    fig.savefig(p, dpi=130)
    plt.close(fig)
    paths.append(str(p))

    # 4) Histograma expansión euclídea
    exp_o = []
    exp_i = []
    for n in sample[:1000]:
        exp_o.append(euclidean_orbit_metrics(orbit(n, 40, False))["expansion"])
        exp_i.append(euclidean_orbit_metrics(orbit(n, 40, True))["expansion"])
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(exp_o, bins=40, alpha=0.55, label="original", color="#1f77b4")
    ax.hist(exp_i, bins=40, alpha=0.55, label="invertida", color="#d62728")
    ax.axvline(1.0, color="k", ls="--", lw=0.8)
    ax.set_xlabel("expansión radio log-fase (final/inicial)")
    ax.set_ylabel("cuentas")
    ax.set_title("Expansión euclídea de órbitas")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    p = fig_dir / "euclidean_expansion_hist.png"
    fig.savefig(p, dpi=130)
    plt.close(fig)
    paths.append(str(p))

    return paths


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-max", type=int, default=20_000)
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "results")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    hw = type('H', (), {'backend':'numpy','device_name':'cpu','cpu_threads':1})()
    print("=" * 70)
    print(" DOS ECUACIONES COLLATZ — Euclídeo + Esférico")
    print("=" * 70)
    print(f"hardware: {hw.backend} | {hw.device_name} | threads={hw.cpu_threads}")
    print(f"n_max = {args.n_max}")

    horizons = [5, 10, 20, 50, 100]
    print("\n[1/5] Batería discreta 0..n_max-1 …")
    disc = discrete_battery(args.n_max, horizons)

    sample = list(range(1, min(args.n_max, 2000)))
    print("[2/5] Métricas euclídeas …")
    eucl = euclidean_battery(sample, steps=40)

    print("[3/5] Métricas esféricas …")
    sph = spherical_battery(sample, steps=40, window=8)

    print("[4/5] Expansión local de los mapas …")
    loc_o = local_expansion_stats(min(args.n_max, 10000), inverted=False)
    loc_i = local_expansion_stats(min(args.n_max, 10000), inverted=True)

    print("[5/5] Figuras …")
    figs = make_figures(args.out.parent, sample)

    # Síntesis
    report = {
        "equations": {
            "original": {"even": "n/2", "odd": "3n+1", "behavior": "atractor 4→2→1"},
            "inverted": {"even": "3n+1", "odd": "n/2", "behavior": "crecimiento ~50% a horizonte corto; muchos ciclos a largo plazo"},
        },
        "hardware": {
            "backend": hw.backend,
            "device": hw.device_name,
            "threads": hw.cpu_threads,
        },
        "discrete": disc,
        "euclidean": eucl,
        "spherical": sph,
        "local_expansion": {"original": loc_o, "inverted": loc_i},
        "figures": figs,
        "conclusions_es": [],
    }

    # conclusiones automáticas
    c = []
    c.append(
        f"ORIGINAL: {100*disc['original']['frac_reach_421']:.2f}% de semillas en 0..{args.n_max-1} "
        f"alcanzan el atractor {{1,2,4}} (pasos med. {disc['original']['median_steps_to_421']})."
    )
    h10 = disc["inverted_horizons"].get("h10", {})
    c.append(
        f"INVERTIDA: tras 10 pasos, {100*h10.get('frac_x_h_gt_x0', 0):.2f}% cumple x_10 > x_0 "
        f"(métrica empírica ~50–52%)."
    )
    c.append(
        f"INVERTIDA ciclos: {100*disc['inverted_cycles']['frac_hit_cycle_within_200_steps']:.1f}% "
        f"de muestra corta cae en ciclo ≤200 pasos → no es escape universal a ∞."
    )
    c.append(
        f"EUCLÍDEO: expansión media original={eucl['original']['mean_expansion']:.3f}, "
        f"invertida={eucl['inverted']['mean_expansion']:.3f} "
        f"(fracción expand>1: orig {eucl['original']['frac_expansion_gt_1']:.2f}, "
        f"inv {eucl['inverted']['frac_expansion_gt_1']:.2f})."
    )
    c.append(
        f"ESFÉRICO: longitud geodésica media original={sph['original']['mean_geodesic_path']:.3f}, "
        f"invertida={sph['inverted']['mean_geodesic_path']:.3f}; "
        f"distancia fin-inicio media orig={sph['original']['mean_end_to_start_geodesic']:.3f}, "
        f"inv={sph['inverted']['mean_end_to_start_geodesic']:.3f}."
    )
    c.append(
        f"EXPANSIÓN LOCAL: original mean_ratio={loc_o['mean_ratio']:.3f} "
        f"(par {loc_o['even_mean_ratio']:.3f} / impar {loc_o['odd_mean_ratio']:.3f}); "
        f"invertida mean_ratio={loc_i['mean_ratio']:.3f} "
        f"(par {loc_i['even_mean_ratio']:.3f} / impar {loc_i['odd_mean_ratio']:.3f})."
    )
    c.append(
        "MATEMÁTICA APLICABLE: ambos mapas son funciones N→N (o Z) medibles; "
        "se embuten isométricamente en R^d (fase) y se proyectan a S^{d-1}. "
        "La diferencia no es 'euclídeo vs esférico' sino la dinámica del mapa: "
        "el original colapsa al ciclo; la invertida explora más la esfera y el plano log "
        "antes de ciclar o crecer a horizonte finito."
    )
    report["conclusions_es"] = c

    out_json = args.out / "two_maps_euclidean_spherical.json"
    out_json.write_text(json.dumps(report, indent=2))

    # markdown legible
    md = args.out / "REPORT.md"
    lines = [
        "# Dos ecuaciones Collatz — Euclídeo y Esférico",
        "",
        "## Ecuaciones",
        "",
        "| Mapa | Par | Impar | Comportamiento |",
        "|------|-----|-------|----------------|",
        "| Original | n/2 | 3n+1 | Atractor 4→2→1 |",
        "| Invertida | 3n+1 | n/2 | Crecimiento ~50% a h corto; ciclos frecuentes |",
        "",
        "## Resultados clave",
        "",
    ]
    for line in c:
        lines.append(f"- {line}")
    lines += [
        "",
        f"## Hardware",
        f"- {hw.backend} / {hw.device_name}",
        "",
        f"JSON completo: `{out_json.name}`",
        "",
        "## Figuras",
        "",
    ]
    for f in figs:
        lines.append(f"- `{Path(f).name}`")
    md.write_text("\n".join(lines))

    print("\n" + "=" * 70)
    print(" CONCLUSIONES")
    print("=" * 70)
    for line in c:
        print(" •", line)
    print(f"\nJSON: {out_json}")
    print(f"MD:   {md}")
    for f in figs:
        print(f"FIG:  {f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
