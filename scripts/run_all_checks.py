#!/usr/bin/env python3
"""Reproduce core numerical checks and write results/summary.json."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.collatz_maps import reaches_421, g_inverted, orbit
from src.markov_chain import verify_unique_successor, parity_transition_matrix
from src.even_dot5_theorem import check_even_half_family
from src.geometry_models import orbit_4d, hyperbolic_step, HypState


def main() -> int:
    n_max = 5000
    hit = sum(1 for n in range(n_max) if reaches_421(n)[0])
    # inverted short-horizon growth
    grow10 = sum(
        1
        for n in range(n_max)
        if (lambda s: s[-1] > s[0])(orbit(n, 10, inverted=True))
    )
    markov = verify_unique_successor(list(range(0, 300)), 40, True)
    parity = parity_transition_matrix(list(range(0, 300)), 25, True)
    half = check_even_half_family(100, 12)

    # 4D product means
    N, steps = 64, 30
    logE_f, logE_g, lvl_f, lvl_g = [], [], [], []
    for s0 in range(1, 81):
        of = orbit_4d(s0, N, steps, mode=0)
        og = orbit_4d(s0, N, steps, mode=1)
        logE_f.append(__import__("math").log1p(of[-1].x_E))
        logE_g.append(__import__("math").log1p(og[-1].x_E))
        lvl_f.append(of[-1].x_H_level)
        lvl_g.append(og[-1].x_H_level)

    # spherical visit 1
    from src.geometry_models import f_mod, g_mod

    def frac_visit_1(inv: bool, N: int = 128, steps: int = 80) -> float:
        c = 0
        for s in range(N):
            x = s
            seen = False
            for _ in range(steps):
                x = g_mod(x, N) if inv else f_mod(x, N)
                if x == 1:
                    seen = True
                    break
            if seen or s == 1:
                c += 1
        return c / N

    summary = {
        "original_reach_421_frac_0_to_n": hit / n_max,
        "n_max": n_max,
        "inverted_frac_x10_gt_x0": grow10 / n_max,
        "markov_deterministic": markov,
        "parity_transitions_inverted": parity,
        "even_dot5": half,
        "spherical_frac_visit_1": {
            "normal_f": frac_visit_1(False),
            "inverted_g": frac_visit_1(True),
        },
        "product_4d": {
            "mean_log1p_E_normal": sum(logE_f) / len(logE_f),
            "mean_log1p_E_inverted": sum(logE_g) / len(logE_g),
            "mean_H_level_normal": sum(lvl_f) / len(lvl_f),
            "mean_H_level_inverted": sum(lvl_g) / len(lvl_g),
        },
        "main_observation": (
            "The inverted map explores more than the normal map across the discrete "
            "geometry models implemented here (Z, Z/NZ, hyperbolic tree, product). "
            "Empirical statement; not a proof of universal divergence under g."
        ),
    }
    out = ROOT / "results" / "summary_reproduced.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print("wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
