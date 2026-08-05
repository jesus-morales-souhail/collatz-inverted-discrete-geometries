"""
Coverage residual on Collatz graphs — not a proof of the conjecture.

Product: map or inverse monoid + space + budget D
  -> measure of covered states and residual R = 1 - covered.

Generators of the inverse tree of f:
  alpha(n) = 2n
  beta(n)  = (n-1)/3  when that is a positive odd integer
"""

from __future__ import annotations

from collections import deque
from typing import Dict, Iterable, List, Set, Tuple

import numpy as np

from .collatz_maps import f_normal, g_inverted


def alpha(n: int) -> int:
    return 2 * int(n)


def beta_ok(n: int) -> bool:
    """Standard inverse branch: (n-1)/3 positive odd integer."""
    n = int(n)
    if (n - 1) % 3 != 0:
        return False
    x = (n - 1) // 3
    return x >= 1 and x % 2 == 1


def beta(n: int) -> int | None:
    if not beta_ok(n):
        return None
    return (int(n) - 1) // 3


def inverse_children(n: int) -> List[int]:
    """Preimages under the usual Collatz graph (forward f)."""
    out = [alpha(n)]
    b = beta(n)
    if b is not None:
        out.append(b)
    return out


def inverse_tree_levels(root: int = 1, max_depth: int = 20) -> Dict[int, Set[int]]:
    """
    BFS layers of the inverse tree from root.
    levels[d] = nodes first reached at depth d (shortest inverse word length).
    """
    levels: Dict[int, Set[int]] = {0: {int(root)}}
    seen: Set[int] = {int(root)}
    frontier = {int(root)}
    for d in range(max_depth):
        nxt: Set[int] = set()
        for n in frontier:
            for c in inverse_children(n):
                if c not in seen:
                    # cap explosion for very large ints in pure BFS without bound
                    if c > 10**18:
                        continue
                    seen.add(c)
                    nxt.add(c)
        levels[d + 1] = nxt
        frontier = nxt
        if not frontier:
            break
    return levels


def inverse_coverage_up_to(
    max_depth: int,
    X_max: int,
    root: int = 1,
) -> dict:
    """
    Grow inverse tree to depth max_depth; measure coverage of {1..X_max}.
    """
    levels = inverse_tree_levels(root=root, max_depth=max_depth)
    covered: Set[int] = set()
    N_d = []
    cum = []
    for d in range(max_depth + 1):
        layer = levels.get(d, set())
        N_d.append(len(layer))
        for n in layer:
            if 1 <= n <= X_max:
                covered.add(n)
        cum.append(len(covered))

    total = int(X_max)
    R = 1.0 - len(covered) / total
    # T_eff between consecutive layers
    T_eff = []
    for d in range(len(N_d) - 1):
        a, b = N_d[d], N_d[d + 1]
        if a > 0 and b > 0:
            T_eff.append(float(np.log(b / a)))
        else:
            T_eff.append(float("nan"))

    return {
        "root": root,
        "max_depth": max_depth,
        "X_max": X_max,
        "N_d": N_d,
        "covered_cumulative_in_1_X": cum,
        "covered_final": len(covered),
        "coverage_frac": len(covered) / total,
        "residual_R": R,
        "T_eff_log_branch": T_eff,
        "mean_T_eff": float(np.nanmean(T_eff)) if T_eff else None,
    }


def forward_coverage(
    seeds: Iterable[int],
    steps: int,
    inverted: bool,
    X_max: int | None = None,
    mod_N: int | None = None,
) -> dict:
    """
    From a set of seeds, iterate f or g for `steps` and measure coverage.
    - on Z: unique values (optionally restricted to 1..X_max for residual)
    - on Z/NZ: fraction of residue classes visited
    """
    fn = g_inverted if inverted else f_normal
    seeds = list(seeds)
    visited: Set[int] = set()
    visited_mod: Set[int] = set() if mod_N else set()

    for s0 in seeds:
        x = max(int(s0), 0)
        for _ in range(steps + 1):
            visited.add(x)
            if mod_N:
                visited_mod.add(x % mod_N)
            x = fn(x)

    out: dict = {
        "map": "g" if inverted else "f",
        "n_seeds": len(seeds),
        "steps": steps,
        "n_unique": len(visited),
    }

    if X_max is not None:
        hit = {n for n in visited if 1 <= n <= X_max}
        out["X_max"] = X_max
        out["covered_in_1_X"] = len(hit)
        out["coverage_frac"] = len(hit) / X_max
        out["residual_R"] = 1.0 - len(hit) / X_max

    if mod_N is not None:
        out["mod_N"] = mod_N
        out["classes_hit"] = len(visited_mod)
        out["coverage_frac_mod"] = len(visited_mod) / mod_N
        out["residual_R_mod"] = 1.0 - len(visited_mod) / mod_N

    return out


def choice_entropy_on_cycle() -> dict:
    """
    On the attractor cycle of f: 4 -> 2 -> 1 -> 4 ...
    each state has exactly one successor: Shannon choice entropy = 0.
    """
    cycle = [4, 2, 1]
    # under f
    succ = {4: 2, 2: 1, 1: 4}
    return {
        "cycle": cycle,
        "successors": succ,
        "branching_factor": 1.0,
        "choice_entropy_bits": 0.0,
        "note": "one legal next state per node on the cycle under f",
    }


def compare_coverage_bundle(
    X_max: int = 5000,
    inv_depth: int = 24,
    forward_steps: int = 80,
    n_seeds: int = 200,
    mod_N: int = 128,
) -> dict:
    seeds = list(range(1, n_seeds + 1))

    inv = inverse_coverage_up_to(max_depth=inv_depth, X_max=X_max, root=1)
    # also report coverage of 1..X only among nodes ever generated (may miss large branches)
    fwd_f = forward_coverage(seeds, forward_steps, inverted=False, X_max=X_max, mod_N=mod_N)
    fwd_g = forward_coverage(seeds, forward_steps, inverted=True, X_max=X_max, mod_N=mod_N)

    # pure modular coverage from all residues as seeds (full base)
    all_res = list(range(mod_N))
    mod_f = forward_coverage(all_res, forward_steps, inverted=False, mod_N=mod_N)
    mod_g = forward_coverage(all_res, forward_steps, inverted=True, mod_N=mod_N)

    cycle = choice_entropy_on_cycle()

    # quality flags for the coverage product (not Collatz proof)
    flags = {
        "cycle_choice_entropy_zero": cycle["choice_entropy_bits"] == 0.0,
        "inverse_tree_grows": inv["N_d"][-1] > inv["N_d"][0] if inv["N_d"] else False,
        "inverse_residual_below_half": inv["residual_R"] < 0.5,
        "g_covers_more_mod_than_f": mod_g["coverage_frac_mod"] > mod_f["coverage_frac_mod"],
        "g_covers_more_Z_than_f": fwd_g["coverage_frac"] > fwd_f["coverage_frac"],
        "mean_T_eff_positive": (inv["mean_T_eff"] or 0) > 0,
    }
    score = sum(1 for v in flags.values() if v)

    return {
        "product": "coverage residual (not a Collatz proof)",
        "choice_entropy_nucleus": cycle,
        "inverse_from_1": inv,
        "forward_Z": {"f": fwd_f, "g": fwd_g},
        "forward_mod_all_residues": {"f": mod_f, "g": mod_g},
        "quality_flags": flags,
        "quality_score": f"{score}/{len(flags)}",
        "reading": {
            "nucleus": "cycle 4-2-1 has choice entropy 0 under f",
            "expansion": "inverse tree N_d and T_eff measure microstate growth from 1",
            "residual": "R = 1 - coverage; lower R = more coverage under budget",
            "not": "does not prove every integer reaches 1",
        },
    }
