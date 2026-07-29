"""
Deterministic Markov chain for inverted (and normal) Collatz.

  X_{n+1} = f(X_n)  with probability 1

  P(X_{n+1} | X_n, ..., X_0) = P(X_{n+1} | X_n)

Successive parities are NOT i.i.d. coin flips.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from typing import Dict, List, Set

import numpy as np

from .collatz_maps import f_normal, g_inverted


def verify_unique_successor(seeds: List[int], steps: int = 40, inverted: bool = True) -> dict:
    fn = g_inverted if inverted else f_normal
    state_to_next: Dict[int, Set[int]] = defaultdict(set)
    for s0 in seeds:
        x = int(s0)
        for _ in range(steps):
            y = fn(x)
            state_to_next[x].add(y)
            x = y
    multi = {s: list(ns) for s, ns in state_to_next.items() if len(ns) > 1}
    return {
        "n_states": len(state_to_next),
        "deterministic_ok": len(multi) == 0,
        "violations": multi,
    }


def parity_transition_matrix(seeds: List[int], steps: int = 25, inverted: bool = True) -> dict:
    fn = g_inverted if inverted else f_normal
    pairs = []
    for s0 in seeds:
        x = int(s0)
        for _ in range(steps):
            y = fn(x)
            pairs.append((x % 2, y % 2))
            x = y
    cnt = Counter(pairs)
    cond = {}
    for p0 in (0, 1):
        sub = {k: v for k, v in cnt.items() if k[0] == p0}
        s = sum(sub.values()) or 1
        cond[str(p0)] = {str(p1): sub.get((p0, p1), 0) / s for p1 in (0, 1)}
    # under inverted, even → 3n+1 is always odd ⇒ P(odd|even)=1
    return {"n_transitions": len(pairs), "P_next_parity_given_now": cond}
