"""
Two Collatz maps on integers.

Normal (attractor 4-2-1, open conjecture on Z+):
    even → n/2,  odd → 3n+1

Inverted (more exploration; NOT always → ∞):
    even → 3n+1, odd → n/2
"""

from __future__ import annotations

from typing import List


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


def orbit(n0: int, steps: int, inverted: bool = False) -> List[int]:
    fn = g_inverted if inverted else f_normal
    x = max(int(n0), 0)
    seq = [x]
    for _ in range(steps):
        x = fn(x)
        seq.append(x)
    return seq


def reaches_421(n0: int, max_steps: int = 10000) -> tuple[bool, int]:
    """Whether normal Collatz hits {1,2,4}. Returns (hit, steps)."""
    x = int(n0)
    if x == 0:
        return True, 0
    for t in range(max_steps):
        if x in (1, 2, 4):
            return True, t
        x = f_normal(x)
    return x in (1, 2, 4), max_steps
