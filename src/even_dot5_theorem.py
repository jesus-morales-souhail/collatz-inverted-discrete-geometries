"""
Theorem (real family, inverted rule with floor parity):

  If x_0 = 2k + 1/2 for integer k ≥ 0, and
    f(x) = 3x+1 when floor(x) is even,
    f(x) = x/2   when floor(x) is odd,

  then floor(x_0) is even, so
    x_1 = 3x_0 + 1 = 2(3k+1) + 1/2,
  again of the form even + 1/2, and x_1 > x_0.
  By induction x_n → +∞.

This does NOT claim every inverted integer orbit diverges.
"""

from __future__ import annotations

import math


def inv_real_floor(x: float) -> float:
    fl = int(math.floor(abs(x) + 1e-15))
    return (3.0 * x + 1.0) if fl % 2 == 0 else x / 2.0


def is_even_half(x: float, tol: float = 1e-9) -> bool:
    """x = 2k + 0.5 for some integer k ≥ 0 (exact on dyadics used here)."""
    a = (x - 0.5) / 2.0
    return a >= -tol and abs(a - round(a)) < tol


def check_even_half_family(max_k: int = 200, steps: int = 30) -> dict:
    """Algebraic induction holds; float checks for small k."""
    ok = True
    samples = []
    for k in range(0, max_k + 1):
        x = 2 * k + 0.5
        # algebraic one step
        y = 3 * x + 1
        if not is_even_half(y) or not (y > x):
            ok = False
        if k <= 5 or k in (10, 50, 100):
            xs = [x]
            z = x
            for _ in range(min(steps, 15)):
                z = inv_real_floor(z)
                xs.append(z)
            samples.append({"k": k, "x0": x, "x_final": xs[-1], "grew": xs[-1] > x})
    return {
        "theorem_one_step_preserves_even_half_and_increases": ok,
        "samples": samples,
        "note": "Full ∞ follows by induction on reals; float may break form at huge magnitude.",
    }
