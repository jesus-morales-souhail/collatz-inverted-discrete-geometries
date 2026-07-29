"""
Discrete models that *imitate* geometries (no continuous parity without discretization).

E: Z+
S: Z/NZ (compact — no infinity)
H: (value, tree level) exponential capacity 2^level
4D product: (x_E, x_S, x_H_val, x_H_level)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List

from .collatz_maps import f_normal, g_inverted


def f_mod(n: int, N: int) -> int:
    return f_normal(n) % N


def g_mod(n: int, N: int) -> int:
    return g_inverted(n) % N


@dataclass
class HypState:
    value: int
    level: int

    def capacity(self) -> int:
        return 2 ** max(self.level, 0)


def hyperbolic_step(state: HypState, inverted: bool) -> HypState:
    fn = g_inverted if inverted else f_normal
    v0, v1 = state.value, fn(state.value)
    if v1 > v0:
        lvl = state.level + 1
    elif v1 < v0:
        lvl = max(0, state.level - 1)
    else:
        lvl = state.level
    return HypState(value=v1, level=lvl)


@dataclass
class State4D:
    x_E: int
    x_S: int
    x_H_val: int
    x_H_level: int
    mode: int  # 0=f, 1=g


def step_4d(X: State4D, N: int) -> State4D:
    inv = X.mode == 1
    fn = g_inverted if inv else f_normal
    fnm = g_mod if inv else f_mod
    e = fn(X.x_E)
    s = fnm(X.x_S, N)
    h = hyperbolic_step(HypState(X.x_H_val, X.x_H_level), inv)
    return State4D(e, s, h.value, h.level, X.mode)


def orbit_4d(n0: int, N: int, steps: int, mode: int) -> List[State4D]:
    X = State4D(n0, n0 % N, n0, 0, mode)
    seq = [X]
    for _ in range(steps):
        X = step_4d(X, N)
        seq.append(X)
    return seq
