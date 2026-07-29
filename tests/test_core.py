"""Minimal tests for publishable claims."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.collatz_maps import f_normal, g_inverted, orbit, reaches_421
from src.markov_chain import verify_unique_successor, parity_transition_matrix
from src.even_dot5_theorem import check_even_half_family, is_even_half
from src.geometry_models import step_4d, State4D, f_mod, g_mod


def test_maps_table():
    assert f_normal(10) == 5
    assert f_normal(7) == 22
    assert g_inverted(10) == 31
    assert g_inverted(7) == 3


def test_original_small_reach_421():
    for n in range(1, 200):
        ok, _ = reaches_421(n, 5000)
        assert ok


def test_markov_deterministic():
    r = verify_unique_successor(list(range(0, 100)), steps=30, inverted=True)
    assert r["deterministic_ok"]


def test_parity_even_goes_odd_under_g():
    m = parity_transition_matrix(list(range(0, 200)), steps=20, inverted=True)
    # P(odd | even) = 1
    assert abs(m["P_next_parity_given_now"]["0"]["1"] - 1.0) < 1e-12


def test_even_half_theorem():
    r = check_even_half_family(50, 10)
    assert r["theorem_one_step_preserves_even_half_and_increases"]
    assert is_even_half(2.5)
    assert is_even_half(0.5)


def test_spherical_no_escape():
    N = 32
    x = 7
    for _ in range(100):
        x = g_mod(x, N)
        assert 0 <= x < N


def test_4d_sphere_bounded():
    X = State4D(10, 10 % 64, 10, 0, mode=1)
    for _ in range(40):
        X = step_4d(X, 64)
        assert 0 <= X.x_S < 64


if __name__ == "__main__":
    test_maps_table()
    test_original_small_reach_421()
    test_markov_deterministic()
    test_parity_even_goes_odd_under_g()
    test_even_half_theorem()
    test_spherical_no_escape()
    test_4d_sphere_bounded()
    print("all tests ok")
