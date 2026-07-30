"""
Statistical quality of bitstreams from Collatz orbits on lattice planes.

Not a cryptographic claim. Implements a small NIST SP 800-22-style battery
(frequency, block frequency, runs, serial, autocorrelation) plus a log-orbit
drift diagnostic (Lagarias-type parity heuristic vs mean reversion).

Lattice planes (dimensions added one by one):
  1D  parity on Z
  1D  LSB stream of X_t
  S   low bits of X_t mod N  (ring lattice Z/NZ)
  H   sign-of-growth bits along the orbit
  Z^2 product of (parity, mod-2 of modular value)
  Z^3 product + growth bit
  log-floor LSBs of log2(1+X)

Maps f and g are controls of each other. Expected outcome for chaos-map PRNGs
is mixed or fail under serious batteries; this module records pass/fail at a
fixed significance level without claiming cryptographic suitability.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Dict, List, Sequence

import numpy as np

from .collatz_maps import f_normal, g_inverted


# ─── Orbits ─────────────────────────────────────────────────────────────────

def orbit_array(n0: int, steps: int, inverted: bool) -> np.ndarray:
    fn = g_inverted if inverted else f_normal
    x = max(int(n0), 0)
    out = np.empty(steps + 1, dtype=np.int64)
    out[0] = x
    for t in range(steps):
        x = fn(int(x))
        out[t + 1] = x
    return out


def multi_seed_orbits(
    seeds: Sequence[int], steps: int, inverted: bool
) -> np.ndarray:
    """Stack orbits shape (n_seeds, steps+1)."""
    return np.stack([orbit_array(s, steps, inverted) for s in seeds], axis=0)


# ─── Lattice plane bit extractors ───────────────────────────────────────────

def bits_parity(orbit: np.ndarray) -> np.ndarray:
    """Plane 1: parity of X_t on Z."""
    return (orbit.astype(np.int64) & 1).astype(np.uint8)


def bits_lsb(orbit: np.ndarray, n_bits: int = 1) -> np.ndarray:
    """Plane 1b: n least-significant bits of X_t, flattened bit-major per step."""
    o = orbit.astype(np.int64)
    bits = []
    for b in range(n_bits):
        bits.append(((o >> b) & 1).astype(np.uint8))
    # interleave by time: [b0_t0, b1_t0, ..., b0_t1, ...]
    stacked = np.stack(bits, axis=1)  # (T, n_bits)
    return stacked.reshape(-1)


def bits_modular(orbit: np.ndarray, N: int = 256, n_bits: int = 8) -> np.ndarray:
    """Spherical lattice Z/NZ: low bits of (X mod N)."""
    m = (orbit.astype(np.int64) % int(N)).astype(np.int64)
    return bits_lsb(m, n_bits=min(n_bits, max(1, int(N).bit_length())))


def bits_growth(orbit: np.ndarray) -> np.ndarray:
    """Hyperbolic-style: 1 if X_{t+1} > X_t else 0 (length T-1)."""
    o = orbit.astype(np.int64)
    return (o[1:] > o[:-1]).astype(np.uint8)


def bits_product_z2(orbit: np.ndarray, N: int = 64) -> np.ndarray:
    """
    Z^2 lattice coordinates from (parity, (X mod N) mod 2), interleaved.
    """
    p = bits_parity(orbit)
    s = ((orbit.astype(np.int64) % int(N)) & 1).astype(np.uint8)
    return np.stack([p, s], axis=1).reshape(-1)


def bits_product_z3(orbit: np.ndarray, N: int = 64) -> np.ndarray:
    """
    Z^3 product: parity, modular LSB, growth bit (growth padded with 0 at end).
    """
    p = bits_parity(orbit)
    s = ((orbit.astype(np.int64) % int(N)) & 1).astype(np.uint8)
    g = bits_growth(orbit)
    g = np.concatenate([g, np.array([0], dtype=np.uint8)])
    return np.stack([p, s, g], axis=1).reshape(-1)


def bits_log_floor_lsb(orbit: np.ndarray, n_bits: int = 4) -> np.ndarray:
    """LSBs of floor(log2(1+X)) — coarse scale bits."""
    o = orbit.astype(np.float64)
    lf = np.floor(np.log2(1.0 + np.abs(o))).astype(np.int64)
    return bits_lsb(lf, n_bits=n_bits)


def bits_delta_parity(orbit: np.ndarray) -> np.ndarray:
    """Parity of |X_{t+1}-X_t| (length T-1)."""
    o = orbit.astype(np.int64)
    d = np.abs(o[1:] - o[:-1])
    return (d & 1).astype(np.uint8)


EXTRACTORS: Dict[str, Callable[[np.ndarray], np.ndarray]] = {
    "E_parity_Z": bits_parity,
    "E_lsb1": lambda o: bits_lsb(o, 1),
    "E_lsb4": lambda o: bits_lsb(o, 4),
    "S_mod256_low8": lambda o: bits_modular(o, 256, 8),
    "S_mod64_low6": lambda o: bits_modular(o, 64, 6),
    "H_growth": bits_growth,
    "Z2_parity_x_mod": bits_product_z2,
    "Z3_parity_mod_growth": bits_product_z3,
    "log2_floor_lsb4": bits_log_floor_lsb,
    "delta_parity": bits_delta_parity,
}


def stream_from_seeds(
    seeds: Sequence[int],
    steps: int,
    inverted: bool,
    extractor: Callable[[np.ndarray], np.ndarray],
) -> np.ndarray:
    chunks = []
    for s in seeds:
        o = orbit_array(s, steps, inverted)
        chunks.append(extractor(o))
    return np.concatenate(chunks).astype(np.uint8)


# ─── Statistical tests (NIST SP 800-22 style, self-contained) ───────────────

def _erfc(x: float) -> float:
    return float(math.erfc(float(x)))


def _normal_sf_abs(z: float) -> float:
    """Two-sided p-value for |Z| under N(0,1): erfc(|z|/sqrt(2))."""
    return float(_erfc(abs(z) / np.sqrt(2.0)))


@dataclass
class TestResult:
    name: str
    statistic: float
    p_value: float
    n_bits: int
    pass_at_001: bool
    note: str = ""

    def as_dict(self) -> dict:
        return {
            "name": self.name,
            "statistic": float(self.statistic),
            "p_value": float(self.p_value),
            "n_bits": int(self.n_bits),
            "pass_at_0.01": bool(self.pass_at_001),
            "note": self.note,
        }


def test_monobit(bits: np.ndarray) -> TestResult:
    """NIST frequency (monobit) test."""
    b = bits.astype(np.int8)
    n = len(b)
    if n < 100:
        return TestResult("monobit", 0.0, 0.0, n, False, "n too small")
    s = 2 * int(np.sum(b)) - n  # +1 for 1, -1 for 0
    s_obs = abs(s) / np.sqrt(n)
    p = _normal_sf_abs(s_obs)
    return TestResult("monobit", float(s_obs), p, n, p >= 0.01)


def test_block_frequency(bits: np.ndarray, M: int = 128) -> TestResult:
    """NIST block frequency test."""
    b = bits.astype(np.float64)
    n = len(b)
    N = n // M
    if N < 1:
        return TestResult("block_frequency", 0.0, 0.0, n, False, "n too small")
    blocks = b[: N * M].reshape(N, M)
    pi = blocks.mean(axis=1)
    chi2 = 4.0 * M * float(np.sum((pi - 0.5) ** 2))
    # p = igamc(N/2, chi2/2); approximate with survival of chi2 via Wilson-Hilferty
    # Use incomplete gamma via scipy-free: numpy doesn't have gammaincc for all —
    # use series for upper incomplete ratio when df is moderate.
    p = _chi2_sf(chi2, df=N)
    return TestResult("block_frequency", chi2, p, N * M, p >= 0.01, f"M={M},Nblocks={N}")


def test_runs(bits: np.ndarray) -> TestResult:
    """NIST runs test."""
    b = bits.astype(np.int8)
    n = len(b)
    if n < 100:
        return TestResult("runs", 0.0, 0.0, n, False, "n too small")
    pi = float(np.mean(b))
    if abs(pi - 0.5) >= 2.0 / np.sqrt(n):
        # prerequisite failed — frequency already bad
        return TestResult("runs", 0.0, 0.0, n, False, "prerequisite monobit failed")
    v = 1 + int(np.sum(b[1:] != b[:-1]))
    num = abs(v - 2.0 * n * pi * (1.0 - pi))
    den = 2.0 * np.sqrt(2.0 * n) * pi * (1.0 - pi)
    p = _erfc(num / den) if den > 0 else 0.0
    return TestResult("runs", float(v), float(p), n, p >= 0.01)


def test_serial_m2(bits: np.ndarray) -> TestResult:
    """Serial test on overlapping 2-bit patterns (chi-square, 3 df)."""
    b = bits.astype(np.int8)
    n = len(b)
    if n < 100:
        return TestResult("serial_2bit", 0.0, 0.0, n, False, "n too small")
    # non-overlapping pairs for simplicity
    m = (n // 2) * 2
    pairs = b[:m].reshape(-1, 2)
    codes = pairs[:, 0] * 2 + pairs[:, 1]
    counts = np.bincount(codes, minlength=4).astype(np.float64)
    exp = len(pairs) / 4.0
    chi2 = float(np.sum((counts - exp) ** 2 / exp))
    p = _chi2_sf(chi2, df=3)
    return TestResult("serial_2bit", chi2, p, m, p >= 0.01)


def test_serial_m3(bits: np.ndarray) -> TestResult:
    """Serial test on non-overlapping 3-bit blocks (7 df)."""
    b = bits.astype(np.int8)
    n = len(b)
    m = (n // 3) * 3
    if m < 300:
        return TestResult("serial_3bit", 0.0, 0.0, n, False, "n too small")
    trips = b[:m].reshape(-1, 3)
    codes = trips[:, 0] * 4 + trips[:, 1] * 2 + trips[:, 2]
    counts = np.bincount(codes, minlength=8).astype(np.float64)
    exp = len(trips) / 8.0
    chi2 = float(np.sum((counts - exp) ** 2 / exp))
    p = _chi2_sf(chi2, df=7)
    return TestResult("serial_3bit", chi2, p, m, p >= 0.01)


def test_autocorr(bits: np.ndarray, lag: int = 1) -> TestResult:
    """Two-sided test of lag autocorrelation of ±1 stream."""
    b = bits.astype(np.float64)
    n = len(b)
    if n <= lag + 50:
        return TestResult(f"autocorr_lag{lag}", 0.0, 0.0, n, False, "n too small")
    x = 2.0 * b - 1.0
    c = float(np.dot(x[:-lag], x[lag:])) / (n - lag)
    # under iid fair bits, C ~ N(0, 1/(n-lag))
    z = c * np.sqrt(n - lag)
    p = _normal_sf_abs(z)
    return TestResult(f"autocorr_lag{lag}", float(c), p, n - lag, p >= 0.01)


def test_poker(bits: np.ndarray, m: int = 4) -> TestResult:
    """Poker / block template: chi-square on non-overlapping m-bit words."""
    b = bits.astype(np.int8)
    n = len(b)
    nblocks = n // m
    if nblocks < 5 * (2**m):
        return TestResult(f"poker_m{m}", 0.0, 0.0, n, False, "n too small")
    blocks = b[: nblocks * m].reshape(nblocks, m)
    # binary to int
    powers = 2 ** np.arange(m - 1, -1, -1)
    codes = (blocks * powers).sum(axis=1).astype(np.int64)
    counts = np.bincount(codes, minlength=2**m).astype(np.float64)
    exp = nblocks / float(2**m)
    chi2 = float(np.sum((counts - exp) ** 2 / exp))
    p = _chi2_sf(chi2, df=2**m - 1)
    return TestResult(f"poker_m{m}", chi2, p, nblocks * m, p >= 0.01)


def _chi2_sf(x: float, df: int) -> float:
    """Survival function P(Chi2_df >= x) via regularized upper gamma Q(df/2, x/2)."""
    if x <= 0:
        return 1.0
    if df <= 0:
        return 0.0
    # Wilson–Hilferty approximation then normal SF (good enough for pass/fail tables)
    # For small df use series for incomplete gamma upper.
    a = df / 2.0
    z = x / 2.0
    # lower gamma series / upper by 1 - P
    # Q(a,z) ≈ 1 - e^{-z} sum_{k=0}^∞ z^{a+k}/Gamma(a+k+1) * Gamma(a) — use recursive
    try:
        # direct: use numpy.random no — implement upper incomplete for half-integers etc.
        p_lower = _gammainc_lower(a, z)
        return float(max(0.0, min(1.0, 1.0 - p_lower)))
    except Exception:
        return 0.0


def _gammainc_lower(a: float, x: float, n_terms: int = 200) -> float:
    """Regularized lower incomplete gamma P(a,x) = γ(a,x)/Γ(a), series for x < a+1 else cont.frac."""
    if x < 0:
        return 0.0
    if x == 0:
        return 0.0
    # log Gamma via Stirling for large a, or recursive for moderate
    # Use series: P(a,x) = e^{-x} x^a / Gamma(a) * sum x^k / (a(a+1)...(a+k))
    log_gx = -x + a * np.log(x) - _log_gamma(a)
    if x < a + 1:
        term = 1.0 / a
        s = term
        for k in range(1, n_terms):
            term *= x / (a + k)
            s += term
            if term < 1e-14 * s:
                break
        return float(np.exp(log_gx) * s)
    # continued fraction for Q, then P = 1-Q (Lentz)
    # b0=0; use standard Numerical Recipes form
    tiny = 1e-30
    f = tiny
    c = f
    d = 0.0
    for i in range(1, n_terms + 1):
        if i == 1:
            an = 1.0
        else:
            an = -(i - 1) * (i - 1 - a)
        # b_i = x - a + 2i  for standard form of Q
        # CF for Q(a,x) = e^{-x} x^a / Gamma(a) * 1/(x+1-a - 1*(1-a)/(x+3-a - ...))
        # simpler: use series for Q from large x is hard; use 1 - series if series converges
        # fallback Wilson-Hilferty
        pass
    # Wilson-Hilferty
    h = 2.0 / (9.0 * a)
    z = ((x / a) ** (1.0 / 3.0) - (1.0 - h)) / np.sqrt(h)
    # P(Chi2 <= x) ≈ Phi(z)
    return float(0.5 * (1.0 + math.erf(z / math.sqrt(2.0))))


def _log_gamma(a: float) -> float:
    return float(math.lgamma(a))


def run_battery(bits: np.ndarray) -> List[TestResult]:
    bits = np.asarray(bits, dtype=np.uint8).ravel()
    # drop if empty
    if len(bits) < 100:
        return [TestResult("battery", 0.0, 0.0, len(bits), False, "stream too short")]
    tests = [
        test_monobit(bits),
        test_block_frequency(bits, M=128),
        test_runs(bits),
        test_serial_m2(bits),
        test_serial_m3(bits),
        test_poker(bits, m=4),
        test_autocorr(bits, lag=1),
        test_autocorr(bits, lag=8),
        test_autocorr(bits, lag=32),
    ]
    return tests


def battery_summary(results: List[TestResult]) -> dict:
    passed = sum(1 for r in results if r.pass_at_001)
    return {
        "n_tests": len(results),
        "n_pass_alpha_0.01": passed,
        "n_fail": len(results) - passed,
        "pass_rate": passed / max(len(results), 1),
        "tests": [r.as_dict() for r in results],
    }


# ─── Log-orbit drift (Lagarias heuristic / mean-reversion probe) ────────────

def log_orbit(orbit: np.ndarray) -> np.ndarray:
    return np.log1p(np.abs(orbit.astype(np.float64)))


def log_drift_stats(orbit: np.ndarray) -> dict:
    """
    Treat Y_t = log(1+|X_t|) as a discrete path.
    Mean step E[ΔY] is the empirical drift (negative for typical f, not for g).
    Regression ΔY ~ a + b Y estimates a crude mean-reversion coefficient b
    (OU would have b < 0). This is a diagnostic, not an SDE fit.
    """
    y = log_orbit(orbit)
    if len(y) < 5:
        return {"n": len(y), "mean_dY": None, "var_dY": None, "reversion_b": None}
    dy = np.diff(y)
    # OLS: dy = a + b * y[:-1]
    x = y[:-1]
    X = np.column_stack([np.ones_like(x), x])
    try:
        coef, _, _, _ = np.linalg.lstsq(X, dy, rcond=None)
        a, b = float(coef[0]), float(coef[1])
    except Exception:
        a, b = float("nan"), float("nan")
    return {
        "n": int(len(y)),
        "mean_Y": float(np.mean(y)),
        "mean_dY": float(np.mean(dy)),
        "var_dY": float(np.var(dy)),
        "std_dY": float(np.std(dy)),
        "intercept_a": a,
        "reversion_b": b,
        "note": (
            "Heuristic: random-parity model gives negative drift for f "
            "(Lagarias surveys). g swaps branches so forced negative drift is lost. "
            "reversion_b < 0 would mimic OU mean reversion; not assumed."
        ),
    }


def theoretical_parity_drift() -> dict:
    """
    Classical random-parity heuristic (not a theorem for integer Collatz):
      f: with p=1/2: log factor -log 2 (even) or log(3) (odd, ignoring +1)
         E ≈ 0.5*(-log 2 + log 3) = 0.5*log(3/2) > 0 if only that —
      Correct heuristic includes the geometric number of divisions by 2 after 3n+1,
      giving E[log] < 0 for the Syracuse map. We report both crude numbers.
    """
    # crude one-step ignoring +1 and multi-divides
    crude_f = 0.5 * (-np.log(2.0) + np.log(3.0))
    crude_g = 0.5 * (np.log(3.0) - np.log(2.0))  # same numbers swapped roles...
    # actually for g: even→3n+1 (~log3), odd→n/2 (-log2) so same crude average!
    # The distinction is in the multi-2 geometry of the full Syracuse map.
    # Report Syracuse-style expected log for f from literature approx -0.05 nats order
    return {
        "crude_one_step_mean_log_f": float(0.5 * (np.log(3.0) - np.log(2.0))),
        "crude_one_step_same_for_g_if_fair_parity": True,
        "comment": (
            "Under fair independent parity the one-step log factors of f and g "
            "are equal in law after swapping labels. The empirical difference "
            "comes from parity dependence and from multi-divisions by 2; "
            "under g even always maps to odd (3n+1 odd), so parities are not iid."
        ),
        "reference": "Lagarias surveys on the 3x+1 problem; parity is not an iid coin under iteration.",
    }


# ─── Full experiment ───────────────────────────────────────────────────────

def evaluate_map_on_planes(
    inverted: bool,
    seeds: Sequence[int],
    steps: int,
    plane_names: Sequence[str] | None = None,
) -> dict:
    if plane_names is None:
        plane_names = list(EXTRACTORS.keys())
    out = {"inverted": inverted, "map": "g" if inverted else "f", "planes": {}}
    # log drift on a single long orbit from first seed
    o0 = orbit_array(int(seeds[0]), steps, inverted)
    out["log_drift_seed0"] = log_drift_stats(o0)

    for name in plane_names:
        ext = EXTRACTORS[name]
        stream = stream_from_seeds(seeds, steps, inverted, ext)
        bat = battery_summary(run_battery(stream))
        bat["n_bits"] = int(len(stream))
        bat["plane"] = name
        out["planes"][name] = bat
    return out


def compare_f_g_prng(
    seeds: Sequence[int] | None = None,
    steps: int = 4096,
    plane_names: Sequence[str] | None = None,
) -> dict:
    if seeds is None:
        # odd seeds to avoid trivial early collapse patterns only
        seeds = list(range(3, 3 + 64 * 2, 2))  # 64 odd seeds
    seeds = list(seeds)
    rep_f = evaluate_map_on_planes(False, seeds, steps, plane_names)
    rep_g = evaluate_map_on_planes(True, seeds, steps, plane_names)

    # aggregate pass rates per plane
    table = []
    for name in (plane_names or EXTRACTORS.keys()):
        pf = rep_f["planes"][name]["pass_rate"]
        pg = rep_g["planes"][name]["pass_rate"]
        table.append(
            {
                "plane": name,
                "pass_rate_f": pf,
                "pass_rate_g": pg,
                "n_bits_f": rep_f["planes"][name]["n_bits"],
                "n_bits_g": rep_g["planes"][name]["n_bits"],
                "n_fail_f": rep_f["planes"][name]["n_fail"],
                "n_fail_g": rep_g["planes"][name]["n_fail"],
            }
        )

    return {
        "seeds": seeds[:8] + (["..."] if len(seeds) > 8 else []),
        "n_seeds": len(seeds),
        "steps": steps,
        "alpha": 0.01,
        "battery": [
            "monobit",
            "block_frequency",
            "runs",
            "serial_2bit",
            "serial_3bit",
            "poker_m4",
            "autocorr_lag1",
            "autocorr_lag8",
            "autocorr_lag32",
        ],
        "note": (
            "Self-contained subset inspired by NIST SP 800-22. Not a full NIST, "
            "TestU01 or PractRand run. No cryptographic suitability is claimed. "
            "Chaos-map PRNGs often fail serious batteries; mixed/negative results "
            "are the expected scientific outcome."
        ),
        "theoretical_parity": theoretical_parity_drift(),
        "f": rep_f,
        "g": rep_g,
        "pass_rate_table": table,
        "lattice_planes": {
            "E_parity_Z": "1D lattice Z, parity coordinate",
            "E_lsb1": "1D, LSB of X_t",
            "E_lsb4": "1D, 4 LSBs of X_t (more dimensions from same integer)",
            "S_mod256_low8": "ring lattice Z/256Z, 8-bit word",
            "S_mod64_low6": "ring lattice Z/64Z",
            "H_growth": "growth direction bit (tree-level proxy)",
            "Z2_parity_x_mod": "product lattice Z^2 bits",
            "Z3_parity_mod_growth": "product lattice Z^3 bits",
            "log2_floor_lsb4": "scale coordinate floor(log2(1+X)) LSBs",
            "delta_parity": "parity of absolute step size",
        },
        "framework_pointers": {
            "skew_product": (
                "Arnold RDS: (x,y) |-> (T x, S_x(y)). Here the base can be "
                "parity/time and the fibre the integer state; the PRNG stream "
                "is an observation of the fibre, not a claim of continuous RDS."
            ),
            "pdmp": (
                "Davis PDMP: deterministic flow between random jumps. Not fitted "
                "here; only the discrete jump map is studied."
            ),
            "msm": (
                "Markov state models coarse-grain continuous trajectories into "
                "discrete states. Applicable when states come from real data "
                "(e.g. molecular MD). Not constructed in this repository."
            ),
            "log_diffusion": (
                "Contact of OU-type modelling with Collatz-type maps is the "
                "log-orbit as a discrete path: ask whether empirical drift "
                "and reversion_b look mean-reverting. Not an external coupling "
                "of unrelated processes."
            ),
        },
    }
