"""
Classical DFT ↔ inverse DFT, then Quantum Fourier Transform (QFT matrix)
on Collatz parity / orbit signals — to score the exploration crack g vs f.

Pipeline
--------
1) Build real 1D signals from orbits of f (normal) and g (inverted):
   - parity bit stream  b_t = X_t mod 2
   - log-growth stream  s_t = log1p(X_t)

2) Classical Fourier:
   S(ω) = DFT[s],  s_hat = iDFT[S]  (reconstruction check)

3) Quantum Fourier Transform (simulated):
   Encode a normalized amplitude vector |ψ⟩ from the signal,
   apply unitary QFT_n on n = log2(N) qubits (N power of 2),
   measure |⟨k|QFT|ψ⟩|² as a probability mass over frequency bins.

4) Crack score:
   Distance between spectral measures of f vs g (TV, L2, KL, peak shift).
   Higher distance ⇒ easier to tell “grieta” in frequency domain.

5) Optional variational weight w(ω; θ) on frequencies; gradient ascent of
   separation score (simple differentiable “variable function”).

This is linear algebra / simulation — not a claim of hardware QFT advantage.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

from .collatz_maps import f_normal, g_inverted


# ═══════════════════════════════════════════════════════════════════════════
# Signals from Collatz
# ═══════════════════════════════════════════════════════════════════════════

def orbit_int(n0: int, steps: int, inverted: bool) -> np.ndarray:
    fn = g_inverted if inverted else f_normal
    x = max(int(n0), 0)
    seq = [x]
    for _ in range(steps - 1):
        x = fn(x)
        seq.append(x)
    return np.asarray(seq, dtype=np.float64)


def parity_signal(n0: int, steps: int, inverted: bool) -> np.ndarray:
    return (orbit_int(n0, steps, inverted) % 2).astype(np.float64)


def log_signal(n0: int, steps: int, inverted: bool) -> np.ndarray:
    return np.log1p(orbit_int(n0, steps, inverted))


def parity_matrix(
    seeds: List[int], steps: int, inverted: bool
) -> np.ndarray:
    """Rows = seeds, cols = time (parity)."""
    return np.stack([parity_signal(s, steps, inverted) for s in seeds], axis=0)


def mean_parity_spectrum(
    seeds: List[int], steps: int, inverted: bool
) -> np.ndarray:
    """Mean power spectrum of parity rows."""
    M = parity_matrix(seeds, steps, inverted)
    # zero-mean each row
    M = M - M.mean(axis=1, keepdims=True)
    F = np.fft.rfft(M, axis=1)
    return np.mean(np.abs(F) ** 2, axis=0)


# ═══════════════════════════════════════════════════════════════════════════
# Classical DFT / iDFT
# ═══════════════════════════════════════════════════════════════════════════

def dft(x: np.ndarray) -> np.ndarray:
    return np.fft.fft(np.asarray(x, dtype=np.complex128))


def idft(X: np.ndarray) -> np.ndarray:
    return np.fft.ifft(np.asarray(X, dtype=np.complex128))


def dft_roundtrip_error(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.complex128)
    y = idft(dft(x))
    return float(np.max(np.abs(y - x)))


def power_spectrum(x: np.ndarray) -> np.ndarray:
    X = dft(x - np.mean(x))
    return np.abs(X) ** 2


# ═══════════════════════════════════════════════════════════════════════════
# Quantum Fourier Transform (unitary matrix simulation)
# ═══════════════════════════════════════════════════════════════════════════

def next_pow2(n: int) -> int:
    p = 1
    while p < n:
        p <<= 1
    return p


def qft_matrix(n_qubits: int) -> np.ndarray:
    """
    Unitary QFT on n qubits, size N=2^n:
      F_{jk} = N^{-1/2} exp(2π i j k / N)
    (standard math QFT; same as DFT unitary form).
    """
    N = 1 << n_qubits
    j = np.arange(N).reshape(-1, 1)
    k = np.arange(N).reshape(1, -1)
    omega = np.exp(2j * np.pi / N)
    F = (omega ** (j * k)) / np.sqrt(N)
    return F.astype(np.complex128)


def encode_amplitudes(x: np.ndarray, n_qubits: Optional[int] = None) -> Tuple[np.ndarray, int]:
    """
    Map a real signal to a unit vector |ψ⟩ of length 2^n (pad / truncate).
    Uses soft-positive embedding: abs + epsilon, then L2 normalize.
    """
    x = np.asarray(x, dtype=np.float64).ravel()
    # shift to non-negative for a simple Born-friendly encoding
    x = x - x.min()
    x = x + 1e-12
    N_need = len(x)
    if n_qubits is None:
        n_qubits = int(np.ceil(np.log2(max(N_need, 2))))
    N = 1 << n_qubits
    amp = np.zeros(N, dtype=np.complex128)
    m = min(N, N_need)
    amp[:m] = np.sqrt(x[:m])  # amplitude ~ sqrt(intensity)
    amp /= np.linalg.norm(amp) + 1e-15
    return amp, n_qubits


def apply_qft(psi: np.ndarray, n_qubits: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    |φ⟩ = QFT |ψ⟩. Returns (phi, probs=|φ|²).
    """
    psi = np.asarray(psi, dtype=np.complex128).ravel()
    N = len(psi)
    if n_qubits is None:
        if N & (N - 1):
            raise ValueError("psi length must be power of 2")
        n_qubits = int(np.log2(N))
    F = qft_matrix(n_qubits)
    if F.shape[0] != N:
        # pad/trim psi
        amp = np.zeros(F.shape[0], dtype=np.complex128)
        m = min(N, F.shape[0])
        amp[:m] = psi[:m]
        amp /= np.linalg.norm(amp) + 1e-15
        psi = amp
    phi = F @ psi
    probs = np.abs(phi) ** 2
    probs = probs / (probs.sum() + 1e-15)
    return phi, probs


def qft_roundtrip_fidelity(psi: np.ndarray, n_qubits: int) -> float:
    """QFT† QFT = I ⇒ fidelity of roundtrip."""
    F = qft_matrix(n_qubits)
    phi = F @ psi
    back = F.conj().T @ phi  # inverse QFT
    return float(np.abs(np.vdot(psi, back)) ** 2)


# ═══════════════════════════════════════════════════════════════════════════
# Crack scores in frequency / QFT domain
# ═══════════════════════════════════════════════════════════════════════════

def total_variation(p: np.ndarray, q: np.ndarray) -> float:
    p = p / (p.sum() + 1e-15)
    q = q / (q.sum() + 1e-15)
    return 0.5 * float(np.sum(np.abs(p - q)))


def l2_distance(p: np.ndarray, q: np.ndarray) -> float:
    return float(np.linalg.norm(p - q))


def spectral_entropy(p: np.ndarray) -> float:
    p = p / (p.sum() + 1e-15)
    p = p[p > 1e-15]
    return float(-np.sum(p * np.log(p)))


@dataclass
class CrackReport:
    dft_tv: float
    dft_l2: float
    qft_tv: float
    qft_l2: float
    qft_entropy_f: float
    qft_entropy_g: float
    reconstruction_err_f: float
    reconstruction_err_g: float
    qft_fidelity_f: float
    qft_fidelity_g: float
    crack_probability: float  # derived score in [0,1]
    notes: str


def crack_probability_from_distances(
    dft_tv: float, qft_tv: float, qft_l2: float
) -> float:
    """
    Map spectral separation to a [0,1] “probability of detecting the crack”
    if one randomly samples a frequency bin / computational basis outcome
    and tests which map it came from (heuristic, not Bayes-optimal).

    Uses a smooth map: 1 - exp(-α · score).
    """
    score = 0.45 * dft_tv + 0.45 * qft_tv + 0.10 * min(qft_l2, 2.0)
    return float(1.0 - np.exp(-2.2 * score))


def analyze_crack(
    seeds: List[int],
    steps: int = 64,
    n_qubits: int = 6,
) -> Tuple[CrackReport, dict]:
    """
    Build mean power spectra (classical) and mean QFT probs for f vs g.
    """
    # classical mean power on parity
    Pf = mean_parity_spectrum(seeds, steps, inverted=False)
    Pg = mean_parity_spectrum(seeds, steps, inverted=True)
    # align lengths (rfft same length if steps same)
    dft_tv = total_variation(Pf, Pg)
    dft_l2 = l2_distance(Pf / (Pf.sum() + 1e-15), Pg / (Pg.sum() + 1e-15))

    # QFT: encode mean log-orbit of each seed, average Born probs
    N = 1 << n_qubits
    probs_f = np.zeros(N)
    probs_g = np.zeros(N)
    rec_f = []
    rec_g = []
    fid_f = []
    fid_g = []

    for s in seeds:
        sf = log_signal(s, steps, inverted=False)
        sg = log_signal(s, steps, inverted=True)
        rec_f.append(dft_roundtrip_error(sf))
        rec_g.append(dft_roundtrip_error(sg))

        af, nq = encode_amplitudes(sf, n_qubits)
        ag, _ = encode_amplitudes(sg, n_qubits)
        _, pf = apply_qft(af, nq)
        _, pg = apply_qft(ag, nq)
        probs_f += pf
        probs_g += pg
        fid_f.append(qft_roundtrip_fidelity(af, nq))
        fid_g.append(qft_roundtrip_fidelity(ag, nq))

    probs_f /= len(seeds)
    probs_g /= len(seeds)
    qft_tv = total_variation(probs_f, probs_g)
    qft_l2 = l2_distance(probs_f, probs_g)
    p_crack = crack_probability_from_distances(dft_tv, qft_tv, qft_l2)

    report = CrackReport(
        dft_tv=dft_tv,
        dft_l2=dft_l2,
        qft_tv=qft_tv,
        qft_l2=qft_l2,
        qft_entropy_f=spectral_entropy(probs_f),
        qft_entropy_g=spectral_entropy(probs_g),
        reconstruction_err_f=float(np.mean(rec_f)),
        reconstruction_err_g=float(np.mean(rec_g)),
        qft_fidelity_f=float(np.mean(fid_f)),
        qft_fidelity_g=float(np.mean(fid_g)),
        crack_probability=p_crack,
        notes=(
            "crack_probability is a spectral separability score in [0,1], "
            "Spectral separability score in [0,1]; not a proof of divergence."
        ),
    )
    extras = {
        "Pf_classical": Pf,
        "Pg_classical": Pg,
        "probs_f_qft": probs_f,
        "probs_g_qft": probs_g,
        "freqs_rfft": np.fft.rfftfreq(steps),
    }
    return report, extras


# ═══════════════════════════════════════════════════════════════════════════
# Variational weight on QFT bins (variable function)
# ═══════════════════════════════════════════════════════════════════════════

def variational_separation(
    probs_f: np.ndarray,
    probs_g: np.ndarray,
    n_iter: int = 300,
    lr: float = 0.25,
) -> dict:
    """
    Learn w_k = softmax(θ_k) maximizing separability
      J(θ) = sum_k w_k |p_g(k) - p_f(k)|
    (absolute spectral contrast — reweights bins where maps differ most).
    """
    pf = probs_f / (probs_f.sum() + 1e-15)
    pg = probs_g / (probs_g.sum() + 1e-15)
    delta = pg - pf
    target = np.abs(delta)
    # init θ proportional to contrast
    theta = np.log(target + 1e-6)
    history = []
    for _ in range(n_iter):
        w = np.exp(theta - theta.max())
        w = w / (w.sum() + 1e-15)
        J = float(np.dot(w, target))
        history.append(J)
        grad = w * (target - J)
        theta = theta + lr * grad
    w = np.exp(theta - theta.max())
    w = w / (w.sum() + 1e-15)
    J = float(np.dot(w, target))
    top = np.argsort(-w * target)[:8]
    return {
        "J_final": J,
        "J_history": history,
        "weights": w,
        "delta": delta,
        "abs_delta": target,
        "top_bins": top.tolist(),
        "top_weights": w[top].tolist(),
        "top_delta": delta[top].tolist(),
        "variable_function": "w(ω;θ)=softmax(θ),  J=E_w[|p_g - p_f|]",
        "uniform_J": float(np.mean(target)),
    }
