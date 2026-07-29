"""
Constitutive flux on a 1D grid, and the same product on Collatz cell fields.

Continuous forms:
  q = -k ∇T          (Fourier)
  J =  σ E = -σ ∇V   (Ohm)

Discrete (spacing Δx, forward difference):
  (∇u)_n ≈ (u_{n+1} - u_n)/Δx = (d * u)_n
  flux_n = -κ (d * u)_n

With constant κ and a circular embedding, the DFT factors the product:
  F{flux}_k = -κ H_k F{u}_k
  H_k = (exp(2π i k/N) - 1)/Δx

An orbit (X_t) supplies scalar fields φ on cells t = 0,...,T-1
(e.g. log(1+|X|), X mod N, level-weighted log). The maps f and g are not
claimed to be linear transport; only the flux skeleton on those fields is shared.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np

from .collatz_maps import f_normal, g_inverted


# ─── Discrete gradient (Fourier / Ohm skeleton) ────────────────────────────

def forward_diff_kernel(dx: float = 1.0) -> np.ndarray:
    """Spatial FIR kernel d with d*u ≈ ∇u (forward)."""
    return np.array([-1.0, 1.0], dtype=np.float64) / float(dx)


def apply_conv_full(u: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    """Linear convolution; result length len(u)+len(kernel)-1, then trim to len(u) (same as causal pad)."""
    u = np.asarray(u, dtype=np.float64)
    y = np.convolve(u, kernel, mode="full")
    # align so y[n] uses u[n], u[n+1] for forward diff: pad u on the right
    # simpler: explicit forward difference length N-1, pad last
    return y


def discrete_gradient(u: np.ndarray, dx: float = 1.0) -> np.ndarray:
    """(∇u)_n for n=0..N-2, length N-1; pad end with 0 for length N."""
    u = np.asarray(u, dtype=np.float64)
    g = (u[1:] - u[:-1]) / dx
    return np.concatenate([g, [0.0]])


def constitutive_flux(u: np.ndarray, kappa: float, dx: float = 1.0) -> np.ndarray:
    """
    flux = -κ ∇u   (same form for Fourier heat and Ohm current).
    """
    return -float(kappa) * discrete_gradient(u, dx)


def dft_product_identity(
    u: np.ndarray, kappa: float, dx: float = 1.0
) -> Dict[str, np.ndarray | float]:
    """
    Check: DFT(flux) ≈ -κ · H(ω) · DFT(u)
    where H is the frequency response of the forward-difference operator.
    """
    u = np.asarray(u, dtype=np.float64)
    N = len(u)
    flux = constitutive_flux(u, kappa, dx)

    U = np.fft.fft(u)
    Q = np.fft.fft(flux)

    # Frequency response of forward difference on circular embedding:
    # for circular: (u_{n+1}-u_n)/dx  ↔  (e^{2πik/N} - 1)/dx
    k = np.arange(N)
    H = (np.exp(2j * np.pi * k / N) - 1.0) / dx
    Q_pred = -kappa * H * U

    # linear (non-circular) gradient pads last zero — compare on first N-1 via circular proxy
    err = float(np.linalg.norm(Q - Q_pred) / (np.linalg.norm(Q) + 1e-15))
    return {
        "U": U,
        "Q": Q,
        "Q_pred": Q_pred,
        "H": H,
        "relative_err_circular_model": err,
        "product_form": "F{flux}_k = -κ · H_k · F{u}_k",
    }


# ─── Collatz fields on the cell line ───────────────────────────────────────

def collatz_orbit(n0: int, steps: int, inverted: bool) -> np.ndarray:
    fn = g_inverted if inverted else f_normal
    x = max(int(n0), 0)
    seq = [float(x)]
    for _ in range(steps - 1):
        x = fn(x)
        seq.append(float(x))
    return np.asarray(seq, dtype=np.float64)


def potential_log(orbit: np.ndarray) -> np.ndarray:
    """Potential-like field on cells (scale-friendly)."""
    return np.log1p(np.abs(orbit))


def branch_conductivity(orbit: np.ndarray, inverted: bool) -> np.ndarray:
    """κ_t from parity: larger on the 3n+1 branch, smaller on n/2."""
    parity = (orbit.astype(np.int64) % 2)
    if inverted:
        kappa = np.where(parity == 0, 3.0, 0.5)  # g: even → 3n+1
    else:
        kappa = np.where(parity == 0, 0.5, 3.0)  # f: odd → 3n+1
    return kappa.astype(np.float64)


def flux_variable_kappa(u: np.ndarray, kappa: np.ndarray, dx: float = 1.0) -> np.ndarray:
    """flux_n = -κ_n · (∇u)_n  (Hadamard product in space)."""
    grad = discrete_gradient(u, dx)
    return -kappa * grad


@dataclass
class PlaneChannels:
    """Scalar fields on the same cell line."""

    cell: np.ndarray
    euclidean: np.ndarray
    spherical: np.ndarray
    hyperbolic: np.ndarray


def build_planes(orbit: np.ndarray, N_mod: int = 64) -> PlaneChannels:
    """
    Fields on cells t=0..T-1:
      E: log(1+|X|);  S: X mod N;  H: log weighted by a growth-level proxy.
    """
    T = len(orbit)
    cell = np.arange(T, dtype=np.float64)
    eucl = potential_log(orbit)
    sph = (orbit.astype(np.int64) % N_mod).astype(np.float64)
    growth = np.concatenate([[0.0], np.sign(np.diff(orbit))])
    level = np.maximum.accumulate(np.cumsum(growth))
    hyp = potential_log(orbit) * (1.0 + 0.25 * level)
    return PlaneChannels(cell=cell, euclidean=eucl, spherical=sph, hyperbolic=hyp)


def plane_constitutive_spectra(
    planes: PlaneChannels,
    kappa: float | np.ndarray,
    dx: float = 1.0,
) -> Dict[str, dict]:
    """
    For each plane field φ, compute flux = -κ ∇φ and DFT product residuals.
    """
    out = {}
    fields = {
        "E_euclidean": planes.euclidean,
        "S_spherical": planes.spherical,
        "H_hyperbolic": planes.hyperbolic,
        "cell_base": planes.cell,
    }
    for name, u in fields.items():
        if np.isscalar(kappa):
            flux = constitutive_flux(u, float(kappa), dx)
            info = dft_product_identity(u, float(kappa), dx)
        else:
            flux = flux_variable_kappa(u, np.asarray(kappa, dtype=np.float64), dx)
            # variable κ: DFT product with constant κ does not hold exactly
            info = {
                "product_form": "flux = -κ_t ⊙ ∇u  (variable κ; full DFT factorisation only if κ constant)",
                "relative_err_circular_model": None,
            }
        U = np.fft.fft(u)
        Q = np.fft.fft(flux)
        out[name] = {
            "potential": u,
            "flux": flux,
            "U_abs": np.abs(U),
            "Q_abs": np.abs(Q),
            "info": info,
        }
    return out


def compare_f_g_constitutive(
    n0: int = 10,
    steps: int = 64,
    kappa: float = 1.0,
    N_mod: int = 64,
) -> dict:
    """
    Same constitutive skeleton on planes from f-orbit vs g-orbit.
    Reports spectral flux energy and whether circular DFT identity holds
    for constant κ.
    """
    of = collatz_orbit(n0, steps, inverted=False)
    og = collatz_orbit(n0, steps, inverted=True)
    pf, pg = build_planes(of, N_mod), build_planes(og, N_mod)

    # constant κ (Ohm/Fourier exact product in frequency for circular model)
    spec_f = plane_constitutive_spectra(pf, kappa)
    spec_g = plane_constitutive_spectra(pg, kappa)

    # branch-dependent κ (Collatz "material" switches with parity)
    kf = branch_conductivity(of, inverted=False)
    kg = branch_conductivity(og, inverted=True)
    spec_f_var = plane_constitutive_spectra(pf, kf)
    spec_g_var = plane_constitutive_spectra(pg, kg)

    def energy(spec, key="Q_abs"):
        return {name: float(np.sum(ch[key] ** 2)) for name, ch in spec.items()}

    return {
        "n0": n0,
        "steps": steps,
        "constant_kappa": kappa,
        "dft_identity_err_f_E": spec_f["E_euclidean"]["info"].get("relative_err_circular_model"),
        "dft_identity_err_g_E": spec_g["E_euclidean"]["info"].get("relative_err_circular_model"),
        "flux_energy_const_f": energy(spec_f),
        "flux_energy_const_g": energy(spec_g),
        "flux_energy_branch_f": energy(spec_f_var),
        "flux_energy_branch_g": energy(spec_g_var),
        "structure": {
            "continuous": "q=-k∇T , J=σE=-σ∇V",
            "discrete": "flux=-κ (d*u)",
            "frequency": "F{flux}_k = -κ H_k F{u}_k  (κ constant, circular)",
            "on_planes": "same product on φ^E, φ^S, φ^H, cell index",
            "collatz": "nonlinear map supplies the field φ; κ may be constant or parity-dependent",
        },
    }
