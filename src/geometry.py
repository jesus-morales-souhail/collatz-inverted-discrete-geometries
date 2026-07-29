"""
Regla 7: Geometría esférica (Cúpula).
Estados normalizados viven sobre la esfera unitaria.
"""

from __future__ import annotations

import numpy as np


def project_to_sphere(x: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """Normaliza el vector de estado a la esfera unitaria."""
    v = np.asarray(x, dtype=float).ravel()
    n = float(np.linalg.norm(v))
    if n < eps:
        # vector nulo → polo arbitrario
        out = np.zeros_like(v)
        out[0] = 1.0
        return out
    return v / n


def spherical_geodesic_distance(a: np.ndarray, b: np.ndarray) -> float:
    """Distancia geodésica en la esfera: arccos(<a,b>)."""
    ua = project_to_sphere(a)
    ub = project_to_sphere(b)
    cos = float(np.clip(np.dot(ua, ub), -1.0, 1.0))
    return float(np.arccos(cos))


def spherical_target_from_direction(direction: int, dim: int) -> np.ndarray:
    """
    Estado objetivo sobre la esfera: sesgado hacia long (+e0) o short (-e0),
    con resto uniforme normalizado.
    """
    t = np.ones(dim, dtype=float)
    t[0] = 2.0 if direction > 0 else 0.25
    if direction < 0 and dim > 1:
        t[1] = 2.0
    return project_to_sphere(t)


def dome_state_ok(
    signal_vector: np.ndarray,
    direction: int,
    max_geodesic: float = 1.2,
) -> tuple[bool, float]:
    """Acepta si el estado actual no está demasiado lejos del objetivo en la cúpula."""
    target = spherical_target_from_direction(direction, len(signal_vector))
    d = spherical_geodesic_distance(signal_vector, target)
    return d <= max_geodesic, d
