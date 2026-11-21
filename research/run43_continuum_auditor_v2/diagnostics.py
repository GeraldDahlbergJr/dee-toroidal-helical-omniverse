"""Roundoff and angular diagnostics for the Run 43 v1 momentum stencil."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run43_continuum_auditor_v1"))
from auditor import G, contraction, metric_vectors, radial_first, radial_second

EPS = np.finfo(np.float64).eps
SAFETY_FACTOR = 64
ORDER_MARGIN = 10


def divergence(audit, w, spectral=False):
    """Evaluate with the input dtype or replace only the angular derivatives."""
    vectors, gradients = metric_vectors(audit)
    spacing = [w.dtype.type(h) for h in (audit.dr, audit.dt, audit.dp)]

    def derivative(i, value, order=1):
        if i == 0:
            return radial_first(value, spacing[0])
        if spectral:
            n = value.shape[i + 1]
            k = np.arange(n // 2 + 1)
            multiplier = (1j * k)**order
            if n % 2 == 0 and order == 1:
                multiplier[-1] = 0
            shape = [1] * value.ndim
            shape[i + 1] = len(k)
            return np.fft.irfft(np.fft.rfft(value, axis=i + 1) *
                                multiplier.reshape(shape), n=n, axis=i + 1)
        return audit.angular_derivative(value, i + 1, spacing[i], 2)

    first = [derivative(i, w) for i in range(3)]
    out = np.zeros_like(w)
    for i in range(3):
        for j in range(i, 3):
            if i == j == 0:
                second = radial_second(w, spacing[0])
            elif spectral and i == j:
                second = derivative(i, w, 2)
            else:
                second = derivative(j, first[i])
            out += contraction(second, vectors[j], vectors[i])
            if i != j:
                out += contraction(second, vectors[i], vectors[j])
                out += contraction(first[i], gradients[j][i], vectors[j])
            if gradients[i][j] is not None:
                out += contraction(first[j], gradients[i][j], vectors[i])
    return out


def operator_scale(audit, psi, w, rhs):
    """Pointwise uncancelled scale using stencil row sums and vector norms."""
    vectors, gradients = metric_vectors(audit)
    lengths = [np.linalg.norm(vector, axis=0) for vector in vectors]
    dr1 = np.full_like(audit.r, 1 / audit.dr)
    dr2 = np.full_like(audit.r, 4 / audit.dr**2)
    dr1[[0, -1]] *= 4
    dr2[[0, -1]] *= 3
    d1 = [dr1, 1 / audit.dt, 1 / audit.dp]
    gain = np.zeros(audit.shape)
    for i in range(3):
        for j in range(3):
            d2 = dr2 if i == j == 0 else d1[i] * d1[j]
            gain += d2 * lengths[j] * lengths[i]
            if gradients[i][j] is not None:
                gain += d1[j] * np.linalg.norm(gradients[i][j], axis=0) * lengths[i]
    return ((8/3) * psi**-10 * gain * np.linalg.norm(w, axis=0).max()
            + np.linalg.norm(rhs, axis=0))


def classify(residual, floor):
    if residual <= floor:
        return "roundoff-limited"
    if residual <= ORDER_MARGIN * floor:
        return "near-floor"
    return "above-floor"


def rms(audit, psi, value, mask):
    return audit.norms(value, psi, mask)["volume_weighted_rms"]


def angular_tail(value, axis):
    """Diagnostic on native samples, weighted by the squared k² multiplier."""
    n = value.shape[axis]
    k = np.arange(n // 2 + 1)
    spectrum = np.fft.rfft(value, axis=axis, norm="forward")
    energy = np.sum(abs(spectrum)**2, axis=tuple(i for i in range(value.ndim) if i != axis))
    multiplicity = np.full(k.shape, 2)
    multiplicity[0] = 1
    if n % 2 == 0:
        multiplicity[-1] = 1
    weighted = multiplicity * k.astype(float)**4 * energy
    amplitude = np.sqrt(np.sum(multiplicity * energy))
    derivative = np.sqrt(weighted.sum())
    noise = 64 * EPS * (n // 2)**2 * amplitude
    active = derivative > noise
    # Compose centered differences: sin(kh)^2 / (kh)^2.
    flagged = np.sinc(2 * k / n)**2 < 0.9
    tail = float(np.sqrt(weighted[flagged].sum()) / derivative) if active else 0.0
    nyquist = float(np.sqrt(energy[-1]) / amplitude) if n % 2 == 0 and amplitude else 0.0
    return {"high_mode_second_derivative_fraction": tail,
            "nyquist_amplitude_fraction": nyquist,
            "classification": "underresolved-risk" if tail > 0.05 or nyquist > 64 * EPS
            else "screen-clear-on-samples"}


def physical_rhs(audit, psi):
    return 8 * np.pi * G * psi**-4 * audit.S
