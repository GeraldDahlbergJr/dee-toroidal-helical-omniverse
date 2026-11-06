"""Certify exact q6 Jacobian invertibility at archived spline coefficients."""
import argparse
from fractions import Fraction
import json
import math

import numpy as np


def as_fraction(value):
    return Fraction(float(value))


def absolute_difference_upper(values, axis, *, periodic):
    if periodic:
        difference = values - np.roll(values, 1, axis=axis)
    else:
        difference = np.diff(values, axis=axis)
    maximum = np.max(np.abs(difference))
    return as_fraction(np.nextafter(maximum, np.inf))


def upper_float(value):
    return math.nextafter(float(value), math.inf)


def lower_float(value):
    return math.nextafter(float(value), -math.inf)


def certify(path):
    with np.load(path) as archive:
        n = int(archive["n"])
        coefficients = archive["coefficients"]
    if n < 8 or coefficients.shape != (4, n + 2, n, n):
        raise ValueError("Expected the open-radial, periodic-angular degree-four space")
    if coefficients.dtype != np.float64 or not np.isfinite(coefficients).all():
        raise ValueError("Expected finite IEEE float64 coefficients")

    pi_lower = Fraction("3.14159")
    pi_upper = Fraction("3.14160")
    poincare_sq = Fraction(".034")
    scalar = coefficients[0]
    psi_lower = 1 + min(Fraction(0), as_fraction(scalar.min()))
    psi_upper = 1 + max(Fraction(0), as_fraction(scalar.max()))
    if psi_lower <= 0:
        raise ValueError("Positive conformal factor could not be certified")

    vector = coefficients[1:]
    radial = np.pad(vector, ((0, 0), (1, 1), (0, 0), (0, 0)))
    derivative_limits = [
        4 * n / Fraction(".20") * absolute_difference_upper(radial, 1, periodic=False),
        4 * n / (2 * pi_lower * Fraction(".05")) * absolute_difference_upper(vector, 2, periodic=True),
        4 * n / (2 * pi_lower * Fraction(".75")) * absolute_difference_upper(vector, 3, periodic=True),
    ]
    # Each limit bounds one chart-orthonormal derivative, componentwise.
    # The sum over all nine entries bounds the Frobenius norm.
    strain_upper = 6 * sum(derivative_limits)

    g = Fraction(".001")
    b_upper = Fraction(".068614375") * (400 + 9 / Fraction(".75")**2)
    potential_upper = 2 * pi_upper * g * (
        5 * Fraction(".0127638") * psi_upper**4 + b_upper
    )
    scalar_lower = 1 - poincare_sq * potential_upper
    vector_lower = Fraction(".969")
    if scalar_lower <= 0:
        raise ValueError("Scalar diagonal block coercivity was not certified")

    cross_scalar = Fraction(".185") * strain_upper / (2 * psi_lower**7)
    cross_vector = poincare_sq * 48 * pi_upper * g * psi_upper**5 * Fraction(".62")
    schur_upper = cross_scalar * cross_vector / (scalar_lower * vector_lower)
    return {
        "n": n,
        "psi_lower": lower_float(psi_lower),
        "psi_upper": upper_float(psi_upper),
        "strain_frobenius_upper": upper_float(strain_upper),
        "scalar_gradient_coercivity_lower": lower_float(scalar_lower),
        "vector_gradient_coercivity_lower": lower_float(vector_lower),
        "scalar_cross_gradient_norm_upper": upper_float(cross_scalar),
        "vector_cross_gradient_norm_upper": upper_float(cross_vector),
        "schur_gradient_norm_upper": upper_float(schur_upper),
        "exact_q6_jacobian_invertible": schur_upper < 1,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archives", nargs="+")
    args = parser.parse_args()
    results = [certify(path) for path in args.archives]
    print(json.dumps(results, indent=2))
    if not all(row["exact_q6_jacobian_invertible"] for row in results):
        raise SystemExit(1)
