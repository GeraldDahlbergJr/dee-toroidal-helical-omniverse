"""
Dynamic Equilibrium Events (DEE)
Toroidal-Helical Geometry Demonstration

Author: Gerald Dahlberg Jr.
ORCID: 0009-0001-0672-7636
DOI: 10.5281/zenodo.21477633

This program generates a parameterized toroidal-helical trajectory
for computational exploration of the DEE framework.

It is a computational demonstration, not experimental validation.
"""

import numpy as np


def toroidal_helical_path(
    major_radius=1.0,
    minor_radius=0.25,
    helical_mode=3,
    samples=1000,
):
    """Generate a three-dimensional toroidal-helical trajectory."""

    theta = np.linspace(0.0, 2.0 * np.pi, samples)

    radius = (
        major_radius
        + minor_radius * np.cos(helical_mode * theta)
    )

    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    z = minor_radius * np.sin(helical_mode * theta)

    return theta, x, y, z


def main():
    theta, x, y, z = toroidal_helical_path()

    print("DEE Toroidal-Helical Geometry")
    print("-----------------------------")
    print(f"Samples: {len(theta)}")
    print(f"x range: {x.min():.4f} to {x.max():.4f}")
    print(f"y range: {y.min():.4f} to {y.max():.4f}")
    print(f"z range: {z.min():.4f} to {z.max():.4f}")


if __name__ == "__main__":
    main()
