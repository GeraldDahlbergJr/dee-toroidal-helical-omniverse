"""
Dynamic Equilibrium Events (DEE)
Time-Dependent Evolution Demonstration

Author: Gerald Dahlberg Jr.
ORCID: 0009-0001-0672-7636
DOI: 10.5281/zenodo.21477633

Computational exploration of a simplified
time-dependent DEE evolution model.

This is a computational demonstration and is
not experimental confirmation of the theory.
"""

import numpy as np


def dee_evolution(
    delta_initial=1.0,
    equilibrium=0.0,
    relaxation_rate=0.35,
    drive_amplitude=0.20,
    drive_frequency=1.0,
    duration=20.0,
    samples=2000,
):
    """Generate a simple driven, relaxing DEE trajectory."""

    time = np.linspace(0.0, duration, samples)
    dt = time[1] - time[0]

    delta = np.zeros(samples)
    delta[0] = delta_initial

    for i in range(1, samples):
        drive = drive_amplitude * np.sin(
            2.0 * np.pi * drive_frequency * time[i - 1]
        )

        rate = (
            -relaxation_rate * (delta[i - 1] - equilibrium)
            + drive
        )

        delta[i] = delta[i - 1] + rate * dt

    return time, delta


def main():
    time, delta = dee_evolution()

    print("DEE Time-Dependent Evolution")
    print("============================")
    print(f"Samples: {len(time)}")
    print(f"Initial delta: {delta[0]:.6f}")
    print(f"Final delta: {delta[-1]:.6f}")
    print(f"Minimum delta: {delta.min():.6f}")
    print(f"Maximum delta: {delta.max():.6f}")


if __name__ == "__main__":
    main()
