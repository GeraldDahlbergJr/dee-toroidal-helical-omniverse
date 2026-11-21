"""Run 43 / W7-X-inspired flux-potential pilot.

A discrete exterior-derivative identity check on a regular toroidal coordinate
grid, NOT a reconstruction of Wendelstein 7-X's VMEC equilibrium.
Uses angular fivefold periodicity, without claiming experimental validation.
"""
import numpy as np


def centered_periodic(a, spacing, axis):
    return (np.roll(a, -1, axis=axis) - np.roll(a, 1, axis=axis)) / (2 * spacing)


def derivatives(a, dr, dt, dp):
    return (
        np.gradient(a, dr, axis=0, edge_order=2),
        centered_periodic(a, dt, 1),
        centered_periodic(a, dp, 2),
    )


def audit(n, epsilon=0.025, m=2, n_field_periods=5):
    nr = n + 1
    radius = np.linspace(0.2, 0.8, nr)[:, None, None]
    theta = (2 * np.pi * np.arange(n) / n)[None, :, None]
    phi = (2 * np.pi * np.arange(n) / n)[None, None, :]
    dr = 0.6 / n
    dt = dp = 2 * np.pi / n
    chi = m * theta - n_field_periods * phi
    # Magnetic potential one-form A = Ar dr + At dtheta + Ap dphi.
    Ar = np.zeros((nr, n, n))
    At = np.broadcast_to(0.3 * radius**2, Ar.shape).copy()
    Ap = np.broadcast_to(-0.2 * radius**2 + epsilon * radius**3 * np.cos(chi), Ar.shape).copy()
    Dr, Dt, Dp = lambda a: derivatives(a, dr, dt, dp)
    ar_r, ar_t, ar_p = Dr(Ar)
    at_r, at_t, at_p = Dr(At)
    ap_r, ap_t, ap_p = Dr(Ap)

    # The two-form F=dA. In curvilinear coordinates these are covariant
    # coordinate components, NOT orthonormal physical B components.
    F_rt = at_r - ar_t
    F_rp = ap_r - ar_p
    F_tp = ap_t - at_p

    # Discrete dF = 0; commuting coordinate differences are essential.
    closure = Dr(F_tp)[0] - Dr(F_rp)[1] + Dr(F_rt)[2]
    closure_max = float(np.max(np.abs(closure)))

    # Independent analytic manufactured fields test derivative accuracy.
    exact_rt = np.broadcast_to(0.6 * radius, Ar.shape)
    exact_rp = np.broadcast_to(-0.4 * radius + 3 * epsilon * radius**2 * np.cos(chi), Ar.shape)
    exact_tp = np.broadcast_to(-m * epsilon * radius**3 * np.sin(chi), Ar.shape)
    error = float(max(
        np.max(np.abs(F_rt - exact_rt)),
        np.max(np.abs(F_rp - exact_rp)),
        np.max(np.abs(F_tp - exact_tp)),
    ))
    return {"n": n, "closure_max": closure_max, "field_error": error}


def main():
    results = [audit(n) for n in (20, 40, 80)]
    for r in results:
        print(f"n={r['n']:3d} |dF|_inf={r['closure_max']:.3e} "
              f"|F-F_exact|_inf={r['field_error']:.3e}")
        assert r["closure_max"] < 1e-10, "Discrete exterior-derivative closure failed"
    for coarse, fine in zip(results, results[1:]):
        order = np.log2(coarse["field_error"] / fine["field_error"])
        print(f"manufactured field order {coarse['n']}->{fine['n']}: {order:.4f}")
        assert order >= 1.5, "Field convergence order below 1.5"


if __name__ == "__main__":
    main()
