"""Isolated Run 43 mixed-flux / SBP-SAT recovery benchmark.

Synthetic 1-D second-order elliptic problem: -u'' + u = f on [0,1].
The first-order auxiliary flux is j = D u.  Boundary values are imposed
weakly with symmetric consistent SAT (Nitsche form).
Not yet integrated into Einstein–DEE; no Run 43 12-check claims.
Requires NumPy.
"""
import numpy as np


def sbp_first_derivative(n):
    """Diagonal-norm SBP(2,1): H D + D.T H = B."""
    h = 1.0 / n
    H = np.eye(n + 1) * h
    H[0, 0] = H[-1, -1] = h / 2.0
    D = np.zeros((n + 1, n + 1))
    D[0, 0:2] = [-1, 1]
    D[-1, -2:] = [-1, 1]
    for i in range(1, n):
        D[i, i - 1:i + 2] = [-0.5, 0, 0.5]
    D /= h
    B = np.zeros_like(D)
    B[0, 0] = -1.0
    B[-1, -1] = 1.0
    return H, D, B


def solve(n, penalty=12.0):
    H, D, B = sbp_first_derivative(n)
    h = 1.0 / n
    x = np.linspace(0.0, 1.0, n + 1)
    exact = 1.0 + x + np.sin(np.pi * x)
    f = np.pi**2 * np.sin(np.pi * x) + exact
    g = np.zeros(n + 1)
    g[0], g[-1] = 1.0, 2.0
    W = np.zeros_like(D)
    W[0, 0] = W[-1, -1] = penalty / h
    # Mixed auxiliary flux j=D u, eliminated in the symmetric weak problem.
    # -B D and -D.T B are the two consistency/adjoint boundary SAT terms.
    # W enforces Dirichlet data; deleting the consistency terms gives only
    # first-order convergence in this manufactured example.
    A = D.T @ H @ D - B @ D - D.T @ B + H + W
    rhs = H @ f - D.T @ B @ g + W @ g
    u = np.linalg.solve(A, rhs)
    j = D @ u
    error = float(np.sqrt((u - exact) @ H @ (u - exact)))
    sbp_defect = float(np.max(np.abs(H @ D + D.T @ H - B)))
    weak_residual = float(np.max(np.abs(A @ u - rhs)))
    min_eigenvalue = float(np.linalg.eigvalsh(A)[0])
    return dict(n=n, h=h, error=error, sbp_defect=sbp_defect,
                weak_residual=weak_residual, min_eigenvalue=min_eigenvalue,
                flux_norm=float(np.linalg.norm(j)))


def main():
    grids = (16, 32, 64, 128)
    rows = [solve(n) for n in grids]
    for row in rows:
        print("n={n:3d} L2(H)={error:.9e} SBP={sbp_defect:.2e} "
              "res={weak_residual:.2e} lambda_min={min_eigenvalue:.4e}".format(**row))
        assert row["sbp_defect"] < 1e-12
        assert row["weak_residual"] < 1e-9
        assert row["min_eigenvalue"] > 0.0
    for previous, current in zip(rows, rows[1:]):
        order = float(np.log2(previous["error"] / current["error"]))
        print(f"order {previous['n']}->{current['n']}: {order:.6f}")
        assert order >= 1.5
    print("PASS: isolated mixed-flux SBP-SAT manufactured benchmark")


if __name__ == "__main__":
    main()
