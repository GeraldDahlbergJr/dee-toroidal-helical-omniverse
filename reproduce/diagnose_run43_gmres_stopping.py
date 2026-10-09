"""Run43 Newton linear-solve diagnostic (isolated; not used by production workflow).

Tests the hypothesis that GMRES's absolute tolerance can report success while
Run43's independent relative-residual gate rejects the step. Also evaluates
RHS normalization, which leaves A*x=b unchanged but changes GMRES's effective
absolute stopping threshold in original units. No production settings altered.
Run: python reproduce/diagnose_run43_gmres_stopping.py
"""
import json
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import gmres

RTOL = 1e-11
ATOL = 1e-13
GATE = 1e-10

def measure(A, b, normalized):
    scale = float(np.linalg.norm(b))
    if scale == 0:
        return {"normalized": normalized, "info": 0, "true_relative_residual": 0.0}
    rhs = b / scale if normalized else b
    x, info = gmres(A, rhs, rtol=RTOL, atol=ATOL, restart=30, maxiter=20)
    if normalized:
        x = x * scale
    rel = float(np.linalg.norm(A @ x - b) / scale)
    return {"normalized": normalized, "info": int(info),
            "rhs_norm": scale, "true_relative_residual": rel,
            "passes_frozen_gate": bool(info == 0 and rel <= GATE)}

def main():
    # Nontrivial, well-conditioned synthetic example; NOT Run43 checkpoint.
    n = 100
    A = diags([-0.2 * np.ones(n-1), 2.0 * np.ones(n),
               -0.2 * np.ones(n-1)], [-1, 0, 1], format="csr")
    b = 1e-8 * np.sin(np.arange(1, n + 1))
    raw = measure(A, b, False)
    scaled = measure(A, b, True)
    print(json.dumps({"case": "synthetic_only", "raw": raw, "normalized": scaled}, indent=2))
    assert scaled["passes_frozen_gate"], "Normalized diagnostic failed frozen gate"
    if raw["passes_frozen_gate"]:
        print("Raw solve also passed; this case does not demonstrate the mismatch.")
    else:
        print("Raw solve failed frozen gate; compare with normalized solve.")

if __name__ == "__main__":
    main()
