"""Reproduce GMRES absolute-tolerance false acceptance on a scaled linear system.

Diagnostic only. No production solver, physics, or Run #43 checkpoint is changed.
Run: python reproduce/test_newton_gmres_scale_invariance.py
"""
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import gmres

n = 80
jac = diags([np.full(n-1, -1.), np.full(n, 4.), np.full(n-1, -1.)], [-1, 0, 1], format="csr")
reference = np.sin(np.arange(n) + 0.1)
base_rhs = jac @ reference
for scale in (1e-7, 1e-9, 1e-11):
    rhs = scale * base_rhs
    rhs_norm = np.linalg.norm(rhs)
    old_step, old_info = gmres(jac, rhs, rtol=1e-11, atol=1e-13, restart=30, maxiter=20)
    normalized_step, normalized_info = gmres(jac, rhs / rhs_norm, rtol=1e-11, atol=0., restart=30, maxiter=20)
    new_step = rhs_norm * normalized_step
    old_rel = np.linalg.norm(jac @ old_step - rhs) / rhs_norm
    new_rel = np.linalg.norm(jac @ new_step - rhs) / rhs_norm
    print(f"scale={scale:.0e} old_info={old_info} old_rel={old_rel:.6g} "
          f"normalized_info={normalized_info} normalized_rel={new_rel:.6g}")
    assert normalized_info == 0 and new_rel <= 1e-10
print("PASS: normalized relative-only solve meets the unchanged 1e-10 gate at all tested scales")
