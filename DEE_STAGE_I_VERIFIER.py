#!/usr/bin/env python3
"""Two bounded exact product identities for the Stage I cutoff attempt.
No PDE estimate, spectral exclusion, physical experiment, or file write.
"""
import sympy as s


def main():
    print('SymPy', s.__version__, flush=True)
    x = s.symbols('x0:3', real=True)
    omega = s.Symbol('Omega', real=True, nonzero=True)
    fun = lambda name: s.Function(name)(*x)
    a = s.Matrix(3, 3, lambda i, j: fun(f'A{min(i,j)}{max(i,j)}'))
    b = s.Matrix([fun(f'B{i}') for i in range(3)])
    c0, u, chi, phi = map(fun, ('C0', 'u', 'chi', 'phi'))
    grad = lambda v: s.Matrix([s.diff(v, xx) for xx in x])
    div = lambda v: sum(s.diff(v[i], x[i]) for i in range(3))

    def q(v):
        return -div(a * grad(v)) - s.I * omega * (
            2 * b.dot(grad(v)) + div(b) * v) - omega**2 * c0 * v

    count = 0

    def check(name, residual):
        nonlocal count
        if s.expand(residual.doit()) != 0:
            raise AssertionError(name)
        count += 1
        print('PASS:', name, flush=True)

    cutoff = -2 * (a * grad(chi)).dot(grad(u)) - div(a * grad(chi)) * u
    cutoff -= 2 * s.I * omega * b.dot(grad(chi)) * u
    check('full anisotropic cutoff commutator', q(chi * u) - chi * q(u) - cutoff)

    conjugated = q(u) + 2 * (a * grad(phi)).dot(grad(u))
    conjugated += (div(a * grad(phi)) - grad(phi).dot(a * grad(phi))
                   + 2 * s.I * omega * b.dot(grad(phi))) * u
    check('real exponential conjugation retains principal tensor',
          s.exp(phi) * q(s.exp(-phi) * u) - conjugated)
    print(f'PASS: Stage I {count} exact local identities; exterior uniqueness NOT ESTABLISHED.', flush=True)


if __name__ == '__main__':
    main()
