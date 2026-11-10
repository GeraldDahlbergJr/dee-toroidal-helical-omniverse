"""Exact-rational cancellation-preserving bivariate polynomial bounds.

Feasibility component only: not a Run 43 residual certificate.
"""
from fractions import Fraction
from math import comb


class Poly2:
    """Power coefficients in normalized cell coordinates in [0,1]^2."""
    def __init__(self, terms=None):
        self.terms = {tuple(k): Fraction(v) for k, v in (terms or {}).items() if v}

    @classmethod
    def constant(cls, value):
        return cls({(0, 0): value})

    def __add__(self, other):
        if not isinstance(other, Poly2):
            other = Poly2.constant(other)
        out = dict(self.terms)
        for index, value in other.terms.items():
            out[index] = out.get(index, Fraction(0)) + value
            if not out[index]:
                del out[index]
        return Poly2(out)

    __radd__ = __add__

    def __neg__(self):
        return Poly2({k: -v for k, v in self.terms.items()})

    def __sub__(self, other):
        return self + (-other if isinstance(other, Poly2) else -Fraction(other))

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if not isinstance(other, Poly2):
            other = Poly2.constant(other)
        out = {}
        for (i, j), a in self.terms.items():
            for (k, l), b in other.terms.items():
                idx = (i + k, j + l)
                out[idx] = out.get(idx, Fraction(0)) + a * b
        return Poly2(out)

    __rmul__ = __mul__

    def derivative(self, ds=0, dt=0):
        p = self
        for axis, times in [(0, ds), (1, dt)]:
            for _ in range(times):
                out = {}
                for indices, value in p.terms.items():
                    if indices[axis]:
                        new_idx = list(indices)
                        new_idx[axis] -= 1
                        out[tuple(new_idx)] = value * indices[axis]
                p = Poly2(out)
        return p

    def bound(self):
        """Exact Bernstein convex-hull enclosure on a unit square."""
        if not self.terms:
            return Fraction(0), Fraction(0)
        ns = max(i for i, _ in self.terms)
        nt = max(j for _, j in self.terms)
        coefficients = []
        for a in range(ns + 1):
            for b in range(nt + 1):
                coefficients.append(sum(
                    c * Fraction(comb(a, i), comb(ns, i))
                    * Fraction(comb(b, j), comb(nt, j))
                    for (i, j), c in self.terms.items() if i <= a and j <= b))
        return min(coefficients), max(coefficients)
