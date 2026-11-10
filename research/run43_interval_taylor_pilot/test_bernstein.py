import unittest
from fractions import Fraction
from .bernstein import Poly2


class BernsteinTests(unittest.TestCase):
    def test_cancellation(self):
        s = Poly2({(1, 0): 1})
        t = Poly2({(0, 1): 1})
        a = (s + t) * (s - t)
        self.assertEqual((a - a).bound(), (0, 0))

    def test_binary_rational(self):
        p = Poly2.constant(0.1)
        self.assertEqual(p.bound(), (Fraction(0.1), Fraction(0.1)))

    def test_derivatives(self):
        s = Poly2({(1, 0): 1})
        p = s * s - s
        self.assertEqual(p.derivative(ds=2).bound(), (2, 2))
        # Bernstein hull encloses true minimum -1/4; need not be tight.
        self.assertEqual(p.bound(), (Fraction(-1, 2), 0))
        self.assertEqual(p.derivative(ds=16).bound(), (0, 0))


if __name__ == "__main__":
    unittest.main()
