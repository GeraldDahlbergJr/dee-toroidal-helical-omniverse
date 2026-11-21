"""Strong continuum momentum of quartic/Fourier reference reconstructions."""
import numpy as np
from numpy.polynomial import Polynomial

# Generator of Cartesian rotation about the cylinder axis.
J = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 0.]])


def rotate_derivative(v):
    return np.einsum("ij,j...->i...", J, v)


def longitude(v, phi):
    x, y, z = v[..., None]
    c, s = np.cos(phi), np.sin(phi)
    return np.array(np.broadcast_arrays(c*x-s*y, s*x+c*y, z))


def conformal_gradient(gradient):
    trace = np.einsum("ii...->...", gradient)
    return gradient + gradient.swapaxes(0, 1) - (2/3)*np.eye(3)[:, :, None, None]*trace


def divergence(w, wr, wrr, wt, wrt, wtt, r, theta, major):
    """Differentiate in a cylindrical orthonormal frame at phi=0."""
    c, s = np.cos(theta), np.sin(theta)
    q = major+r*c
    jw = rotate_derivative(w)
    # Columns are directional derivatives along cylindrical rho, phi and z.
    def gradient(a, b, spin):
        return np.stack((c*a-s*b, spin, s*a+c*b), axis=1)

    L = conformal_gradient(gradient(wr, wt/r, jw/q))
    Lr = conformal_gradient(gradient(
        wrr, wrt/r-wt/r**2, rotate_derivative(wr)/q-c*jw/q**2))
    Lt = conformal_gradient(gradient(
        wrt-wt/r, wtt/r+wr, rotate_derivative(wt)/q+r*s*jw/q**2))
    Lphi = np.einsum("ij,jk...->ik...", J, L) - np.einsum("ij...,jk->ik...", L, J)
    return c*Lr[:, 0]+s*Lr[:, 2]+(-s*Lt[:, 0]+c*Lt[:, 2])/r+Lphi[:, 1]/q


def traces(w, n, indices, r, theta, major, width):
    """Return left and right strong divergence traces at interior element joins."""
    radial, angles = w.shape[1]-1, w.shape[2]
    k = np.fft.fftfreq(angles, d=1/angles)
    first = 1j*k
    first[angles//2] = 0  # real Nyquist cosine has zero first derivative at nodes
    spectrum = np.fft.fft(w, axis=-1)
    wt = np.fft.ifft(spectrum*first, axis=-1).real
    wtt = np.fft.ifft(-spectrum*k**2, axis=-1).real
    angular = slice(None, None, angles//n)
    nodes = indices*(radial//n)
    grid = np.linspace(0, 1, 5)
    basis = [Polynomial.fromroots(np.delete(grid, j)) for j in range(5)]
    basis = [p/p(x) for p, x in zip(basis, grid)]
    scale = radial/(4*width)
    for endpoint in (1, 0):
        element = nodes[:, None]+np.arange(5)-4*endpoint
        d1 = scale*np.array([p.deriv()(endpoint) for p in basis])
        d2 = scale**2*np.array([p.deriv(2)(endpoint) for p in basis])
        def radial_derivative(field, coeff):
            samples = field[:, element, angular]
            # Subtract the endpoint value to reduce cancellation of constants.
            return np.einsum("j,irjt->irt", coeff, samples-field[:, nodes, None, angular])

        yield divergence(
            w[:, nodes, angular], radial_derivative(w, d1), radial_derivative(w, d2),
            wt[:, nodes, angular], radial_derivative(wt, d1), wtt[:, nodes, angular],
            r, theta, major)
