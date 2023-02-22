import numpy as np

from directional_central_difference import central_difference_matrix


def boundary_conditions(g, a, M, x, y):
    """

    Adjust rhs wrt boundary

    :param g: g=(g0, g1, g2, g3) tuple of functions defined on the boundary
    :param x: array from 0 to 1 of size M+1
    :param y: array from 0 to 2 of size M+2
    :return: adjustments to the rhs for elements close to the boundary
    """
    g0, g1, g2, g3 = g
    rhs = np.zeros((M-1, M-1))

    # Adjust for derivative along d1
    rhs[:, 0] += a*g1(y[1:-1])  # Left
    rhs[:, -1] += a*g2(y[1:-1])  # Right

    # Adjust for derivative along d2
    rhs[0, :] += g0(x[:-2])  # Bottom
    rhs[-1, :] += g3(x[2:])  # Top
    rhs[1:, 0] += g1(y[1:-2])  # Left, skip first entry to avoid overlap with bottom
    rhs[:-1, -1] += g2(y[2:-1])  # Right, skip last entry to avoid overlap with top

    return rhs.reshape((M-1)**2)


def solve(a, d, g, f, M):
    """
    :param d: tuple of vectors (d1, d2) specifying the direction of the heat flows
    :param g: g=(g0, g1, g2, g3) tuple of functions defined on the boundary
    :param f: rhs of the analytic equation
    :return: meshgrid (X, Y) and numeric solution u on the interior
    """
    x, h = np.linspace(0, 1, M+1, retstep=True)
    y, k = np.linspace(0, 2, M+1, retstep=True)

    X, Y = np.meshgrid(x[1:-1], y[1:-1])

    d1, d2 = d
    A = - a*central_difference_matrix(M, d1)/h**2-central_difference_matrix(M, d2)/h**2
    rhs = f(X, Y).reshape((M-1)**2)
    rhs += boundary_conditions(g, a, M, x, y)/h**2

    u = np.linalg.solve(A, rhs).reshape((M-1, M-1))

    return X, Y, u
