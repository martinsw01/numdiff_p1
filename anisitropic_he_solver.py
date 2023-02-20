import matplotlib.pyplot as plt
import numpy as np

from directional_central_difference import central_difference_matrix


def boundary_conditions(g, a, M, x, y):
    """

    Add g to rhs

    :param g: g=(g0, g1, g2, g3) tuple of functions defined on the boundary
    :param x: array from 0 to 1 of size M+1
    :param y: array from 0 to 2 of size M+2
    :return:
    """
    g0, g1, g2, g3 = g
    rhs = np.zeros((M-1, M-1))

    rhs[0, 1:-1] = g0(x[1:M-2])  # Bottom
    rhs[:, 0] = a*g1(y[1:M]) + g1(y[:M-1])  # Left
    rhs[:, -1] = a*g2(y[1:M]) + g2(y[2:])  # Right
    rhs[-1, 1:-1] = g3(x[3:-1])  # Top

    # rhs[:M-1] += g0(x[:M-1])
    # rhs[::M-1] += g1(y[:M-1])
    # rhs[M-2::M-1] += g2(y[2:])
    # rhs[(M-1)*(M-2):] += g3(x[2:])

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

    A = - a*central_difference_matrix(M, d[0])/h**2-central_difference_matrix(M, d[1])/h**2
    rhs = f(X, Y).reshape((M-1)**2)
    rhs += boundary_conditions(g, a, M, x, y)/h**2

    u = np.linalg.solve(A, rhs).reshape((M-1, M-1))

    return X, Y, u
