import numpy as np

from anisitropic_he_solver import solve
from test_functions import Test3, boundaries
from viz import plot_solution, subplots_3d


def main():
    a = 1
    r = 2
    d1 = (1, 0)
    d2 = (1, 1)  # Adjusted after the grid, taking r into account (is actually (1, r))
    d = d1, d2

    M = 50

    T, f = Test3(r, a)
    g = boundaries(T)
    X, Y, u = solve(a, d, g, f, r, M)

    _, (ax1, ax2) = subplots_3d(ncols=2)
    plot_solution(X, Y, u, ax1)
    plot_solution(X, Y, T(X, Y), ax2).show()


def main_fatten_bndry():
    a = 1
    r = np.sqrt(2)
    d1 = (1, 0)
    d2 = (1, 1)  # Adjusted after the grid, taking r into account (is actually (1, r))
    d = d1, d2

    M = 50
    N = int(2*M/r) + 1

    T, f = Test3(r, a)
    g = boundaries(T)
    X, Y, u = solve(a, d, g, f, r, M, N)

    _, (ax1, ax2) = subplots_3d(ncols=2)
    plot_solution(X, Y, u, ax1)
    plot_solution(X, Y, T(X, Y), ax2).show()


if __name__ == '__main__':
    main_fatten_bndry()
