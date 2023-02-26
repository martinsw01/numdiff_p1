import numpy as np

from isotropic_he_solver import solve
from test_functions import Test5, general_boundries
from viz import plot_solution, subplots_3d


def main():
    M = 50

    T, f = Test5()
    g = general_boundries(T)
    X, Y, u = solve(g, f, M)

    _, (ax1, ax2) = subplots_3d(ncols=2)
    plot_solution(X[1:-1, 1:-1], Y[1:-1, 1:-1], u[1:-1, 1:-1], ax1)
    plot_solution(X[1:-1, 1:-1], Y[1:-1, 1:-1], T(X, Y)[1:-1, 1:-1], ax2).show()


if __name__ == '__main__':
    main()
