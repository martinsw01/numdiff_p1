import numpy as np

from isotropic_he_solver import solve
from test_functions import Test5, Test6, general_boundries
from utilities import pack_interior_nodes, unpack_interior_nodes
from viz import plot_solution, subplots_3d


def main():
    M = 50

    T, f = Test6()
    g = general_boundries(T)
    X, Y, u = solve(g, f, M, modified_scheme=True)

    x_grid = X[0, :]
    h = x_grid[1] - x_grid[0]

    # TODO: create a separate function perfoming this tast vvvvvv
    # Find the number of interior points at each y
    endpoint_at_row = np.zeros(M-1, dtype=np.int64)
    for i, x in enumerate(x_grid[1:-1]):
        n = (np.sqrt(1-x)-np.sqrt(1-x)%h)/h
        endpoint_at_row[i] = endpoint_at_row[i-1]
        endpoint_at_row[i] += n if n*h != np.sqrt(1-x) else n-1

    _, (ax1, ax2) = subplots_3d(ncols=2)

    exact_sol = unpack_interior_nodes(pack_interior_nodes(T(X, Y), endpoint_at_row), endpoint_at_row)

    plot_solution(X[1:-1, 1:-1], Y[1:-1, 1:-1], u[1:-1, 1:-1], ax1)
    plot_solution(X[1:-1, 1:-1], Y[1:-1, 1:-1], exact_sol[1:-1, 1:-1], ax2).show()


if __name__ == '__main__':
    main()
