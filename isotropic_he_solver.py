import numpy as np
from directional_central_difference import central_difference_matrix_irregular_bndry

def boundary_conditions(g, a, end_pts, x_grid, y_grid):
    """

    Adjust rhs wrt boundary

    :param g: g=(g0, g1, g2) tuple of functions defined on the boundary
    :param endpoints: array containing rightmost interior endpoint index 
    :param x_grid: array from 0 to 1 of size M+1
    :param y_grid: array from 0 to 1 of size M+2
    :return: adjustments to the rhs for elements close to the boundary
    """
    g0, g1, g2 = g

    N = end_pts[-1]
    end_pts = np.insert(end_pts, 0, 0)
    rhs = np.zeros(N)

    # Create rhs by iterating over sections/blocks
    for i, n in enumerate(end_pts[:-1]):
        block = np.zeros(np.diff(end_pts)[i])

        # Bottom boundry
        if i == 0:
            block += [g0(x) for x in x_grid[1:-1]]

        # Left and right boundry
        block[0] += g1(y_grid[i+1])
        block[-1] += g2(x_grid[len(block) + 1])

        # Top boundry
        # Iterate over the rightmost points with boundry above
        m = np.diff(end_pts)[i] - np.diff(end_pts)[i+1]
        for i, xi in enumerate(x_grid[len(block) - m: len(block)]):
            block[-m + i] = g2(xi)
        
        rhs[n:end_pts[i+1]] = block

    return rhs


def solve(a, d, g, f, M):
    """
    :param d: tuple of vectors (d1, d2) specifying the direction of the heat flows
    :param g: g=(g0, g1, g2, g3) tuple of functions defined on the boundary
    :param f: rhs of the analytic equation
    :return: meshgrid (X, Y) and numeric solution u on the interior
    """
    x, h = np.linspace(0, 1, M+1, retstep=True)
    y, k = np.copy(x, h)

    X, Y = np.meshgrid(x[1:-1], y[1:-1])

    # Find the number of interior points at each y
    endpoint_at_row = np.zeros(M-1, dtype=np.int64)
    for i, x in enumerate(x[1:-1]):
        n = (np.sqrt(1-x) - np.sqrt(1-x)%h)/h
        endpoint_at_row[i] = endpoint_at_row[i-1]
        endpoint_at_row[i] += n if n*h != np.sqrt(1-x) else n - 1

    A = central_difference_matrix_irregular_bndry(endpoint_at_row)

    # d1, d2 = d
    # A = - a*central_difference_matrix(M, d1)/h**2-central_difference_matrix(M, d2)/h**2
    # rhs = f(X, Y).reshape((M-1)**2)
    # rhs += boundary_conditions(g, a, M, x, y)/h**2

    # u = np.linalg.solve(A, rhs).reshape((M-1, M-1))

    u = A # PLACEHOLDER
    return X, Y, u
    
