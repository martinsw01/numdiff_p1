import numpy as np
from directional_central_difference import central_difference_matrix_irregular_bndry
from utilities import get_projection_point, pack_interior_nodes, unpack_interior_nodes


def boundry_conditions(g, end_pts, x_grid, y_grid):
    """

    Adjust rhs wrt boundary

    :param g: g=(g0, g1, g2) tuple of functions defined on the boundary
    :param endpoints: array containing rightmost interior endpoint index 
    :param x_grid: array from 0 to 1 of size M+1
    :param y_grid: array from 0 to 1 of size M+1
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
            block += g0(x_grid[1:-1])

        # Left and right boundry
        block[0] += g1(y_grid[i+1])
        block[-1] += g2(x_grid[len(block)+1])

        # Top boundry
        # Iterate over the rightmost points with boundry above
        m = np.diff(end_pts)[i]-np.diff(end_pts)[i+1] if i < len(end_pts)-2 else np.diff(end_pts)[-1] # number of righ.pts without interior nodes above
        for j, xi in enumerate(x_grid[len(block)-m: len(block)]):
            block[-m+j] += g2(get_projection_point((xi, y_grid[i+1]))[0])

        rhs[n:end_pts[i+1]] = block

    return rhs


def get_rhs(f, g, end_pts, x_grid, y_grid):
    '''
    Create the rhs b-vector containing boundry conditions and rhs analytic.
    '''
    N = end_pts[-1]
    rhs = np.zeros(N)

    h = x_grid[1] - x_grid[0]

    rhs += -boundry_conditions(g, end_pts, x_grid, y_grid)/h**2

    X, Y = np.meshgrid(x_grid, y_grid)
    rhs += pack_interior_nodes(f(X, Y), end_pts)

    return rhs


def solve(g, f, M):
    """
    :param g: g=(g0, g1, g2, g3) tuple of functions defined on the boundary
    :param f: rhs of the analytic equation
    :return: meshgrid (X, Y) and numeric solution u on the interior
    """
    x_grid, h = np.linspace(0, 1, M+1, retstep=True)
    y_grid = np.linspace(0, 1, M+1)
    X, Y = np.meshgrid(x_grid, y_grid)

    # TODO: create a separate function perfoming this tast vvvvvv
    # Find the number of interior points at each y
    endpoint_at_row = np.zeros(M-1, dtype=np.int64)
    for i, x in enumerate(x_grid[1:-1]):
        n = (np.sqrt(1-x)-np.sqrt(1-x)%h)/h
        endpoint_at_row[i] = endpoint_at_row[i-1]
        endpoint_at_row[i] += n if n*h != np.sqrt(1-x) else n-1

    A = central_difference_matrix_irregular_bndry(endpoint_at_row)/h**2
    rhs = get_rhs(f, g, endpoint_at_row, x_grid, y_grid)

    u = unpack_interior_nodes(np.linalg.solve(A, rhs), endpoint_at_row)

    return X, Y, u
