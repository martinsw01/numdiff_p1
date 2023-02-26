import matplotlib.pyplot as plt
import numpy as np


def central_difference_matrix(M, N, d):
    """
    :param M: number of grid points along the x-axis
    :param N: number of grid points along the y-axis
    :param d: direction
    :return: Central difference matrix approximating twice derivative along d
    """

    x, y = d

    offset = (M-1) * y + x

    A = -2 * np.eye((N-1)*(M-1))

    off_diag = np.ones((N-1)*(M-1)-offset)
    if x:
        off_diag[M-2::M-1] = 0
    A += np.diag(off_diag, k=offset) + np.diag(off_diag, k=-offset)

    return A

def central_difference_matrix_irregular_bndry(end_pts, **kwargs):
    """
    Valid kwargs:
        :x_grid: M+1 array linspace from 0 to 1 with M steps
        :y_grid: M+1 array linspace from 0 to 1 with M steps
        :modified_scheme: Bool, create a matrix for the modified disretization-scheme
    :param endpoints: array containing rightmost interior endpoint index
    :return: Central difference matrix approximating twice derivative along x and y
                using the iregular domain
    """
    if "modified_scheme" in kwargs and not ("x_grid" and "y_grid" in kwargs):
        raise ValueError("Please add an x- and y-grid to use modified scheme")
    elif kwargs.get("modified_scheme"):
        mod_scheme = True
        x_grid = kwargs.get("x_grid")
        y_grid = kwargs.get("y_grid")
        h = x_grid[1]-x_grid[0]

    def B(n, lower_diag=None, diag=None):
        '''
        Create a tridiag(lower_diag, diag, 1) block.
        :param n: dimention of block
        :param lower_diag: n-1 array with values for lower diag
        :param diag: n array with values for diag
        '''
        upper_diag = np.repeat(1.0, n-1)
        if lower_diag is None:
            lower_diag = np.copy(upper_diag)
        if diag is None:
            diag = np.repeat(-4.0, n)
        return np.diag(lower_diag, -1) + np.diag(diag, 0) + np.diag(upper_diag, 1)
    def I(m, n, vals):
        '''Create an m*n block with vals (can be 1) along the diagonal and 0 else'''
        A = np.zeros(shape=(m, n)) 
        np.fill_diagonal(A, vals)
        return A
    

    N = end_pts[-1]
    A = np.zeros(shape=(N, N))
    end_pts = np.insert(end_pts, 0, 0) # Prepend 0

    # Iterate "blocks" down A
    for i, _ in enumerate(end_pts[:-1]):
        range_i = np.arange(end_pts[i], end_pts[i+1])
        n = np.diff(end_pts)[i] # number of interior points at row i
        m = n - np.diff(end_pts)[i+1] if i < len(end_pts)-2 else 0 # number of righ.pts without interior nodes above

        # Left-offset block
        if i > 0:
            if mod_scheme:
                diag = np.ones(n)
                
                # Iterate over rightmost point without interior nodes above
                for j, xi in enumerate(x_grid[n-m:n]):
                    diag[-m+j] = innward(eta1(xi, y_grid[i+1], h))
                left_block = I(np.diff(end_pts)[i], np.diff(end_pts)[i-1], diag)
            else:
                left_block = I(np.diff(end_pts)[i], np.diff(end_pts)[i-1], 1)
            A[
                range_i,
                end_pts[i-1]:end_pts[i]
            ] = left_block

        # Right-offset block
        if i < len(end_pts)-2:
            A[
                range_i,
                end_pts[i+1]:end_pts[i+2]
            ] = I(np.diff(end_pts)[i], np.diff(end_pts)[i+1], 1)

        # Center/Diagonal block
        if mod_scheme:
            diag = np.repeat(-4.0, n)
            lower_diag = np.repeat(1.0, n-1)

            # Correct for the right boundry
            eta0_i = eta0(x_grid[np.diff(end_pts)[i]], y_grid[i+1], h)
            diag[-1]        += 2 - at(eta0_i)
            lower_diag[-1]   = innward(eta0_i)

            # Iterate over the rightmost points without interior nodes above, and
            # correct for the upper boundry
            # m = n - np.diff(end_pts)[i+1] if i < len(end_pts)-2 else 0
            for j, xi in enumerate(x_grid[n-m:n]):
                diag[-m+j] += 2 - innward(eta1(xi, y_grid[i+1], h))

            center_block = B(n, lower_diag, diag)
        else:
            center_block = B(n)
        A[
            range_i,
            end_pts[i]:end_pts[i+1]
        ] = center_block

    return A

#
#   Functions for modified disretization scheme
#
def outward(eta: float) -> float:
    return 2 / (eta*(eta + 1))
def innward(eta: float) -> float:
    return 2 / (eta + 1)
def at(eta: float) -> float:
    return 2 / eta

def eta0(x: float, y: float, h: float) -> float:
    '''Finds horizontal eta given a spatial point x, y. And stepsize h'''
    return (np.sqrt(1-y) - x) / h
def eta1(x: float, y: float, h: float) -> float:
    '''Finds vertical eta given a spatial point x, y. And stepsize h'''
    return ((1-x**2) - y) / h

#
#   Tests
#
def test_irregular_matrix(**kwargs):
    '''
    Test implementation.
    kwargs:
        modified_scheme: bool
    '''
    # Test implementation
    M = 7
    x_grid, h = np.linspace(0, 1, M+1, retstep=True)
    y_grid = np.copy(x_grid)
    if kwargs.get("modified_scheme"):
        kwargs["x_grid"] = x_grid
        kwargs["y_grid"] = y_grid

    # Find the number of interior points at each y
    endpoint_at_row = np.zeros(M-1, dtype=np.int64)
    for i, x in enumerate(x_grid[1:-1]):
        n = (np.sqrt(1-x) - np.sqrt(1-x)%h)/h
        endpoint_at_row[i] = endpoint_at_row[i-1]
        endpoint_at_row[i] += n if n*h != np.sqrt(1-x) else n - 1

    A = central_difference_matrix_irregular_bndry(endpoint_at_row, **kwargs)

    A = np.flip(A, axis=0)
    plt.pcolormesh(A)
    plt.show()

if __name__ == '__main__':
    # for d in [(1, 0), (0, 1), (1, 1)]:
    #     plt.subplots()
    #     plt.pcolormesh(np.flip(central_difference_matrix(5, d), axis=0))
    #     plt.title(d)
    # plt.show()

    test_irregular_matrix(modified_scheme=True)

