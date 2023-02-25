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

def central_difference_matrix_irregular_bndry(end_pts):
    """
    :param endpoints: array containing rightmost interior endpoint index
    :return: Central difference matrix approximating twice derivative along x and y
                using the iregular domain
    """
    def B(n):
        '''Create a tridiag(1, -4, 1) block'''
        off_diag = np.repeat(1, n-1)
        diag = np.repeat(-4, n)
        return np.diag(off_diag, -1) + np.diag(diag, 0) + np.diag(off_diag, 1)
    def I(m, n):
        '''Create an m*n block with 1's along the diagonal and 0 else'''
        A = np.zeros(shape=(m, n)) 
        np.fill_diagonal(A, 1)
        return A
    

    N = end_pts[-1]
    A = np.zeros(shape=(N, N))
    end_pts = np.insert(end_pts, 0, 0) # Prepend 0

    # Iterate "blocks" down A
    for i, n in enumerate(end_pts[:-1]):
        range_i = np.arange(end_pts[i], end_pts[i+1])

        # Left-offset block
        if i > 0:
            A[
                range_i,
                end_pts[i-1]:end_pts[i]
            ] = I(np.diff(end_pts)[i], np.diff(end_pts)[i-1])

        # Right-offset block
        if i < len(end_pts)-2:
            A[
                range_i,
                end_pts[i+1]:end_pts[i+2]
            ] = I(np.diff(end_pts)[i], np.diff(end_pts)[i+1])

        # Center/Diagonal block
        A[
            range_i,
            end_pts[i]:end_pts[i+1]
        ] = B(np.diff(end_pts)[i])

    return A

def test_irregular_matrix():
    # Test implementation
    M = 7
    x, h = np.linspace(0, 1, M+1, retstep=True)

    # Find the number of interior points at each y
    endpoint_at_row = np.zeros(M-1, dtype=np.int64)
    for i, x in enumerate(x[1:-1]):
        n = (np.sqrt(1-x) - np.sqrt(1-x)%h)/h
        endpoint_at_row[i] = endpoint_at_row[i-1]
        endpoint_at_row[i] += n if n*h != np.sqrt(1-x) else n - 1

    # endpoint_at_row = [5, 9, 11]
    A = central_difference_matrix_irregular_bndry(endpoint_at_row)

    A = np.flip(A, axis=0)
    plt.pcolormesh(A)
    plt.show()

if __name__ == '__main__':
    for d in [(1, 0), (0, 1), (1, 1)]:
        plt.subplots()
        problem2_modified_discretization
        plt.pcolormesh(np.flip(central_difference_matrix(5, d), axis=0))
        plt.title(d)
    plt.show()

    test_irregular_matrix()
