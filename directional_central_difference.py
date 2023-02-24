import matplotlib.pyplot as plt
import numpy as np


def central_difference_matrix(M, d):
    """
    :param M: number of grid points along the axes
    :param d: direction
    :return: Central difference matrix approximating twice derivative along d
    """

    x, y = d

    offset = (M-1) * y + x

    N = (M - 1) ** 2
    A = -2 * np.eye(N)

    off_diag = np.ones(N - offset)
    off_diag[M-2::M-1] = 0
    A += np.diag(off_diag, k=offset) + np.diag(off_diag, k=-offset)

    return A

def central_difference_matrix_irregular_bndry(end_pts):
    """
    :param endpoints: array containing rightmost endpoint interior index
    :return: Central difference matrix approximating twice derivative along x and y
    """
    def B(N):
        off_diag = np.repeat(1, N-1)
        diag = np.repeat(-4, N)
        return np.diag(off_diag, -1) + np.diag(diag, 0) + np.diag(off_diag, 1)
    
    N = end_pts[-1]
    A = np.zeros(shape=(N, N))
    for n, i in enumerate(end_pts):
        # Left-offset block
        A[end_pts[i-1]:end_pts[i]-1, end_pts[i-1]:end_pts[i]-1] # dim = N_i * N_(i-1)

    return A

if __name__ == '__main__':
    # for d in [(1, 0), (0, 1), (1, 1)]:
    #     plt.subplots()
    #     plt.pcolormesh(central_difference_matrix(5, d))
    #     plt.title(d)
    # plt.show()
    
    # Test implementation
    M = 8
    x, h = np.linspace(0, 1, M+1, retstep=True)

    # Find the number of interior points at each y
    endpoint_at_row = np.zeros(M-1, dtype=np.int64)
    for i, x in enumerate(x[1:-1]):
        n = (np.sqrt(1-x) - np.sqrt(1-x)%h)/h
        endpoint_at_row[i] = endpoint_at_row[i-1]
        endpoint_at_row[i] += n if n*h != np.sqrt(1-x) else n - 1

    print(endpoint_at_row)

    A = central_difference_matrix_irregular_bndry(endpoint_at_row)    
    plt.pcolormesh(A)
    plt.show()
