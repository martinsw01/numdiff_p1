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


if __name__ == '__main__':
    plt.pcolormesh(central_difference_matrix(10, (5,3)))
    plt.show()
