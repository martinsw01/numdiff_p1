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


if __name__ == '__main__':
    for d in [(1, 0), (0, 1), (1, 1)]:
        plt.subplots()
        plt.pcolormesh(central_difference_matrix(5, 5, d))
        plt.title(d)
    plt.show()
