import matplotlib.pyplot as plt
import numpy as np


def central_difference_matrix(M, k):

    N = (M - 1) ** 2
    A = -2 * np.eye(N)

    off_diag = np.ones(N - k)
    off_diag[M-2::M-1] = 0
    A += np.diag(off_diag, k=k) + np.diag(off_diag, k=-k)

    return A


if __name__ == '__main__':
    plt.pcolormesh(central_difference_matrix(5, 5))
    plt.show()
