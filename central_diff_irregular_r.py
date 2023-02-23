import matplotlib.pyplot as plt
import numpy as np


def central_difference_matrix(M, N, eta, d):
    """
    Create central difference matrix for system with irregular grid, s.t. r does not divide 2M.
    Using the method of undetermined coefficients, we get different coefficients for the last M-1 coefficients.

    :param M: number of grid points along the x-axis
    :param N: number of grid points along the y-axis
    :param eta: 0 < eta ≤ 1, s.t. (eta+N)*k=2, where k=rh is the step size along the y-axis
    :param d: direction
    :return: Central difference matrix approximating twice derivative along d
    """

    x, y = d

    offset = (M-1) * y + x

    A = -2 * np.eye((M - 1) * (N-1))

    A[-1, -1] = - 2/eta

    upper_off_diag = np.ones((M-1)*(N-1) - offset)
    lower_off_diag = np.ones((M-1)*(N-1) -offset)

    upper_off_diag[1-M:] = 2/(eta+1)
    if x:
        upper_off_diag[M-2::M-1] = 0
        lower_off_diag[M-2::M-1] = 0

    A += np.diag(upper_off_diag, k=offset) + np.diag(lower_off_diag, k=-offset)

    return A


if __name__ == '__main__':
    for d in [(1, 1)]:
        plt.subplots()
        plt.pcolormesh(central_difference_matrix(7, 5, 0.5, d))
        plt.title(d)
    plt.show()
