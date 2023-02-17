import numpy as np


def central_difference_matrix(d, M):
    """
    Write a function that returns the matrix of the central difference
    
    :param d: direction (x, y), tuple of integers
    :param M: horizontal and vertical grid points
    :return:
    """
    dummy_test = 0

    x, y = d

    def coeff(i, j):
        if i == j:
            return - 2
        elif i == j + y * (M - 1) - x:
            return 1
        elif i == j - y * (M - 1) + x:
            return 1
        else:
            return 0

    A = [[coeff(i, j)
          for i in range((M - 1) ** 2)]
         for j in range((M - 1) ** 2)]

    return np.array(A)


def test_along_x_axis():
    M = 10

    def tridiag(c, a, b, M):
        e = np.ones(M)
        return np.diag(a * e) + np.diag(b * e[:-1], 1) + np.diag(c * e[:-1], -1)

    A = tridiag(1, -2, 1, (M - 1) ** 2)

    assert (central_difference_matrix((1, 0), M) == A).all()


if __name__ == '__main__':
    test_along_x_axis()
