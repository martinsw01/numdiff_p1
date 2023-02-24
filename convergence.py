import numpy as np

from anisitropic_he_solver import solve
from test_functions import boundaries, Test3
from viz import loglogplot_error


def main():
    a = 1
    r = 2
    d1 = (1, 0)
    d2 = (1, 1)
    d = d1, d2

    T, f = Test3(r, a)
    g = boundaries(T)

    M = np.logspace(1, 2, 10, dtype=int)
    h = 1/M
    e = np.empty_like(h)

    for i in range(len(M)):
        X, Y, u = solve(a, d, g, f, M[i])
        h[i] = 1/M[i]
        e[i] = np.max(np.abs(u-T(X, Y)))

    loglogplot_error(h, e).show()


if __name__ == '__main__':
    main()
