import numpy as np

from anisitropic_he_solver import solve
from test_functions import *
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
        X, Y, u = solve(a, d, g, f, r, M[i])
        e[i] = np.max(np.abs(u-T(X, Y)))

    loglogplot_error(h, e).show()


def main_fatten_bdry(test_functions):
    a = 1
    r = np.sqrt(2)
    d1 = (1, 0)
    d2 = (1, 1)
    d = d1, d2

    T, f = test_functions(r, a)
    g = boundaries(T)

    M = np.logspace(1, 1.5, 20, dtype=int)
    N = np.array(2*M/r, dtype=int) + 1
    h = 1/M
    e = np.empty_like(h)

    for i in range(len(M)):
        X, Y, u = solve(a, d, g, f, r, M[i], N[i])
        e[i] = np.max(np.abs(u-T(X, Y)))

    loglogplot_error(h, e).show()


if __name__ == '__main__':
    main_fatten_bdry(Test2)
    main_fatten_bdry(Test3)

