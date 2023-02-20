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

    M = [int(2**i) for i in np.linspace(3, 6.4, 10)]
    h = np.empty_like(M, dtype=float)
    e = np.empty_like(M, dtype=float)

    for i in range(len(M)):
        X, Y, u = solve(a, d, g, f, M[i])
        h[i] = 1/M[i]
        e[i] = np.max(np.abs(u-T(X, Y)))

    loglogplot_error(h, e).show()
    print(np.polyfit(np.log(h), np.log(e), 1)[0])


if __name__ == '__main__':
    main()
