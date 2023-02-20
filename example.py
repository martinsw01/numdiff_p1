from anisitropic_he_solver import solve
from viz import plot_solution, subplots_3d


def T(x, y):
    return x*(1-x)*y*(2-y)


def g0(x):  # Bottom
    return T(x, 0)


def g1(y):  # Left
    return T(0, y)


def g2(y):  # Right
    return T(1, y)


def g3(x):  # Top
    return T(x, 2)


def rhs(r, a):
    return lambda x, y: -(2*r**2*(x**2-x)+2*r*(x*y-x*(-y+2)-y*(-x+1)+(-x+1)*(-y+2))+2*(a+1)*(y**2-2*y))


def main():
    a = 1
    d1 = (1, 0)
    d2 = (1, 1)  # Adjusted after the grid, taking r into account (is actually (1, r))
    d = d1, d2

    g = g0, g1, g2, g3
    M = 50

    f = rhs(2, a)
    X, Y, u = solve(a, d, g, f, M)

    _, (ax1, ax2) = subplots_3d(ncols=2)
    plot_solution(X, Y, u, ax1)
    plot_solution(X, Y, T(X, Y), ax2).show()


if __name__ == '__main__':
    main()
