import numpy as np


class TestFunction:
    def __init__(self, r, a):
        self.r = r
        self.a = a

    def __iter__(self):  # Make it possible to destructure into T and f
        return iter([self.T, self.f])

    def T(self, x, y):
        pass

    def f(self, x, y):
        pass


def boundaries(T):
    def g0(x):  # Bottom
        return T(x, 0)

    def g1(y):  # Left
        return T(0, y)

    def g2(y):  # Right
        return T(1, y)

    def g3(x):  # Top
        return T(x, 2)

    return g0, g1, g2, g3


class Test1(TestFunction):
    def T(self, x, y):
        return np.exp(-x-1/2*y)

    def f(self, x, y):
        return -(self.r + self.a+1 + 0.25*self.r**2) * self.T(x, y)


class Test2(TestFunction):
    def T(self, x, y):
        return 0.5 * x * (1-x)

    def f(self, x, y):
        return np.full_like(x, self.a + 1)


class Test3(TestFunction):
    def T(self, x, y):
        return 0.5 * y * (2-y)

    def f(self, x, y):
        return np.full_like(x, self.r**2)

