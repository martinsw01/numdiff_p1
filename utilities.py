import numpy as np
import matplotlib.pyplot as plt
from sympy.solvers import solve
from sympy import Symbol

'''
this file holdes general utility-functions used in the project in genereal
'''

def isreal(num):
    '''
    :param num: number which we want to check is real or not
    '''
    if num.imag != 0:
        result = False
    else:
        result = True

def get_projection_point(p: tuple):
    '''
    :param p: A tuple reåresenting the point from which we want to project onto 1-x^2
    '''
    r = Symbol('r')
    sol = solve(p[0] + r*(1 - 2*p[1]) - 2*r**3, r)
    r = sol[0]
    return (r, 1-r**2)

def unpack_interior_nodes(vec, end_pts):
    '''
    Unpacks a 1D vector to a 2D meshgrid
    :param vec: 1d vector contianing values
    :param end_pts: array containint indexes of each rows rightmost interior point
    :return: 2d array of vec's values (0 where vec is not defined)
    '''
    N = end_pts[-1]
    if len(vec) != N:
        raise ValueError(f"Incorrect shape: {len(vec)}. Expected: {N}")

    grid = np.zeros(shape=(len(end_pts)+2, len(end_pts)+2))

    end_pts = np.insert(end_pts, 0, 0) # Prepend 0
    for i, n in enumerate(end_pts[:-1]):
        grid[i+1, 1:np.diff(end_pts)[i]+1] = vec[n:end_pts[i+1]]


    return grid

def pack_interior_nodes(grid, end_pts):
    '''
    Packs a 2d grid into a 1d vector, only keeping end_pts[i] nodes
    at row i. Note that grid[0] corresponds to x = 0 and so on
    :param grid: 2d array containing values on a grid
    :param end_pts: array containint indexes of each rows rightmost interior point
    :return: 1d array of interior nodes on the grid
    '''

    N = end_pts[-1]
    expected_shape = (len(end_pts) + 2, len(end_pts) + 2)
    if grid.shape != expected_shape:
        raise ValueError(f"Incorrect shape: {grid.shape}. Expected: ({expected_shape})")
    
    # We are only interrested in interior points
    grid = grid[1:-1, 1:]

    end_pts = np.insert(end_pts, 0, 0) # Prepend 0
    vec = np.zeros(N)
    print(end_pts)
    for i, row in enumerate(grid):
        vec[end_pts[i]:end_pts[i+1]] = row[:np.diff(end_pts)[i]]

    return vec

def test_pack():
    M = 10
    x, h = np.linspace(0, 1, M+1, retstep=True)
    y = np.copy(x)

    X, Y = np.meshgrid(x, y)
    Z = X**2 + Y**2


    # TODO: create a separate function perfoming this tast vvvvvv
    # Find the number of interior points at each y
    end_pts = np.zeros(M-1, dtype=np.int64)
    for i, x in enumerate(x[1:-1]):
        n = (np.sqrt(1-x) - np.sqrt(1-x)%h)/h
        end_pts [i] = end_pts[i-1]
        end_pts [i] += n if n*h != np.sqrt(1-x) else n - 1


    plt.pcolormesh(Z)

    Z = pack_interior_nodes(Z, end_pts)
    Z = unpack_interior_nodes(Z, end_pts)
    plt.pcolormesh(Z)

    plt.show()


if __name__=="__main__":
    test_pack()

