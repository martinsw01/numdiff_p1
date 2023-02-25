import numpy as np
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
    pass
