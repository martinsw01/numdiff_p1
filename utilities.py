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

print(get_projection_point((0.6, 1)))

