import numpy as np
from sympy.solvers import solve
from sympy import Symbol

'''
this file holdes general utility-functions used in the project in genereal
'''

def get_projection_point(p: tuple):
    '''
    :param p: A tuple reåresenting the point from which we want to project onto 1-x^2
    '''
    r = Symbol('r')
    sol = solve(p[0] + r(1 - 2*p[1]) - 2*r**3, r)
    if isinstance(sol[0], complex) or sol[0]<0: # picks the solution that is non-negative and non-complex
        r = sol[1]
    else: 
        r = sol[0]
    return (r, 1 - r**2) # returns the projection-point in 1-x^2


