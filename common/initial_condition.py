from math import sqrt, exp, pi
from parameters import A, D

def val(X, t):
    x = X[0]+0.5
    y = X[1]
    r = sqrt(x*x + y*y)
    if r<0.25:
        return A*(exp((-r**2)/(4*D*t))/(4*pi*D*t))
    else:
        return 0.0

