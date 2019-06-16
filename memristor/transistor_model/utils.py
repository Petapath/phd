import sys
import numpy as np


# this is equialent of C++ bind 
#
def partial(func, *args, **kwargs):
    def f(*args_rest, **kwargs_rest):
        kw = kwargs.copy()
        kw.update(kwargs_rest)
        return func(*(args + args_rest), **kw) 
    return f


# we use this to convert discrete domains into
# continous functions using the partial defined above
# tt is the time domain array of the function 
# vv is the corresponding range array of the function 
# t  is the point were we are interpolating 
#
def interpolator(tt,vv,t):
    return np.interp(t,tt,vv)


# takes a func(t) defined over a domain period 
# and makes it periodic with regards to t
#
def make_periodic(func,period,t):
    return func(t-(t//period)*period)


#
# 
# 
