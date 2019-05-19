import numpy as np
import sympy as sp
from scipy.integrate import odeint
from scipy.optimize import root
from scipy.stats import linregress
%matplotlib inline
import matplotlib.pyplot as plt
import sys


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
#
def interpolator(tt,vv,t):
    return np.interp(t,tt,vv)


class TransistorModel:
    
    # 
    # public interface 
    #

    # constructor
    #
    def __init__(self, Cion, Rion):
        
        self.Cion=7.2e-6
        self.Rion=6.7e4
        
    # single transistor model
    # the tt defines the domain on which 
    # to evaluate the vfunc
    #
    def apply_1T_model(self,Qstart,tt,vfunc):

        v_array=np.array([vfunc(t) for t in tt])
        q_swipe=self.evolve_charge(Qstart,tt,vfunc)

        return (self.j1_1t(q_swipe/self.Cion,v_array),v_array,q_swipe)     
    
    
    #
    # private
    #


    # defines a differential equation for the 
    # evolution of change vs applied voltage over time 
    #
    # returns a difference term for a point in time
    def memdiode(self,q,t,vfunc):
        R = self.Rion
        C = self.Cion

        dqdt = (vfunc(t) - 2*q/C)/R

        return dqdt
    
    
    def evolve_charge(self,Qstart,tt,vfunc):

        # solve and reshape the result to an array from a 1D vector
        q=odeint(self.memdiode,Qstart,tt,args=(vfunc,))
        q=q.reshape(len(tt))

        return q    
    

    # C=Q/V
    #
    # V=Q/C
    #
    # V kg·m2·s−3·A−1   (or J/C)
    #
    # single transistor model 
    #
    def j1_1t(self,Vc,V):

        v1=V-Vc
        
        j = self.js1 * ( np.exp((self.e*v1)/(self.m1*self.kb*self.T),dtype=np.float64) - \
                         np.exp((self.e*(v1-V))/(self.m1*self.kb*self.T),dtype=np.float64) ) - \
                         self.jphoto

        return j  


    ###################################################################
    # constants 

    e   = 1.60217662e-19 # A * s
    kb  = 1.38064852e-23 # J/K
        
    # physical parameters
        
    mss=1.93
    fc=0.7
    m1=mss*(1-fc/2)

    js1 = 6.1e-10        # A/m  -10
    js2 = 6.1e-7        # A/m    used to be -7  , -3 gave us the flat Vn on 0
    jphoto=0.          # our generation photo current
    T=300


