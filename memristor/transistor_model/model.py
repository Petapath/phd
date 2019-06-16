import numpy as np
import sympy as sp
from scipy.integrate import odeint
from scipy.optimize import root



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
    
    # two transistor model
    # the tt defines the domain on which 
    # to evaluate the vfunc
    #
    def apply_2T_model(self,Qstart,v0,tt,vfunc):

        v_array=np.array([vfunc(t) for t in tt])
        q_swipe=self.evolve_charge(Qstart,tt,vfunc)

        (v,j)=self.solve_for_j(q_swipe,v_array,self.jphoto,self.Cion,v0)

        return (j,v,q_swipe)     
    
    
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


    
    ###################################################################
    #
    # one transistor model
    # 

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
    #
    # two transistor model
    # 

    # C=Q/V
    #
    # V=Q/C
    #
    # V kg·m2·s−3·A−1   (or J/C)
    #
    # two transistor transistor model 
    #
    def j1_2t(self,Vn,v1,V):
                
        j = self.js1 * ( np.exp((self.e*(v1-Vn))/(self.m1*self.kb*self.T),dtype=np.float64) - \
                         np.exp((self.e*(v1-V))/(self.m1*self.kb*self.T),dtype=np.float64) )

        return j

    # C=Q/V
    #
    # V=Q/C
    #
    # V kg·m2·s−3·A−1   (or J/C)
    #
    # two transistor transistor model 
    #
    def j2_2t(self,Vn,v2,V):
                
        j = self.js2 * ( np.exp((self.e*v2)/(self.m1*self.kb*self.T)) - \
                         np.exp((self.e*(v2-Vn))/(self.m1*self.kb*self.T)) )

        return j

    # j2-j1-jphoto=0
    #
    def kirchoff(self,Vn,v1,v2,V,jphoto):
        return self.j2_2t(Vn,v2,V)-self.j1_2t(Vn,v1,V)+jphoto

    # j is equivalent to simply j2
    # first, we need to solve: j2=j1+jphoto for the unknown Vn
    # then substitue Vn into expression for j2
    #
    def solve_for_j(self,q,V,jphoto,C,v0):

        timesteps=q.size
        
        # q array stores charges on the caps (q same accross all caps)
        # V array has input voltages corresponding to charges above 

        v2=q/C
        v1=V-q/C       
        
        v_n=np.array([root(self.kirchoff,x0=v0,args=(v1[i],v2[i],V[i],jphoto)).x 
                           for i in range(timesteps)]) 
        
        v_n=v_n.reshape(timesteps);
        
        j_array=self.j2_2t(v_n,v2,V) 
            
        return (v_n,j_array)



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


