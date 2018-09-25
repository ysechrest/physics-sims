# -*- coding: utf-8 -*-
"""
Created on Sat Oct  8 22:17:39 2016

@author: randale
"""

import numpy as np
class particle:
    
    def __init__(self):
        self.mass = 1.0
        self.charge = 1.0
        self.r = np.array([0.0,0.0,0.0])
        self.v = np.array([0.0,0.0,0.0])
        
    def update_state(self,dt):
        F = lorentz_force(self.charge,self.v,e_field,b_field)
        a = F/self.mass
        dv = a*dt/2.0
        dr = v*dt/2.0
        
        r_
        
        
def lorentz_force(t,y,*f_args):
    q = f_args[0]
    m = f_args[1]
    e_field = f_args[2]
    b_field = f_args[3]
    r = [y[0],y[1],y[2]]
    v = [y[3],y[4],y[5]]    
    F = q*e_field + np.cross(v,b_field)
    R = [F[0]/m,F[1]/m,F[2]/m,v[0],v[1],v[2]]
    return R