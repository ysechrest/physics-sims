# -*- coding: utf-8 -*-
"""
Created on Sat May 26 07:38:24 2018

@author: randale
"""

import scipy.signal as spsig
import numpy as np
import scipy.integrate as spint

def simulate():
    Nt = 80000
    dt = 1.0e-4
    
    T = np.ones(Nt)*293.0
    r = np.ones(Nt)*1.0e-3
    t = np.arange(Nt)*dt
    
    solver = spint.ode(f)
    solver.set_integrator('dopri5')
    solver.set_initial_value(np.array([r[0],T[0]]),0.0)
    
    ii=1
    while solver.successful() and ii<Nt:
        res = solver.integrate(solver.t+dt)
        r[ii] = res[0]        
        T[ii] = res[1]
        ii+=1
    
    return t,r,T

def simulate2():
    Nt = 39000
    dt = 1.0e-3
    
    T = np.ones(Nt)*293.0
    r = np.ones(Nt)*1.0e-3
    t = np.arange(Nt)*dt
    
    solver = spint.ode(f2)
    solver.set_integrator('dopri5')
    solver.set_initial_value(np.array([r[0],T[0]]),0.0)
    
    ii=1
    while solver.successful() and ii<Nt:
        res = solver.integrate(solver.t+dt)
        r[ii] = res[0]        
        T[ii] = res[1]
        ii+=1
    
    return t,r,T

def f(t,y):
    csb = 5.67e-12
    cp = 0.39
    rho = 7.14
    ems = 0.15
    alpha = 0.87
    SLIE = 9.6e-3
    Tenv = 293.0
    I_euv = I_euv_of_t(t,SLIE)
    
    pel_abs = 0.14
    mask_ref= 0.65
    
    r = y[0]
    T = y[1]
    
    Qeuv = ((1-pel_abs)*I_euv+mask_ref*(1-pel_abs)*I_euv)*alpha
    Qrad = -ems*csb*(T**4.0-Tenv**4.0) 
    Tdot = (0.25*Qeuv + Qrad)/(rho*cp*0.333333*r)
    rdot = -G_of_T(T)/rho 
    
    return np.array([rdot,Tdot])
    
def f2(t,y):
    csb = 5.67e-12
    cp = 0.39
    rho = 7.14
    ems = 0.15
    alpha = 0.87
    SLIE = 9.6e-3
    Tenv = 293.0
    I_euv = (5.0*SLIE/(0.240*0.60))
    
    pel_abs = 0.14
    mask_ref= 0.65
    
    r = y[0]
    T = y[1]
    
    Qeuv = ((1-pel_abs)*I_euv+mask_ref*(1-pel_abs)*I_euv)*alpha
    Qrad = -ems*csb*(T**4.0-Tenv**4.0) 
    Tdot = (0.25*Qeuv + Qrad)/(rho*cp*0.333333*r)
    rdot = -G_of_T(T)/rho 
    
    return np.array([rdot,Tdot])

def I_euv_of_t(t,SLIE):
    slit_width = 0.008    
    v = 0.240
    t_revisit = 0.6
    
    I0 = (SLIE*5.0)/slit_width
    
    return (I0/2.0)*(spsig.square(2*np.pi*(t/t_revisit-0.01),duty=(slit_width/(v*t_revisit)))+1)
    
def G_of_T(T):
    A = 6.102
    B = -6776.0
    P_v = (760.0)*10.0**(A+B/T)
    #G = 1.0e4 #g 1/cm^2 1/yr
    #G = G*(1/365.)*(1/24.)*(1/60.)*(1/60.)
    G = 5.04e3*P_v*(65.4/T)**0.5
    G = G*(1/(24.0*60.0*60.0))
    return G