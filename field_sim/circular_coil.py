# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 09:42:22 2016

@author: randale
"""

import numpy as np
import scipy.constants as spconsts
import matplotlib.pyplot as pplot
import scipy.integrate as spint
import pdb

class Circular_Coil_Element:
    
    def __init__(self,I_coil,r_coil,r_vec_axis):
        self.current = I_coil
        self.r_vec_axis = r_vec_axis
        self.radius = r_coil
        
    def calc_field(self,r_vec):        
        r_vec_prime = r_vec - self.r_vec_axis
        z = r_vec_prime[2]
        rho_prime = np.sqrt(r_vec_prime[0]**2+r_vec_prime[1]**2)
        a = self.radius        
        
        nphi = 100
        dphi = 2*np.pi/float(nphi)
        phis = np.arange(nphi)*dphi
        B_xp = 0.0
        B_yp = 0.0
        B_z = 0.0
        for phi in phis:
            I_phi_x = z*np.cos(phi)/np.sqrt(rho_prime**2.0+a**2.0+z**2.0-2.0*rho_prime*a*np.cos(phi))**3.0
            I_phi_y = z*np.sin(phi)/np.sqrt(rho_prime**2.0+a**2.0+z**2.0-2.0*rho_prime*a*np.cos(phi))**3.0
            I_phi_z = (a-rho_prime*np.cos(phi))/np.sqrt(rho_prime**2.0+a**2.0+z**2.0-2.0*rho_prime*a*np.cos(phi))**3.0
            B_xp += I_phi_x*dphi
            B_yp += I_phi_y*dphi
            B_z += I_phi_z*dphi
            
        theta = np.arctan2(r_vec_prime[1],r_vec_prime[0])
        B_x = B_xp*np.cos(theta)-B_yp*np.sin(theta)
        B_y = B_xp*np.sin(theta)+B_yp*np.cos(theta)

        return (spconsts.mu_0*self.current*self.radius/(4*np.pi))*np.array([B_x,B_y,B_z])
     
    def plot_coil(self):
        r_axis = np.sqrt(self.r_vec_axis[0]**2+self.r_vec_axis[1]**2)
        r_out= self.radius + r_axis
        r_in= r_axis - self.radius
        z = self.r_vec_axis[2]
        pplot.plot([r_out],[z],'s',color='r')
        pplot.plot([r_in],[z],'s',color='b')
        
class Inf_Line_Element:
    def __init__(self,I_coil,r_coil,r_vec_axis):
        self.current = I_coil
        self.r_vec_axis = r_vec_axis
        self.radius = r_coil
        
    def calc_field(self,r_vec):        
        r_vec_prime = r_vec - self.r_vec_axis
        z = r_vec_prime[2]
        rho_prime = np.sqrt(r_vec_prime[0]**2+r_vec_prime[1]**2)
        a = self.radius        
        
        if rho_prime <= a:        
            B_phi = spconsts.mu_0*(self.current*(rho_prime/a)**2)/(2.0*np.pi*rho_prime)
        else:
            B_phi = spconsts.mu_0*self.current/(2.0*np.pi*rho_prime)

        theta = np.arctan2(r_vec_prime[1],r_vec_prime[0])
        B_x = -B_phi*np.sin(theta)
        B_y = B_phi*np.cos(theta)

        return np.array([B_x,B_y,0.0])
     
    def plot_coil(self):
        r_axis = np.sqrt(self.r_vec_axis[0]**2+self.r_vec_axis[1]**2)
        r_out= self.radius + r_axis
        r_in= r_axis - self.radius
        z = self.r_vec_axis[2]
        pplot.plot([r_out],[z],'s',color='r')
        #pplot.plot([r_in],[z],'s',color='b')        
        
def trace_ode(t,r_vec,coils,sign):
    b_field = np.zeros(3)
    for coil in coils:
        b_field += coil.calc_field(r_vec)
    
    mag_b = np.sqrt(np.sum(b_field*b_field))
    b_hat = b_field/mag_b
    return sign*b_hat
        
def field_trace(r_vec_0,dt,length,coils,dir_sign=1.0):
    odesolver = spint.ode(trace_ode)
    odesolver.set_integrator('dopri5',atol=0.01)
    
    odesolver.set_initial_value(r_vec_0,0.0)
    odesolver.set_f_params(coils,dir_sign)
    dl_cum = 0.0
    line = [r_vec_0]
    r_vec_ii_minus = r_vec_0
    while odesolver.successful() and dl_cum < length:
        r_vec_ii = odesolver.integrate(odesolver.t+dt)
        dl_cum += np.sqrt(np.sum((r_vec_ii-r_vec_ii_minus)**2))
        r_vec_ii_minus = r_vec_ii
        line.append(r_vec_ii)
    
    return np.array(line)
    
def field_line_plot(z0,rs,dt,length,coils,dir_sign=1.0):
    
    #pplot.figure()
    for r_ii in rs:
        r_vec_0 = [r_ii,0,z0]
        line = field_trace(r_vec_0,dt,length,coils,dir_sign=dir_sign)     
        pplot.plot(line[:,0],line[:,2])
    