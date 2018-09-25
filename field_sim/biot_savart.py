# -*- coding: utf-8 -*-
"""
Created on Sat Sep 15 09:54:37 2018

@author: randale
"""
import numpy as np
import scipy.constants as spconsts
import matplotlib.pyplot as pplot
import copy

def calc_bfield(dI_vec,rs,ro,mode='vectorize'):
    Ix,Iy,Iz = dI_vec
    xs,ys,zs = rs
    xo,yo,zo = ro
    
    if mode == 'vectorize':
        ns = len(xs)
        no = len(xo)
        Ix = np.repeat(Ix[:,np.newaxis],no,axis=1)
        Iy = np.repeat(Iy[:,np.newaxis],no,axis=1)
        Iz = np.repeat(Iz[:,np.newaxis],no,axis=1)
        xs = np.repeat(xs[:,np.newaxis],no,axis=1)
        ys = np.repeat(ys[:,np.newaxis],no,axis=1)
        zs = np.repeat(zs[:,np.newaxis],no,axis=1) 
        xo = np.repeat(xo[np.newaxis,:],ns,axis=0)
        yo = np.repeat(yo[np.newaxis,:],ns,axis=0)
        zo = np.repeat(zo[np.newaxis,:],ns,axis=0)
    
    xp = xs-xo
    yp = ys-yo
    zp = zs-zo
    
    rp_mag = np.sqrt(xp**2.0+yp**2.0+zp**2.0)
    
    Bx = (spconsts.mu_0/(4.*spconsts.pi))*(Iy*zp - Iz*yp)/(rp_mag)**3.0
    By = (spconsts.mu_0/(4.*spconsts.pi))*(Iz*xp - Ix*zp)/(rp_mag)**3.0
    Bz = (spconsts.mu_0/(4.*spconsts.pi))*(Ix*yp - Iy*xp)/(rp_mag)**3.0

    Bx = np.sum(Bx,axis=0)
    By = np.sum(By,axis=0)
    Bz = np.sum(Bz,axis=0)

    return (Bx,By,Bz)
    
def field_trace(dI,rs,ro_0,nsteps=100,dstep=0.001):
    
    ro = [copy.copy(ro_0)]
    Bx = []
    By = []
    Bz = []
    for ii in np.arange(nsteps):
        r0 = ro[ii]        
        Bx_ii,By_ii,Bz_ii = calc_bfield(dI,rs,r0,mode='vectorize')
        Bx.append(Bx_ii)
        By.append(By_ii)
        Bz.append(Bz_ii)
        
        bmag = np.sqrt(Bx_ii**2.0+By_ii**2.0+Bz_ii**2.0)
        bx_hat = Bx_ii/bmag
        by_hat = By_ii/bmag
        bz_hat = Bz_ii/bmag
        
        rk =  (r0[0] + bx_hat*dstep/2.,
               r0[1] + by_hat*dstep/2.,
               r0[2] + bz_hat*dstep/2.)
               
        ro.append((r0[0] + bx_hat*dstep/6.,
                   r0[1] + by_hat*dstep/6.,
                   r0[2] + bz_hat*dstep/6.))
                
        Bx_ii,By_ii,Bz_ii = calc_bfield(dI,rs,rk,mode='vectorize')
        bmag = np.sqrt(Bx_ii**2.0+By_ii**2.0+Bz_ii**2.0)
        bx_hat = Bx_ii/bmag
        by_hat = By_ii/bmag
        bz_hat = Bz_ii/bmag
        
        rk  = (r0[0] + bx_hat*dstep/2.,
               r0[1] + by_hat*dstep/2.,
               r0[2] + bz_hat*dstep/2.)   
               
        ro[ii+1] = (ro[ii+1][0] + bx_hat*dstep/3.,
                    ro[ii+1][1] + by_hat*dstep/3.,
                    ro[ii+1][2] + bz_hat*dstep/3.)
        
        Bx_ii,By_ii,Bz_ii = calc_bfield(dI,rs,rk,mode='vectorize')
        bmag = np.sqrt(Bx_ii**2.0+By_ii**2.0+Bz_ii**2.0)
        bx_hat = Bx_ii/bmag
        by_hat = By_ii/bmag
        bz_hat = Bz_ii/bmag
        
        rk  = (r0[0] + bx_hat*dstep,
               r0[1] + by_hat*dstep,
               r0[2] + bz_hat*dstep)   
               
        ro[ii+1] = (ro[ii+1][0] + bx_hat*dstep/3.,
                    ro[ii+1][1] + by_hat*dstep/3.,
                    ro[ii+1][2] + bz_hat*dstep/3.)
        
        Bx_ii,By_ii,Bz_ii = calc_bfield(dI,rs,rk,mode='vectorize')
        bmag = np.sqrt(Bx_ii**2.0+By_ii**2.0+Bz_ii**2.0)
        bx_hat = Bx_ii/bmag
        by_hat = By_ii/bmag
        bz_hat = Bz_ii/bmag

        ro[ii+1] = (ro[ii+1][0] + bx_hat*dstep/6.,
                    ro[ii+1][1] + by_hat*dstep/6.,
                    ro[ii+1][2] + bz_hat*dstep/6.)
                   
    return (Bx,By,Bz),ro
        

def circular_coil(r,I,nphi=360):
    phi = (2*np.pi/360.)*np.linspace(0,360,nphi)
    Ix = -I*np.sin(phi)
    Iy = I*np.cos(phi)
    Iz = Ix*0.0
    
    x = r*np.cos(phi)
    y = r*np.sin(phi)
    z = x*0.0
    
    return (Ix,Iy,Iz),(x,y,z)