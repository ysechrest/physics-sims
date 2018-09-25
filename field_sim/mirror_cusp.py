# -*- coding: utf-8 -*-
"""
Created on Sun Sep 25 22:31:32 2016

@author: randale
"""
from field_sim import circular_coil
import numpy as np
import matplotlib.pyplot as pplot
import pdb

def coil_set_02():
    coils = []
    I_c = 1.0
    r_c = 2.0
    #pplot.figure()
    r_i = r_c/2.0
    r_o = r_c
    I_i = I_c
    I_o = I_c*2.0
    if True:
        c_ii = circular_coil.Circular_Coil_Element(-I_i,r_i,np.array([0.0,0.0,r_i/2.0]))
        c_ii.plot_coil()
        coils.append(c_ii)
        c_ii = circular_coil.Circular_Coil_Element(I_o,r_o,np.array([0.0,0.0,r_o/2.0]))
        c_ii.plot_coil()
        coils.append(c_ii)
        c_ii = circular_coil.Circular_Coil_Element(I_i,r_i,np.array([0.0,0.0,-r_i/2.0]))
        c_ii.plot_coil()
        coils.append(c_ii)   
        c_ii = circular_coil.Circular_Coil_Element(-I_o,r_o,np.array([0.0,0.0,-r_o/2.0]))
        c_ii.plot_coil()
        coils.append(c_ii)
 

    
    if True:
        I_m = 7.0*I_c   
        r_m = 3.0*r_c
        c_ii = circular_coil.Circular_Coil_Element(I_m,r_m,np.array([0.0,0.0,2.0*r_m]))
        c_ii.plot_coil()
        coils.append(c_ii)
        c_ii = circular_coil.Circular_Coil_Element(I_m,r_m,np.array([0.0,0.0,-2.0*r_m]))
        c_ii.plot_coil()   
        coils.append(c_ii)
     

    
    return coils
    
def coil_set_01():
    coils = []
    I_c = 1.0
    r_c = 1.0
    #pplot.figure()
    if True:
        c_ii = circular_coil.Circular_Coil_Element(I_c,r_c,np.array([0.0,0.0,-1.5*r_c]))
        c_ii.plot_coil()
        coils.append(c_ii)
        c_ii = circular_coil.Circular_Coil_Element(I_c,r_c,np.array([0.0,0.0,1.5*r_c]))
        c_ii.plot_coil()
        coils.append(c_ii)
        c_ii = circular_coil.Circular_Coil_Element(-0.66*I_c,r_c,np.array([0.0,0.0,0.5*r_c]))
        c_ii.plot_coil()
        coils.append(c_ii)   
        c_ii = circular_coil.Circular_Coil_Element(-0.66*I_c,r_c,np.array([0.0,0.0,-0.5*r_c]))
        c_ii.plot_coil()
        coils.append(c_ii)
        c_ii = circular_coil.Circular_Coil_Element(-0.66*I_c,r_c*0.75,np.array([0.0,0.0,0.0]))
        c_ii.plot_coil()
        coils.append(c_ii)

    
    if True:
        I_m = 5.0*I_c   
        r_m = 3.0*r_c
        c_ii = circular_coil.Circular_Coil_Element(I_m,r_m,np.array([0.0,0.0,2.0*r_m]))
        c_ii.plot_coil()
        coils.append(c_ii)
        c_ii = circular_coil.Circular_Coil_Element(I_m,r_m,np.array([0.0,0.0,-2.0*r_m]))
        c_ii.plot_coil()   
        coils.append(c_ii)
        

    if False:
        I_d = I_c
        r_d = 1.5*r_c
        c_ii = circular_coil.Circular_Coil_Element(2.0*I_d,r_d,np.array([0.0,0.0,0.0]))
        c_ii.plot_coil()
        coils.append(c_ii)
        c_ii = circular_coil.Circular_Coil_Element(-I_d,r_d,np.array([0.0,0.0,1.5*r_c]))
        c_ii.plot_coil()
        coils.append(c_ii)
        c_ii = circular_coil.Circular_Coil_Element(-I_d,r_d,np.array([0.0,0.0,-1.5*r_c]))
        c_ii.plot_coil()
        coils.append(c_ii)
    
    if False:
        I_d = 3.0*I_c
        r_d = 1.5*r_c
        c_ii = circular_coil.Circular_Coil_Element(I_d,r_d,np.array([0.0,0.0,r_c]))
        c_ii.plot_coil()
        coils.append(c_ii)
        c_ii = circular_coil.Circular_Coil_Element(I_d,r_d,np.array([0.0,0.0,-r_c]))
        c_ii.plot_coil()
        coils.append(c_ii)
        c_ii = circular_coil.Circular_Coil_Element(-I_d,r_d,np.array([0.0,0.0,1.5*r_c]))
        c_ii.plot_coil()
        coils.append(c_ii)
        c_ii = circular_coil.Circular_Coil_Element(-I_d,r_d,np.array([0.0,0.0,-1.5*r_c]))
        c_ii.plot_coil()
        coils.append(c_ii)

    
    return coils

def define_dipole_cusp_coil_set():
    coils = []
    I_c = 1.0
    r_c = 1.0
    c_ii = circular_coil.Circular_Coil_Element(-I_c,r_c,np.array([0.0,0.0,2.5*r_c/2.0]))
    c_ii.plot_coil()
    coils.append(c_ii)
    c_ii = circular_coil.Circular_Coil_Element(I_c,r_c,np.array([0.0,0.0,4*r_c/2.0]))
    c_ii.plot_coil()   
    coils.append(c_ii)    
    c_ii = circular_coil.Circular_Coil_Element(-I_c,r_c,np.array([0.0,0.0,-3*r_c/2.0]))
    c_ii.plot_coil()
    coils.append(c_ii)
    c_ii = circular_coil.Circular_Coil_Element(I_c,r_c,np.array([0.0,0.0,-4*r_c/2.0]))
    c_ii.plot_coil()   
    coils.append(c_ii)
    
    c_ii = circular_coil.Circular_Coil_Element(-1.1*I_c,r_c,np.array([0.0,0.0,r_c/2.0]))
    c_ii.plot_coil()
    coils.append(c_ii)
    c_ii = circular_coil.Circular_Coil_Element(-1.1*I_c,r_c,np.array([0.0,0.0,-r_c/2.0]))
    c_ii.plot_coil()
    coils.append(c_ii)
    
    return coils
    
def define_dipole_cusp_coil_set2():
    coils = []
    I_c = 1.0
    r_c = 1.0
    c_ii = circular_coil.Circular_Coil_Element(-I_c,r_c,np.array([0.0,0.0,r_c/2.0]))
    c_ii.plot_coil()
    coils.append(c_ii)
    c_ii = circular_coil.Circular_Coil_Element(I_c,r_c,np.array([0.0,0.0,3*r_c/2.0]))
    c_ii.plot_coil()   
    coils.append(c_ii)    
    c_ii = circular_coil.Circular_Coil_Element(-I_c,r_c,np.array([0.0,0.0,-r_c/2.0]))
    c_ii.plot_coil()
    coils.append(c_ii)
    c_ii = circular_coil.Circular_Coil_Element(I_c,r_c,np.array([0.0,0.0,-3*r_c/2.0]))
    c_ii.plot_coil()   
    coils.append(c_ii)
    
    c_ii = circular_coil.Circular_Coil_Element(1.66*I_c,3.0*r_c,np.array([0.0,0.0,1.5*r_c]))
    c_ii.plot_coil()   
    coils.append(c_ii)
    c_ii = circular_coil.Circular_Coil_Element(1.66*I_c,3.0*r_c,np.array([0.0,0.0,-1.5*r_c]))
    c_ii.plot_coil()
    coils.append(c_ii)
    c_ii = circular_coil.Circular_Coil_Element(1.66*I_c,5.0*r_c,np.array([0.0,0.0,2.5*r_c]))
    c_ii.plot_coil()
    coils.append(c_ii)
    c_ii = circular_coil.Circular_Coil_Element(1.66*I_c,5.0*r_c,np.array([0.0,0.0,-2.5*r_c]))
    c_ii.plot_coil()
    coils.append(c_ii) 
    
    return coils
    
def define_cusp_mirror_coils():
    coils = []
    I_c = 1.0
    r_c = 1.0
    pplot.figure()
    if True:
        c_ii = circular_coil.Circular_Coil_Element(-I_c,r_c,np.array([0.0,0.0,-r_c/2.0]))
        c_ii.plot_coil()
        coils.append(c_ii)
        c_ii = circular_coil.Circular_Coil_Element(-I_c,r_c,np.array([0.0,0.0,r_c/2.0]))
        c_ii.plot_coil()
        coils.append(c_ii)   
        c_ii = circular_coil.Circular_Coil_Element(-I_c,r_c,np.array([0.0,0.0,-3*r_c/2.0]))
        c_ii.plot_coil()
        coils.append(c_ii)
        c_ii = circular_coil.Circular_Coil_Element(-I_c,r_c,np.array([0.0,0.0,3*r_c/2.0]))
        c_ii.plot_coil()
        coils.append(c_ii)   
        c_ii = circular_coil.Circular_Coil_Element(I_c,r_c,np.array([0.0,0.0,5.0*r_c/2.0]))
        c_ii.plot_coil()
        coils.append(c_ii)
        c_ii = circular_coil.Circular_Coil_Element(I_c,r_c,np.array([0.0,0.0,-5.0*r_c/2.0]))
        c_ii.plot_coil()   
        coils.append(c_ii)
    
    if True:
        I_m = 7.0*I_c   
        r_m = 3.0*r_c
        c_ii = circular_coil.Circular_Coil_Element(I_m,r_m,np.array([0.0,0.0,1.5*r_m]))
        c_ii.plot_coil()
        coils.append(c_ii)
        c_ii = circular_coil.Circular_Coil_Element(I_m,r_m,np.array([0.0,0.0,-1.5*r_m]))
        c_ii.plot_coil()   
        coils.append(c_ii)
        
        #c_ii = circular_coil.Circular_Coil_Element(I_m,r_m,np.array([0.0,0.0,1*r_m/3.0]))
        #c_ii.plot_coil()
        #coils.append(c_ii)
        #c_ii = circular_coil.Circular_Coil_Element(I_m,r_m,np.array([0.0,0.0,-1*r_m/3.0]))
        #c_ii.plot_coil()   
        #coils.append(c_ii)    
        #c_ii = circular_coil.Circular_Coil_Element(I_m,r_m,np.array([0.0,0.0,r_m]))
        #c_ii.plot_coil()
        #coils.append(c_ii)
        #c_ii = circular_coil.Circular_Coil_Element(I_m,r_m,np.array([0.0,0.0,-r_m]))
        #c_ii.plot_coil()
        #coils.append(c_ii)
        #c_ii = circular_coil.Circular_Coil_Element(2.5*I_m,r_m,np.array([0.0,0.0,5.0*r_m/3.0]))
        #c_ii.plot_coil()
        #coils.append(c_ii)
        #c_ii = circular_coil.Circular_Coil_Element(2.5*I_m,r_m,np.array([0.0,0.0,-5.0*r_m/3.0]))
        #c_ii.plot_coil()
        #coils.append(c_ii)

    if True:
        I_d = 6.0*I_c
        r_d = 2.0*r_c
        c_ii = circular_coil.Circular_Coil_Element(I_d,r_d,np.array([0.0,0.0,0.0]))
        c_ii.plot_coil()
        coils.append(c_ii)
        c_ii = circular_coil.Circular_Coil_Element(-2*I_d/5.0,r_d,np.array([0.0,0.0,2.0*r_d]))
        c_ii.plot_coil()
        coils.append(c_ii)
        c_ii = circular_coil.Circular_Coil_Element(-2*I_d/5.0,r_d,np.array([0.0,0.0,-2.0*r_d]))
        c_ii.plot_coil()
        coils.append(c_ii)

    #c_ii = circular_coil.Circular_Coil_Element(6*I_c,2.0*r_c,np.array([0.0,0.0,4.0*r_c]))
    #c_ii.plot_coil()   
    #coils.append(c_ii)
    #c_ii = circular_coil.Circular_Coil_Element(6*I_c,2.0*r_c,np.array([0.0,0.0,-4.0*r_c]))
    #c_ii.plot_coil()
    #coils.append(c_ii)
    
    
    #c_ii = circular_coil.Circular_Coil_Element(1.66*I_c,5.0*r_c,np.array([0.0,0.0,2.5*r_c]))
    #c_ii.plot_coil()
    #coils.append(c_ii)
    #c_ii = circular_coil.Circular_Coil_Element(1.66*I_c,5.0*r_c,np.array([0.0,0.0,-2.5*r_c]))
    #c_ii.plot_coil()
    #coils.append(c_ii) 
    
    return coils
    
def define_cusp_mirror_coils2():
    coils = []
    I_c = 1.0
    r_c = 1.0
    c_ii = circular_coil.Circular_Coil_Element(-I_c,2.0*r_c,np.array([0.0,0.0,r_c]))
    c_ii.plot_coil()
    coils.append(c_ii)
    c_ii = circular_coil.Circular_Coil_Element(I_c,2.0*r_c,np.array([0.0,0.0,3*r_c]))
    c_ii.plot_coil()   
    coils.append(c_ii)    
    c_ii = circular_coil.Circular_Coil_Element(-I_c,2.0*r_c,np.array([0.0,0.0,-r_c]))
    c_ii.plot_coil()
    coils.append(c_ii)
    c_ii = circular_coil.Circular_Coil_Element(I_c,2.0*r_c,np.array([0.0,0.0,-3*r_c]))
    c_ii.plot_coil()
    coils.append(c_ii)
 
    #c_ii = circular_coil.Circular_Coil_Element(-I_c,r_c,np.array([0.0,0.0,r_c]))
    #c_ii.plot_coil()
    #coils.append(c_ii)
    #c_ii = circular_coil.Circular_Coil_Element(I_c,r_c,np.array([0.0,0.0,0.0]))
    #c_ii.plot_coil()   
    #coils.append(c_ii)    
    #c_ii = circular_coil.Circular_Coil_Element(-I_c,r_c,np.array([0.0,0.0,-r_c]))
    #c_ii.plot_coil()
    #coils.append(c_ii)

    #c_ii = circular_coil.Circular_Coil_Element(6*I_c,2.0*r_c,np.array([0.0,0.0,4.0*r_c]))
    #c_ii.plot_coil()   
    #coils.append(c_ii)
    #c_ii = circular_coil.Circular_Coil_Element(6*I_c,2.0*r_c,np.array([0.0,0.0,-4.0*r_c]))
    #c_ii.plot_coil()
    #coils.append(c_ii)
    
    
    #c_ii = circular_coil.Circular_Coil_Element(1.66*I_c,5.0*r_c,np.array([0.0,0.0,2.5*r_c]))
    #c_ii.plot_coil()
    #coils.append(c_ii)
    #c_ii = circular_coil.Circular_Coil_Element(1.66*I_c,5.0*r_c,np.array([0.0,0.0,-2.5*r_c]))
    #c_ii.plot_coil()
    #coils.append(c_ii) 
    
    return coils
   
def cusp_field_trace():

    dt = 0.1
    nx = 10
    xmax = 3.0
    #xs = (2*np.arange(nx)/np.float(nx-1)-1)*xmax
    xs = np.array([2.25,1.75,1.25,0.75,0.25,-2.25,-1.75,-1.25,-0.75,-0.25])
    z0 = 1.0
    length = 2.0    
    
    p1 = pplot.figure()
    #coils = define_cusp_exp_coil_set()
    #coils = define_dipole_cusp_coil_set()
    coils = coil_set_02()
    #for coil in coils:
    #    coil.plot_coil()
    
    for x_ii in xs:
        r_vec_0 = np.array([0.1,x_ii,z0])
        line1 = circular_coil.field_trace(r_vec_0,dt,length,coils,dir_sign=-1.0)
        r_vec_0 = np.array([0.1,x_ii,z0])
        line2 = circular_coil.field_trace(r_vec_0,dt,length,coils,dir_sign=1.0)
        p1_ax = p1.gca()        
        p1_ax.plot(line1[:,1],line1[:,2],'k')
        p1_ax.plot(line2[:,1],line2[:,2],'k')
 
def field_plot_xz(B_field=None):
    nx = 31
    nz = 31
    x_max = 4.0
    z_max = 4.0
    xs = (2*np.arange(nx)/float(nx-1)-1)*x_max
    zs = (2*np.arange(nz)/float(nz-1)-1)*z_max
    if B_field is None:
        coils = define_dipole_cusp_coil_set()
        B_field = np.zeros((nx,nz,3))
        for ii in np.arange(len(xs)):
            for jj in np.arange(len(zs)):
                x_ii = xs[ii]
                z_ii = zs[jj]
                r_vec = [x_ii,0.0,z_ii]
                for coil in coils:
                    B_field[ii,jj,:] += coil.calc_field(r_vec)

    mag_B = np.sqrt(np.sum(B_field*B_field,axis=2))
    max_B = np.max(mag_B)
    bbox = [np.min(xs),np.max(xs),np.min(zs),np.max(zs)]
    pplot.figure()
    pplot.imshow(np.transpose(B_field[:,:,0]),cmap='RdBu',origin='lower',extent=bbox,interpolation='bicubic',vmax=max_B,vmin=-1.0*max_B)
    pplot.title('Bx')    
    pplot.figure()
    pplot.title('By')
    pplot.imshow(np.transpose(B_field[:,:,1]),cmap='RdBu',origin='lower',extent=bbox,interpolation='bicubic',vmax=max_B,vmin=-1.0*max_B)
    pplot.figure()
    pplot.title('Bz')
    pplot.imshow(np.transpose(B_field[:,:,2]),cmap='RdBu',origin='lower',extent=bbox,interpolation='bicubic',vmax=max_B,vmin=-1.0*max_B)
    pplot.figure()
    vmax = np.max(np.log(mag_B))
    vmin = np.min(np.log(mag_B))
    print np.min(mag_B)
    print mag_B[10,10]
    pplot.imshow(np.transpose(np.log(mag_B)),cmap='YlGnBu',origin='lower',extent=bbox,interpolation='bicubic',vmax=vmax,vmin=vmin)
    levels = np.array([0,0.001,0.01,0.1,0.33,0.66,1.0])*max_B
    oc = pplot.contour(np.transpose(mag_B),levels,extent=bbox,origin='lower',cmap='gnuplot2')

    return B_field

def field_plot_xy():
    nx = 21
    nz = 21
    x_max = 0.9
    z_max = 2.0
    xs = (2*np.arange(nx)/float(nx-1)-1)*x_max
    zs = (2*np.arange(nz)/float(nz-1)-1)*z_max
    #coils = define_cusp_exp_coil_set()
    coils = define_dipole_cusp_coil_set()
    B_field = np.zeros((nx,nz,3))
    for ii in np.arange(len(xs)):
        for jj in np.arange(len(zs)):
            x_ii = xs[ii]
            z_ii = zs[jj]
            r_vec = [x_ii,z_ii,0.0]
            for coil in coils:
                B_field[ii,jj,:] += coil.calc_field(r_vec)

    mag_B = np.sqrt(np.sum(B_field*B_field,axis=2))
    max_B = np.max(mag_B)
    Br = np.sqrt(B_field[:,:,0]**2+ B_field[:,:,1]**2)
    bbox = [np.min(xs),np.max(xs),np.min(zs),np.max(zs)]
    pplot.figure()
    pplot.imshow(np.transpose(B_field[:,:,0]),cmap='RdBu',origin='lower',extent=bbox,interpolation='bicubic',vmax=max_B,vmin=-1.0*max_B)
    pplot.title('Bx')    
    pplot.figure()
    pplot.title('By')
    pplot.imshow(np.transpose(B_field[:,:,1]),cmap='RdBu',origin='lower',extent=bbox,interpolation='bicubic',vmax=max_B,vmin=-1.0*max_B)
    pplot.figure()
    pplot.title('Bz')
    pplot.imshow(np.transpose(B_field[:,:,2]),cmap='RdBu',origin='lower',extent=bbox,interpolation='bicubic',vmax=max_B,vmin=-1.0*max_B)
    pplot.figure()
    pplot.title('Br')
    pplot.imshow(np.transpose(Br),cmap='RdBu',origin='lower',extent=bbox,interpolation='bicubic',vmax=max_B,vmin=-1.0*max_B)
    pplot.figure()
    pplot.imshow(np.transpose(mag_B),cmap='RdBu',origin='lower',extent=bbox,interpolation='bicubic',vmax=max_B,vmin=-1.0*max_B)
    pplot.figure()    
    pplot.plot(mag_B[nx/2,:])    
    levels = np.array([0,0.1,0.33,0.66,1.0])*max_B
    oc = pplot.contour(np.transpose(mag_B),levels,extent=bbox,origin='lower',cmap='gnuplot2')

def field_plot_y(y=0.0,coils=None):
    nx = 101
    x_max = 4.0
    xs = (2*np.arange(nx)/float(nx-1)-1)*x_max
    #coils = define_cusp_exp_coil_set()
    #coils = define_dipole_cusp_coil_set()
    #coils = define_cusp_mirror_coils()
    if coils is None:
        coils = coil_set_01()
    B_field = np.zeros((nx,3))
    for ii in np.arange(len(xs)):
            x_ii = xs[ii]
            r_vec = [0.0,y,x_ii]
            for coil in coils:
                B_field[ii,:] += coil.calc_field(r_vec)

    mag_B = np.sqrt(np.sum(B_field*B_field,axis=1))
    max_B = np.max(mag_B)
    pplot.figure()
    pplot.plot(xs,B_field[:,0])
    pplot.plot(xs,B_field[:,1])
    pplot.plot(xs,B_field[:,2])
    pplot.plot(xs,mag_B)
    return (xs,B_field)
    
def field_plot_r(phi,coils):
    nx = 101
    x_max = 3.0
    rs = (2*np.arange(nx)/float(nx-1)-1)*x_max
    xs = rs*np.cos(phi*np.pi/180.0)
    ys = rs*np.sin(phi*np.pi/180.0)
    #coils = define_cusp_exp_coil_set()
    B_field = np.zeros((nx,3))
    for ii in np.arange(len(xs)):
            x_ii = xs[ii]
            y_ii = ys[ii]
            r_vec = [x_ii,y_ii,0.0]
            for coil in coils:
                B_field[ii,:] += coil.calc_field(r_vec)

    mag_B = np.sqrt(np.sum(B_field*B_field,axis=1))
    max_B = np.max(mag_B)
    pplot.figure()
    pplot.plot(rs,B_field[:,0])
    pplot.plot(rs,B_field[:,1])
    pplot.plot(rs,B_field[:,2])
    pplot.plot(rs,mag_B)
    
def gradb_dot_b_plot():
    nx = 21
    ny = 21
    r_max = 0.6
    dx = 2*r_max/np.float(nx)
    dy = 2*r_max/np.float(ny)
    dz = 2*r_max/np.float(nx) 
    xs = (2*np.arange(nx)/np.float(nx-1)-1)*r_max
    ys = (2*np.arange(ny)/np.float(ny-1)-1)*r_max
    coils = define_cusp_exp_coil_set()
    B_field = np.zeros((nx,ny,2,3))
    for ii in np.arange(len(xs)):
        for jj in np.arange(len(ys)):
            for kk in np.arange(2):
                x_ii = xs[ii]
                y_ii = ys[jj]
                z_ii = dz*kk + dz/2.0
                r_vec = [x_ii,y_ii,z_ii]
                for coil in coils:
                    B_field[ii,jj,kk,:] += coil.calc_field(r_vec)

    mag_B = np.sqrt(np.sum(B_field*B_field,axis=3))
    b_hat = B_field/np.repeat(mag_B[:,:,:,np.newaxis],3,axis=3)
    dBdx = np.mean((np.roll(mag_B,-1,axis=0) - np.roll(mag_B,1,axis=0))/(2*dx),axis=2)
    dBdy = np.mean((np.roll(mag_B,-1,axis=1) - np.roll(mag_B,1,axis=1))/(2*dy),axis=2)
    dBdz = (mag_B[:,:,0] - mag_B[:,:,1])/dz

    gradbdotb = dBdx*(b_hat[ii,jj,0,0] + b_hat[ii,jj,1,0])/2.0 + \
        dBdy*(b_hat[ii,jj,0,1] + b_hat[ii,jj,1,1])/2.0 + \
        dBdz*(b_hat[ii,jj,0,2] + b_hat[ii,jj,1,2])/2.0
    #pdb.set_trace()    
    
    max_val = np.max(np.abs(gradbdotb))
    bbox = [np.min(xs),np.max(xs),np.min(ys),np.max(ys)]
    pplot.figure()
    pplot.imshow(gradbdotb,cmap='RdBu',origin='lower',extent=bbox,interpolation='bicubic',vmax=max_val,vmin=-1.0*max_val)

    max_val = np.max(np.abs(dBdx))    
    bbox = [np.min(xs),np.max(xs),np.min(ys),np.max(ys)]
    pplot.figure()
    pplot.imshow(dBdx,cmap='RdBu',origin='lower',extent=bbox,interpolation='bicubic',vmax=max_val,vmin=-1.0*max_val)

    max_val = np.max(np.abs(dBdy))    
    bbox = [np.min(xs),np.max(xs),np.min(ys),np.max(ys)]
    pplot.figure()
    pplot.imshow(dBdy,cmap='RdBu',origin='lower',extent=bbox,interpolation='bicubic',vmax=max_val,vmin=-1.0*max_val)

    max_val = np.max(np.abs(dBdz))    
    bbox = [np.min(xs),np.max(xs),np.min(ys),np.max(ys)]
    pplot.figure()
    pplot.imshow(dBdz,cmap='RdBu',origin='lower',extent=bbox,interpolation='bicubic',vmax=max_val,vmin=-1.0*max_val)
   