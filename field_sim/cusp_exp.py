# -*- coding: utf-8 -*-
"""
Created on Sun Sep 25 22:31:32 2016

@author: randale
"""
from field_sim import circular_coil
import numpy as np
import matplotlib.pyplot as pplot
import pdb

def define_hexapole_coil_set(side_len,current=1.0):
    coils = [] 
    for ii in np.arange(6):
        theta = (ii*60.0)*np.pi/180.0
        sign = (-1.0)**ii
        r_vec = np.array([side_len*np.cos(theta),side_len*np.sin(theta),0.0])
        c_ii = circular_coil.Inf_Line_Element(sign*current,0.01,r_vec)
        coils.append(c_ii)
    return coils

def define_dipole_coil_set(sep=1.0,current=1.0):
    coils = []
    #coils.append(circular_coil.Inf_Line_Element(current,0.01,np.array([sep,sep/2.0,0.0])))
    #coils.append(circular_coil.Inf_Line_Element(-current,0.01,np.array([-sep,sep/2.0,0.0])))
    #coils.append(circular_coil.Inf_Line_Element(current,0.01,np.array([sep,-sep/2.0,0.0])))
    #coils.append(circular_coil.Inf_Line_Element(-current,0.01,np.array([-sep,-sep/2.0,0.0])))
    coils.append(circular_coil.Inf_Line_Element(current,0.01,np.array([sep,0.0,0.0])))
    coils.append(circular_coil.Inf_Line_Element(-current,0.01,np.array([-sep,0.0,0.0])))
    return coils
    
def define_cusp_exp_coil_set():

    I_mc1 = 1.0
    I_mc2 = 8.0
    I_mc3 = 4.0
    I_hpc = 1.0
    r_mc1 = 0.33
    z_mc1 = r_mc1
    r_mc2 = 2*r_mc1/3.0
    z_mc2 = 2*r_mc1
    r_mc3 = 0.66
    z_mc3 = r_mc1/2.0
    xside_hpc = 5*2*r_mc1*np.tan(30.0*np.pi/180.0)
    #xside_hpc = 1.75*r_mc1*np.sin(30.0*np.pi/180.0)    
    
    coils = []
    coils.append(circular_coil.Circular_Coil_Element(I_mc1,r_mc1,np.array([0.0,0.0,z_mc1])))
    coils.append(circular_coil.Circular_Coil_Element(-I_mc1,r_mc1,np.array([0.0,0.0,-z_mc1])))
    coils.append(circular_coil.Circular_Coil_Element(I_mc2,r_mc2,np.array([0.0,0.0,z_mc2])))
    coils.append(circular_coil.Circular_Coil_Element(-I_mc2,r_mc2,np.array([0.0,0.0,-z_mc2])))
    coils.append(circular_coil.Circular_Coil_Element(I_mc3,r_mc3,np.array([0.0,0.0,z_mc3])))
    coils.append(circular_coil.Circular_Coil_Element(-I_mc3,r_mc3,np.array([0.0,0.0,-z_mc3])))
    
    I_mc4 = 1.0
    r_mc4 = 1.0
    z_mc4 = z_mc1
    I_mc5 = 4.0
    r_mc5 = 1.33
    z_mc5 = z_mc3
    coils.append(circular_coil.Circular_Coil_Element(I_mc4,r_mc4,np.array([0.0,0.0,z_mc4])))
    coils.append(circular_coil.Circular_Coil_Element(-I_mc4,r_mc4,np.array([0.0,0.0,-z_mc4])))
    coils.append(circular_coil.Circular_Coil_Element(I_mc5,r_mc5,np.array([0.0,0.0,z_mc5])))
    coils.append(circular_coil.Circular_Coil_Element(-I_mc5,r_mc5,np.array([0.0,0.0,-z_mc5])))    
    
    #coils = coils + define_hexapole_coil_set(xside_hpc,I_hpc)
    return coils
    
def define_dipole_cusp_coil_set():

   coils = []
   coils = coils + define_dipole_coil_set(sep=1.0, current=-2.0)
   
   
   I_c = 1.0
   r_c = 0.66
   
   c_ii = circular_coil.Circular_Coil_Element(-I_c,r_c,np.array([0.0,1.0,r_c/4.0]))
   c_ii.plot_coil()
   coils.append(c_ii)
   c_ii = circular_coil.Circular_Coil_Element(I_c,r_c,np.array([0.0,1.0,-r_c/4.0]))
   c_ii.plot_coil()   
   coils.append(c_ii)    

   c_ii = circular_coil.Circular_Coil_Element(I_c,r_c,np.array([0.0,-1.0,r_c/4.0]))
   c_ii.plot_coil()
   coils.append(c_ii)
   c_ii = circular_coil.Circular_Coil_Element(-I_c,r_c,np.array([0.0,-1.0,-r_c/4.0]))
   c_ii.plot_coil()   
   coils.append(c_ii)
   
   return coils
   
def cusp_field_trace():

    dt = 0.1   
    nx = 10
    xmax = 1.50
    #xs = (2*np.arange(nx)/np.float(nx-1)-1)*xmax
    xs = np.array([1.50,1.33,1.16,0.75,0.66,0.5,-1.50,-1.33,-1.16,-0.75,-0.66,-0.5])
    z0 = 0.1
    length = 2.0    
    
    p1 = pplot.figure()
    #coils = define_cusp_exp_coil_set()
    coils = define_dipole_cusp_coil_set()
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
 

def field_plot_xz():
    nx = 21
    nz = 21
    x_max = 1.0
    z_max = 1.0
    xs = (2*np.arange(nx)/float(nx-1)-1)*x_max
    zs = (2*np.arange(nz)/float(nz-1)-1)*z_max
    coils = define_cusp_exp_coil_set()
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
    pplot.imshow(np.transpose(mag_B),cmap='YlGnBu',origin='lower',extent=bbox,interpolation='bicubic',vmax=max_B,vmin=0.0)
    #levels = np.array([1,3,10,30,50,75,100])*1.0e-4
    #oc = pplot.contour(np.transpose(mag_B),levels,extent=bbox,origin='lower',cmap='gnuplot2')


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
    #levels = np.array([1,3,10,30,50,75,100])*1.0e-4
    #oc = pplot.contour(np.transpose(mag_B),levels,extent=bbox,origin='lower',cmap='gnuplot2')

def field_plot_y():
    nx = 101
    x_max = 1.0
    xs = (2*np.arange(nx)/float(nx-1)-1)*x_max
    coils = define_cusp_exp_coil_set()
    B_field = np.zeros((nx,3))
    for ii in np.arange(len(xs)):
            x_ii = xs[ii]
            r_vec = [0.0,x_ii,0.0]
            for coil in coils:
                B_field[ii,:] += coil.calc_field(r_vec)

    mag_B = np.sqrt(np.sum(B_field*B_field,axis=1))
    max_B = np.max(mag_B)
    pplot.figure()
    pplot.plot(xs,B_field[:,0])
    pplot.plot(xs,B_field[:,1])
    pplot.plot(xs,B_field[:,2])
    pplot.plot(xs,mag_B)
    
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
   