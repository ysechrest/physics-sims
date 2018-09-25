# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 12:36:00 2016

@author: randale
"""

from field_sim import circular_coil
import scipy.constants as spconsts
import numpy as np
import matplotlib.pyplot as pplot
import pdb
def test_coil_element():
    c1 = circular_coil.Circular_Coil_Element(1.0,.1,np.array([0,0,0.05]))
    c2 = circular_coil.Circular_Coil_Element(1.0,.1,np.array([0,0,-0.05]))
    x = np.arange(11)/100.0
    z = np.arange(21)/100.0-0.1

    B_field = np.zeros((11,21,3))
    for ii in np.arange(len(x)):
        for jj in np.arange(len(z)):
            B_field[ii,jj,:] = c1.calc_field(np.array([x[ii],0,z[jj]]))
            B_field[ii,jj,:] += c2.calc_field(np.array([x[ii],0,z[jj]]))
    mag_B = np.sqrt(np.sum(B_field*B_field,axis=2))
    pplot.imshow(mag_B,cmap='YlOrRd')
    pplot.figure()
    pplot.plot(z,mag_B[0,:])
    pdb.set_trace()
    
def define_coil_group(current,r,z,w,h,nr,nz,make_plot=False):
    dr = w/float(nr)
    dz = h/float(nz)
    r_min = r-w/2.0
    z_min = z-h/2.0

    coils = []
    for ii in np.arange(nr):
        for jj in np.arange(nz):
            r_ii = r_min+(ii+0.5)*dr
            z_jj = z_min+(jj+0.5)*dz
            c_ij = circular_coil.Circular_Coil_Element(current,r_ii,np.array([0.0,0.0,z_jj]))
            coils.append(c_ij)
            if make_plot:
                pplot.plot([r_ii],[z_jj],'s',color='#b87333')
                pplot.plot([-r_ii],[z_jj],'s',color='#b87333')
            
    return coils

def bipolar_mirror(I,r,I2,r2,dz,w,h,nr,nz):
    c1 = define_coil_group(I,r,-dz,w,h,nr,nz)
    c2 = define_coil_group(I2,r2,0,w,h,nr,nz)
    c3 = define_coil_group(I,r,dz,w,h,nr,nz)
    
    #c1 = define_coil_group(I,r,-dz,w,h,nr,nz)
    #c2 = define_coil_group(I2,r2,dz,w,h,nr,nz)    
    
    coils = []
    coils.extend(c1)
    coils.extend(c2)
    coils.extend(c3)
    
    return coils
    
def helmholtz_model(I,r,z):
    coils = []
    coils.append(circular_coil.Circular_Coil_Element(I,r,np.array([0.0,0.0,z])))
    coils.append(circular_coil.Circular_Coil_Element(I,r,np.array([0.0,0.0,-z])))
    
    nsamps = 128
    r_max = 4*r
    rs = np.arange(nsamps)/float(nsamps-1)*r_max
    bx = np.zeros(nsamps)
    by = np.zeros(nsamps)
    bz = np.zeros(nsamps)
    for ii in np.arange(nsamps):
        for coil in coils:
            b_field = coil.calc_field(np.array([rs[ii],0.0,0.0]))
            bx[ii] += b_field[0]
            by[ii] += b_field[1]
            bz[ii] += b_field[2]
    
    pplot.figure()        
    #pplot.plot(rs,bx)
    #pplot.plot(rs,by)
    pplot.plot(rs,bz)
    
    nsamps = 128
    r_max = 4*r
    rs = np.arange(nsamps)/float(nsamps-1)*r_max
    bx = np.zeros(nsamps)
    by = np.zeros(nsamps)
    bz = np.zeros(nsamps)
    for ii in np.arange(nsamps):
        for coil in coils:
            b_field = coil.calc_field(np.array([0.0,0.0,rs[ii]]))
            bx[ii] += b_field[0]
            by[ii] += b_field[1]
            bz[ii] += b_field[2]
       
    #pplot.plot(rs,bx)
    #pplot.plot(rs,by)
    pplot.plot(rs,bz)

def solenoid_model(current=0.5,r=0.127,z=0.0,w=0.001,h=0.100,nr=1.0,nz=100.0):
    pplot.figure()    
    coils = []
    coils = coils + define_coil_group(current,r,z,w,h,nr,nz,make_plot=True)
    
    nsamps = 128
    z_max = .305
    b_field = np.zeros((nsamps,3))
    for ii in np.arange(nsamps):
        r_vec_ii = np.array([0.0,0.0,ii*z_max/(1.0*nsamps-1.0)])
        for coil in coils:
            b_field[ii,:] += coil.calc_field(r_vec_ii)
            
    pplot.figure()
    pplot.plot(z_max/(1.0*nsamps-1.0)*np.arange(nsamps),b_field[:,2])
            
def pf5_model(current):
    coils = []
    coils = coils + define_coil_group(current,2.0118,0.6489,0.1359,0.0685,6,2)
    coils = coils + define_coil_group(current,2.0118,0.5751,0.1359,0.0685,6,2)
    coils = coils + define_coil_group(current,2.0118,-0.6489,0.1359,0.0685,6,2)
    coils = coils + define_coil_group(current,2.0118,-0.5751,0.1359,0.0685,6,2)
    return coils    

def pf5_field_vs_radius():
    nsamps = 128
    r_max = 4.0
    rs = np.arange(nsamps)/float(nsamps-1)*r_max
    bx = np.zeros(nsamps)
    by = np.zeros(nsamps)
    bz = np.zeros(nsamps)
    coils = pf5_model(1000.0)
    for ii in np.arange(nsamps):
        for coil in coils:
            b_field = coil.calc_field(np.array([rs[ii],0.0,0.0]))
            bx[ii] += b_field[0]
            by[ii] += b_field[1]
            bz[ii] += b_field[2]
            
    pplot.plot(rs,bx)
    pplot.plot(rs,by)
    pplot.plot(rs,bz)
    pplot.plot(rs,np.sqrt(bx**2+by**2+bz**2))
    pplot.plot((2.0118+1.118+.07)*np.array([1,1]),1.1*np.max(np.abs(bz))*np.array([-1,1]),'k--')
    pplot.xlabel('radius [m]')
    pplot.ylabel('field [Tesla]')
    
def pf5_field():
    nr = 41
    nz = 13
    r_max = 4.0
    z_max = 1.2
    rs = (2*np.arange(nr)/float(nr-1)-1.0)*r_max
    zs = (2*np.arange(nz)/float(nz-1)-1.0)*z_max
    mag_B = np.zeros((nr,nz))
    coils = pf5_model(1000.0)
    for ii in np.arange(nr):
        for jj in np.arange(nz):
            for coil in coils:
                b_field = coil.calc_field(np.array([rs[ii],0.0,zs[jj]]))
                mag_B[ii,jj] += np.sqrt(np.sum(b_field*b_field))
            
    bbox = [np.min(rs),np.max(rs),np.min(zs),np.max(zs)]
    pplot.imshow(np.transpose(np.log(mag_B)),cmap='YlOrRd',origin='lower',extent=bbox,interpolation='bicubic')
    levels = np.array([20,40,60,80,100,120,140,160,180,200,220,240,260,280,300])*1.0e-4
    oc = pplot.contour(np.transpose(mag_B),levels,extent=bbox,origin='lower',cmap='terrain')
    for coil in coils:
       coil.plot_coil()

def define_hexapole_coil_set(side_len,current=1.0):
    coils = []
    #for ii in np.arange(6):
    #    theta = (ii*60+30)*np.pi/180.0
    #    #sign = (-1.0)**ii
    #    sign = 1.0
    #    r_vec = np.array([side_len*np.cos(theta),side_len*np.sin(theta),0.0])
    #    c_ii = circular_coil.Inf_Line_Element(sign*current,0.01,r_vec)
    #    coils.append(c_ii)

    for ii in np.arange(4):
        theta = (ii*90+45)*np.pi/180.0
        sign = (-1.0)**ii
        #sign = 1.0
        r_vec = np.array([side_len*np.cos(theta),side_len*np.sin(theta),0.0])
        c_ii = circular_coil.Inf_Line_Element(sign*current,0.01,r_vec)
        coils.append(c_ii)        
        
    return coils
     
def hexapole_field():
    nx = 41
    ny = 41
    x_max = 2.0
    y_max = 2.0
    xs = (2*np.arange(nx)/float(nx-1)-1)*x_max
    ys = (2*np.arange(ny)/float(ny-1)-1)*y_max
    #coils = define_hexapole_coil_set(0.5)
    coils = define_hexapole_coil_set(0.75,current=1.0e4)
    #coils = coils + define_hexapole_coil_set(0.5,current=-0.5e4)
    B_field = np.zeros((nx,ny,3))
    for ii in np.arange(len(xs)):
        for jj in np.arange(len(ys)):
            x_ii = xs[ii]
            y_ii = ys[jj]
            r_vec = [x_ii,y_ii,0.0]
            for coil in coils:
                B_field[ii,jj,:] += coil.calc_field(r_vec)

    mag_B = np.sqrt(np.sum(B_field*B_field,axis=2))
    max_B = np.max(mag_B)
    bbox = [np.min(xs),np.max(xs),np.min(ys),np.max(ys)]
    pplot.figure()
    pplot.imshow(np.transpose(B_field[:,:,0]),cmap='RdBu',origin='lower',extent=bbox,interpolation='bicubic',vmax=max_B,vmin=-1.0*max_B)
    pplot.figure()
    pplot.imshow(np.transpose(B_field[:,:,1]),cmap='RdBu',origin='lower',extent=bbox,interpolation='bicubic',vmax=max_B,vmin=-1.0*max_B)
    pplot.figure()
    pplot.imshow(np.transpose(mag_B),cmap='YlGnBu',origin='lower',extent=bbox,interpolation='bicubic',vmax=max_B,vmin=0.0)
    levels = np.array([1,3,10,30,50,75,100])*1.0e-4
    oc = pplot.contour(np.transpose(mag_B),levels,extent=bbox,origin='lower',cmap='gnuplot2')

    #max_B = np.max(np.sum((B_field*B_field)[17:26,17:26,:],axis=2))
    max_B = np.log(max_B)    
    xx = np.repeat(xs[:,np.newaxis],len(ys),axis=1)
    yy = np.transpose(np.repeat(ys[:,np.newaxis],len(xs),axis=1))
    rs = np.sqrt(xx**2.0+yy**2.0)
    thetas = np.arctan2(yy,xx)
    Br = B_field[:,:,0]*np.cos(thetas)+B_field[:,:,1]*np.sin(thetas)
    Btheta = -B_field[:,:,0]*np.sin(thetas)+B_field[:,:,1]*np.cos(thetas)
    pplot.figure()
    maxbr = np.max(np.abs(Br))
    pplot.imshow(np.transpose(Br),cmap='RdBu',origin='lower',extent=bbox,interpolation='bicubic',vmin=-1.0*maxbr,vmax=1.0*maxbr)
    pplot.colorbar()
    pplot.figure()
    maxbr = np.max(np.abs(Btheta))
    pplot.imshow(np.transpose(Btheta),cmap='RdBu',origin='lower',extent=bbox,interpolation='bicubic',vmin=-1.0*maxbr,vmax=1.0*maxbr)
    pplot.colorbar()
    levels = np.array([0,10,100])*1.0e-4
    oc = pplot.contour(np.transpose(Btheta),levels,extent=bbox,origin='lower',cmap='gnuplot2')


def hexapole_field_vs_r(theta=0.0,fig_id=None,col_idx=None):
    import matplotlib.cm as cm
    nx = 61
    x_max = 1.0
    rs = (np.arange(nx)/float(nx-1))*x_max
    xs = rs*np.cos(theta*np.pi/180.0)
    ys = rs*np.sin(theta*np.pi/180.0)
    coils = define_hexapole_coil_set(0.75,current=1.0e4)
    coils = coils + define_hexapole_coil_set(0.5,current=-0.5e4)
    B_field = np.zeros((nx,3))
    for ii in np.arange(len(xs)):
            x_ii = xs[ii]
            y_ii = ys[ii]
            r_vec = np.array([x_ii,y_ii,0.0])
            for coil in coils:
                B_field[ii,:] += coil.calc_field(r_vec)
                
    mag_B = np.sqrt(np.sum(B_field*B_field,axis=1))
    #pplot.figure(fig_id)
    #pplot.plot(xs,B_field[:,0])
    #pplot.plot(xs,B_field[:,1])
    #pplot.plot(xs,mag_B)
    
    Br = B_field[:,0]*np.cos(theta*np.pi/180.)+B_field[:,1]*np.sin(theta*np.pi/180.)
    Btheta = -B_field[:,0]*np.sin(theta*np.pi/180.)+B_field[:,1]*np.cos(theta*np.pi/180.)
       
    pplot.figure(fig_id)   
    if not(col_idx is None):   
        cmap = cm.get_cmap(name='spectral')
        color= cmap(col_idx)
        pplot.plot(rs,Br,color=color)
        pplot.plot(rs,Btheta,'--',color=color)
    else:
        pplot.plot(rs,Br)
        pplot.plot(rs,Btheta,'--')
    
def coaxial_coilset():
    I = 1.0
    r = 1.0
    coils = []
    coils.append(circular_coil.Circular_Coil_Element(I,r,np.array([0.0,0.0,0.1])))
    coils.append(circular_coil.Circular_Coil_Element(I,r,np.array([0.0,0.0,-0.1])))
    coils.append(circular_coil.Circular_Coil_Element(-2*I,2*r,np.array([0.0,0.0,0.1])))
    coils.append(circular_coil.Circular_Coil_Element(-2*I,2*r,np.array([0.0,0.0,-0.1])))
  
    nx = 101
    x_max = 3.0
    xs = (2*np.arange(nx)/float(nx-1)-1)*x_max
    B_field = np.zeros((nx,3))
    for ii in np.arange(len(xs)):
            x_ii = xs[ii]
            r_vec = [x_ii,0.0,0.0]
            for coil in coils:
                B_field[ii,:] += coil.calc_field(r_vec)

    mag_B = np.sqrt(np.sum(B_field*B_field,axis=1))
    #pplot.plot(xs,B_field[:,0])
    #pplot.plot(xs,B_field[:,1])
    pplot.plot(xs,B_field[:,2])
    
    if True:
        nx = 31
        nz = 31
        x_max = 3.0
        z_max = 3.0
        xs = (2*np.arange(nx)/float(nx-1)-1)*x_max
        zs = (2*np.arange(nz)/float(nz-1)-1)*z_max
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
        pplot.imshow(np.transpose(mag_B),cmap='YlGnBu',origin='lower',extent=bbox,interpolation='bicubic',vmax=max_B,vmin=0.0)
        #levels = np.array([1,3,10,30,50,75,100])*1.0e-4
        #oc = pplot.contour(np.transpose(mag_B),levels,extent=bbox,origin='lower',cmap='gnuplot2')
        
def plot_B_on_axis(coils):
    #I_coil = 1.0
    #r_coil = 0.09525/2.0
    #r_vec_axis = np.array([0,0,0])
    #coil = circular_coil.Circular_Coil_Element(I_coil,r_coil,r_vec_axis)
    #z_max = 0.2667
    z_max = 0.160
    n = 101
    z = z_max*(np.arange(n)/np.float(n-1)-1)
    Bz = []
    for z_ii in z:
        B_field = 0.0
        for coil in coils:
            B_field += coil.calc_field(np.array([0.0,0.0,z_ii]))
        Bz.append(B_field[2])
    
    Bz = np.array(Bz)
    pplot.plot(z,1.0e4*Bz)
    
def small_mirror_test(): 
    coils = bipolar_mirror(50000,0.1,-5000,0.1,0.2,0.02,0.02,2,2)
    pplot.figure()
    for coil in coils:
        coil.plot_coil()
    
    circular_coil.field_line_plot(0,np.arange(0.01,0.2,0.01),0.005,0.1,coils,dir_sign=-1.0)
    circular_coil.field_line_plot(0,np.arange(0.01,0.2,0.01),0.005,0.1,coils,dir_sign=1.0)
    pplot.figure()
    plot_B_on_axis(coils)

def plot_B_vs_r(coils,rmax,dr,z0):
    rs = np.arange(0,rmax+dr,dr)
    Bz = []
    Br = []
    Bphi = []
    for r_ii in rs:
        r_vec = np.array([r_ii,0.0,z0])
        B_ii = np.array([0.0,0.0,0.0])
        for coil in coils:
            B_ii += coil.calc_field(r_vec)
        
        Bz.append(B_ii[2])
        Bphi.append(B_ii[1])
        Br.append(B_ii[0])
    
    Bz = np.array(Bz)
    Br = np.array(Br)
    Bphi = np.array(Bphi)    

    pdb.set_trace()    
    
    pplot.figure()
    pplot.plot(rs,Bz)
    pplot.plot(rs,Bphi)
    pplot.plot(rs,Br)