# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 12:07:47 2015

@author: randale
"""
from scipy.misc import factorial
from scipy.special import hyp2f1
import scipy.constants as scon
from math import sqrt
from math import pi
import numpy as np
import quantum.qm_ang_mom as qm
import matplotlib.pyplot as pplot
def reduced_matrix_ele_radial(n,l,nprime):
    f1 = (-1)**(nprime-l)/(4.0*factorial((2*l-1)))
    f2 = sqrt((factorial(n+l)*factorial(nprime+l-1))/(factorial(n-l-1)*factorial(nprime-l)))
    f3 = ((4*n*nprime)**(l+1)*(n-nprime)**(n+nprime-2*l-2))/(float(n+nprime)**(n+nprime))
    F1 = hyp2f1(-n+l+1,-nprime+l,2*l,-(4*n*nprime)/float(n-nprime)**2)
    F2 = hyp2f1(-n+l-1,-nprime+l,2*l,-(4*n*nprime)/float(n-nprime)**2)
    R = f1*f2*f3*(F1-F2*((n-nprime)/float(n+nprime))**2)
    return R
    
def balmer_dipole_r_integral(n,l,l_up):
    if l == 0 and l_up == 1:
        Rsq_nl = 2.0**(17)*(n**7)*(n**2-1)*(n-2)**(2*n-6)/\
            (n+2)**(2*n+6)
    elif l == 1 and l_up == 2:
        Rsq_nl = 2.0**(19)*(n**9)*(n**2-1)*(n-2)**(2*n-7)/\
            (3*(n+2)**(2*n+7))
    elif l == 1 and l_up == 0:
        Rsq_nl = 2.0**(15)*(n**9)*(n-2)**(2*n-6)/\
            (3*(n+2)**(2*n+6))
    else:
        Rsq_nl = 0
    return sqrt(Rsq_nl)
    
#def dipole_matrix_element(n1,j1,m1,n2,j2,m2):
#    R_nl = reduced_matrix_ele_radial(n1,j1,n2)
#    r_x = sqrt(2*pi/3.0)*(ylm_coupling(j1,m1,1,1,j2,m2)-ylm_coupling(j1,m1,1,-1,j2,m2))
#    r_y = sqrt(2*pi/3.0)*(ylm_coupling(j1,m1,1,1,j2,m2)+ylm_coupling(j1,m1,1,-1,j2,m2))*1J
#    r_z = sqrt(4*pi/3.0)*ylm_coupling(j1,m1,1,0,j2,m2)
  
def dipole_matrix(n1,n2):
    #sign prolbem on rz?
    basis_n1 = qm.get_basis_j1j2((0,n1-1),(-0.5,0.5))
    basis_n2 = qm.get_basis_j1j2((0,n2-1),(-0.5,0.5))
    r_x = np.zeros((len(basis_n2),len(basis_n1)),dtype=np.complex)
    r_y = r_x.copy()
    r_z = r_x.copy()
    ii=0
    for b1_i in basis_n1:
        jj=0
        for b2_j in basis_n2:
            j1 = b1_i[0]
            m1 = b1_i[1]
            ms1 = b1_i[3]
            j2 = b2_j[0]
            m2 = b2_j[1]
            ms2 = b2_j[3]
            if np.abs(j2-j1) <= 1 and np.abs(m1-m2) < 2 and j1+j2 > 0:
                if j2 > j1:
                    R_nl = -1.0*reduced_matrix_ele_radial(n2,j2,n1)
                else:
                    R_nl = -1.0*reduced_matrix_ele_radial(n1,j1,n2)
                r_x[jj,ii] = sqrt(2*pi/3.0)*(-1)**(-m2+1)*(qm.ylm_coupling(j1,m1,1,1,j2,-m2)-qm.ylm_coupling(j1,m1,1,-1,j2,-m2))*R_nl
                r_y[jj,ii] = sqrt(2*pi/3.0)*(-1)**(-m2)*(qm.ylm_coupling(j1,m1,1,1,j2,-m2)+qm.ylm_coupling(j1,m1,1,-1,j2,-m2))*1J*R_nl
                r_z[jj,ii] = sqrt(4*pi/3.0)*(-1)**(-m2)*qm.ylm_coupling(j1,m1,1,0,j2,-m2)*R_nl
            jj=jj+1
        ii=ii+1
    return (r_x,r_y,r_z)
    
def stark_hybrid_dipole_matrix(n1,n2,Ex,Ey,Ez):
    delE_1,v_1 = qm.get_stark_states(n1,Ex,Ey,Ez)
    delE_2,v_2 = qm.get_stark_states(n2,Ex,Ey,Ez)
    dn_x,dn_y,dn_z = dipole_matrix(n2,n1)
    A = np.dot(np.conjugate(np.transpose(v_1)),dn_x)
    dk_x = np.dot(A,v_2)
    A = np.dot(np.conjugate(np.transpose(v_1)),dn_y)
    dk_y = np.dot(A,v_2)
    A = np.dot(dn_z,v_2)
    dk_z = np.dot(np.conjugate(np.transpose(v_1)),A)
    
    delE = dk_z.copy()*0.0
    for ii in np.arange(0,delE.shape[1]):
        for jj in np.arange(0,delE.shape[0]):
            delE[jj,ii] = delE_2[ii] - delE_1[jj]
            
    amp = np.sqrt(np.absolute(dk_x)**2+ \
        np.absolute(dk_y)**2+np.absolute(dk_z)**2)
    u,v = bin_values(np.real(delE.ravel()),np.real(amp.ravel()),17.0)
    pplot.figure()
    pplot.bar(v,u)
            
    return dk_x,dk_y,dk_z,delE
    
def mse_hybrid_dipole_matrix(n1,n2,B,V):
    H1,dW_0 = qm.motional_stark_matrix(n1,[0,0,0],B,V)
    H2,dW_0 = qm.motional_stark_matrix(n2,[0,0,0],B,V)
    delE_1,v_1 = np.linalg.eig(H1)
    delE_2,v_2 = np.linalg.eig(H2)
    dn_x,dn_y,dn_z = dipole_matrix(n2,n1)
    A = np.dot(np.conjugate(np.transpose(v_1)),dn_x)
    dk_x = np.dot(A,v_2)
    A = np.dot(np.conjugate(np.transpose(v_1)),dn_y)
    dk_y = np.dot(A,v_2)
    A = np.dot(np.conjugate(np.transpose(v_1)),dn_z)
    dk_z = np.dot(A,v_2)
    
    delE_0 = scon.value('Rydberg constant times hc in J')*(1/4.0-1/9.0)
    
    delE = dk_z.copy()*0.0
    for ii in np.arange(0,delE.shape[1]):
        for jj in np.arange(0,delE.shape[0]):
            delE[jj,ii] = delE_2[ii] - delE_1[jj]
            
    #amp = np.sqrt(np.absolute(dk_x)**2+ \
    #    np.absolute(dk_y)**2+np.absolute(dk_z)**2)
    #u,v = bin_values(np.real(scon.h*scon.c/(dW_0*delE.ravel()+delE_0)),np.real(amp.ravel()),17.0)
    #u,v = bin_values(np.real(delE.ravel()),np.real(amp.ravel()),17.0)
    #pplot.figure()
    #pplot.bar(v,u)
            
    return dk_x,dk_y,dk_z,delE
    
def mse_simulate_spectrum(B,V):
    dk_x,dk_y,dk_z,delE = mse_hybrid_dipole_matrix(2,3,B,V)
    #amp = np.sqrt(np.absolute(dk_x)**2+ \
    #    np.absolute(dk_y)**2+np.absolute(dk_z)**2)
    modB = np.sqrt(B[0]**2+B[1]**2+B[2]**2)
    modV = np.sqrt(V[0]**2+V[1]**2+V[2]**2)
    dW_0 = scon.value('Bohr magneton')*modB
    dE = 1.1*np.max(np.real(delE))*dW_0/scon.h*1.0e-9
    #dE =5*dW_0/(15.0*scon.h)*1.0e-9
    x = np.arange(0,2*dE+.1,.1) - dE
    s_pip = x.copy()*0.0
    s_pim = x.copy()*0.0
    s_sig = x.copy()*0.0
    pplot.figure()
    for ii in np.arange(0,dk_x.shape[0]):
        for jj in np.arange(0,dk_x.shape[1]):
            dE_ij = np.real(delE[ii,jj])
            dnu = (dE_ij*dW_0/scon.h)*1.0e-9
            a_ij = np.absolute(dk_x[ii,jj]+dk_y[ii,jj]*1j)**2
            s_pim_ij = a_ij*np.exp(-4*np.log(2)*((x-dnu)/13.06)**2)
            a_ij = np.absolute(dk_x[ii,jj]-dk_y[ii,jj]*1j)**2
            s_pip_ij = a_ij*np.exp(-4*np.log(2)*((x-dnu)/13.06)**2)
            a_ij = np.absolute(dk_z[ii,jj])**2
            s_sig_ij = a_ij*np.exp(-4*np.log(2)*((x-dnu)/13.06)**2)
            s_pip = s_pip + s_pip_ij
            s_pim = s_pim + s_pim_ij
            s_sig = s_sig + s_sig_ij
    pplot.plot(x,-s_pip,'r')
    pplot.plot(x,-s_pim,'b')
    pplot.plot(x,s_sig,'g')
    #pplot.figure()
    #pplot.plot(delE.ravel()*dW_0/scon.h*1.0e-9,delE.ravel()*0.0,'ro')    
    
def mse_spectrum_view(alpha,beta,gamma,modB,modV,width):
    
    Bp = np.array([0,np.sin(gamma),np.cos(gamma)])*modB
    Vp = np.array([modV,0,0])
    R_x = np.array([[1,0,0],\
        [0,np.cos(gamma),-np.sin(gamma)],\
        [0,np.sin(gamma),np.cos(gamma)]])
    s_hat_p = [-np.cos(alpha)*np.cos(beta),\
        np.sin(beta),np.sin(alpha)*np.cos(beta)]
    s_hat = np.dot(R_x,s_hat_p)
    V = np.dot(R_x,Vp)
    B = np.dot(R_x,Bp)
    #print(s_hat_p)
    #print(s_hat)
    #print(B)
    
    #pplot.figure()
    #pplot.plot([0,B[2]],[0,B[1]])
    #pplot.plot([0,V[2]],[0,V[1]])
    #pplot.plot([0,s_hat[2]],[0,s_hat[1]])
    #ax = pplot.gca()
    #ax.set_xlim([-1,1])
    #ax.set_ylim([-1,1])
    
    e2_hat = np.cross(s_hat,B/modB)
    e2_hat = e2_hat/np.linalg.norm(e2_hat)
    e1_hat = np.cross(s_hat,e2_hat)
    e1_hat = e1_hat/np.linalg.norm(e1_hat)
    
    print(e1_hat)
    print(e2_hat)
    
    dk_x,dk_y,dk_z,delE = mse_hybrid_dipole_matrix(2,3,B,V)
    I_e1 = np.absolute(dk_x*e1_hat[0]+dk_y*e1_hat[1]+dk_z*e1_hat[2])**2
    I_e2 = np.absolute(dk_x*e2_hat[0]+dk_y*e2_hat[1]+dk_z*e2_hat[2])**2
    #I_tot = I_e1 + I_e2
    
    dW_0 = scon.value('Bohr magneton')*modB
    dE = 1.1*np.max(np.real(delE))*dW_0/scon.h*1.0e-9
    #dE =5*dW_0/(15.0*scon.h)*1.0e-9
    x = np.arange(0,2*dE+.1,.1) - dE
    s_e1 = x.copy()*0.0
    s_e2 = x.copy()*0.0
    pplot.figure()
    for ii in np.arange(0,dk_x.shape[0]):
        for jj in np.arange(0,dk_x.shape[1]):
            
            dE_ij = np.real(delE[ii,jj])
            dnu = (dE_ij*dW_0/scon.h)*1.0e-9
            a_ij = I_e1[ii,jj]
            s_e1_ij = a_ij*np.exp(-4*np.log(2)*((x-dnu)/width)**2)
            a_ij = I_e2[ii,jj]
            s_e2_ij = a_ij*np.exp(-4*np.log(2)*((x-dnu)/width)**2)
            s_e1 = s_e1 + s_e1_ij
            s_e2 = s_e2 + s_e2_ij
    pplot.plot(x,s_e1,'r')
    pplot.plot(x,-s_e2,'b')
    pplot.plot(x,s_e1+s_e2,'g')
    #pplot.figure()
    #pplot.plot(delE.ravel()*dW_0/scon.h*1.0e-9,delE.ravel()*0.0,'ro')
    
def stark_zeeman_comp(n):
    b_hat = np.array([0.0,0.0,1.0])
    Vmod = np.sqrt(2*3.0e4*scon.e/scon.m_p)
    V = Vmod*np.array([1.0,0.0,0.0])
    delE_0 = -1.0*scon.value('Rydberg constant times hc in eV')*(1/float(n)**2)    
    
    for bmod in (np.arange(1,30)/10.0):
        B = b_hat*bmod
        H1,dW_0 = qm.motional_stark_matrix(n,[0,0,0],B,V)
        H1 = H1*dW_0
        delE_1,v_1 = np.linalg.eig(H1)
        E_L = np.cross(V,B)
        H2 = qm.stark_matrix(n,E_L[0],E_L[1],E_L[2])
        H2 = H2*scon.e*scon.value('Bohr radius')
        delE_2,v_2 = np.linalg.eig(H2)  
        #E_sz.append(delE_1/scon.e)
        #E_s.append(delE_2/scon.e)
        pplot.plot(bmod+delE_1*0.0,delE_0 + delE_1/scon.e,'ro')
        #pplot.plot(bmod+delE_1*0.0,delE_0 + delE_1/scon.e - dW_0/scon.e ,'go')
        #pplot.plot(bmod+delE_1*0.0,delE_0 + delE_1/scon.e + dW_0/scon.e ,'go')
        pplot.plot(bmod+delE_2*0.0,delE_0 + delE_2/scon.e,'bo')
        #pplot.plot(bmod+delE_1*0.0, delE_1/(scon.e*scon.value('Bohr radius')*E_L),'ro')
    
def bin_values(arr,weight,nbins):
    minv = np.min(arr)
    maxv = np.max(arr)
    dv = (maxv-minv)/float(nbins-1)
    print(minv,maxv,dv)
    v = np.arange(minv-dv/2.0,maxv+dv/2.0+dv,dv)
    u = np.zeros(v.shape)
    for ii in np.arange(0,len(arr)):
        a_ii = arr[ii]
        idx = np.floor((a_ii-np.min(v))/dv)
        u[idx] = u[idx] + weight[ii]
    #v = v + dv/2.0
    return u,v
    
    
def print_matrix(A):
    real_str = ''
    im_str = ''
    for jj in np.arange(0,A.shape[0]):
        for ii in np.arange(0,A.shape[1]):
            real_str = real_str + ' {:+1.3f}'.format(np.real(A)[jj,ii])
            im_str = im_str + ' {:+1.3f}'.format(np.imag(A)[jj,ii])
        real_str = real_str + '\n'
        im_str = im_str + '\n'
    print(real_str)
    print(im_str)
            
def test_ylm_coupling(n):
    basis_n = qm.get_basis_j1j2((0,n-1),(-0.5,0.5))
    r_p = np.zeros((len(basis_n),len(basis_n)),dtype=np.complex)
    r_m = r_p.copy()
    r_0 = r_p.copy()
    ii=0
    for b1_i in basis_n:
        jj=0
        for b2_j in basis_n:
            j1 = b1_i[0]
            m1 = b1_i[1]
            j2 = b2_j[0]
            m2 = b2_j[1]
            r_p[jj,ii] = qm.ylm_coupling(j2,m2,1,1,j1,m1)
            r_m[jj,ii] = qm.ylm_coupling(j2,m2,1,-1,j1,m1)
            r_0[jj,ii] = qm.ylm_coupling(j2,m2,1,0,j1,m1)
            jj=jj+1
        ii=ii+1
    return (r_p,r_m,r_0)