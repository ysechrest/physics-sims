# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 17:24:27 2015

@author: randale
"""

import numpy as np
from scipy.misc import factorial
import scipy.constants as scon
import math

def get_basis_j1j2(j1_range,j2_range):

    basis = []
    
    for j1_i in np.arange(j1_range[0],j1_range[1]+1,1):
        for j2_i in np.arange(j2_range[0],j2_range[1]+1,1):
            for m1 in np.arange(-j1_i,j1_i+1,1.0):
                for m2 in np.arange(-j2_i,j2_i+1.0,1.0):
                    b_i = (j1_i,m1,j2_i,m2)
                    basis.append(b_i)
                
    return basis
    
def get_basis_jm(j1_range,j2_range):
    
    basis = []
    
    for j1_i in np.arange(j1_range[0],j1_range[1]+1,1):
        for j2_i in np.arange(j2_range[0],j2_range[1]+1,1):
            j_range = (math.fabs(j1_i-j2_i),j1_i+j2_i)
            for j_i in np.arange(j_range[0],j_range[1]+1,1):
                for m in np.arange(-j_i,j_i+1,1.0):
                    b_i = (j1_i,j2_i,j_i,m)
                    basis.append(b_i)
            
    return basis
    
def cg_coeff(j1,m1,j2,m2,j,m):
    if m != m1+m2:
        return 0
    elif (math.fabs(j1-j2) > j) or (j1 + j2 < j):
        return 0
    else:
        n1 = (2*j+1)*factorial(j+j1-j2)*factorial(j-j1+j2)*factorial(j1+j2-j)
        d1 = factorial(j1+j2+j+1)
        f1 = math.sqrt(n1/float(d1))
        f2 = math.sqrt(factorial(j+m)*factorial(j-m)*factorial(j1+m1)*\
            factorial(j1-m1)*factorial(j2+m2)*factorial(j2-m2))
        sum = 0
        k_max = min((j1+j2-j,j1-m1,j2+m2))
        k_min = max((-j+j2-m1,-j+j1+m2,0))
        #print k_min
        #print k_max
        for k in np.arange(k_min,k_max+1,1):
            s1 = (j1+j2-j-k)
            s2 = (j1-m1-k)
            s3 = (j2+m2-k)
            s4 = (j-j2+m1+k)
            s5 = (j-j1-m2+k)
            Dk = factorial(k)*factorial(s1)*factorial(s2)*\
                factorial(s3)*factorial(s4)*factorial(s5)
            s_k = ((-1)**k)/float(Dk)
            sum += s_k

        return f1*f2*sum
        
def wigner_3j(j1,m1,j2,m2,j3,m3):
    return (-1)**(j1-j2-m3)*cg_coeff(j1,m1,j2,m2,j3,-m3)
    
def ylm_coupling(j1,m1,j2,m2,j3,m3):
    return math.sqrt((2*j1+1)*(2*j2+1)/(4.0*math.pi*(2*j3+1)))*\
        wigner_3j(j1,m1,j2,m2,j3,m3)*wigner_3j(j1,0,j2,0,j3,0)

def cg_matrix(basis_j1j2,basis_jtot):
    cg_ij = np.zeros((len(basis_j1j2),len(basis_jtot)))
    ii=0
    for b1_i in basis_j1j2:
        jj=0
        for b2_i in basis_jtot:
            j1 = b1_i[0]
            j2 = b1_i[2]
            m1 = b1_i[1]
            m2 = b1_i[3]
            j = b2_i[2]
            m = b2_i[3]
            jj1 = b2_i[0]
            jj2 = b2_i[1]
            if (jj1 != j1) or (jj2 != j2):
                cg_ij[ii,jj] = 0
            else:
                cg_ij[ii,jj] = cg_coeff(j1,m1,j2,m2,j,m)
            jj+=1
        ii+=1
    
    for ii in np.arange(len(basis_j1j2)):
        v_ii = cg_ij[...,ii]
        v_ii_norm = v_ii/math.sqrt(sum(v_ii*v_ii))
        cg_ij[...,ii] = v_ii_norm
    
    return cg_ij
    
def stark_matrix(n,Ex,Ey,Ez):
    #***Atomic Units Used***
    # Ex,Ey,Ez should be [V/m]*(e*a_0/E_hartree)
    basis_j1j2 = get_basis_j1j2((0,n-1),(-0.5,0.5))
    h_stark = np.zeros((len(basis_j1j2),len(basis_j1j2)),dtype=np.complex)
    ii=0
    for b1_i in basis_j1j2:
        jj=0
        for b2_i in basis_j1j2:
            j1 = b1_i[0]
            m1 = b1_i[1]
            ms1 = b1_i[3]
            j2 = b2_i[0]
            m2 = b2_i[1]
            ms2 = b2_i[3]
            if ms1 == ms2:
                if j2 == j1 + 1:
                    R_nl = (3/2.0)*n*math.sqrt(n**2-j2**2)
                elif j1 == j2 + 1:
                    R_nl = (3/2.0)*n*math.sqrt(n**2-j1**2)
                else:
                    R_nl = 0.0
                h_stark[ii,jj] = math.sqrt(4*math.pi/3.0)*(\
                    (-(Ex-Ey*1J)/math.sqrt(2.0))*(-1)**(-m1)*ylm_coupling(j2,m2,1,1,j1,-m1) +\
                    ((Ex+Ey*1J)/math.sqrt(2.0))*(-1)**(-m1)*ylm_coupling(j2,m2,1,-1,j1,-m1) +\
                    (Ez+0J)*(-1)**(-m1)*ylm_coupling(j2,m2,1,0,j1,-m1))*R_nl
            else:
                h_stark[ii,jj]=0.0
            jj+=1
        ii+=1
    #h_stark = h_stark/2.0
        
    return h_stark
    
def get_stark_states(n,Ex,Ey,Ez):
    h_stark = stark_matrix(n,Ex,Ey,Ez)
    w,v = np.linalg.eig(h_stark)
    
    #print(w)
    #print(v[:,0])
    return w,v

def zeeman_matrix(n,dW_0):
    basis_j1j2 = get_basis_j1j2((0,n-1),(-0.5,0.5))
    h_zeeman = np.zeros((len(basis_j1j2),len(basis_j1j2)),dtype=np.complex)
    ii=0
    for b1_i in basis_j1j2:
        jj=0
        for b2_i in basis_j1j2:
            j1 = b1_i[0]
            m1 = b1_i[1]
            ms1 = b1_i[3]
            j2 = b2_i[0]
            m2 = b2_i[1]
            ms2 = b2_i[3]
            if m1 != m2 or j1 !=j2 or ms1 != ms2:
                h_zeeman[ii,jj] = 0
            else:
                h_zeeman[ii,jj] = dW_0*(m1 + 2.0*ms1)
            jj+=1
        ii+=1
        
    return h_zeeman
    
def motional_stark_matrix(n,E,B,V):
    E_L = np.cross(V,B)
    E_tot = E + E_L
    modE = np.linalg.norm(E_tot)
    modB = np.linalg.norm(B)
    dW_0 = scon.value('Bohr magneton')*modB #Energy normalization factor
    beta = (scon.value('Bohr magneton')/(scon.e * scon.value('Bohr radius')))*(modB/modE)
    b_hat = B/modB
    e_hat = E_tot/modE
    
    H_s =  stark_matrix(n,e_hat[0],e_hat[1],e_hat[2])
    H_z = zeeman_matrix(n,1.0)
    H_mse = (H_z + (1/beta)*H_s)
    return H_mse,dW_0
    

def fs_matrix(n,basis_jtot):
    fs_ij = np.zeros((len(basis_jtot),len(basis_jtot)))
    for ii in range(len(basis_jtot)):
        b_i = basis_jtot[ii]
        j = b_i[2]
        fs_ij[ii,ii] = -13.6/float(n**2)*((1/(n*127.0)**2)*(n/(j+0.5)-0.75))
        
    return fs_ij

#def cg_matrix_tostring(cg_ij,basis_j1j2,basis_jtot):
#    