# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 17:24:27 2015

@author: randale
"""

import numpy as np
from scipy.misc import factorial
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
        n1 = (2j+1)*factorial(j+j1-j2)*factorial(j-j1+j2)*factorial(j1+j2-j)
        d1 = factorial(j1+j2+j+1)
        f1 = math.sqrt(n1/float(d1))
        f2 = math.sqrt(factorial(j+m)*factorial(j-m)*factorial(j1+m1)*\
            factorial(j1-m1)*factorial(j2+m2)*factorial(j2-m2))
        sum = 0
        k_max = min((j1+j2-j,j1-m1,j2+m2))
        k_min = max((j2-j-m1,j1+m2-j,0))
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
    
def fs_matrix(n,basis_jtot):
    fs_ij = np.zeros((len(basis_jtot),len(basis_jtot)))
    for ii in range(len(basis_jtot)):
        b_i = basis_jtot[ii]
        j = b_i[2]
        fs_ij[ii,ii] = -13.6/float(n**2)*((1/(n*137.0)**2)*(n/(j+0.5)-0.75))
        
    return fs_ij

#def cg_matrix_tostring(cg_ij,basis_j1j2,basis_jtot):
#    