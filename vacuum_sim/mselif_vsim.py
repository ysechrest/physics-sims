# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 11:25:39 2015

@author: randale
"""

import numpy as np
import matplotlib.pyplot as pplot
import scipy.signal as spsig

def dpc_dt(pc,ps):
    return (C*(ps-pc) - S*pc)/(V_tot)
    
def dps_dt(pc,ps,q):
    return (q - C*(ps-pc))/V_source
    
def main():
    global V_chamber
    global V_flight
    global V_source
    global V_tot
    global C
    global S

    V_chamber = 4094*1.64e-5 #[m^3]
    V_flight = 254*1.64e-5 #[m^3]
    V_source = 24*1.64e-5 #[m^3]
    V_tot = V_chamber+V_flight
    C = 0.018 #[m^3/s]
    S = 2700*.001 #[m^3/s]  

    t = 0.6
    dt = .0005
    nsamps = t/dt
    p_chamber = np.zeros((nsamps))
    p_chamber[0] = 1.0e-6/0.0075
    p_source = np.zeros((nsamps))
    p_source[0] = 1.0e-6/0.0075
    q_gas = np.zeros((2*nsamps-1))+.12
    #q_gas[000:200] = 0.2
    #b,a = spsig.butter(8,0.04)
    #q_gas = spsig.filtfilt(b,a,q_gas,padlen=150)
    t = np.arange(nsamps)*dt    
    
    for ii in np.arange(1,nsamps):
        pc0 = p_chamber[ii-1]
        ps0 = p_source[ii-1]
        q0 = q_gas[2*(ii-1)]
        
        k1_c = dpc_dt(pc0,ps0)
        k1_s = dps_dt(pc0,ps0,q0)
        
        k2_c = dpc_dt(pc0+dt/2.0*k1_c,ps0+dt/2.0*k1_s)
        k2_s = dps_dt(pc0+dt/2.0*k1_c,ps0+dt/2.0*k1_s,q_gas[2*ii-1])
        
        k3_c = dpc_dt(pc0+dt/2.0*k2_c,ps0+dt/2.0*k2_s)
        k3_s = dps_dt(pc0+dt/2.0*k2_c,ps0+dt/2.0*k2_s,q_gas[2*ii-1])
        
        k4_c = dpc_dt(pc0+dt*k3_c,ps0+dt*k3_s)
        k4_s = dps_dt(pc0+dt*k3_c,ps0+dt*k3_s,q_gas[2*ii])
        
        p_chamber[ii] = pc0 + (dt/6.0)*(k1_c+2*k2_c+2*k3_c+k4_c)
        p_source[ii] = ps0 + (dt/6.0)*(k1_s+2*k2_s+2*k3_s+k4_s)
        
    l1 = pplot.plot(t,p_chamber*7.5e-3*1.0e6,label="Pc [microTorr]")
    l2 = pplot.plot(t,p_source*7.5e-3*1.0e3,label="Ps [mTorr]")
    l3 = pplot.plot(t,q_gas[::2]*0.0075*1000*100,label="Q gas x 100 [Torr l/s]")
    pplot.xlabel("time [sec]")
    pplot.legend()