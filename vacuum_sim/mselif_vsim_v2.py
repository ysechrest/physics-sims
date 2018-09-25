# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 17:05:45 2016

@author: randale
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 11:25:39 2015

@author: randale
"""

import numpy as np
import matplotlib.pyplot as pplot
import scipy.signal as spsig

def dpc_dt(pc,pnc):
    return (C1*(pnc-pc) - S*pc)/(V_tot)
    
def dps_dt(pnc,ps,q):
    return (q - C*(ps-pnc))/V_source
    
def dpnc_dt(ps,pc,pnc):
    return(C*(ps-pnc) + q_nc - C1*(pnc-pc))/V_nc
    
def main():
    global V_chamber
    global V_flight
    global V_source
    global V_tot
    global V_nc
    global C
    global C1
    global S
    global q_nc

    V_chamber = 4094*1.64e-5 #[m^3]
    V_flight = 254*1.64e-5 #[m^3]
    V_source = 24*1.64e-5 #[m^3]
    V_tot = V_chamber+V_flight
    V_nc = np.pi*(2.54)**2*33.0*1.0e-6
    #V_nc = V_nc
    #C1 = 11.6*np.pi*(2.54)**2*1.0e-3
    C1 = 0.1505 #1.6 in diameter[m^3/s]
    C = 0.0125 #1.2 cm diameter[m^3/s]
    S = 2200*.001 #2200*0.001(for air, hydrogen=?)[m^3/s]
    #S = 510*.001 #2200*0.001(for air, hydrogen=?)[m^3/s] 

    t = 0.6
    dt = .0005
    nsamps = t/dt
    p_chamber = np.zeros((nsamps))
    p_chamber[0] = 1.0e-6/0.0075
    p_source = np.zeros((nsamps))
    p_source[0] = 1.0e-6/0.0075
    p_neutral = np.zeros((nsamps))
    p_neutral[0] = 1.0e-6/0.0075
    q_gas = np.zeros((2*nsamps-1))+0.75e-1
    q_nc = 0.0
    #q_gas[000:200] = 0.2
    #b,a = spsig.butter(8,0.04)
    #q_gas = spsig.filtfilt(b,a,q_gas,padlen=150)
    t = np.arange(nsamps)*dt    
    
    for ii in np.arange(1,nsamps):
        pc0 = p_chamber[ii-1]
        ps0 = p_source[ii-1]
        pnc0 = p_neutral[ii-1]
        q0 = q_gas[2*(ii-1)]
        
        k1_c = dpc_dt(pc0,pnc0)
        k1_s = dps_dt(pnc0,ps0,q0)
        k1_nc = dpnc_dt(ps0,pc0,pnc0)
        
        k2_c = dpc_dt(pc0+dt/2.0*k1_c,pnc0+dt/2.0*k1_nc)
        k2_s = dps_dt(pnc0+dt/2.0*k1_nc,ps0+dt/2.0*k1_s,q_gas[2*ii-1])
        k2_nc = dpnc_dt(ps0+dt/2.0*k1_s,pc0+dt/2.0*k1_c,pnc0+dt/2.0*k1_nc)
        
        k3_c = dpc_dt(pc0+dt/2.0*k2_c,pnc0+dt/2.0*k2_nc)
        k3_s = dps_dt(pnc0+dt/2.0*k2_nc,ps0+dt/2.0*k2_s,q_gas[2*ii-1])
        k3_nc = dpnc_dt(ps0+dt/2.0*k2_s,pc0+dt/2.0*k2_c,pnc0+dt/2.0*k2_nc)        
        
        k4_c = dpc_dt(pc0+dt*k3_c,pnc0+dt*k3_nc)
        k4_s = dps_dt(pnc0+dt*k3_nc,ps0+dt*k3_s,q_gas[2*ii])
        k4_nc = dpnc_dt(ps0+dt*k3_s,pc0+dt*k3_c,pnc0+dt*k3_nc)        
        
        p_chamber[ii] = pc0 + (dt/6.0)*(k1_c+2*k2_c+2*k3_c+k4_c)
        p_source[ii] = ps0 + (dt/6.0)*(k1_s+2*k2_s+2*k3_s+k4_s)
        p_neutral[ii] = pnc0 + (dt/6.0)*(k1_nc+2*k2_nc+2*k3_nc+k4_nc)        
        
    l1 = pplot.plot(t,p_chamber*7.5e-3*1.0e6/10.0,label="Pc [1.0e-5 Torr]")
    l4 = pplot.plot(t,p_neutral*7.5e-3*1.0e3,label='Pnc [mTorr]')
    l2 = pplot.plot(t,p_source*7.5e-3*1.0e3,label="Ps [mTorr]")
    l3 = pplot.plot(t,q_gas[::2]*0.0075*1000*100,label="Q gas x 100 [Torr l/s]")
    pplot.xlabel("time [sec]")
    pplot.legend(loc=4)