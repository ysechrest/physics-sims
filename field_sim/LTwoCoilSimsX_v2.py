# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 20:56:29 2015

@author: mgalante
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
           
def my_range(start,end,step):
    while start <= end:
        yield start
        start += step

def cross_mg(v1,v2):
    x = (v1[1,:]*v2[2,:])-(v1[2,:]*v2[1,:])
    y = (v1[2,:]*v2[0,:])-(v1[0,:]*v1[2,:])
    z = (v1[0,:]*v2[1,:])-(v1[1,:]*v2[0,:])
    
    return x,y,z

def biot_savart(vdl,va,vpos):
        
    import numpy as np
    
    #vpos = np.transpose(vpos)
    c = 3.0e10
    i = np.shape(va)
    i = i[1]
    ones = np.ones(i)
    vr = np.outer(vpos,ones)-va
    r3 = ((vr[0,:]**2+vr[1,:]**2+vr[2,:]**2)**(1.5))
    #r3 = (np.dot(vr,vr))**(1.5)
    top = np.cross(vdl,vr,axisa=0,axisb=0,axisc=0)
    bottom = np.outer(np.ones(3),r3)
    inside = top/bottom
    B = 1/c * np.sum(inside,axis=1)
    
    return B
    
def mag(vector):
    import numpy as np    
    
    v1 = np.sqrt(np.sum(vector*vector))
    return v1

def field_magnitude_vec(coil_path,rvec):
    import numpy as np    
    
    sizevec = np.shape(rvec)
    length = sizevec[0]
    Bvec = np.zeros([length,3])

 
    for i in range(0,length):
        Bvec[i,:] = biot_savart(coil_path.vdl,coil_path.va,rvec[i,:])
        
        
    return Bvec
    
def field_magnitude(coil_path, rvec):
    import numpy as np    
    
    sizevec = np.shape(rvec)
    length = sizevec[0]
    Bvec = np.zeros([length,3])
    Bmag = np.zeros(length)
 
    for i in range(0,length):
        Bvec[i,:] = biot_savart(coil_path.vdl,coil_path.va,rvec[i,:])
        sign = np.sign(Bvec[i,2])
        Bmag[i] = sign * mag(Bvec[i,:])
    
    return Bmag
    
def field_magnitude_sign(coil_path,rvec):
    import numpy as np
    
    sizevec = np.shape(rvec)
    length = sizevec[0]
    Bsign = np.zeros([length,3])
    
    for i in range(0,length):
        Bvec = biot_savart(coil_path.vdl,coil_path.va,rvec[i,:])
        Bsign[i,:] = np.sign(Bvec)
        
    return Bsign
    
def field_magnitudetd(coil_path,xvec,zvec):
    import numpy as np
     
    sizexvec = np.shape(xvec)
    xlength = sizexvec[0]
    sizezvec = np.shape(zvec)
    zlength = sizezvec[0]
    Bvec = np.zeros([3,xlength,zlength])
    Bmag = np.zeros([xlength,zlength])
    for i in range(xlength):
        for j in range(zlength):
            Bvec[:,i,j] = biot_savart(coil_path.vdl,coil_path.va,[xvec[i],0.0,zvec[j]])
            Bmag[i,j] = mag(Bvec[:,i,j])
        
    return Bmag
    
def field_tracer(coil_path,init_pos,bounds,axes,direc,ax=None):
    
    import numpy as np    
     
    plot_points = np.zeros((int(1e6),2))
    incr = 0
    verbosity = 0
    i = axes[0]
    j = axes[1]
    
    error_target = 1e-2
    
    error_max = error_target
    error_min = error_target/2
    
    error = 1.
    dt = 1e2
    ddt = 1.1
    
    maxstep = 5000
    step = 0
    
    r_0 = np.zeros(3)
    r_1 = np.zeros(3)
    
    r_0 = init_pos
    nodeath = 1
    if ax is None:
        ax=plt.gca()
    color = 'r'
    l1, = ax.plot([],[],color=color,ls='-')
    
    while(r_0[0] >= bounds[0,0] and r_0[0] <= bounds[0,1] and \
        r_0[1] >= bounds[1,0] and r_0[1] <= bounds[1,1] and \
        r_0[2] >= bounds[2,0] and r_0[2] <= bounds[2,1] and nodeath):
            
        B_0 = direc * biot_savart(coil_path.vdl,coil_path.va,r_0)
        
        r_1 = r_0 + dt*B_0
        
        B_1 = direc * biot_savart(coil_path.vdl,coil_path.va,r_1)
        r_1 = r_0 + 0.5 * dt * (B_1 + B_0)
        
        error = np.sqrt(np.sum((B_1-B_0)**2)/np.sum(B_0**2))
        
        while(error > error_max) or (error < error_min):
            r_1 = r_0 + dt * B_0
            B_1 = direc * biot_savart(coil_path.vdl,coil_path.va,r_1)
            
            r_1 = r_0 + 0.5*dt*(B_1 + B_0)
            
            B1B0 = (B_1-B_0)*(B_1-B_0)
            B02 = B_0*B_0
            
            Btop = np.sum(B1B0)
            Bbottom = np.sum(B02)
            error = np.sqrt(Btop/Bbottom)
        
            if(error < error_min):
                dt = dt * ddt
                if(verbosity): print 'dt=',dt,'   error=',error, ' (adapting: coarser step)'
        
            if(error > error_max):
                dt = dt/ddt
                if(verbosity): print 'dt=',dt,'   error=',error, ' (adapting: refining step)'
            
        if(step >= maxstep):nodeath=0
        if(mag(np.array(init_pos)-np.array(r_0))<= mag(np.array(r_1)-np.array(r_0))*2. and (step > 10)): nodeath=0
        
        #plt.plot([r_0[i],r_1[i]],[r_0[j],r_1[j]],'k')        
        plot_points[incr] = np.array([r_0[i],r_0[j]])
        plot_points[incr+1] = np.array([r_1[i],r_1[j]])
        incr += 2        
        r_0=r_1
        step = step+1
        #print step
        
    plt.plot(plot_points[:incr,0],plot_points[:incr,1],color=color,ls='-')
    plt.xlim(bounds[2,0],bounds[2,1])
    plt.ylim(bounds[1,0],bounds[1,1])
    
    #return {'plot_points':plot_points,'incr':incr}


def rotate_ang(vec,theta,axes):
    import numpy as np
    
    i = axes[0]
    j = axes[1]
    
    rotatedvec = np.zeros(3)
    rotatedvec[i] = vec[i] * np.cos(theta) + vec[j] * np.sin(theta)
    rotatedvec[j] = vec[j] * np.cos(theta) - vec[i] * np.sin(theta)
    
    return rotatedvec
    
def field_tracer_plot(coil_path,bounds,init_pos,axes,maxB,ax=None):
    
    import numpy as np
    
    calibration = maxB
    maxlines = 10.
    i = axes[0]
    j = axes[1]
    
    Bvec = np.zeros([3,maxlines+1])    
    for stepdir in my_range(-1,1,2):
        ipos = init_pos
        linenum = 0
        
        while( (ipos[0] >= bounds[0,0]) and (ipos[0] <= bounds[0,1]) and \
            (ipos[1] >= bounds[1,0]) and (ipos[1] <= bounds[1,1]) and \
            (ipos[2] >= bounds[2,0]) and (ipos[2] <= bounds[2,1]) and (linenum <= maxlines) ):
            field_tracer(coil_path,ipos,bounds,axes,1)
            field_tracer(coil_path,ipos,bounds,axes,-1)
            
            B = biot_savart(coil_path.vdl,coil_path.va,ipos)
            #print B.shape
            movevec = stepdir * rotate_ang(B/mag(B),np.pi/2.,axes)/mag(B)*calibration
            Bvec[:,linenum] = B          
            movevec_ij = np.array([0.,0.,0.])
            movevec_ij[i] = movevec[i]
            movevec_ij[j] = movevec[j]
            
            ipos = ipos + movevec_ij
            linenum = linenum +1
            print linenum
            
    return Bvec

def draw_box(inp,color,ax=None):

    
    x0 = inp[0][0]
    x1 = inp[0][1]
    y0 = inp[1][0]
    y1 = inp[1][1]    
    if ax is None:    
        plt.plot([x0,x1,x1,x0,x0,x1,x1,x0],[y0,y0,y1,y1,y0,y1,y0,y1],color)
    else:
        plt.plot([x0,x1,x1,x0,x0,x1,x1,x0],[y0,y0,y1,y1,y0,y1,y0,y1],color)

def draw_coils(coil_path,ax=None):
            
    for a_feature in coil_path.features:
        z = a_feature.z_loc
        OD = a_feature.OD
        ID = a_feature.ID
        w = a_feature.width/2.
        
        col = 'k'        
        
        draw_box([[z-w,z+w],[-OD/2,OD/2]],col,ax=ax)
        draw_box([[z-w,z+w],[-ID/2,ID/2]],col,ax=ax)
        
        

def draw_vv(ax=None):
    
    cm_per_in = 2.54    
    
    OD = 23.5 * cm_per_in
    w = 13.0/2. * cm_per_in
    z = 0.

    inp = [[z-w,z+w],[-OD/2.,OD/2.]]
    x0 = inp[0][0]
    x1 = inp[0][1]
    y0 = inp[1][0]
    y1 = inp[1][1]
    
    plt.plot([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],'b',linewidth=2)
    
def draw_port_center(OD,z,w,s,ax=None):
    
    cm_per_in = 2.54
    
    if s==1:    
        z = z * cm_per_in
    else:
        z = -1* z * cm_per_in 
    OD = OD * cm_per_in
    w = w * cm_per_in

    if s==1:
        inp = [[z+w,z],[-OD/2.,OD/2]]
    else:
        inp = [[z-w,z],[-OD/2.,OD/2]]
    
    
    x0 = inp[0][0]
    x1 = inp[0][1]
    y0 = inp[1][0]
    y1 = inp[1][1]    
    
    plt.plot([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],'b',linewidth=2)    
    
def draw_port_ocenter(xx,z,w,h,u,s,ax=None):
    
    cm_per_in = 2.54
    
    if u==1:
        xx = xx*cm_per_in
    else:
        xx = -1*xx*cm_per_in        
        
    if s==1:
        z = z * cm_per_in
    else:
        z = -1 * z * cm_per_in
        
    w = w*cm_per_in

    if s==1 and u==1:    
        inp = [[z,z+w],[xx,xx+h]]
    if s==1 and u==-1:
        inp = [[z,z+w],[xx,xx-h]]

    if s==-1 and u==1:    
        inp = [[z,z-w],[xx,xx+h]]
    if s==-1 and u==-1:
        inp = [[z,z-w],[xx,xx-h]]
        
        
    x0 = inp[0][0]
    x1 = inp[0][1]
    y0 = inp[1][0]
    y1 = inp[1][1]

    plt.plot([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],'b',linewidth=2)
    
    
class loop_coil(object):
    def _init__(self,current,radius,z_loc):
        
        import numpy as np
    
        res = 25.
        dtheta = 2. * np.pi/res
    
        theta = np.arange(res)*dtheta
    
        dl = current * radius * dtheta
        vdl = np.array([dl*(-np.sin(theta)),dl*np.cos(theta),np.zeros(res)])
    
        a = radius
        va = np.array([a*np.cos(theta),a*np.sin(theta),np.zeros(res)+z_loc])
        
        OD = radius*2.
        ID = radius*2.
        width = 0.
        layers = 1.
        turns = 1.
        i = 0
        init = 0
        
        self.z_loc = z_loc
        self.OD = OD
        self.ID = ID
        self.width = width
        self.layers = layers
        self.turns = turns
        self.current = current
    
        self.vdl = vdl
        self.va = va
        self.i = i
        self.init = init
        
        
class pancake_coil(object):
    def __init__(self,features):
        
        import numpy as np
        
        res = 4.
        dtheta = 2 * np.pi / res
        theta = np.arange(res) * dtheta
        
        #features = loop_coil(current,radius,z_loc)
        
        ID = features.ID
        OD = features.OD
        width = features.width
        turns = features.turns
        layers = features.layers
        current = features.current
        
        wire_thickness = width/layers
        ID = ID + wire_thickness
        OD = OD - wire_thickness
        
        r_spaces = turns/layers - 1
        r_spacing = (OD/2. - ID/2.)/r_spaces            
        
        z_spaces = layers-1
        z_spacing = wire_thickness
        
        radius = ID/2.
        z_loc = features.z_loc - z_spacing * z_spaces/2.
        dl = current * layers * radius * dtheta
        vdl = np.array([dl*(-np.sin(theta)),dl*np.cos(theta),np.zeros(res)])
        
        a = radius
        va = np.array([a*np.cos(theta),a*np.sin(theta),np.zeros(res)+features.z_loc])
        
        for z_space in my_range(0,z_spaces,1):
            z_loc = features.z_loc - z_spacing * z_spaces/2. + z_space*z_spacing
            for r_space in my_range(1,r_spaces,1):
                radius = ID/2. + r_space * r_spacing
                dl = current*radius * dtheta
                nextTermV = np.array([dl*(-np.sin(theta)),dl*np.cos(theta),np.zeros(res)])
                vdl = np.concatenate([vdl,nextTermV],axis=1)
                a = radius
                nextTermA = np.array([a*np.cos(theta),a*np.sin(theta),np.zeros(res)+z_loc])
                va = np.concatenate([va,nextTermA],axis=1)
                
        self.vdl = vdl
        self.va = va
        self.features = features
        self.i = 0
        self.init = 0
        

#Define the L2 coil dimensions.  L2 coils are used to expand the field lines as they exist the chamber.  
#These coils are not yet installed on the machine.

class L2_coil(object):
    def __init__(self,current,z_loc):
        
        cm_per_in = 2.54
        
        ID = 10.0 * cm_per_in
        OD = 19.938 * cm_per_in
        width = 0.75 * cm_per_in
        turns = 33.
        layers = 2.
        
        self.z_loc = z_loc
        self.OD = OD
        self.ID = ID
        self.width = width
        self.layers = layers
        self.turns = turns
        self.current = current
        
def L2_coil_eval(current,z_loc):
    st = L2_coil(current,z_loc)
    pc = pancake_coil(st)
    
    return pc

#Define the S1 Proto coil dimensions.  These are the coils currently installed on the machine.
        
class s1_proto_coil(object):
    def __init__(self,current,z_loc):
        
        cm_per_in = 2.54
        
        ID = 30.0 * cm_per_in
        OD = 37.0 * cm_per_in
        width = 0.825 * cm_per_in
        turns = 18.0
        layers = 2.
        
        self.z_loc = z_loc
        self.OD = OD
        self.ID = ID
        self.width = width
        self.layers = layers
        self.turns = turns
        self.current = current
        
def s1_proto_eval(current,z_loc):
    st = s1_proto_coil(current,z_loc)
    pc = pancake_coil(st)

    return pc    

#Define the dimensions for the "small" coils that are used as a Helmholtz pair around the designed proposed helicon antenna.
#These coils do not exist yet.
    
class small_coil(object):
    def __init__(self,current,z_loc):
        
        cm_per_in = 2.54
        
        ID = 6.0 * cm_per_in
        OD = 10.0 * cm_per_in
        width = 0.75 * cm_per_in
        turns = 18.0
        layers = 2.
        
        self.z_loc = z_loc
        self.OD = OD
        self.ID = ID
        self.width = width
        self.layers = layers
        self.turns = turns
        self.current = current
        
def small_coil_eval(current,z_loc):
    st = small_coil(current,z_loc)
    pc = pancake_coil(st)
    
    return pc
        

class new_coil_mg(object):
    def __init__(self,current,z_loc):

        cm_per_in = 2.54
        
        ID = 36.0 * cm_per_in
        OD = 40.0 * cm_per_in
        width = 0.825 * cm_per_in
        turns = 18.
        layers = 2.
        
        self.z_loc = z_loc
        self.OD = OD
        self.ID = ID
        self.width = width
        self.layers = layers
        self.turns = turns
        self.current = current
        
def new_coil_eval(current,z_loc):
    st = new_coil_mg(current,z_loc)
    pc = pancake_coil(st)
    
    return pc

class inner_coil_mg(object):
    def __init__(self,current,z_loc):
        
        cm_per_in = 2.54
        
        ID = 25. * cm_per_in
        OD = 29. * cm_per_in
        width = 0.825*cm_per_in
        turns = 18.0
        layers = 2.0


        self.z_loc = z_loc
        self.OD = OD
        self.ID = ID
        self.width = width
        self.layers = layers
        self.turns = turns
        self.current = current
        
def inner_coil_eval(current,z_loc):
    
    st = inner_coil_mg(current,z_loc)
    pc = pancake_coil(st)
    
    return pc
    
class Coil_Path(object):
    def __init__(self):
            self.vdl = None#new_path.vdl
            self.va = None#new_path.va
            self.features = []#new_path.features
            self.three_ones = np.ones(3)            
            self.i = 0
            self.i_ones = None
            self.va_shape = None
            self.vdl_shape = None
            
    def addpath(self,new_path):
        try:
            self.vdl = np.concatenate([self.vdl,new_path.vdl],axis=1)
            self.va = np.concatenate([self.va,new_path.va],axis=1)
        except ValueError:
            self.vdl = new_path.vdl
            self.va = new_path.va
            
            
        self.features.append(new_path.features)
        self.i += 1
        self.va_shape = np.shape(self.vdl)
        self.i_ones = np.ones(self.va_shape[1])
            
            

coil_path = Coil_Path()

cm_per_in = 2.54
statamp_per_amp = 3.0e9
cur = 500.
cur1 = 500.
L2cur = -400.

L2cur_f = 100.
cur_sm = 300

totalcoilcur = cur * statamp_per_amp#250
totalcoilcur1 = cur1 * statamp_per_amp#500
l2totalcoilcur = L2cur * statamp_per_amp#-400
l2totalcoil_f = L2cur_f * statamp_per_amp#300
totalcoil_sm = cur_sm * statamp_per_amp#1000

outercur = 0.0#200. * statamp_per_amp
innercur = 0.0#-400. * statamp_per_amp


startz = 4.
number = 6.

for z_in in my_range(startz,startz+(number-1)*0.9,0.9):   
    z_cm = z_in * cm_per_in
    s1 = s1_proto_eval(totalcoilcur1,z_cm)        
    coil_path.addpath(s1)

startz = 4
number = 6.

for z_in in my_range(startz,startz+(number-1)*0.9,0.9):    
    z_cm = z_in * cm_per_in
    s1 = s1_proto_eval(totalcoilcur,-z_cm) 
    coil_path.addpath(s1)
    
startz = 8.5
number = 1.

for z_in in my_range(startz,startz+(number-1)*0.75,0.75):
    z_cm = z_in*cm_per_in
    l2_fa = L2_coil_eval(l2totalcoil_f,-z_cm)
    coil_path.addpath(l2_fa)

    
startz = -14
number = 1.

for z_in in my_range(startz,startz+(number-1)*0.75,0.75):
    z_cm = z_in*cm_per_in
    sm = small_coil_eval(totalcoil_sm,z_cm)
    coil_path.addpath(sm)

startz = -14 + 6/2.
number = 1.

for z_in in my_range(startz,startz+(number-1)*0.75,0.75):
    z_cm = z_in*cm_per_in
    sm = small_coil_eval(totalcoil_sm,z_cm)
    coil_path.addpath(sm)
    
startz = 10.0#8.75#4.36
number = 8.0

for z_in in my_range(startz,startz+(number-1)*0.75,0.75):
    z_cm = z_in*cm_per_in
    l2 = L2_coil_eval(l2totalcoilcur,z_cm)
    coil_path.addpath(l2)


bounds = np.array([[-20.,20.],[-100.,100.],[-100.,100.]])
 
res = 250.
rseries = np.arange(res)/(res-1) * (bounds[2,1]-bounds[2,0])+bounds[2,0]
zeroseries = np.zeros(res)
    
choiceseries = np.zeros(150) + 4.0624*2.54
rvecz_arr = np.array([zeroseries,zeroseries,rseries])
rvecz = np.transpose(rvecz_arr)
#rvecx = np.transpose([rseries],[zeroseries],[zeroseries])
rvecx_arr = np.array([rseries,zeroseries,zeroseries])    
rvecx = np.transpose(rvecx_arr)

fieldseriesz = field_magnitude(coil_path,rvecz)
fieldserieszvec = field_magnitude_vec(coil_path,rvecz)
zsign = field_magnitude_sign(coil_path,rvecz)

fieldseriesxy = field_magnitude(coil_path,rvecx)
xysign = field_magnitude_sign(coil_path,rvecx)
arbseriesx = np.ones(res)*22.5
rvecarbx_arr = np.array([arbseriesx,zeroseries,rseries])
rvecarbx = np.transpose(rvecarbx_arr)
#fieldseriesarbx = field_magnitude(coil_path,rvecarbx)


arbseriesz = np.ones(res)*22.5
rvecarbz_arr = np.array([rseries,zeroseries,arbseriesz])
rvecarbz = np.transpose(rvecarbz_arr)
#fieldseriesarbz = field_magnitude(coil_path,rvecarbz)

b_arr = np.zeros([res,res])
#aseries1 = np.arange(res)/(res-1)*(40-(-20))+(-40)
#aseriessh = np.arange(res)/(res-1)*(10-(-10))+(-10)
#aseries = np.arange(res)/(res-1)*(30-(-30))+(-30)

aseries1 = np.arange(res)/(res-1)*(100-(-100))+(-100)
aseriessh = np.arange(res)/(res-1)*(100-(-100))+(-100)
aseries = np.arange(res)/(res-1)*(100-(-100))+(-100)
for ii in range(0,np.int(res)):
    seriesz = aseries1
    #for j in range(0,100):
    #if aseries1[ii] < -6.5*cm_per_in:
     #   seriesx = aseriessh#np.arange(100)/(100-1) * (bounds[2,1]-bounds[2,0])+bounds[2,0]
    #else:
    seriesx = np.ones(res)*aseries[ii]
    vecser = np.array([seriesx,zeroseries,seriesz])        
    vecs = np.transpose(vecser)
    b_arr[ii,:] = field_magnitude(coil_path,vecs)

b_arrt = np.transpose(b_arr)

b_vec = np.zeros([res,res,3])
for ii in range(0,np.int(res)):
    seriesz = np.ones(res)*aseries1[ii]
    seriesx = aseries
    vecser = np.array([seriesx,zeroseries,seriesz])
    vecs = np.transpose(vecser)
    b_vec[ii,:,:] = field_magnitude_vec(coil_path,vecs)
    
bx_arr = b_vec[:,:,0]
bz_arr = b_vec[:,:,2]

np.savetxt('bz_arr_500_10.csv',bz_arr,delimiter=',')
np.savetxt('bx_arr_500_10.csv',bx_arr,delimiter=',')


plt.close('all')       

plt.figure(4)
plt.plot(rseries,fieldseriesz)
plt.xlabel('z (cm) at x = 0 cm')
plt.ylabel('abs(B) (G)')

ws = 0
we = rseries.size
fieldserieschop = fieldseriesz[ws:we]
rserieschop = rseries[ws:we]
    
plt.figure(5)
plt.plot(rseries,fieldseriesxy)
plt.xlabel('x (cm) at z = 0 cm')
plt.ylabel('abs(B) (G)')


#with PdfPages('drawings/no_coils/separation_config_v2.pdf') as pdf:
#indent here to save pdf
fig = plt.figure()
ax3 = fig.add_subplot(111)
ax3.plot(bounds[2,:],bounds[1,:],visible=0)
ax3.set_xlabel('z (cm)')
ax3.set_ylabel('x (cm)')
draw_coils(coil_path)
draw_vv()


maxB = np.max(np.abs(fieldseriesz))
bv = field_tracer_plot(coil_path,bounds,np.array([0,0,0]),[2,1],maxB,ax=None)
    #pdf.savefig()
    #plt.close()


#this section of the code plots the stream lines using python's streamplot function.  I've commented it out for now.
"""  
plt.figure()
z0,z1=aseries1[0],aseries1[-1]
r0,r1=aseries[0],aseries[-1]
z = np.linspace(z0,z1,250)
r = np.linspace(r0,r1,250)
z,r= np.meshgrid(z,r)
bzt = np.transpose(bz_arr)
bxt = np.transpose(bx_arr)



plt.streamplot(z,r,bzt,bxt ,color='r',density=1.2)
    
    
#draw_coils(coil_path)
#draw_vv()
    
    plt.xlabel('z (cm)')
    plt.ylabel('x (cm)')
        plt.xlim(-20,20)
        plt.ylim(0,30)
        #plt.axhline(y=30,xmin=-100,xmax=100)

    #plt.axhline(y=15.6*2.54,xmin=-100,xmax=100)
    #plt.axvline(x=15.6*2.54)
    #plt.axhline(y=-5)
    #lt.axhline(y=5)
    #pp= PdfPages('sahhib.pdf')
   
    pdf.savefig()
    plt.close()
"""