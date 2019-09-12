from samson_const import *
import matplotlib as mpl
from readsnap_cr import readsnapcr
import Sasha_functions as SF
import graphics_library as GL
import gas_temperature as GT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import samson_functions as SSF
from matplotlib import rcParams
from pylab import *
import cameron_functions as CAMF
import crtestfunction as CRTF


def multapply(testarr,con1,con2):
    testarr0=testarr[con1]
    outarr=testarr0[con2]
    return outarr


def gastrack(runtodo,snaptrack,snapstart,Tcut_t=1.0,highTcut_t=1e10,\
        vcut_t=0.0,vhcut_t=1e10,withinr_t=20.0,zup_t=25.0,zdown_t=20.0,userad=0):
        data_t = CRTF.gaswindphase(runtodo,snaptrack,Tcut=Tcut_t,highTcut=highTcut_t,\
        vcut=vcut_t,vhcut=vhcut_t,withinr=withinr_t,zup=zup_t,zdown=zdown_t,userad=userad)
        data_t = CRTF.gaswindphase(runtodo,snaptrack,Tcut=Tcut_t,highTcut=highTcut_t,\
        vcut=vcut_t,vhcut=vhcut_t,withinr=withinr_t,zup=zup_t,zdown=zdown_t,userad=userad)
        Gid_t = data_t['Gid']
        partZ_t = data_t['partZ']
        #u, c = np.unique(Gid_t, return_counts=True)
        vcut = -1e10 #if we track gas particles, we should consider all velocities, temp, etc
        withinr=1e10
        zup=1e10
        zdown=0.001
        Tcut = 1.0
        highTcut =1e15
        vhcut = 1e15
        data = CRTF.gaswindphase(runtodo,snapstart,Tcut=Tcut,highTcut=highTcut,\
        vcut=vcut,vhcut=vhcut,withinr=withinr,zup=zup,zdown=zdown,userad=userad)
        partX_n = data['partX']
        partY_n = data['partY']
        partZ_n = data['partZ']
        partR_n = data['partR']
        vx_n = data['vx']
        vy_n = data['vy']
        vz_n = data['vz']
        vr_n = data['vr']
        TrueTemp_n = data['TrueTemp']
        converted_rho_n = data['convertedrho']
        Gmass_n = data['Gmass']
        vmax = data['vmax']
        vmin = data['vmin']
        Gid_n = data['Gid']
        Gid_tu, unindex_t = np.unique(Gid_t, return_index=True)
        Gid_nu, repindex_n = np.unique(Gid_n, return_counts=True)
        unindex_n = np.in1d(Gid_n,Gid_nu[repindex_n<2]) #consider only unique id, drop all repeated
        Gidnr=Gid_n[unindex_n]
        idint0=np.in1d(Gidnr,Gid_tu)
        partR=multapply(partR_n,unindex_n,idint0)
        vr=multapply(vr_n,unindex_n,idint0)
        partX=multapply(partX_n,unindex_n,idint0)
        partY=multapply(partY_n,unindex_n,idint0)
        partZ=multapply(partZ_n,unindex_n,idint0)
        vx=multapply(vx_n,unindex_n,idint0)
        vy=multapply(vy_n,unindex_n,idint0)
        vz=multapply(vz_n,unindex_n,idint0)
        TrueTemp=multapply(TrueTemp_n,unindex_n,idint0)
        converted_rho=multapply(converted_rho_n,unindex_n,idint0)
        Gmass=multapply(Gmass_n,unindex_n,idint0)
        Gid=multapply(Gid_n,unindex_n,idint0)
        return {'vmax':vmax,'vmin':vmin,'Gid':Gid,'TrueTemp':TrueTemp,'convertedrho':converted_rho,\
 'vx':vx, 'vy':vy, 'vz':vz, 'vr':vr, 'partX':partX,'partY':partY, 'partZ':partZ, 'partR':partR, 'Gmass':Gmass}