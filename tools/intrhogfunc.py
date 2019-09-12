from struct import *
from samson_const import *
from pathloc import *
import sys
import os
import pylab
import numpy as np
import math
from distcalcs import *
from zip_em_all import *
import time
from scipy.interpolate import interp1d
from datetime import date
import scipy.stats as stats
import scipy.optimize as optimize
import scipy.special as special
import errno
import numpy.linalg as la
from mpl_toolkits.mplot3d import Axes3D
from lum_mag_conversions import luminosity_to_magnitude
import Sasha_functions as SF
from readsnap_samson import *
from readsnap_cr import readsnapcr
import matplotlib.pyplot as plt
import colormaps as cmaps
import crtestfunction as CRTF
import readsnipshot as RSS
import samson_functions as SSF
import cameron_functions as CAMF
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp


def intdpfunc(dp,zeropoint,zlist):
    dz = np.absolute(zlist[1]-zlist[0])*kpc_in_cm
    poscut = zlist>zeropoint
    negdp = dp[~poscut]*dz
    posdp = dp[poscut]*dz
    negintdp = np.cumsum(negdp)
    revarr = np.cumsum(posdp[::-1])
    posintdp = revarr[::-1]
    intdp = np.concatenate([negintdp,posintdp])
    return np.absolute(intdp)


def intrhogfunc(rhog,gzlist,zlist):
    dz = np.absolute(zlist[1]-zlist[0])*kpc_in_cm
    poscut =  gzlist>0
    negrhog = rhog[~poscut]*dz
    posrhog = rhog[poscut]*dz
    negintrhog = np.cumsum(negrhog)
    revarr = np.cumsum(posrhog[::-1])
    posintrhog = revarr[::-1]
    intrhog = np.concatenate([negintrhog,posintrhog])
    return np.absolute(intrhog)


def intrhogfrompredata(predata,withinr,maxlength,\
                        usehalfz=0,cutcold=0,vertical=1,\
                       horizontal=0,withoutr=-1.0,dr=1.0,outHI=0):
    data = SSF.calrhofrompredata(predata,withinr,maxlength,\
                                usehalfz=usehalfz,cutcold=cutcold,vertical=vertical,\
                              horizontal=horizontal,withoutr=withoutr,dr=dr,outHI=outHI)
    zlist=data['zlist']; gzlist = data['gzlist']; rhog=data['rhog'];
    intrhog = intrhogfunc(rhog,gzlist,zlist)
    return intrhog

def intrhogtorlist(predata,rlist,maxlength):
    intrhogl = np.array([])
    for i in range(len(rlist)-1):
        withinr=rlist[i+1]; withoutr=rlist[i];
        intrhog = intrhogfrompredata(predata,withinr,maxlength,vertical=1,withoutr=withoutr)
        intrhogl=np.append(intrhogl,np.amax(intrhog))
    intrhogl=np.append(intrhogl,0.0)
    return intrhogl
        
        

        