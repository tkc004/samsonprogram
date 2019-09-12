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


def gaus(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))

def calsigmaHIlos(G,xlist,ylist,maxlength):
    Gm = G['m']; Gnh = G['nh']; KernalLengths = G['h']; density = G['rho']; Gmetal = G['z']; #all in code unit
    camdict = CAMF.calnH(Gm,Gnh,KernalLengths,density,Gmetal)
    np.amax(camdict['fHI'])
    np.amin(camdict['fHI'])
    fHI = camdict['fHI']
    fH2 = camdict['fH2']
    Gp = G['p'];
    Gx = Gp[:,0]; Gy = Gp[:,1]; Gz = Gp[:,2];
    Gv = G['v']; 
    Gvx = Gv[:,0]; Gvy = Gv[:,1]; Gvz = Gv[:,2];
    sigmaHIlosgrid = np.array([])
    xpos = np.array([])
    ypos = np.array([])
    areagrid = np.array([])
    for i in range(len(xlist)-1):
        for j in range(len(ylist)-1):
            xmax = xlist[i+1]
            xmin = xlist[i]
            ymax = ylist[j+1]
            ymin = ylist[j]
            cutz = (Gz<maxlength/2.0)*(Gz>-maxlength/2.0)
            cutx = (Gx<xmax)*(Gx>xmin)
            cuty = (Gy<ymax)*(Gy>ymin)
            cut = cutx*cuty*cutz
            hist, bin_edges = np.histogram(Gvz[cut], density=True, bins=np.linspace(-50,50,num=100),\
                                           weights=Gm[cut]*fHI[cut])
            x = 0.5*(bin_edges[1:]+bin_edges[:-1]); y = hist;
            sigma=10; mean=0;
            try:
                popt,pcov = curve_fit(gaus,x,y,p0=[0.0175,mean,sigma])
                sigmaHIlos = popt[2]
            except RuntimeError:
                sigmaHIlos = 0.0
            sigmaHIlosgrid = np.append(sigmaHIlosgrid,sigmaHIlos)
            areagrid = np.append(areagrid,np.absolute((xmax-xmin)*(ymax-ymin)))
            HImgrid = np.append(areagrid,np.sum(Gm[cut]*fHI[cut]))
            xpos = np.append(xpos,(xmax+xmin)/2.0)
            ypos = np.append(ypos,(ymax+ymin)/2.0)            
    pos = np.column_stack((xpos,ypos))
    outdata = {'sigmaHIlosgrid':sigmaHIlosgrid, 'areagrid':areagrid, 'pos':pos, 'HImgrid':HImgrid}
    return outdata
    
    
    
    
    
    
    
    
    
    
    
    