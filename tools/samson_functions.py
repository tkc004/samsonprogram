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
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='inferno', cmap=cmaps.inferno)
plt.register_cmap(name='plasma', cmap=cmaps.plasma)

#from tasz import hubble_param

pi  = np.pi
sin = np.sin
cos = np.cos

def funclinear(x, a, b):
        return a+b*x

def ellipse(u,v):
    x = rx*cos(u)*cos(v)
    y = ry*sin(u)*cos(v)
    z = rz*sin(v)
    return x,y,z


def mvee(points, tol = 0.001):
    """
    Finds the ellipse equation in "center form"
    (x-c).T * A * (x-c) = 1
    """
    N, d = points.shape
    Q = np.column_stack((points, np.ones(N))).T
    err = tol+1.0
    u = np.ones(N)/N
    while err > tol:
        # assert u.sum() == 1 # invariant
        X = np.dot(np.dot(Q, np.diag(u)), Q.T)
        M = np.diag(np.dot(np.dot(Q.T, la.inv(X)), Q))
        jdx = np.argmax(M)
        step_size = (M[jdx]-d-1.0)/((d+1)*(M[jdx]-1.0))
        new_u = (1-step_size)*u
        new_u[jdx] += step_size
        err = la.norm(new_u-u)
        u = new_u
    c = np.dot(u,points)
    A = la.inv(np.dot(np.dot(points.T, np.diag(u)), points)
               - np.multiply.outer(c,c))/d
    return A, c

def timeconvert(ascale, firever=2):
        readtimelist=readtime(firever=firever)
        snap2list=readtimelist['snaplist']
        time2list=readtimelist['timelist']
        a2list=readtimelist['alist']
        timeyr = np.interp(ascale,a2list,time2list)*1e9
        return timeyr

def calLfrompar(Gx,Gy,Gz,Gvx,Gvy,Gvz,Gm,rin=5.):
        Gr = np.sqrt(Gx*Gx+Gy*Gy+Gz*Gz)
        cutr = Gr < rin #kpc
        Gmcutr = Gm[cutr]; 
        Gxcutr = Gx[cutr]; Gycutr = Gy[cutr]; Gzcutr = Gz[cutr];
        Gvxcutr = Gvx[cutr]; Gvycutr = Gvy[cutr]; Gvzcutr = Gvz[cutr];
        print 'Gvx,Gvy,Gvz', np.average(Gvxcutr), np.average(Gvycutr), np.average(Gvzcutr)
        Lang = [0.,0.,0.]
        for i in range(len(Gxcutr)):
                Lang += Gmcutr[i]*np.cross([Gxcutr[i],Gycutr[i],Gzcutr[i]],[Gvxcutr[i],Gvycutr[i],Gvzcutr[i]])
        return Lang


def calLfromrun(runtodo,Nsnap,ptype,rin=5.):
    G = readsnapfromrun(runtodo,Nsnap,ptype,rotface=0,datasup=0)
    Lang = calLfrompar(G['p'][:,0],G['p'][:,1],G['p'][:,2],G['v'][:,0],G['v'][:,1],G['v'][:,2],G['m'],rin=5.)
    return Lang

def calangLfromrun(runtodo,Nsnap,ptype,rin=5.):
    Lang = calLfromrun(runtodo,Nsnap,ptype,rin=5.)
    Lx = Lang[0]; Ly = Lang[1]; Lz = Lang[2];
    theta, phi = SF.rotateL_to_z(Lx,Ly,Lz,Lx,Ly,Lz,outangle=1)
    return theta, phi
    


def readsnapfromrun(runtodo,Nsnap,ptype,rotface=0,datasup=0):
        info=outdirname(runtodo, Nsnap)
        M1speed=info['M1speed']
        rundir=info['rundir']
        runtitle=info['runtitle']
        slabel=info['slabel']
        snlabel=info['snlabel']
        dclabel=info['dclabel']
        resolabel=info['resolabel']
        the_snapdir=info['the_snapdir']
        the_prefix = info['the_prefix']
        the_suffix = info['the_suffix']
        Nsnapstring=info['Nsnapstring']
        Fcal=info['Fcal']
        iavesfr=info['iavesfr']
        timestep=info['timestep']
        h0=info['h0']
        haveB=info['haveB']
        havecr=info['havecr']
        newlabel=info['newlabel']
        cosmo=info['cosmo']
        halostr=info['halostr']
        firever=info['firever']
        maindir=info['maindir']
        usepep=info['usepep']
        snumadd=info['snumadd'] 
        G = readsnapwcen(the_snapdir, Nsnapstring, ptype, snapshot_name=the_prefix, extension=the_suffix,\
         havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
         datasup=datasup,runtodo=runtodo,rundir=rundir,halostr=halostr,firever=firever)
        return G


def readcenfromAHF(the_snapdir, Nsnapstring, ptype, snapshot_name='',\
         extension='', havecr=0,h0=1,cosmo=1,usepep=0,rundir='',halostr='',\
         maindir='',firever=2,snumadd=0):
        header = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=snapshot_name,\
         extension=extension, havecr=havecr,h0=h0,cosmological=cosmo, header_only=1)
        Nsnap = int(Nsnapstring)
        ascale = header['time']
        print 'this time', ascale
        thisred = 1./ascale-1.
        hubble = header['hubble']
        omega_matter = header['omega_matter']
        print 'hubble', hubble
        if cosmo==1:
                if usepep==1:
                        halosingle = SF.read_halo_history_pep(rundir, Nsnap, singlesnap=1, firever=firever,halonostr=halostr, hubble=hubble, comoving=0, maindir=maindir)
                        afactor=ascale
                else:
                        halosingle = SF.read_halo_history(rundir, halonostr=halostr,hubble=hubble, comoving=0, maindir=maindir, singlesnap=1, atime=ascale,snumadd=snumadd)
                        afactor=1.0
                xcen = halosingle['x']*afactor
                ycen = halosingle['y']*afactor
                zcen = halosingle['z']*afactor
                xvcen = halosingle['xv']
                yvcen = halosingle['yv']
                zvcen = halosingle['zv']
                Rvirnow = halosingle['R']
                MgAHF = halosingle['Mg']
        else:
                xcen=0
                ycen=0
                zcen=0
                xvcen=0
                yvcen=0
                zvcen=0
        if cosmo==1:
            return {'halosingle':halosingle,'xcen':xcen, 'ycen':ycen, 'zcen':zcen, 'xvcen':xvcen, 'yvcen':yvcen, 'zvcen':zvcen}
        else:
            return {'xcen':xcen, 'ycen':ycen, 'zcen':zcen, 'xvcen':xvcen, 'yvcen':yvcen, 'zvcen':zvcen}
        
        
    
def readcenfromrun(runtodo,Nsnap,ptype,snapshot_name='snapshot', extension='.hdf5'):
        info=outdirname(runtodo, Nsnap)
        M1speed=info['M1speed']
        rundir=info['rundir']
        runtitle=info['runtitle']
        slabel=info['slabel']
        snlabel=info['snlabel']
        dclabel=info['dclabel']
        resolabel=info['resolabel']
        the_snapdir=info['the_snapdir']
        the_prefix = info['the_prefix']
        the_suffix = info['the_suffix']
        Nsnapstring=info['Nsnapstring']
        Fcal=info['Fcal']
        iavesfr=info['iavesfr']
        timestep=info['timestep']
        h0=info['h0']
        haveB=info['haveB']
        havecr=info['havecr']
        newlabel=info['newlabel']
        cosmo=info['cosmo']
        halostr=info['halostr']
        firever=info['firever']
        maindir=info['maindir']
        usepep=info['usepep']
        snumadd=info['snumadd']        
        cendata=readcenfromAHF(the_snapdir, Nsnapstring, ptype, snapshot_name=snapshot_name,\
         extension=extension, havecr=havecr,h0=h0,cosmo=cosmo,usepep=usepep,rundir=rundir,halostr=halostr,\
         maindir=maindir,firever=firever,snumadd=snumadd)
        return cendata

def readsnapwcen(the_snapdir, Nsnapstring, ptype, snapshot_name='snapshot', extension='.hdf5',\
         havecr=0,h0=0,cosmo=0, usepep=0, maindir='',snumadd=0,rotface=0, datasup=0,runtodo='',\
         rundir='',halostr='',firever=2):           
        cendata = readcenfromAHF(the_snapdir, Nsnapstring, 0, snapshot_name=snapshot_name,\
         extension=extension, havecr=havecr,h0=h0,cosmo=cosmo,rundir=rundir,halostr=halostr,\
         usepep=usepep,maindir=maindir,firever=firever,snumadd=snumadd)
        xcen = cendata['xcen']; ycen = cendata['ycen']; zcen = cendata['zcen'];
        xvcen = cendata['xvcen']; yvcen = cendata['yvcen']; zvcen = cendata['zvcen'];
        Nsnap = int(Nsnapstring)
        G = readsnapcr(the_snapdir, Nsnapstring, ptype, snapshot_name=snapshot_name, extension=extension, havecr=havecr,h0=h0,cosmological=cosmo)
        Gpos = G['p']
        Gvel = G['v']
        Gm = G['m']
        if ptype==0:
                Grho = G['rho']
                Gu = G['u']
                Gm = G['m']
                if havecr>0:
                        cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                else:
                        cregy = np.zeros(len(Gm))
                Neb = G['ne']
        try:
                GB = G['B']
        except KeyError:
                GB = []
        Gx = Gpos[:,0]-xcen
        Gy = Gpos[:,1]-ycen
        Gz = Gpos[:,2]-zcen
        Gvx = Gvel[:,0]-xvcen   #km/s
        Gvy = Gvel[:,1]-yvcen
        Gvz = Gvel[:,2]-zvcen
        print 'xvcen, yvcen, zvcen', xvcen, yvcen, zvcen   
        if datasup==1:
                from crtestfunction import findcennew
                xcen, ycen, zcen = findcennew(runtodo,Nsnap,withinr=5.,dir='x',datasup=datasup,Gx=Gx,Gy=Gy,Gz=Gz,Gm=Gm)
                Gx = Gx-xcen
                Gy = Gy-ycen
                Gz = Gz-zcen
        if rotface==1:
                Lang = calLfrompar(Gx,Gy,Gz,Gvx,Gvy,Gvz,Gm,rin=5.)
                Gx, Gy, Gz = SF.rotateL_to_z(Gx,Gy,Gz,Lang[0],Lang[1],Lang[2])
                Gvx, Gvy, Gvz = SF.rotateL_to_z(Gvx,Gvy,Gvz,Lang[0],Lang[1],Lang[2])
        G['p'][:,0]=Gx; G['p'][:,1]=Gy; G['p'][:,2]=Gz;
        G['v'][:,0]=Gvx; G['v'][:,1]=Gvy; G['v'][:,2]=Gvz;
        return G


# This routine uses relative positions of the particles to locate the center. The basic mechanism is to group particles into rectangular grids, and then fit the high density grid with ellipsoidal. 

#DMrelpos is the positions of the particles and nobin^3 is the number of grid. 

#nopafactor* den is the density cutoff (den is the average density of the particles). 

# !It is better to keep nobin <100 or the code will cost a lot of time and the center will be locating the highest density clump but not based on overall shape 

# ! It is not possible to fit with too many or too few grids (if nopa is too small, the number of grid will be too large and you should tune the nopafactor


def ellipsoidal_centering(DMrelposX,DMrelposY,DMrelposZ,nobin,nopafactor):
        DMco=np.array([DMrelposX,DMrelposY,DMrelposZ])
        DMcor=DMco.T

        DMparno=len(DMrelposX)
        H, edges = np.histogramdd(DMcor, bins = (nobin, nobin, nobin))

        NominS=[]
        ellX=[]
        ellY=[]
        ellZ=[]

        den=float(DMparno)/(pi*4./3.)/float(nobin)/float(nobin)/float(nobin)

        nopa=float(den)*nopafactor

        Dpos=[]
        Dposx=[]
        Dposy=[]
        Dposz=[]
        totalno=0
        for i in range(nobin):
                for j in range(nobin):
                        for k in range(nobin):
                                if (H[i,j,k]>nopa):
                                        Dposx=np.append(Dposx,edges[0][i])
                                        Dposy=np.append(Dposy,edges[1][j])
                                        Dposz=np.append(Dposz,edges[2][k])
                                        totalno+=H[i,j,k]
        if len(Dposx)<4:
                print 'Warning: Density threshold is too high'
                return -1, -1, -1 
        if len(Dposx)>1000:
                print 'Warning: Density threshold is too low and the no of grids is too large'
                return -1, -1, -1

        print 'fitting ell'
        points = np.vstack(((Dposx),(Dposy),(Dposz))).T
        A, centroid = mvee(points)
        DXi=centroid[0]
        DYi=centroid[1]
        DZi=centroid[2]
        return DXi, DYi, DZi

def outdirname(runtodo, Nsnap=500):
        dclabel=''
        resolabel=''
        Fcal=0.1
        iavesfr = 1.0
        snlabel=''
        slabel=''
        runtitle=''
        subdir='hdf5'
        timestep=''
        maindir='scratch'
        cosmo=0
        color='k'
        usesnaplist=0
        withinRv=1
        h0=0.702
        halostr='00'
        usepep=0
        beginno=100
        finalno=600
        snapsep=1
        firever=2
        initsnap=0
        exceptcool=0
        kappa=0
        havemetal=1
        haveB=0
        havecr=0
        Rvir=0
        M1speed=0
        the_prefix='snapshot'
        the_suffix='.hdf5'
        correctIa=0
        stron=0
        galcen=0
        newlabel = ''
        strlabel = ''
        crlabel=''
        Sheaform=0
        snumadd=0
        multifile='n'
        rundir='none'
        highres=0
        Msini=0.0
        Mhini=0.0
        newlabel=''


        if (runtodo=='fm10q'):
                rundir='m10q/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m10q'
                firever=2
                maindir='oasis/extra'
                multifile='y'
                fileno=8
                cosmo=1
                usepep=1

        if (runtodo=='fm12m'):
                rundir='m12m_ref13/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='fm12m'
                firever=2
                multifile='y'
                fileno=4
                maindir='oasis/extra'
                usepep=1
                cosmo=1

        if (runtodo=='fm12b'):
                rundir='m12b_ref12/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m12b'
                firever=2
                maindir='oasis/extra'
                usepep=1
                cosmo=1

        if (runtodo=='m11v4cr_b_70'):
                rundir='/m11v/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                firever=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11v4'
                halostr='04'
                halocolor='k'
                highres=2
                usepep=1
                crlabel='cr70'
                cosmo=1
                newlabel=r'MHD $\kappa$=3e28'

        if (runtodo=='m11v2cr_b_70'):
                rundir='/m11v/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                firever=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11v2'
                halostr='02'
                halocolor='k'
                highres=2
                usepep=1
                crlabel='cr70'
                cosmo=1
                newlabel=r'MHD $\kappa$=3e28'

        if (runtodo=='m11v3cr_b_70'):
                rundir='/m11v/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                firever=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11v'
                halostr='03'
                halocolor='k'
                highres=2
                usepep=1
                crlabel='cr70'
                cosmo=1
                newlabel=r'MHD $\kappa$=3e28'

        if (runtodo=='m11v1cr_b_70'):
                rundir='/m11v/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                firever=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11v1'
                halostr='01'
                halocolor='k'
                highres=2
                usepep=1
                crlabel='cr70'
                cosmo=1
                newlabel=r'MHD $\kappa$=3e28'

        if (runtodo=='m11vcr_b_70'):
                rundir='/m11v/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                firever=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11v'
                halostr='00'
                halocolor='k'
                highres=2
                usepep=1
                crlabel='cr70'
                cosmo=1
                newlabel=r'MHD $\kappa$=3e28'

        if (runtodo=='f573'):
                rundir='FIRE_2_0_or_h573_criden1000_noaddm_sggs/'
        #       halostr='00'
                halostr='00'
                subdir='output'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                halocolor='g'
                labelname='m10z'
                icolor=0.01
                firever=2
                maindir='oasis'
                cosmo=1
                dclabel='m10z hydro'
                highres=0

        if (runtodo=='f553'):
                rundir='FIRE_2_0_or_h553_criden1000_noaddm_sggs/'
                halostr='00'
                subdir='output'
                xmax_of_box=40.0
                beginno=100
                finalno=600
                halocolor='r'
                labelname='m11a'
                icolor=0.34
                firever=2
                maindir='oasis'
                usepep=0
                cosmo=1
                dclabel='m11a hydro'
                highres=0


        if (runtodo=='fm10v'):
                rundir='m10v/'
                halostr='00'
                beginno=600
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                multifile='y'
                labelname='m10v'
                firever=2
                maindir='oasis/extra'
                usepep=0
                fileno=8
                cosmo=1
                dclabel='m10v hydro'
                highres=0


        if (runtodo=='fm11v'):
                rundir='m11v/'
                halostr='00'
                beginno=600
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m11v'
                firever=2
                maindir='oasis/extra'
                usepep=1
                cosmo=1
                highres=0

        if (runtodo=='fm11v1'):
                rundir='m11v1/'
                halostr='01'
                beginno=600
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m11v1'
                firever=2
                maindir='oasis/extra'
                usepep=1
                cosmo=1
                highres=0

        if (runtodo=='fm11v2'):
                rundir='m11v2/'
                halostr='02'
                beginno=600
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m11v2'
                firever=2
                maindir='oasis/extra'
                usepep=1
                cosmo=1
                highres=0

        if (runtodo=='fm11v3'):
                rundir='m11v3/'
                halostr='03'
                beginno=600
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m11v3'
                firever=2
                maindir='oasis/extra'
                usepep=1
                cosmo=1
                highres=0

        if (runtodo=='fm11v4'):
                rundir='m11v4/'
                halostr='01'
                beginno=600
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m11v4'
                firever=2
                maindir='oasis/extra'
                usepep=1
                cosmo=1
                highres=0

        if (runtodo=='fm12f'):
                rundir='m12f_ref13/'
                halostr='00'
                beginno=600
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m12fhr'
                firever=2
                multifile='y'
                fileno=4
                maindir='oasis/extra'
                usepep=1
                cosmo=1
                highres=0

        if (runtodo=='fm12fmd'):
                rundir='m12f_res7100/'
                halostr='00'
                beginno=600
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m12f'
                firever=2
                multifile='y'
                fileno=4
                maindir='oasis/metal_diffusion'
                usepep=1
                cosmo=1
                Sheaform=1
                highres=4



        if (runtodo=='fm11d'):
                rundir='m11d_res7000/'
                halostr='01'
                beginno=600
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m11e'
                firever=2
                h0=0.68
                maindir='oasis'
                cosmo=1
                highres=0

        if (runtodo=='fm11dmd'):
                rundir='m11d_res7100/'
                halostr='01'
                beginno=600
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m11emd'
                firever=2
                h0=0.68
                maindir='oasis/metal_diffusion'
                cosmo=1
                highres=4

        if (runtodo=='f383'):
                rundir='FIRE_2_0_or_h383_criden1000_noaddm_sggs/'
                halostr='00'
                subdir='output'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                halocolor='y'
                icolor=0.85
                labelname='m11c'
                firever=2
                maindir='oasis'
                cosmo=1
                usepep=0
                highres=0

        if (runtodo=='fm11q'):
                rundir='m11q/'
                halostr='00'
                beginno=600
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m11q'
                firever=2
                maindir='oasis/extra'
                multifile='y'
                fileno=4
                usepep=1
                cosmo=1
                highres=0

        if (runtodo=='fm11qmd'):
                rundir='m11q_res7100/'
                halostr='00'
                beginno=0
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='fm11qmd'
                firever=2
                maindir='oasis/metal_diffusion'
                multifile='y'
                fileno=4
                usepep=1
                Sheaform=1
                cosmo=1
                highres=4


        if (runtodo=='fm11qhrmd'):
                rundir='m11q_res880/'
                halostr='00'
                beginno=0
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='fm11qhr'
                firever=2
                maindir='oasis/metal_diffusion'
                multifile='y'
                fileno=4
                usepep=1
                cosmo=1
                highres=4


        if (runtodo=='f1146'):
                rundir='FIRE_2_0_or_h1146_criden1000_noaddm_sggs/'
                halostr='01'
                beginno=600
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m10y'
                firever=2
                maindir='oasis'
                cosmo=1
                highres=0

        if (runtodo=='f61'):
                rundir='FIRE_2_0_or_h61/'
                halostr='00'
                beginno=600
                finalno=600
                subdir='output'
                xmax_of_box=60.0
                halocolor='b'
                labelname='m11f'
                firever=2
                maindir='oasis'
                multifile='y'
                icolor=0.93
                cosmo=1
                highres=0


        if (runtodo=='f46'):
                rundir='FIRE_2_0_or_h46/'
                halostr='00'
                beginno=600
                finalno=600
                subdir='output'
                xmax_of_box=60.0
                halocolor='b'
                labelname='m11g'
                firever=2
                maindir='oasis'
                multifile='y'
                icolor=0.93
                cosmo=1
                dclabel='m11g hydro'
                crlabel='noCR'
                highres=0


        if (runtodo=='fm12i'):
                rundir='m12i_res7000/'
                halostr='00'
                beginno=600
                finalno=600
                subdir='output'
                xmax_of_box=60.0
                halocolor='b'
                labelname='m12i'
                firever=2
                maindir='oasis'
                multifile='y'
                icolor=0.93
                cosmo=1
                crlabel='noCR'
                h0=0.68
                dclabel='m12i hydro'
                snumadd=1
                highres=0



        if (runtodo=='mw_cr_lr_dc28_1_23_17_test16c_bridges'):
                rundir='mw_cr_lr_dc28_1_23_17_test16c_bridges'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='b'
        if (runtodo=='mw_cr_lr_dc28_1_23_17_test16chole_bridges'):
                rundir='mw_cr_lr_dc28_1_23_17_test16chole_bridges'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='g'
        if (runtodo=='mw_cr_lr_dc28_1_23_17_test16cout_bridges'):
                rundir='mw_cr_lr_dc28_1_23_17_test16cout_bridges'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='r'
        if (runtodo=='mw_cr_lr_dc28_1_23_17_testhole'):
                rundir='mw_cr_lr_dc28_1_23_17_testhole'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc29_1_23_17_test6'):
                rundir='mw_cr_lr_dc29_1_23_17_test6'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='b'
        if (runtodo=='mw_cr_lr_dc29_6_6'):
                rundir='mw_cr_lr_dc29_6_6'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='b'

        if (runtodo=='mw_cr_lr_dc29_1_23_17_test6chole_bridges'):
                rundir='mw_cr_lr_dc29_1_23_17_test6chole_bridges'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='g'
        if (runtodo=='mw_cr_lr_dc28_1_23_17_test6'):
                rundir='mw_cr_lr_dc28_1_23_17_test6'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_3_31_M1'):
                rundir='mw_cr_lr_dc28_3_31_M1'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_4_19_M1'):
                rundir='mw_cr_lr_dc28_4_19_M1'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_4_19_M1_equaltimestep'):
                rundir='mw_cr_lr_dc28_4_19_M1_equaltimestep'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc29_4_19_M1'):
                rundir='mw_cr_lr_dc29_4_19_M1'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_5_10_M1'):
                rundir='mw_cr_lr_dc28_5_10_M1'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_5_27_M1_equaltimestep'):
                rundir='mw_cr_lr_dc28_5_27_M1_equaltimestep'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_5_18_equaltimestep'):
                rundir='mw_cr_lr_dc28_5_18_equaltimestep'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc26_5_27_M1'):
                rundir='mw_cr_lr_dc26_5_27_M1'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e26'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc27_5_27_ets'):
                rundir='mw_cr_lr_dc27_5_27_ets'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc27_hole_5_27_ets'):
                rundir='mw_cr_lr_dc27_hole_5_27_ets'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc29_5_27_M1'):
                rundir='mw_cr_lr_dc29_5_27_M1'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_5_27_M1'):
                rundir='mw_cr_lr_dc28_5_27_M1'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_5_27_M1_c5000'):
                rundir='mw_cr_lr_dc28_5_27_M1_c5000'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_5_27_M1_c1000'):
                rundir='mw_cr_lr_dc28_5_27_M1_c1000'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_5_27_M1_c100'):
                rundir='mw_cr_lr_dc28_5_27_M1_c100'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_5_27_M1_relaxfluxlimit'):
                rundir='mw_cr_lr_dc28_5_27_M1_relaxfluxlimit'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_5_18_M1'):
                rundir='mw_cr_lr_dc28_5_18_M1'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_5_18'):
                rundir='mw_cr_lr_dc28_5_18'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_5_10_M1_modad'):
                rundir='mw_cr_lr_dc28_5_10_M1_modad'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_5_10'):
                rundir='mw_cr_lr_dc28_5_10'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_5_10_equaltimestep'):
                rundir='mw_cr_lr_dc28_5_10_equaltimestep'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_5_10_M1_equaltimestep'):
                rundir='mw_cr_lr_dc28_5_10_M1_equaltimestep'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_5_18_M1_equaltimestep'):
                rundir='mw_cr_lr_dc28_5_18_M1_equaltimestep'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_5_10_M1_equaltimestep_modad'):
                rundir='mw_cr_lr_dc28_5_10_M1_equaltimestep_modad'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_3_31_M1_hole'):
                rundir='mw_cr_lr_dc28_3_31_M1_hole'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_3_31_M1_equaltimestep'):
                rundir='mw_cr_lr_dc28_3_31_M1_equaltimestep'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_5_27_equaltimestep'):
                rundir='mw_cr_lr_dc28_5_27_equaltimestep'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc27_6_6_ets'):
                rundir='mw_cr_lr_dc27_6_6_ets'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_6_6_ets'):
                rundir='mw_cr_lr_dc28_6_6_ets'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc27_hole_6_6_ets'):
                rundir='mw_cr_lr_dc27_hole_6_6_ets'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_hole_6_6_ets'):
                rundir='mw_cr_lr_dc28_hole_6_6_ets'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_6_15_ets'):
                rundir='mw_cr_lr_dc28_6_15_ets'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 ETS'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='g'
        if (runtodo=='mw_cr_lr_dc28_6_15_ets_freeze'):
                rundir='mw_cr_lr_dc28_6_15_ets_freeze'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 ETS'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='g'
        if (runtodo=='mw_cr_lr_dc28_6_15_ets_freeze_nocooling'):
                rundir='mw_cr_lr_dc28_6_15_ets_freeze_nocooling'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 ETS'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='g'
                exceptcool=1
        if (runtodo=='mw_cr_lr_dc29_6_15_ets'):
                rundir='mw_cr_lr_dc29_6_15_ets'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29 ETS'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='g'

        if (runtodo=='mw_cr_lr_dc28_6_15_ets_fhfcr'):
                rundir='mw_cr_lr_dc28_6_15_ets_fhfcr'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 ETS freeze hydro'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='g'

        if (runtodo=='mw_cr_lr_dc29_6_15_ets_fhfcr'):
                rundir='mw_cr_lr_dc29_6_15_ets_fhfcr'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29 ETS freeze hydro'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='g'






        if (runtodo=='mw_cr_lr_dc28_6_15_ets_fhfcr_t1'):
                rundir='mw_cr_lr_dc28_6_15_ets_fhfcr_t1'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 ETS freeze hydro'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'

        if (runtodo=='mw_cr_lr_nodiff_6_15'):
                rundir='mw_cr_lr_nodiff_6_15'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='0'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='g'
        if (runtodo=='mw_cr_lr_dc27_6_15_ets'):
                rundir='mw_cr_lr_dc27_6_15_ets'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='g'
        if (runtodo=='mw_cr_lr_dc28_hole_6_15_ets'):
                rundir='mw_cr_lr_dc28_hole_6_15_ets'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 with hole'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='b'
        if (runtodo=='mw_cr_lr_dc28_6_15_M1'):
                rundir='mw_cr_lr_dc28_6_15_M1'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='r'



        if (runtodo=='mw_cr_lr_dc28_6_15_M1_ca'):
                rundir='mw_cr_lr_dc28_6_15_M1_ca'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='r'
        if (runtodo=='mw_cr_lr_dc27_6_15_M1_nofluxterm'):
                rundir='mw_cr_lr_dc27_6_15_M1_nofluxterm'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27 M1 c500 no fluxterm'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'

        if (runtodo=='mw_cr_lr_dc27_2_10_M1'):
                rundir='mw_cr_lr_dc27_2_10_M1'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27 M1 c500 2/10'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'

        if (runtodo=='M1dc26'):
                rundir='CRdifftest/M1dc26'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e26 M1 c500'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'

        if (runtodo=='dc26'):
                rundir='CRdifftest/dc26'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e26'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'


        if (runtodo=='dc27_2_21'):
                rundir='CRdifftest/dc27_2_21'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27'
                resolabel='lr'
                Fcal=0.1
                kappa=10
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'

        if (runtodo=='dc27g'):
                rundir='bw/CRdifftest/dc27g'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='Diff'
                resolabel='lr'
                Fcal=0.1
                kappa=10
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'



        if (runtodo=='dc28g'):
                rundir='bw/CRdifftest/dc28g'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='Diff'
                resolabel='lr'
                Fcal=0.1
                kappa=100
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'


        if (runtodo=='M1dc27_2_21'):
                rundir='CRdifftest/M1dc27_2_21'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27 M1 c500'
                resolabel='lr'
                Fcal=0.1
                kappa=10
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'

        if (runtodo=='M1dc27_2_21_smallts'):
                rundir='CRdifftest/M1dc27_2_21_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'$\tilde{c}=500$'
                resolabel='lr'
                Fcal=0.1
                kappa=10
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'

        if (runtodo=='M1dc27_r128_2_21'):
                rundir='CRdifftest/M1dc27_r128_2_21'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27 M1 c500 r128'
                resolabel='lr'
                Fcal=0.1
                kappa=10
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'


        if (runtodo=='M1dc27_2_21_c100'):
                rundir='CRdifftest/M1dc27_2_21_c100'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27 M1 c100'
                resolabel='lr'
                Fcal=0.1
                kappa=10
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'


        if (runtodo=='M1dc27_2_21_c100_smallts'):
                rundir='CRdifftest/M1dc27_2_21_c100_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'$\tilde{c}=100$'
                resolabel='lr'
                Fcal=0.1
                kappa=10
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'

        if (runtodo=='M1dc28_2_21_c1000'):
                rundir='CRdifftest/M1dc28_2_21_c1000'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1 c1000'
                resolabel='lr'
                Fcal=0.1
                kappa=100
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'

        if (runtodo=='M1dc28_2_21_c1000_smallts'):
                rundir='CRdifftest/M1dc28_2_21_c1000_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1 c1000 small timestep'
                resolabel='lr'
                Fcal=0.1
                kappa=100
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'


        if (runtodo=='M1dc28_2_21_c2000_smallts'):
                rundir='CRdifftest/M1dc28_2_21_c2000_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1 c2000 small timestep'
                resolabel='lr'
                Fcal=0.1
                kappa=100
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'





        if (runtodo=='M1dc28_2_21_c4000'):
                rundir='CRdifftest/M1dc28_2_21_c4000'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1 c4000'
                resolabel='lr'
                Fcal=0.1
                kappa=100
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'

        if (runtodo=='M1dc28_2_21_smallts'):
                rundir='CRdifftest/M1dc28_2_21_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1 c500 small timestep'
                resolabel='lr'
                Fcal=0.1
                kappa=100
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'





        if (runtodo=='M1dc27'):
                rundir='CRdifftest/M1dc27'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27 M1 c500'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'


        if (runtodo=='M1dc27g_2_21_c100'):
                rundir='CRdifftest/M1dc27g_2_21_c100'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27 M1 c100 gaussian'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'
                kappa=10


        if (runtodo=='M1dc27g_2_21_smallts'):
                rundir='bw/CRdifftest/M1dc27g_2_21_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='M1=500'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'
                kappa=10


        if (runtodo=='M1dc27g_2_21_c100_smallts'):
                rundir='bw/CRdifftest/M1dc27g_2_21_c100_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='M1=100'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'
                kappa=10

        if (runtodo=='M1dc27g_2_21_c1000_smallts'):
                rundir='bw/CRdifftest/M1dc27g_2_21_c1000_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='M1=1000'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'
                kappa=10

        if (runtodo=='M1strg_2_21_c500_smallts'):
                rundir='bw/CRdifftest/M1strg_2_21_c500_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='M1=500'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'
                kappa=0


        if (runtodo=='M1dc28g_2_21_smallts'):
                rundir='CRdifftest/M1dc28g_2_21_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1 c500 gaussian'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'
                kappa=100


        if (runtodo=='M1dc29g_2_21_smallts'):
                rundir='CRdifftest/M1dc29g_2_21_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29 M1 c500 gaussian'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'
                kappa=1000



        if (runtodo=='M1dc28g_2_21_c100_smallts'):
                rundir='CRdifftest/M1dc28g_2_21_c100_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1 c100 gaussian'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'
                kappa=100


        if (runtodo=='M1dc29g_2_21_smallts_cutcr'):
                rundir='CRdifftest/M1dc29g_2_21_smallts_cutcr'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29 M1 c500 gaussian'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'
                kappa=1000

        if (runtodo=='M1dc28g_2_21_smallts_cutcr'):
                rundir='CRdifftest/M1dc28g_2_21_smallts_cutcr'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1 c500 gaussian'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'
                kappa=100



        if (runtodo=='M1dc29g_2_21_smallts'):
                rundir='CRdifftest/M1dc29g_2_21_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='M1=500'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'
                kappa=1000

        if (runtodo=='M1dc29g_2_21_c1000_smallts'):
                rundir='CRdifftest/M1dc29g_2_21_c1000_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='M1=1000'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'
                kappa=1000

        if (runtodo=='M1dc29g_2_21_c2000_smallts'):
                rundir='CRdifftest/M1dc29g_2_21_c2000_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='M1=2000'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'
                kappa=1000


        if (runtodo=='M1dc29g_2_21_smallts_gassphere'):
                rundir='CRdifftest/M1dc29g_2_21_smallts_gassphere'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29 M1 c500 gaussian w cooling'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'
                kappa=1000
                havemetal=0

        if (runtodo=='M1dc29g_2_21_c1000_smallts_gassphere'):
                rundir='CRdifftest/M1dc29g_2_21_c1000_smallts_gassphere'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29 M1 c1000 gaussian w cooling'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'
                kappa=1000
                havemetal=0


        if (runtodo=='M1dc28g_2_21_smallts_gassphere'):
                rundir='CRdifftest/M1dc28g_2_21_smallts_gassphere'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1 c500 gaussian w cooling'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'
                kappa=100
                havemetal=0



        if (runtodo=='M1dc28g_2_21_c100_smallts_gassphere'):
                rundir='CRdifftest/M1dc28g_2_21_c100_smallts_gassphere'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='M1=100'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'
                kappa=100
                havemetal=0

        if (runtodo=='M1dc28g_2_21_c1000_smallts_gassphere'):
                rundir='CRdifftest/M1dc28g_2_21_c1000_smallts_gassphere'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='M1=1000'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'
                kappa=100
                havemetal=0



        if (runtodo=='dc27'):
                rundir='CRdifftest/dc27'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='1.5e27'
                resolabel='lr'
                kappa=5
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'


        if (runtodo=='dc27_slowc'):
                rundir='CRdifftest/dc27_slowc'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='1.5e27 cmax=100'
                resolabel='lr'
                kappa=5
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'

        if (runtodo=='dc28'):
                rundir='CRdifftest/dc28'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='1.5e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'
                kappa=50

        if (runtodo=='dc28_2_21'):
                rundir='CRdifftest/dc28_2_21'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                kappa=100
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'

        if (runtodo=='M1dc28'):
                rundir='CRdifftest/M1dc28'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1 c500'
                resolabel='lr'
                Fcal=0.1
                kappa=100
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'

        if (runtodo=='M1dc28_2_21'):
                rundir='CRdifftest/M1dc28_2_21'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1 c500'
                kappa=100
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'


        if (runtodo=='mw_cr_lr_dc28_6_15_M1_nofluxterm'):
                rundir='mw_cr_lr_dc28_6_15_M1_nofluxterm'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1 c500 no fluxterm'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'
        if (runtodo=='mw_cr_lr_dc29_6_15_M1_nofluxterm'):
                rundir='mw_cr_lr_dc29_6_15_M1_nofluxterm'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29 M1 c500 no fluxterm'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='r'
                color='k'
        if (runtodo=='mw_cr_lr_dc29_6_15_c5000_M1_nofluxterm'):
                rundir='mw_cr_lr_dc29_6_15_c5000_M1_nofluxterm'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29 M1 c5000 no fluxterm'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='r'
        if (runtodo=='mw_cr_lr_dc29_6_15_c2000_M1_nofluxterm'):
                rundir='mw_cr_lr_dc29_6_15_c2000_M1_nofluxterm'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29 M1 c2000 no fluxterm'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='b'
        if (runtodo=='mw_cr_lr_dc28_6_15_M1_c5000'):
                rundir='mw_cr_lr_dc28_6_15_M1_c5000'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='r'
        if (runtodo=='mw_cr_lr_dc28_6_15_M1_c5000_nofluxterm'):
                rundir='mw_cr_lr_dc28_6_15_M1_c5000_nofluxterm'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1 c5000 no fluxterm'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='r'

        if (runtodo=='mw_cr_lr_dc28_6_15_M1_freezehydro_freecr'):
                rundir='mw_cr_lr_dc28_6_15_M1_freezehydro_freecr'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1 freeze hydro'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'

        if (runtodo=='mw_cr_lr_dc28_6_15_M1_fhfcr_ets'):
                rundir='mw_cr_lr_dc28_6_15_M1_fhfcr_ets'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1 freeze hydro'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'


        if (runtodo=='mw_cr_lr_dc27_6_15_M1_fhfcr_ets'):
                rundir='mw_cr_lr_dc27_6_15_M1_fhfcr_ets'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27 M1 freeze hydro'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'

        if (runtodo=='mw_cr_lr_dc28_6_15_M1_fhfcr_ets_relaxHll'):
                rundir='mw_cr_lr_dc28_6_15_M1_fhfcr_ets_relaxHll'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1 freeze hydro relax'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'


        if (runtodo=='mw_cr_lr_dc29_6_15_M1_fhfcr_ets'):
                rundir='mw_cr_lr_dc29_6_15_M1_fhfcr_ets'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29 M1 freeze hydro'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'

        if (runtodo=='mw_cr_lr_dc27_6_15_ets_fhfcr'):
                rundir='mw_cr_lr_dc27_6_15_ets_fhfcr'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27 freeze hydro'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'



        if (runtodo=='mw_cr_lr_dc27_2_10_M1_fhfcr_ets'):
                rundir='mw_cr_lr_dc27_2_10_M1_fhfcr_ets'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27 M1 freeze hydro'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'


        if (runtodo=='mw_cr_lr_dc27_2_19_M1_fhfcr_ets'):
                rundir='mw_cr_lr_dc27_2_19_M1_fhfcr_ets'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27 M1 freeze hydro'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'


        if (runtodo=='mw_cr_lr_dc27_2_21_M1'):
                rundir='mw_cr_lr_dc27_2_21_M1'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27 M1'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=10
                color='k'

        if (runtodo=='mw_cr_lr_dc27_2_21_M1_smallts'):
                rundir='mw_cr_lr_dc27_2_21_M1_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC27}$'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=10
                color='b'
                Rvir=234
                M1speed=500


        if (runtodo=='mw_cr_lr_dc27_2_21_M1_smallts_ets'):
                rundir='mw_cr_lr_dc27_2_21_M1_smallts_ets'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27 M1 small time step ETS'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=10
                color='k'



        if (runtodo=='mw_cr_lr_dc27_2_21_M1_ssts'):
                rundir='mw_cr_lr_dc27_2_21_M1_ssts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27 M1 very small time step'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=10
                color='k'


        if (runtodo=='mw_cr_lr_dc27_2_21_M1_fhfcr_ets'):
                rundir='mw_cr_lr_dc27_2_21_M1_fhfcr_ets'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27 M1 freeze hydro'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=10
                color='k'

        if (runtodo=='mw_cr_lr_dc27_2_21_M1_fhfcr_ets_c100'):
                rundir='mw_cr_lr_dc27_2_21_M1_fhfcr_ets_c100'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27 M1 c100 freeze hydro'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=10
                color='k'


        if (runtodo=='mw_cr_lr_dc27_2_21_M1_fhfcr_ets_smallts'):
                rundir='mw_cr_lr_dc27_2_21_M1_fhfcr_ets_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='kappa=3e27 M1 freeze hydro small time step'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=10
                color='k'


        if (runtodo=='mw_cr_lr_dc27_2_21_M1_fhfcr_ets_ssts'):
                rundir='mw_cr_lr_dc27_2_21_M1_fhfcr_ets_ssts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='kappa=3e27 M1 freeze hydro very small time step'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=10
                color='k'


        if (runtodo=='mw_cr_lr_dc28_2_21_M1_fhfcr_ets'):
                rundir='mw_cr_lr_dc28_2_21_M1_fhfcr_ets'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1 freeze hydro'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=10
                color='k'

        if (runtodo=='mw_cr_lr_dc28_2_21_M1'):
                rundir='mw_cr_lr_dc28_2_21_M1'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='kappa=3e28 M1 full hydro'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=10
                color='k'
                Rvir=234

        if (runtodo=='mw_cr_lr_dc28_2_21_M1_smallts'):
                rundir='mw_cr_lr_dc28_2_21_M1_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='kappa=3e28 M1 c500 full hydro Courant'
                newlabel=r'$\kappa$=3e28'
                resolabel='lr'
                Fcal=0.1
                M1speed=500
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=100
                color='r'
                Rvir=234

        if (runtodo=='mw_cr_lr_dc28_2_21_M1_nomaxCR'):
                rundir='mw_cr_lr_dc28_2_21_M1_nomaxCR'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='kappa=3e28 M1 c500 no max CR'
                resolabel='lr'
                Fcal=0.1
                M1speed=500
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=100
                color='b'
                Rvir=234



        if (runtodo=='mw_cr_lr_dc28_2_21_M1_smallts_testsn'):
                rundir='mw_cr_lr_dc28_2_21_M1_smallts_testsn'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1 c500 full hydro Courant testsn'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=10
                color='g'
                Rvir=234


        if (runtodo=='mw_cr_lr_dc29_2_21_M1'):
                rundir='mw_cr_lr_dc29_2_21_M1'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29 M1 full hydro'
                M1speed=500
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=10
                color='y'
                Rvir=234



        if (runtodo=='mw_cr_lr_dc29_2_21_M1_nomaxCR'):
                rundir='mw_cr_lr_dc29_2_21_M1_nomaxCR'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29 M1 nomaxCR'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=10
                color='g'
                Rvir=234

        if (runtodo=='mw_cr_lr_dc29_2_21_M1_c1000'):
                rundir='mw_cr_lr_dc29_2_21_M1_c1000'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29 M1 c1000 full hydro'
                M1speed=1000
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=10
                color='k'
                Rvir=234

        if (runtodo=='mw_cr_lr_dc29_2_21_M1_c1000_smallts'):
                rundir='mw_cr_lr_dc29_2_21_M1_c1000_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29 M1 c1000 full hydro Courant'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=10
                color='y'
                Rvir=234
                M1speed=1000

        if (runtodo=='bwmwlrmhd'):
                rundir='mw_lr_2_21_mhd'
                slabel='no CR'
                havecr=0
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf BnoCR}$'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=0
                #color='k'
                #TK test
                color=cmaps.plasma(0.01)
                Rvir=234
                haveB=1


        if (runtodo=='bwmwmrmhd'):
                rundir='mw_mr_2_21_mhd'
                slabel='no CR'
                havecr=0
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf BnoCR}$'
                newlabel = 'MHD          no CR'
                resolabel='mr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=0
                #color='k'
                #TK test
                color=cmaps.plasma(0.01)
                Rvir=234
                haveB=1
                highres=1


        if (runtodo=='bwmwlrstr'):
                rundir='mw_cr_lr_2_21_purestream'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf STREAM}$'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=0
                #color='c'
                #TK test
                color=cmaps.plasma(0.93)
                Rvir=234
                M1speed=500
                haveB=1
                stron=1


        if (runtodo=='bwmwlrstrc1000'):
                rundir='mw_cr_lr_2_21_purestream_c1000'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf STREAM}$'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=0
                color='c'
                Rvir=234
                M1speed=1000
                haveB=1
                stron=1


        if (runtodo=='bwmwmrstr'):
                rundir='mw_cr_mr_2_21_purestream'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf STREAM}$'
                newlabel='MHD          Streaming'
                resolabel='mr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=0
                #color='c'
                #TK test
                color='b'
                Rvir=234
                M1speed=500
                haveB=1
                stron=1

                
        if (runtodo=='bwmwmrstreve'):
                rundir='mw_cr_mr_2_21_str_ll_va_eve'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf STREAM}$'
                newlabel='MHD          Streaming'
                resolabel='mr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=0
                #color='c'
                #TK test
                color=cmaps.plasma(0.93)
                Rvir=234
                M1speed=500
                haveB=1
                stron=1


        if (runtodo=='bwmwmrstrts'):
                rundir='mw_cr_mr_2_21_str_ll_va_eve_testts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf STREAM}$'
                newlabel='MHD          Streaming'
                resolabel='mr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=0
                #color='c'
                #TK test
                color=cmaps.plasma(0.93)
                Rvir=234
                M1speed=500
                haveB=1
                stron=1


        if (runtodo=='bwmwmrstrll4va'):
                rundir='mw_cr_mr_2_21_str_ll_4va'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf STREAM}$'
                newlabel='MHD          Streaming'
                resolabel='mr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=0
                color='c'
                Rvir=234
                M1speed=500
                haveB=1
                stron=1


        if (runtodo=='bwmwlrstrll'):
                rundir='mw_cr_lr_2_21_str_lowloss'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf STREAMLL}$'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=0
                color='purple'
                Rvir=234
                M1speed=500
                haveB=1

        if (runtodo=='mwlrstrllva'):
                rundir='mw_cr_lr_2_21_str_ll_va'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf STREAMLLva}$'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=0
                color='navy'
                Rvir=234
                M1speed=500
                haveB=1

        if (runtodo=='mwlrstrll4va'):
                rundir='mw_cr_lr_2_21_str_ll_4va'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf STREAMLL4va}$'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=0
                color='brown'
                Rvir=234
                M1speed=500
                haveB=1


        if (runtodo=='bwmwlrdc28mhd'):
                rundir='mw_cr_lr_dc28_2_21_M1_mhd_c1000'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf BDC28}$'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=10
                #color='r'
                #TK test
                color=cmaps.plasma(0.7)
                Rvir=234
                M1speed=1000
                haveB=1



        if (runtodo=='bwmwlrdc28mhdref'):
                rundir='mw_cr_lr_dc28_2_21_M1_mhd_c1000_referee'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf BDC28}$'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/refereeruns/'
                kappa=10
                color='k'
                Rvir=234
                M1speed=1000
                haveB=1


        if (runtodo=='bwmwmrdc28mhd'):
                rundir='mw_cr_mr_dc28_2_21_M1_mhd_c1000'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf BDC28}$'
                newlabel=r'MHD $\kappa$=3e28'
                resolabel='mr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=10
                #color='r'
                #TK test
                color=cmaps.plasma(0.7)
                Rvir=234
                M1speed=1000
                haveB=1


        if (runtodo=='bwmwlrdc28str'):
                rundir='mw_cr_lr_dc28_2_21_M1_mhd_stream_c1000'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf BDC28STR}$'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=10
                #color='m'
                #TK test
                color=cmaps.plasma(0.85)
                Rvir=234
                M1speed=1000
                haveB=1
                stron=1

        if (runtodo=='bwmwlrdc29'):
                rundir='mw_cr_lr_dc29_2_21_M1_c2000'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC29}$'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=10
                #color='y'
                #TK test
                color=cmaps.plasma(0.85)
                Rvir=234
                M1speed=2000

        if (runtodo=='bwmwmrdc29'):
                rundir='mw_cr_mr_dc29_2_21_M1_c2000'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC29}$'
                newlabel=r'$\kappa$=3e29'       
                resolabel='mr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=10
                #color='y'
                #TK test
                color=cmaps.plasma(0.85)
                Rvir=234
                M1speed=2000
                highres=3



        if (runtodo=='bwmwlrdc29cutcr'):
                rundir='mw_cr_lr_dc29_2_21_M1_c2000_cutcr'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC29cutcr}$'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=10
                color='k'
                Rvir=234
                M1speed=2000


        if (runtodo=='mw_cr_lr_dc29_2_21_M1_c2000_smallts'):
                rundir='mw_cr_lr_dc29_2_21_M1_c2000_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC29}$'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=10
                color='y'
                Rvir=234
                M1speed=2000

        if (runtodo=='mw_cr_lr_dc29_2_21_M1_c4000_smallts'):
                rundir='mw_cr_lr_dc29_2_21_M1_c4000_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='kappa=3e29 M1 c4000 full hydro Courant'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=10
                color='y'
                Rvir=234
                M1speed=4000


        if (runtodo=='mw_cr_lr_dc29_2_21_M1_smallts'):
                rundir='mw_cr_lr_dc29_2_21_M1_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29 M1 c500 full hydro Courant'
                newlabel=r'$\kappa$=3e29'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=10
                color='y'
                Rvir=234
                M1speed=500

        if (runtodo=='mw_cr_lr_dc0_2_21'):
                rundir='mw_cr_lr_dc0_2_21'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf ADV}$'
                resolabel='lr'
                M1speed=0
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=0
                color='g'
                Rvir=234


        if (runtodo=='mw_cr_lr_dc28_2_21_M1_mhd_c1000'):
                rundir='mw_cr_lr_dc28_2_21_M1_mhd_c1000'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='kappa=3e28 M1 MHD c1000'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=100
                color='r'
                haveB=1
                Rvir=234
                M1speed=1000

        if (runtodo=='mw_cr_lr_dc28_2_21_M1_mhd'):
                rundir='mw_cr_lr_dc28_2_21_M1_mhd'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='kappa=3e28 M1 MHD'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=100
                color='r'
                haveB=1
                Rvir=234

        if (runtodo=='mw_cr_lr_2_21_purestream'):
                rundir='mw_cr_lr_2_21_purestream'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='MHD streaming'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=0
                color='c'
                haveB=1
                Rvir=234


        if (runtodo=='bwmwlrdc28str'):
                rundir='mw_cr_lr_dc28_2_21_M1_mhd_stream_c1000'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf BDC28STR}$'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=100
                color='m'
                haveB=1
                M1speed=1000
                Rvir=234
                stron=1


        if (runtodo=='bwmwmrdc28str'):
                rundir='mw_cr_mr_dc28_2_21_str_ll_va_eve'
                #rundir='mw_cr_mr_dc28_2_21_M1_mhd_stream_c1000'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf BDC28STR}$'
                newlabel=r'MHD $\kappa$=3e28       Streaming'
                resolabel='mr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=100
                #color='m'
                #TK test
                color=cmaps.plasma(0.85)
                haveB=1
                M1speed=1000
                Rvir=234
                stron=1


        if (runtodo=='mw_cr_lr_dc28_2_21_M1_mhd_stream'):
                rundir='mw_cr_lr_dc28_2_21_M1_mhd_stream'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='kappa=3e28 M1 MHD streaming'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=100
                color='m'
                haveB=1
                Rvir=234

        if (runtodo=='bwmwlrdc27nc'):
                rundir='mw_cr_lr_dc27_2_21_M1_c500_nc'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC27 no cooling}$'
                resolabel='lr'
                M1speed=500
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=10
                color='g'
                Rvir=234



        if (runtodo=='bwmwlrdc27ehh'):
                rundir='mw_cr_lr_dc27_2_21_M1_c500_ehh'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC27 hole}$'
                resolabel='lr'
                M1speed=500
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=10
                color='b'
                Rvir=234

        if (runtodo=='bwmwlrdc27'):
                rundir='mw_cr_lr_dc27_2_21_M1_c500'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC27}$'
                resolabel='lr'
                M1speed=500
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=10
                #color='b'
                #TK test
                color=cmaps.plasma(0.5)
                Rvir=234

        if (runtodo=='bwmwlrdc27ds'):
                rundir='mw_cr_lr_dc27_2_21'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC27 DS}$'
                resolabel='lr'
                M1speed=-1
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=10
                color='b'
                Rvir=234


        if (runtodo=='bwmwllrdc28'):
                rundir='mw_cr_llr_dc28_2_21_M1_c1000'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC28}$'
                resolabel='llr'
                M1speed=1000
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=100
                color='r'
                Rvir=234




        if (runtodo=='bwmwlrdc28ehh'):
                rundir='mw_cr_lr_dc28_2_21_M1_c1000_ehh'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC28 hole}$'
                resolabel='lr'
                M1speed=1000
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=100
                color='r'
                Rvir=234

        if (runtodo=='bwmwlrdc28'):
                rundir='mw_cr_lr_dc28_2_21_M1_c1000'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC28}$'
                resolabel='lr'
                M1speed=1000
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=100
                #color='r'
                #TK test
                color=cmaps.plasma(0.7)
                Rvir=234



        if (runtodo=='bwmwlrdc28ref'):
                rundir='mw_cr_lr_dc28_2_21_M1_mhd_c1000_referee'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC28}$'
                resolabel='lr'
                M1speed=1000
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/refereeruns'
                kappa=100
                color='r'
                Rvir=234



        if (runtodo=='bwmwmrdc0'):
                rundir='mw_cr_mr_dc0_2_21'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf ADV}$'
                newlabel='Advection'
                resolabel='mr'
                M1speed=0
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=0
                #color='g'
                #TK test
                color=cmaps.plasma(0.34)
                Rvir=234

        if (runtodo=='bwmwmrdc27'):
                rundir='mw_cr_mr_dc27_2_21_M1_c500'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC27}$'
                newlabel = r'$\kappa$=3e27'
                resolabel='mr'
                M1speed=500
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=10
                #color='b'
                #TK test
                color=cmaps.plasma(0.5)
                Rvir=234


        if (runtodo=='bwmwmrdc28'):
                rundir='mw_cr_mr_dc28_2_21_M1_c1000'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC28}$'
                newlabel=r'$\kappa$=3e28'
                resolabel='mr'
                M1speed=1000
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=100
                #color='r'
                #TK test
                color=cmaps.plasma(0.7)
                Rvir=234
                highres=2




        if (runtodo=='bwmwlrdc28cutcr'):
                rundir='mw_cr_lr_dc28_2_21_M1_c1000_cutcr'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC28cutcr}$'
                resolabel='lr'
                M1speed=1000
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=100
                color='b'
                Rvir=234


        if (runtodo=='mw_cr_lr_dc28_2_21_M1_c1000_smallts'):
                rundir='mw_cr_lr_dc28_2_21_M1_c1000_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC28}$'
                resolabel='lr'
                M1speed=1000
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=100
                color='r'
                Rvir=234

        if (runtodo=='mw_cr_lr_dc28_2_21_M1_c2000_smallts'):
                rundir='mw_cr_lr_dc28_2_21_M1_c2000_smallts'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='kappa=3e28 M1 c2000 full hydro Courant'
                resolabel='lr'
                M1speed=2000
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=100
                color='r'
                Rvir=234

        if (runtodo=='m82_lr_2_21'):
                rundir='m82_lr_2_21'
                slabel='CR'
                havecr=0
                runtitle='M82'
                snlabel=r'$f_{mec}=1$'
                dclabel=''
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=0
                color='k'

        if (runtodo=='bwmwlrehh'):
                rundir='mw_lr_2_21_ehh'
                slabel='no CR'
                havecr=0
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf noCR hole}$'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=0
                color='k'
                Rvir=234

        if (runtodo=='bwmwmr'):
                rundir='mw_mr_2_21'
                slabel='no CR'
                havecr=0
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf noCR}$'
                newlabel='Hydro          no CR'
                resolabel='mr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=0
                #color='k'
                #TK test
                color=cmaps.plasma(0.01)
                Rvir=234
                highres=0


        if (runtodo=='bwmwlr'):
                rundir='mw_lr_2_21'
                slabel='no CR'
                havecr=0
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf noCR}$'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis/bw/mw'
                kappa=0
                color='k'
                Rvir=234



        if (runtodo=='mw_lr_2_21_IC'):
                rundir='mw_lr_2_21_IC'
                slabel='no CR'
                havecr=0
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR IC'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=0
                color='k'
                Rvir=234

        if (runtodo=='mw_mr_2_21'):
                rundir='mw_mr_2_21'
                slabel='no CR'
                havecr=0
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf mr noCR}$'
                resolabel='mr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=0
                color='k'
                Rvir=234


        if (runtodo=='mw_mr_2_21_IC'):
                rundir='mw_mr_2_21_IC'
                slabel='no CR'
                havecr=0
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='mr no CR IC'
                resolabel='mr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=0
                color='y'
                Rvir=234

        if (runtodo=='mw_cr_mr_dc28_2_21_M1'):
                rundir='mw_cr_mr_dc28_2_21_M1'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='mr kappa=3e28 cM1=1000'
                resolabel='mr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=0
                color='m'
                Rvir=234


        if (runtodo=='mw_lr_2_21_mhd'):
                rundir='mw_lr_2_21_mhd'
                slabel='no CR'
                havecr=0
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR MHD'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=0
                color='k'
                haveB=1
                Rvir=234

        if (runtodo=='mw_lr_2_21_mhd_IC'):
                rundir='mw_lr_2_21_mhd_IC'
                slabel='no CR'
                havecr=0
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR MHD IC'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=0
                color='m'
                haveB=1
                Rvir=234


        if (runtodo=='mw_lr_2_21_strongB_IC'):
                rundir='mw_lr_2_21_strongB_IC'
                slabel='no CR'
                havecr=0
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR strong MHD IC'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=0
                color='m'
                haveB=1
                Rvir=234



        if (runtodo=='mw_lr_2_21_weakB'):
                rundir='mw_lr_2_21_waekB'
                slabel='no CR'
                havecr=0
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='MHD'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=0
                color='g'
                haveB=1
                Rvir=234




        if (runtodo=='mw_cr_lr_dc28_2_21_M1_allowimbalance'):
                rundir='mw_cr_lr_dc28_2_21_M1_allowimbalance'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1 full hydro'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=10
                color='k'


        if (runtodo=='mw_cr_lr_dc27_2_21_M1_smallts_ets'):
                rundir='mw_cr_lr_dc27_2_21_M1_smallts_ets'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27 M1 full hydro Courant ETS'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=10
                color='k'


        if (runtodo=='mw_cr_lr_dc27_2_21_M1_allowimbalance'):
                rundir='mw_cr_lr_dc27_2_21_M1_allowimbalance'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27 M1 full hydro'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                kappa=10
                color='k'


        if (runtodo=='mw_cr_lr_dc28_2_10_M1_fhfcr_ets'):
                rundir='mw_cr_lr_dc28_2_10_M1_fhfcr_ets'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1 freeze hydro'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'

        if (runtodo=='mw_cr_lr_dc29_2_10_M1_fhfcr_ets'):
                rundir='mw_cr_lr_dc29_2_10_M1_fhfcr_ets'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29 M1 freeze hydro'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'


        if (runtodo=='mw_cr_lr_dc28_6_15_c2000_M1_freezehydro_freecr'):
                rundir='mw_cr_lr_dc28_6_15_c2000_M1_freezehydro_freecr'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1 freeze hydro'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='b'

        if (runtodo=='mw_cr_lr_dc29_6_15_M1_freezehydro_freecr'):
                rundir='mw_cr_lr_dc29_6_15_M1_freezehydro_freecr'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29 M1 freeze hydro'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'

        if (runtodo=='mw_cr_lr_dc29_6_15_c2000_M1_freezehydro_freecr'):
                rundir='mw_cr_lr_dc29_6_15_c2000_M1_freezehydro_freecr'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29 M1 freeze hydro'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='b'



        if (runtodo=='mw_cr_lr_dc29_6_15_ets_limitspeed'):
                rundir='mw_cr_lr_dc29_6_15_ets_limitspeed'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29 ETS speed limit=500'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='r'
                initsnap=150 


        if (runtodo=='mw_cr_lr_dc29_6_15_ets_extremelimit'):
                rundir='mw_cr_lr_dc29_6_15_ets_extremelimit'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29 ETS speed limit=5'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='b'
                initsnap=150


        if (runtodo=='mw_cr_lr_dc29_6_15_ets_nolimitspeed'):
                rundir='mw_cr_lr_dc29_6_15_ets_nolimitspeed'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29 ETS no speed limit'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='k'
                initsnap=150

        if (runtodo=='mw_cr_lr_dc28_hole_6_15_M1'):
                rundir='mw_cr_lr_dc28_hole_6_15_M1'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 with hole M1'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc28_hole_6_6'):
                rundir='mw_cr_lr_dc28_hole_6_6'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
        if (runtodo=='mw_cr_lr_dc29_6_15_c2000_M1_nofluxterm_relaxHll'):
                rundir='mw_cr_lr_dc29_6_15_c2000_M1_nofluxterm_relaxHll'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29 M1 c2000 no fluxterm relaxHll'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='y'
        if (runtodo=='mw_cr_lr_dc29_6_15_c5000_M1_nofluxterm_relaxHll'):
                rundir='mw_cr_lr_dc29_6_15_c5000_M1_nofluxterm_relaxHll'
                slabel='CR'
                havecr=6
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29 M1 c5000 no fluxterm relaxHll'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='y'
        if (runtodo=='mw_effeos'):
                rundir='mw_effeos'
                slabel='no CR'
                havecr=0
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='0'
                resolabel='llr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'

        if (runtodo=='mw_effeos_cr'):
                rundir='mw_effeos_cr'
                slabel='CR'
                havecr=1
                runtitle='MW'
                snlabel=r'$f_{mec}=1$'
                dclabel='0'
                resolabel='llr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='oasis'
                color='b'


        if (runtodo=='bwsmclrmhd'):
                rundir='smc_lr_2_21_mhd'
                slabel='CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf BnoCR}$'
                newlabel = 'MHD          no CR'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/smc'
                haveB=1
                #color='k'
                #TK test
                color=cmaps.plasma(0.01)
                Rvir=63
                M1speed=0
                highres=1


        if (runtodo=='bwsmclr'):
                rundir='smc_lr_2_21'
                slabel='CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf noCR}$'
                newlabel='Hydro          no CR'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/smc'
                #color='k'
                #TK test
                color=cmaps.plasma(0.01)
                Rvir=63
                M1speed=0
                highres=0





        if (runtodo=='bwsmcmr'):
                rundir='smc_mr_2_21'
                slabel='CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf noCR}$'
                resolabel='mr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/smc'
                color='k'
                Rvir=63
                M1speed=0

        if (runtodo=='bwsmclrstrpts'):
                rundir='smc_cr_lr_2_21_str_ll_va_evepts'
                slabel='CR'
                havecr=6
                haveB=1
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf STREAM LL}$'
                strlabel=r'$v_{\rm st}=v_{\rm A+c_s}$;       $\Gamma_{\rm st}\propto v_{\rm A}$'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/smc'
                color='purple'
                Rvir=63
                M1speed=500
                stron=1
                
                
        if (runtodo=='bwsmclrstrll'):
                rundir='smc_cr_lr_2_21_str_lowloss'
                slabel='CR'
                havecr=6
                haveB=1
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf STREAM LL}$'
                strlabel=r'$v_{\rm st}=v_{\rm A+c_s}$;       $\Gamma_{\rm st}\propto v_{\rm A}$'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/smc'
                color='purple'
                Rvir=63
                M1speed=500
                stron=1

        if (runtodo=='smclrstrllva'):
                rundir='smc_cr_lr_2_21_str_ll_va'
                slabel='CR'
                havecr=6
                haveB=1
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf STREAM LL vA}$'
                strlabel=r'$v_{\rm st}=v_{\rm A}$;       $\Gamma_{\rm st}\propto v_{\rm A}$'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/smc'
                color='navy'
                Rvir=63
                M1speed=500
                stron=1


        if (runtodo=='smclrstrnl4va'):
                rundir='smc_cr_lr_2_21_str_nl_4va'
                slabel='CR'
                havecr=6
                haveB=1
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf STREAM NL 4vA}$'
                strlabel=r'$v_{\rm st}=4v_{\rm A}$;       $\Gamma_{\rm st}\propto v_{\rm A}$'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/smc'
                color='brown'
                Rvir=63
                M1speed=1000
                stron=1

        if (runtodo=='smclrstrll4va'):
                rundir='smc_cr_lr_2_21_str_ll_4va'
                slabel='CR'
                havecr=6
                haveB=1
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf STREAM LL 4vA}$'
                strlabel=r'$v_{\rm st}=4v_{\rm A}$;       $\Gamma_{\rm st}\propto v_{\rm A}$'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/smc'
                color='brown'
                Rvir=63
                M1speed=1000
                stron=1

        if (runtodo=='bwsmclrstrold'):
                rundir='smc_cr_lr_2_21_purestream'
                slabel='CR'
                havecr=6
                haveB=1
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf STREAM}$'
                newlabel='MHD          Streaming'
                strlabel=r'$v_{\rm st}=v_{\rm A+c_s}$;       $\Gamma_{\rm st}\propto v_{\rm str}$'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/smc'
                #color='c'
                #TK test
                color='b'
                Rvir=63
                M1speed=500
                stron=1

                
        if (runtodo=='bwsmclrstr4vAeve'):
                rundir='smc_cr_lr_2_21_str_ll_4va_eve'
                slabel='CR'
                havecr=6
                haveB=1
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf STREAM}$'
                newlabel='MHD          Streaming'
                strlabel=r'$v_{\rm st}=4 v_{\rm A}$;       $\Gamma_{\rm st}\propto v_{\rm A}$'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/smc'
                #color='c'
                #TK test
                color='k'
                Rvir=63
                M1speed=500
                stron=1
                
        if (runtodo=='bwsmclrstr'):
                rundir='smc_cr_lr_2_21_str_ll_va_eve'
                slabel='CR'
                havecr=6
                haveB=1
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf STREAM}$'
                newlabel='MHD          Streaming'
                strlabel=r'$v_{\rm st}=v_{\rm A}$;       $\Gamma_{\rm st}\propto v_{\rm A}$'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/smc'
                #color='c'
                #TK test
                color=cmaps.plasma(0.93)
                Rvir=63
                M1speed=500
                stron=1

        if (runtodo=='bwsmclrdc0'):
                rundir='smc_cr_lr_dc0_2_21'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf ADV}$'
                newlabel='Advection'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/smc'
                #color='g'
                #TK test
                color=cmaps.plasma(0.34)
                Rvir=63
                M1speed=0


        if (runtodo=='bwsmclrdc27'):
                rundir='smc_cr_lr_dc27_2_21_M1_c500'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC27}$'
                newlabel = r'$\kappa$=3e27'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/smc'
                #color='b'
                #TK test
                color=cmaps.plasma(0.5)
                Rvir=63
                M1speed=500

        if (runtodo=='bwsmcmrdc27'):
                rundir='smc_cr_mr_dc27_2_21_M1_c500'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC27}$'
                resolabel='mr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/smc'
                color='b'
                Rvir=63
                M1speed=500


        if (runtodo=='bwsmclrdc28'):
                rundir='smc_cr_lr_dc28_2_21_M1_c1000'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC28}$'
                newlabel=r'$\kappa$=3e28'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/smc'
                #color='r'
                #TK test
                color=cmaps.plasma(0.7)
                Rvir=63
                M1speed=1000
                highres=2

        if (runtodo=='bwsmcmrdc28'):
                rundir='smc_cr_mr_dc28_2_21_M1_c1000'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC28}$'
                resolabel='mr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/smc'
                color='r'
                Rvir=63
                M1speed=1000

        if (runtodo=='bwsmclrdc28cutcr'):
                rundir='smc_cr_lr_dc28_2_21_M1_c1000_cutcr'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC28cutcr}$'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/smc'
                color='b'
                Rvir=63
                M1speed=1000


        if (runtodo=='bwsmclrdc29'):
                rundir='smc_cr_lr_dc29_2_21_M1_c2000'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC29}$'
                newlabel=r'$\kappa$=3e29'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/smc'
                #color='y'
                #TK test
                color=cmaps.plasma(0.85)
                Rvir=63
                M1speed=2000
                highres=3

        if (runtodo=='bwsmclrdc29cutcr'):
                rundir='smc_cr_lr_dc29_2_21_M1_c2000_cutcr'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC29cutcr}$'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/smc'
                color='k'
                Rvir=63
                M1speed=2000




        if (runtodo=='bwsmclrdc27str'):
                rundir='smc_cr_lr_dc27_2_21_M1_mhd_stream_c500'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC27STR}$'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/smc'
                color='k'
                Rvir=63
                M1speed=500
                haveB=1
                stron=1



        if (runtodo=='bwsmclrdc28strold'):
                rundir='smc_cr_lr_dc28_2_21_M1_mhd_stream_c1000'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf BDC28STR}$'
                newlabel=r'MHD $\kappa$=3e28       Streaming'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/smc'
                #color='m'
                #TK test
                color=cmaps.plasma(0.85)
                Rvir=63
                M1speed=1000
                haveB=1
                stron=1
                
        if (runtodo=='bwsmclrdc28str'):
                rundir='smc_cr_lr_dc28_2_21_str_ll_va_eve'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf BDC28STR}$'
                newlabel=r'MHD $\kappa$=3e28       Streaming'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/smc'
                #color='m'
                #TK test
                color=cmaps.plasma(0.85)
                Rvir=63
                M1speed=1000
                haveB=1
                stron=1


        if (runtodo=='bwsmclrdc28mhd'):
                rundir='smc_cr_lr_dc28_2_21_M1_mhd_c1000_original'
                #rundir='smc_cr_lr_dc28_2_21_M1_mhd_c1000'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf BDC28}$'
                newlabel=r'MHD $\kappa$=3e28'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                #maindir='oasis/refereeruns'
                maindir='oasis/bw/smc'
                #color='r'
                #TK test
                color=cmaps.plasma(0.7)
                Rvir=63
                M1speed=1000
                haveB=1


        if (runtodo=='bwsmclrdc28mhdref'):
                rundir='smc_cr_lr_dc28_2_21_M1_mhd_c1000_referee'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf BDC28}$'
                newlabel=r'MHD $\kappa$=3e28'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/refereeruns'
                color='k'
                Rvir=63
                M1speed=1000
                haveB=1



        if (runtodo=='bwsmclrori'):
                rundir='smc_cr_lr_dc28_2_21_M1_mhd_c1000_original'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf BDC28}$'
                newlabel=r'MHD $\kappa$=3e28'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/refereeruns'
                color='r'
                Rvir=63
                M1speed=1000
                haveB=1



        if (runtodo=='bwsmclrdc29str'):
                rundir='smc_cr_lr_dc29_2_21_M1_mhd_stream_c2000'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC29STR}$'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/smc'
                #TK test
                #color='y'
                color=cmaps.plasma(0)
                Rvir=63
                M1speed=2000
                haveB=1
                stron=1



        if (runtodo=='smc_cr_lr_dc28_2_21_M1_smallts'):
                rundir='smc_cr_lr_dc28_2_21_M1_smallts'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='kappa=3e28 M1 c500'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='r'
                Rvir=63
                M1speed=500

        if (runtodo=='smc_cr_lr_dc28_2_21_M1_c1000_smallts'):
                rundir='smc_cr_lr_dc28_2_21_M1_c1000_smallts'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='kappa=3e28 M1 c1000'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='r'
                Rvir=63
                M1speed=1000


        if (runtodo=='smc_cr_lr_dc28_2_21_M1_c2000_smallts'):
                rundir='smc_cr_lr_dc28_2_21_M1_c2000_smallts'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='kappa=3e28 M1 c2000'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='r'
                Rvir=63
                M1speed=2000


        if (runtodo=='smc_cr_lr_dc29_2_21_M1_c1000_smallts'):
                rundir='smc_cr_lr_dc29_2_21_M1_c1000_smallts'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='kappa=3e29 M1 c1000'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='g'
                Rvir=63
                M1speed=1000

        if (runtodo=='smc_cr_lr_dc29_2_21_M1_c2000_smallts'):
                rundir='smc_cr_lr_dc29_2_21_M1_c2000_smallts'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='kappa=3e29 M1 c2000'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='g'
                Rvir=63
                M1speed=2000

        if (runtodo=='smc_cr_lr_dc29_2_21_M1_c4000_smallts'):
                rundir='smc_cr_lr_dc29_2_21_M1_c4000_smallts'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='kappa=3e29 M1 c4000'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='g'
                Rvir=63
                M1speed=4000



        if (runtodo=='smc_cr_lr_dc28_2_21_M1_mhd'):
                rundir='smc_cr_lr_dc28_2_21_M1_mhd'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='kappa=3e28 M1 c500 MHD'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='r'
                haveB=1
                Rvir=63

        if (runtodo=='smc_cr_lr_dc28_2_21_M1_mhd_stream'):
                rundir='smc_cr_lr_dc28_2_21_M1_mhd_stream'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='kappa=3e28 M1 c500 MHD streaming'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='m'
                haveB=1
                Rvir=63


        if (runtodo=='smc_cr_lr_dc28_2_21_M1_mhd_stream_c1000'):
                rundir='smc_cr_lr_dc28_2_21_M1_mhd_stream_c1000'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='kappa=3e28 M1 c1000 MHD streaming'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw'
                color='m'
                haveB=1
                Rvir=63

        if (runtodo=='smc_cr_lr_2_21_purestream'):
                rundir='smc_cr_lr_2_21_purestream'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='MHD streaming'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw'
                color='c'
                haveB=1
                Rvir=63



        if (runtodo=='smc_cr_lr_dc28_2_21_M1_smallts_cutcr'):
                rundir='smc_cr_lr_dc28_2_21_M1_smallts_cutcr'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1 c500 cutcr'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='r'
                Rvir=63
                M1speed=500


        if (runtodo=='smc_cr_lr_dc27_2_21_M1_smallts'):
                rundir='smc_cr_lr_dc27_2_21_M1_smallts'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='kappa=3e27 M1 c500'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='b'
                Rvir=63
                M1speed=500


        if (runtodo=='smc_cr_lr_dc26_6_15'):
                rundir='smc_cr_lr_dc26_6_15'
                slabel='CR'
                havecr=6
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e26'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'


        if (runtodo=='sbc_lr_6_15'):
                rundir='sbc_lr_6_15'
                slabel='no CR'
                havecr=0
                runtitle='SBC'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'


        if (runtodo=='smc_lr_6_15_sn300'):
                rundir='smc_lr_6_15_sn300'
                slabel='no CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='0'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='k'

        if (runtodo=='smc_lr_2_21_ehalo_IC'):
                rundir='smc_lr_2_21_ehalo_IC'
                slabel='no CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='0'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='k'

        if (runtodo=='smc_lr_2_21_thickdisk'):
                rundir='smc_lr_2_21_thickdisk'
                slabel='no CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR thickdisk'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='k'
                Rvir=63


        if (runtodo=='smc_lr_2_21_vtd'):
                rundir='smc_lr_2_21_vtd'
                slabel='no CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR very thickdisk'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='k'
                Rvir=63



        if (runtodo=='smc_lr_2_21_vtd_t1'):
                rundir='smc_lr_2_21_vtd_t1'
                slabel='no CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR vtdt1'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='k'
                Rvir=63




        if (runtodo=='smc_lr_2_21_su'):
                rundir='smc_lr_2_21_su'
                slabel='no CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR Su IC'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='k'
                Rvir=63


        if (runtodo=='smc_lr_2_21_thickdisk_mhd'):
                rundir='smc_lr_2_21_thickdisk_mhd'
                slabel='no CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR thickdisk MHD init B=1e-8G'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='k'
                Rvir=63
                haveB=1




        if (runtodo=='smc_lr_2_21_mhd'):
                rundir='smc_lr_2_21_mhd'
                slabel='no CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR MHD init B=1e-8G'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='k'
                Rvir=63
                haveB=1


        if (runtodo=='smc_mr_2_21_thickdisk'):
                rundir='smc_mr_2_21_thickdisk'
                slabel='no CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR mr thick disk'
                resolabel='mr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='y'
                Rvir=63
                haveB=1


        if (runtodo=='smc_mr_2_21_su'):
                rundir='smc_mr_2_21_su'
                slabel='no CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR mr'
                resolabel='mr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='y'
                Rvir=63
                haveB=1

        if (runtodo=='smc_mr_2_21_vtd_t1'):
                rundir='smc_mr_2_21_vtd_t1'
                slabel='no CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR mr vtd t1'
                resolabel='mr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='y'
                Rvir=63


        if (runtodo=='smc_mr_2_21_vtd'):
                rundir='smc_mr_2_21_vtd'
                slabel='no CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR mr vtd'
                resolabel='mr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='y'
                Rvir=63
                haveB=1





        if (runtodo=='smc_mr_2_21_mhd'):
                rundir='smc_mr_2_21_mhd'
                slabel='no CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR mr MHD init B=5e-10G'
                resolabel='mr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='y'
                Rvir=63
                haveB=1

        if (runtodo=='smc_mr_2_21_vtd_t1_mhd'):
                rundir='smc_mr_2_21_vtd_t1_mhd'
                slabel='no CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR mr MHD vtdt1'
                resolabel='mr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='y'
                Rvir=63
                haveB=1



        if (runtodo=='smc_lr_2_21_thickdisk_mhd'):
                rundir='smc_lr_2_21_thickdisk_mhd'
                slabel='no CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR MHD thick disk'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='y'
                Rvir=100
                haveB=1

        if (runtodo=='smc_lr_2_21_vtd_t1_mhd'):
                rundir='smc_lr_2_21_vtd_t1_mhd'
                slabel='no CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR MHD vtd'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='k'
                Rvir=100
                haveB=1



        if (runtodo=='smc_lr_2_21_mhd_su'):
                rundir='smc_lr_2_21_mhd_su'
                slabel='no CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR MHD Su IC'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='y'
                Rvir=100
                haveB=1

        if (runtodo=='smc_lr_2_21_weakB_p10G_su'):
                rundir='smc_lr_2_21_weakB_p10G_su'
                slabel='no CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR weak B=1e-10G Su IC'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='y'
                Rvir=100
                haveB=1


        if (runtodo=='smc_lr_2_21_weakB'):
                rundir='smc_lr_2_21_weakB'
                slabel='no CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR weak init B=1e-9G'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='m'
                Rvir=100
                haveB=1



        if (runtodo=='smc_lr_2_21_weakB_p12G'):
                rundir='smc_lr_2_21_weakB_p12G'
                slabel='no CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR weak init B=1e-12G'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='y'
                Rvir=100
                haveB=1



        if (runtodo=='smc_lr_2_21'):
                rundir='smc_lr_2_21'
                slabel='no CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='k'
                Rvir=100


        if (runtodo=='sbc_lr_2_21_su'):
                rundir='sbc_lr_2_21_su'
                slabel='no CR'
                havecr=0
                runtitle='SBC'
                snlabel=r'$f_{mec}=1$'
                dclabel='0'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='k'


        if (runtodo=='sbc_lr_2_21_tk'):
                rundir='sbc_lr_2_21_tk'
                slabel='no CR'
                havecr=0
                runtitle='SBC'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='k'
                Rvir=140


        if (runtodo=='sbc_lr_2_21_extend'):
                rundir='sbc_lr_2_21_extend'
                slabel='no CR'
                havecr=0
                runtitle='SBC'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/sbc'
                color='k'
                Rvir=140


        if (runtodo=='sbc_lr_2_21_extend2'):
                rundir='sbc_lr_2_21_extend2'
                slabel='no CR'
                havecr=0
                runtitle='SBC'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/sbc'
                color='k'
                Rvir=140


        if (runtodo=='sbc_lr_2_21_extend3'):
                rundir='sbc_lr_2_21_extend3'
                slabel='no CR'
                havecr=0
                runtitle='SBC'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/sbc'
                color='k'
                Rvir=140

        if (runtodo=='bwsbclr'):
                rundir='sbc_lr_2_21_extend4'
                slabel='no CR'
                havecr=0
                runtitle='SBC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf noCR}$'
                newlabel='Hydro no CR'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/sbc'
                #color='k'
                #TK test
                color=cmaps.plasma(0.01)
                Rvir=140
                highres=0

        if (runtodo=='bwsbclrmhd'):
                rundir='sbc_lr_2_21_extend4_mhd'
                slabel='no CR'
                havecr=0
                haveB=1
                runtitle='SBC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf noCR}$'
                newlabel='MHD no CR'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/sbc'
                color='k'
                Rvir=140
                highres=1


        if (runtodo=='bwsbclrdc0'):
                rundir='sbc_cr_lr_dc0_2_21_extend4'
                slabel='CR'
                havecr=6
                runtitle='SBC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf ADV}$'
                newlabel='Advection'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/sbc'
                #color='g'
                #TK test
                color=cmaps.plasma(0.34)
                Rvir=140
                M1speed=0


        if (runtodo=='bwsbclrdc27'):
                rundir='sbc_cr_lr_dc27_2_21_M1_extend4'
                slabel='CR'
                havecr=6
                runtitle='SBC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC27}$'
                newlabel = r'$\kappa$=3e27'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/sbc'
                #color='b'
                #TK test
                color=cmaps.plasma(0.5)
                Rvir=140
                M1speed=500

        if (runtodo=='bwsbclrdc28'):
                rundir='sbc_cr_lr_dc28_2_21_M1_c1000_extend4'
                slabel='CR'
                havecr=6
                runtitle='SBC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC28}$'
                newlabel = r'$\kappa$=3e28'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/sbc'
                #color='r'
                #TK test
                color=cmaps.plasma(0.7)
                Rvir=140
                M1speed=1000
                highres=2


        if (runtodo=='bwsbclrdc29'):
                rundir='sbc_cr_lr_dc29_2_21_M1_c2000_extend4'
                slabel='CR'
                havecr=6
                runtitle='SBC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf DC29}$'
                newlabel = r'$\kappa$=3e29'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/sbc'
                #color='y'
                #TK test
                color=cmaps.plasma(0.85)
                Rvir=140
                M1speed=2000
                highres=3


        if (runtodo=='bwsbclrstrold'):
                rundir='sbc_cr_lr_2_21_purestream'
                slabel='CR'
                havecr=6
                haveB=1
                runtitle='SBC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf BDC28}$'
                newlabel = r'MHD     Streaming'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/sbc'
                #color='c'
                #TK test
                color=cmaps.plasma(0.93)
                Rvir=140
                M1speed=500
                stron=1


        if (runtodo=='bwsbclrstr'):
                rundir='sbc_cr_lr_2_21_str_ll_va_eve'
                #rundir='sbc_cr_lr_2_21_purestream'
                slabel='CR'
                havecr=6
                haveB=1
                runtitle='SBC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf BDC28}$'
                newlabel = r'MHD     Streaming'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/sbc'
                #color='c'
                #TK test
                color=cmaps.plasma(0.93)
                Rvir=140
                M1speed=500
                stron=1


        if (runtodo=='bwsbclrdc28mhd'):
                rundir='sbc_cr_lr_dc28_2_21_M1_mhd_c1000'
                slabel='CR'
                havecr=6
                haveB=1
                runtitle='SBC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf BDC28}$'
                newlabel = r'MHD $\kappa$=3e28'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/sbc'
                #color='r'
                #TK test
                color=cmaps.plasma(0.7)
                Rvir=140
                M1speed=1000


        if (runtodo=='bwsbclrdc28str'):
                rundir='sbc_cr_lr_dc28_2_21_M1_mhd_stream_c1000'
                slabel='CR'
                havecr=6
                haveB=1
                runtitle='SBC'
                snlabel=r'$f_{mec}=1$'
                dclabel=r'${\bf BDC28STR}$'
                newlabel = r'MHD $\kappa$=3e28    Streaming'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis/bw/sbc'
                #color='m'
                #TK test
                color=cmaps.plasma(0.85)
                Rvir=140
                M1speed=1000
                stron=1




        if (runtodo=='sbc_cr_lr_dc28_2_21_tk'):
                rundir='sbc_cr_lr_dc28_2_21_tk'
                slabel='CR'
                havecr=6
                runtitle='SBC'
                snlabel=r'$f_{mec}=1$'
                dclabel='kappa=3e28'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='r'
                Rvir=140


        if (runtodo=='sbc_cr_lr_dc28_2_21_tk_initCR1'):
                rundir='sbc_cr_lr_dc28_2_21_tk_initCR1'
                slabel='CR'
                havecr=6
                runtitle='SBC'
                snlabel=r'$f_{mec}=1$'
                dclabel='kappa=3e28 init CR=1eV/cm3'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='b'
                Rvir=140


        if (runtodo=='sbc_cr_lr_dc28_2_21_tk_initCR10'):
                rundir='sbc_cr_lr_dc28_2_21_tk_initCR10'
                slabel='CR'
                havecr=6
                runtitle='SBC'
                snlabel=r'$f_{mec}=1$'
                dclabel='kappa=3e28 init CR=10eV/cm3'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='g'
                Rvir=140


        if (runtodo=='sbc_mr_2_21_su'):
                rundir='sbc_mr_2_21_su'
                slabel='no CR'
                havecr=0
                runtitle='SBC'
                snlabel=r'$f_{mec}=1$'
                dclabel='0'
                resolabel='mr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'
                color='g'




        if (runtodo=='smc_booth_5_27'):
                rundir='smc_booth_5_27'
                slabel='no CR'
                havecr=0
                runtitle='SMC'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.01
                iavesfr = 0.1
                subdir='/output/'
                maindir='oasis'


        if (runtodo=='FIRE_2_0_or_h573_criden1000_noaddm_sggs'):
                rundir='FIRE_2_0_or_h573_criden1000_noaddm_sggs'
                slabel='no CR'
                havecr=0
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                maindir='oasis'

        if (runtodo=='FIRE_2_0_h573_cr_6_15_ets'):
                rundir='FIRE_2_0_h573_cr_6_15_ets'
                slabel='CR'
                havecr=6
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                maindir='/oasis/'
                usesnaplist=1
                color='g'

        if (runtodo=='FIRE_2_0_h573_cr_6_15_M1'):
                rundir='FIRE_2_0_h573_cr_6_15_M1'
                slabel='CR'
                havecr=6
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                maindir='/oasis/'
                usesnaplist=1
                color='g'

        if (runtodo=='FIRE_2_0_h573_cr_6_15_M1_c2000'):
                rundir='FIRE_2_0_h573_cr_6_15_M1_c2000'
                slabel='CR'
                havecr=6
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                maindir='/oasis/'
                usesnaplist=1
                color='g'


        if (runtodo=='FIRE_2_0_h573_cr_6_15_ets_dc27'):
                rundir='FIRE_2_0_h573_cr_6_15_ets_dc27'
                slabel='CR'
                havecr=6
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                maindir='/oasis/'
                usesnaplist=1
                color='b'

        if (runtodo=='m12_tmp'):
                rundir='m12_tmp'
                slabel='CR'
                havecr=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                maindir='/scratch/'
                usesnaplist=1
                color='b'
                h0=0.702
                halostr='00'


        if (runtodo=='m10q_M1CR'):
                rundir='m10q_M1CR'
                slabel='CR'
                havecr=1
                runtitle='m10q'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                maindir='/scratch/'
                usesnaplist=1
                color='b'
                withinRv = 1
                h0=0.702
                halostr='00'

        if (runtodo=='m12i_res7000'):
                rundir='m12i_res7000'
                slabel='no CR'
                havecr=0
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                maindir='/oasis/'
                usesnaplist=1
                color='b'
                h0=0.702




        if (runtodo=='m12i_res7000_output'):
                rundir='m12i_res7000_output'
                slabel='no CR'
                havecr=0
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/'
                cosmo=1
                maindir='/ctrapp_scratch/'
                usesnaplist=1
                color='b'
                h0=0.702


        if (runtodo=='FIRE_2_0_h573_cr_6_15_ets_dc0'):
                rundir='FIRE_2_0_h573_cr_6_15_ets_dc0'
                slabel='CR'
                havecr=6
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='0'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                maindir='/oasis/'
                usesnaplist=1
                color='r'
                h0=0.702

        if (runtodo=='FIRE_2_0_h573_cr_6_15_M1_z0'):
                rundir='FIRE_2_0_h573_cr_6_15_M1_z0'
                slabel='CR'
                havecr=6
                runtitle='573'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                maindir='/oasis/'
                usesnaplist=1
                color='r'
                h0=0.702
                withinRv = 1
                halostr='00'


        if (runtodo=='FIRE_2_0_h573_2_21'):
                rundir='FIRE_2_0_h573_2_21'
                slabel='no CR'
                havecr=0
                runtitle='573'
                snlabel=r'$f_{mec}=1$'
                dclabel='no CR'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                maindir='/oasis/'
                usesnaplist=1
                color='k'
                h0=0.702
                withinRv = 0
                halostr='00'


        if (runtodo=='f476'):
                rundir='FIRE_2_0_or_h476_noaddm/'
                havecr=0
                slabel='no CR'
                dclabel='Frozen FIRE-2 no MHD'
                runtitle='573'
                snlabel=r'$f_{mec}=1$'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                maindir='/oasis/'
                color='c'
                firever=2
                withinRv = 1
                h0=0.702
                halostr='00'
                cosmo=1



        if (runtodo=='FIRE_2_0_h573_cr_2_21_M1_z0'):
                rundir='FIRE_2_0_h573_cr_2_21_M1_z0'
                slabel='CR'
                havecr=6
                runtitle='573'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                maindir='/oasis/'
                usesnaplist=1
                color='r'
                h0=0.702
                withinRv = 0
                halostr='00'

        if (runtodo=='FIRE_2_0_h573_cr_2_21_M1_z0_dc27'):
                rundir='FIRE_2_0_h573_cr_2_21_M1_z0_dc27'
                slabel='CR'
                havecr=6
                runtitle='573'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e27 M1'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                maindir='/oasis/'
                usesnaplist=1
                color='b'
                h0=0.702
                withinRv = 0
                halostr='00'



        if (runtodo=='FIRE_2_0_h573_cr_2_21_M1_z0_cd100'):
                rundir='FIRE_2_0_h573_cr_2_21_M1_z0_cd100'
                slabel='CR'
                havecr=6
                runtitle='573'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1 critical density 100'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                maindir='/oasis/'
                usesnaplist=1
                color='k'
                h0=0.702
                withinRv = 0
                halostr='00'


        if (runtodo=='FIRE_2_0_h573_cr_2_21_M1_z0_defaultts'):
                rundir='FIRE_2_0_h573_cr_2_21_M1_z0_defaultts'
                slabel='CR'
                havecr=6
                runtitle='573'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                maindir='/oasis/'
                usesnaplist=1
                color='r'
                h0=0.702
                withinRv = 0
                halostr='00'


        if (runtodo=='FIRE_2_0_h573_cr_2_21_M1_z0_testts'):
                rundir='FIRE_2_0_h573_cr_2_21_M1_z0_testts'
                slabel='CR'
                havecr=6
                runtitle='573'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                maindir='/oasis/'
                usesnaplist=1
                color='r'
                h0=0.702
                withinRv = 0
                halostr='00'

        if (runtodo=='FIRE_2_0_h573_cr_2_21_M1_z0_ts2'):
                rundir='FIRE_2_0_h573_cr_2_21_M1_z0_ts2'
                slabel='CR'
                havecr=6
                runtitle='573'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28 M1'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                maindir='/oasis/'
                usesnaplist=1
                color='r'
                h0=0.702
                withinRv = 0
                halostr='00'


        if (runtodo=='FIRE_2_0_h573_6_15'):
                rundir='FIRE_2_0_h573_6_15'
                slabel='no CR'
                havecr=0
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='0'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                maindir='/oasis/'
                usesnaplist=1
                color='k'
                h0=0.702
                withinRv=1


        if (runtodo=='FIRE_2_0_or_h573_criden1000_noaddm_sggs'):
                rundir='FIRE_2_0_or_h573_criden1000_noaddm_sggs'
                slabel='no CR'
                havecr=0
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='0'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                maindir='/scratch/'
                usesnaplist=1
                color='k'
                h0=0.702

        if (runtodo=='FIRE_2_0_h573_CRtest2'):
                rundir='FIRE_2_0_h573_CRtest2'
                slabel='CR'
                havecr=6
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                maindir='/oasis/'
                usesnaplist=1
                h0=0.702
                withinRv=1


        if (runtodo=='m11q_mass7000_MHDCR_K100_M1_tkFIX'):
                rundir='m11q_mass7000_MHDCR_K100_M1_tkFIX'
                slabel='CR'
                havecr=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                halostr='00'
                maindir='/scratch/'
                usesnaplist=0
                h0=0.702
                usepep=1
                withinRv=1
                firever=2


        if (runtodo=='m10qmhdcv'):
                rundir='/m10q_mass250_MHDCR_M1_tkFIX/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                havecr=0
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m10q'
                crlabel='mhdcv'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='y'
                highres=1
                newlabel = 'MHD          no CR'


        if (runtodo=='m10qcr_b_70'):
                rundir='/m10q_mass250_MHDCR_M1_tkFIX/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                havecr=1
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m10q cr'
                crlabel='cr70'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='03'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='y'
                highres=2
                newlabel=r'MHD $\kappa$=3e28'


        if (runtodo=='m09mhdcv'):
                rundir='/m09/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                havecr=0
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m09 mhd'
                crlabel='mhdcv'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='m'
                newlabel = 'MHD          no CR'


        if (runtodo=='m09cr_b_70'):
                rundir='/m09/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                havecr=1
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m09 cr'
                crlabel='cr70'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='m'

        if (runtodo=='m10vmhdcv'):
                rundir='/m10v/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                galcen=0
                havecr=0
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m10v MHD'
                crlabel='mhdcv'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='04'
                usesnaplist=0
                h0=0.702
                snapsep=10
                usepep=0
                withinRv=1
                firever=2
                color='brown'
                highres=1
                newlabel = 'MHD          no CR'


        if (runtodo=='m10vcr_b_70'):
                rundir='/m10v/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                galcen=0
                havecr=1
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m10v cr'
                crlabel='cr70'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='04'
                usesnaplist=0
                h0=0.702
                snapsep=10
                usepep=0
                withinRv=1
                firever=2
                color='brown'
                highres=2
                newlabel=r'MHD $\kappa$=3e28'


        if (runtodo=='m10vcr_700'):
                rundir='/m10v/cr_700/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                galcen=0
                havecr=1
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m10v cr dc29'
                crlabel='cr700'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='01'
                usesnaplist=0
                h0=0.702
                snapsep=10
                usepep=0
                withinRv=1
                firever=2
                color='brown'
                highres=3
                newlabel=r'MHD $\kappa$=3e29'


        if (runtodo=='m11bcr_b_70'):
                rundir='/m11b/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                havecr=1
                haveB=1
                galcen=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m11b cr'
                crlabel='cr70'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='r'
                beginno=520
                snapsep=10
                highres=2
                newlabel=r'MHD $\kappa$=3e28'


        if (runtodo=='m11bcr_700'):
                rundir='/m11b/cr_700/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                havecr=1
                haveB=1
                galcen=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m11b cr dc29'
                crlabel='cr700'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='r'
                beginno=520
                snapsep=10
                highres=3
                newlabel=r'MHD $\kappa$=3e29'

        if (runtodo=='m11dcr_70va'):
                rundir='/m11d/cr_70_streamVA/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                beginno=500
                havecr=1
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m11d cr vA'
                crlabel='cr70'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.68
                usepep=1
                withinRv=1
                firever=2
                color='b'
                snapsep=10


        if (runtodo=='m11dcr_b_70'):
                rundir='/m11d/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                havecr=1
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m11d cr'
                crlabel='cr70'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.68
                usepep=0
                withinRv=1
                firever=2
                color='b'
                highres=2
                newlabel=r'MHD $\kappa$=3e28'


        if (runtodo=='m11dcr_700'):
                rundir='/m11d/cr_700/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                havecr=1
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m11d cr dc29'
                crlabel='cr700'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='01'
                usesnaplist=0
                h0=0.68
                usepep=0
                withinRv=1
                firever=2
                color='b'
                beginno=520
                snapsep=10
                highres=3
                newlabel=r'MHD $\kappa$=3e29'

        if (runtodo=='m11hcr_700'):
                rundir='/m11h/cr_700/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                havecr=1
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m11h cr'
                crlabel='cr700'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.68
                usepep=0
                withinRv=1
                firever=2
                color='c'
                highres=3
                newlabel=r'MHD $\kappa$=3e29'
                
        if (runtodo=='m11hcr_b_70'):
                rundir='/m11h/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                havecr=1
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m11h cr'
                crlabel='cr70'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.68
                usepep=0
                withinRv=1
                firever=2
                color='c'
                highres=2
                newlabel=r'MHD $\kappa$=3e28'


        if (runtodo=='m11fcr_700'):
                rundir='/m11f/cr_700/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                havecr=1
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m11f cr 700'
                crlabel='cr700'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='brown'
                highres=3
                newlabel=r'MHD $\kappa$=3e29'


        if (runtodo=='m11fcr_b_70'):
                rundir='/m11f/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                havecr=1
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m11f cr'
                crlabel='cr70'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='brown'
                highres=2
                newlabel=r'MHD $\kappa$=3e28'


        if (runtodo=='m11gcr_700'):
                rundir='/m11g/cr_700/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                havecr=1
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m11g cr 700'
                crlabel='cr700'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='orange'
                highres=3
                newlabel=r'MHD $\kappa$=3e29'
                

        if (runtodo=='m11gcr_b_70'):
                rundir='/m11g/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                havecr=1
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m11g cr'
                crlabel='cr70'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='orange'
                highres=2
                newlabel=r'MHD $\kappa$=3e28'


        if (runtodo=='m11bmhdcv'):
                rundir='/m11b/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='no CR'
                havecr=0
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m11b mhd'
                crlabel='mhdcv'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='r'
                highres=1
                newlabel = 'MHD          no CR'


        if (runtodo=='m11fmhdcv'):
                rundir='/m11f/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='no CR'
                havecr=0
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m11f mhd'
                crlabel='mhdcv'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='r'
                highres=1
                newlabel = 'MHD          no CR'

        if (runtodo=='m11gmhdcv'):
                rundir='/m11g/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='no CR'
                havecr=0
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m11g mhd'
                crlabel='mhdcv'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='r'
                highres=1
                newlabel = 'MHD          no CR'

        if (runtodo=='m11dmhdcv'):
                rundir='/m11d/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='no CR'
                havecr=0
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m11d mhd'
                crlabel='mhdcv'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.68
                usepep=0
                withinRv=1
                firever=2
                color='r'
                highres=1
                newlabel = 'MHD          no CR'

        if (runtodo=='m11hmhdcv'):
                rundir='/m11h/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='no CR'
                havecr=0
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m11h mhd'
                crlabel='mhdcv'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.68
                usepep=0
                withinRv=1
                firever=2
                color='r'
                highres=1
                newlabel = 'MHD          no CR'


        if (runtodo=='m12mmhdcv'):
                rundir='/m12m_mass56000/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='no CR'
                galcen=0
                havecr=0
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m12m mhd'
                multifile='y'
                crlabel=''
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='g'
                snapsep=10
                highres=1
                newlabel = 'MHD          no CR'

        if (runtodo=='m12mcr_b_70'):
                rundir='/m12m_mass56000/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                galcen=1
                havecr=1
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m12m cr'
                crlabel=''
                multifile='y'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='g'
                snapsep=10
                newlabel=r'MHD $\kappa$=3e28'

        if (runtodo=='m12mcr_70va'):
                rundir='/m12m_mass56000/cr_70_streamVA/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                galcen=1
                havecr=1
                haveB=1
                multifile='y'
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m12m cr vA'
                crlabel=''
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='g'
                snapsep=10



        if (runtodo=='m12mcr_700'):
                rundir='/m12m_mass56000/cr_700/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                havecr=1
                haveB=1
                galcen=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m12m cr dc29'
                crlabel='cr700'
                resolabel='lr'
                multifile='y'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='g'
                beginno=100
                snapsep=10
                highres=3
                newlabel=r'MHD $\kappa$=3e29'
                

        if (runtodo=='m12icr_b_70'):
                rundir='/m12i_mass56000/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                galcen=1
                havecr=1
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m12i cr'
                crlabel='cr70'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='g'
                snapsep=10
                highres=2
                newlabel=r'MHD $\kappa$=3e28'


        if (runtodo=='m12fcr_b_70'):
                rundir='/m12f_mass56000/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                galcen=1
                havecr=1
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m12f cr'
                crlabel='cr70'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='g'
                snapsep=10
                newlabel=r'MHD $\kappa$=3e28'


        if (runtodo=='m12icr_70va'):
                rundir='/m12i_mass56000/cr_70_streamVA/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                galcen=1
                havecr=1
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m12i cr vA'
                crlabel='cr70'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='g'
                snapsep=10


        if (runtodo=='m12icr_700'):
                rundir='/m12i_mass56000/cr_700/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                havecr=1
                haveB=1
                galcen=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m12i cr dc29'
                crlabel='cr700'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='g'
                beginno=520
                snapsep=10
                newlabel=r'MHD $\kappa$=3e29'


        if (runtodo=='m12fcr_700'):
                rundir='/m12f_mass56000/cr_700/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                havecr=1
                haveB=1
                galcen=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m12f cr dc29'
                crlabel='cr700'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='g'
                beginno=520
                snapsep=10
                highres=3
                newlabel=r'MHD $\kappa$=3e29'


        if (runtodo=='m12imhdcv'):
                rundir='/m12i_mass56000/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                havecr=0
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m12i mhd'
                crlabel='mhdcv'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='k'
                newlabel = 'MHD          no CR'


        if (runtodo=='m12fmhdcv'):
                rundir='/m12f_mass56000/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                slabel='CR'
                havecr=0
                haveB=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='m12f mhd'
                crlabel='mhdcv'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                cosmo=1
                halostr='00'
                usesnaplist=0
                h0=0.702
                usepep=0
                withinRv=1
                firever=2
                color='k'
                highres=1
                newlabel = 'MHD          no CR'
                


        if (runtodo=='m11q_mass7000_MHDCR_K1000_M1_tkFIX'):
                rundir='m11q_mass7000_MHDCR_K1000_M1_tkFIX'
                slabel='CR'
                havecr=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                halostr='00'
                maindir='/scratch/'
                usesnaplist=0
                h0=0.702
                color='r'
                usepep=0
                withinRv=1
                firever=2

        if (runtodo=='m10q_mass250_MHDCR_K100_M1_tkFIX'):
                rundir='m10q_mass250_MHDCR_K100_M1_tkFIX'
                slabel='CR'
                havecr=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                halostr='00'
                subdir='/output/'
                cosmo=1
                maindir='/scratch/'
                usesnaplist=0
                h0=0.702
                withinRv=1


        if (runtodo=='m10q_mass250_MHDCR_K1000_M1_tkFIX'):
                rundir='m10q_mass250_MHDCR_K1000_M1_tkFIX'
                slabel='CR'
                havecr=1
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e29'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                halostr='00'
                subdir='/output/'
                cosmo=1
                maindir='/scratch/'
                usesnaplist=0
                h0=0.702
                color='r'
                withinRv=1

        if (runtodo=='FIRE_2_0_h573_CRtest2_equaltimestep'):
                rundir='FIRE_2_0_h573_CRtest2_equaltimestep'
                slabel='CR'
                havecr=6
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                maindir='/oasis/'
                usesnaplist=1
                h0=0.702

        if (runtodo=='FIRE_2_0_h553_CRtest2'):
                rundir='FIRE_2_0_h553_CRtest2'
                slabel='CR'
                havecr=6
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                maindir='/oasis/'
                h0=0.702


        if (runtodo=='FIRE_2_0_h553_CRtest2'):
                rundir='FIRE_2_0_h553_CRtest2'
                slabel='CR'
                havecr=6
                runtitle='COSMO'
                snlabel=r'$f_{mec}=1$'
                dclabel='3e28'
                resolabel='lr'
                Fcal=0.1
                iavesfr = 1.0
                subdir='/output/'
                cosmo=1
                maindir='/oasis/'
                h0=0.702
        
        print 'rundir', rundir
        if rundir=='none':
                cosmo=1
                cosmodata=cosmichalo(runtodo)
                rundir=cosmodata['rundir']
                Sheaform=cosmodata['Sheaform']
                crlabel=cosmodata['crlabel']
                suffixadd=cosmodata['suffixadd']
                snapsep=cosmodata['snapsep']
                icolor=cosmodata['icolor'] 
                snumadd=cosmodata['snumadd']
                h0=cosmodata['h0']
                fileno=cosmodata['fileno']
                rundir=cosmodata['rundir']
                subdir=cosmodata['subdir']
                halostr=cosmodata['halostr']
                beginno=cosmodata['beginno']
                finalno=cosmodata['finalno']
                multifile=cosmodata['multifile']
                halocolor=cosmodata['halocolor'] 
                labelname=cosmodata['labelname']
                xmax_of_box=cosmodata['xmax_of_box']
                firever=cosmodata['firever']
                usepep=cosmodata['usepep']
                maindir=cosmodata['maindir']
                highres=cosmodata['highres']
                print 'rundir', rundir
                

        if runtitle=='SMC':
                Mhini = 2.0e10
                Msini = 2.04e8
        if runtitle=='SBC':
                Mhini = 2.1e11
                Msini = 7.1e9
        if runtitle=='MW':
                Mhini = 1.5e12
                Msini = 6.2e10

        Nsnapstring = str(Nsnap)
        if Nsnap<10:
                Nsnapstring = '00'+str(Nsnap)
        elif Nsnap<100:
                Nsnapstring = '0'+str(Nsnap)
        the_snapdir = '/home/tkc004/'+maindir+'/'+rundir+'/'+subdir
        return {'snumadd':snumadd, 'crlabel':crlabel, 'strlabel':strlabel, 'newlabel':newlabel,\
 'galcen':galcen, 'snapsep':snapsep, 'stron':stron,'correctIa':correctIa, 'multifile':multifile,\
'the_prefix':the_prefix,'the_suffix':the_suffix, 'M1speed':M1speed,'Rvir':Rvir,\
'haveB':haveB,'havemetal':havemetal,'kappa':kappa,'exceptcool':exceptcool,\
'initsnap':initsnap, 'firever':firever, 'beginno':beginno,'finalno':finalno,\
'usepep':usepep,'halostr':halostr,'h0':h0,'withinRv':withinRv,'usesnaplist':usesnaplist,\
'color':color,'maindir':maindir,'subdir':subdir,'cosmo':cosmo, 'Nsnapstring':Nsnapstring, 'rundir':rundir,\
'highres':highres,'Mhini':Mhini,'Msini':Msini,\
'runtitle':runtitle,'slabel':slabel,'snlabel':snlabel,'dclabel':dclabel,'resolabel':resolabel,\
'the_snapdir':the_snapdir,'havecr':havecr,'Fcal':Fcal,'iavesfr':iavesfr,'timestep':timestep}




def find_enclosed_radius(Srhalf,Slight,Slightneed,ri,si,ncountmax,relerr):
        rlhalf=ri
        slhalf=si
        ncount=0
        while np.absolute((slhalf-Slightneed)/Slightneed)>relerr:
                slhalf = np.sum(Slight[Srhalf<rlhalf])
                rlhalf *= Slightneed/slhalf
                ncount += 1
#                print 'rlhalf, slhalf', rlhalf, slhalf
                if ncount>ncountmax:
                        break
        rhalf = rlhalf
        shalf = slhalf
        return rhalf, shalf

def readhalos(dirname,halostr,hubble=0.702):
        halofile=open(dirname+'/halos/halo_000'+halostr+'.dat','r')
        halofile.readline()
        halofile.readline()
        dars = halofile.readlines()
        halofile.close()
        zlist=[]
        idlist=[]
        mvirlist=[]
        fMhireslist=[]
        mgaslist=[]
        mstarlist=[]
        rvirlist=[]
        haloXl=[]
        haloYl=[]
        haloZl=[]
        for line in dars:
                xsd = line.split()
                zlist.append(float(xsd[0]))
                idlist.append(int(xsd[1]))
                mvirlist.append(float(xsd[4])/hubble)
                haloXl.append(float(xsd[6])/hubble)
                haloYl.append(float(xsd[7])/hubble)
                haloZl.append(float(xsd[8])/hubble)
                fMhireslist.append(float(xsd[38]))
                mgaslist.append(float(xsd[54])/hubble)
                mstarlist.append(float(xsd[74])/hubble)
                rvirlist.append(float(xsd[12])/hubble)
        zlist=np.array(zlist)
        a_scale = 1./(1.+zlist)
        idlist=np.array(idlist)
        mvirlist=np.array(mvirlist)
        fMhireslist=np.array(fMhireslist)
        #print 'a_scale', a_scale
        #print 'haloXl', haloXl
        haloX=np.array(haloXl)*a_scale
        haloY=np.array(haloYl)*a_scale
        haloZ=np.array(haloZl)*a_scale
        #print 'haloX', haloX
        mgaslist=np.array(mgaslist)
        mstarlist=np.array(mstarlist)
        rvirlist=np.array(rvirlist)
        return {'k':1,'mv':mvirlist,'fM':fMhireslist,'haloX':haloX,'haloY':haloY,'haloZ':haloZ,'mg':mgaslist,'ms':mstarlist,'rv':rvirlist,'z':zlist,'id':idlist};       


def reff_ell(H, edges, haloX, haloY, haloZ, nobin, den):
        NominS=[]
        ellX=[]
        ellY=[]
        ellZ=[]
        Spos=[]
        Sposx=[]
        Sposy=[]
        Sposz=[]
        errlist=[]
        SXi=haloX
        SYi=haloY
        SZi=haloZ
        nit=0
        massrat=0
        upno=1000
        upmr=0.45
        dmr=0.55
        nopa=float(den*200)
        NominS=np.append(NominS,nopa)
        nopau=float(den*400)
        nopad=float(den*10)
        inell=0
        while (nit<1000 and (massrat<dmr or massrat>upmr)):
                nit+=1
                if inell==1:
                        nopa=(nopau+nopad)/2.
                NominS=np.append(NominS,nopa)

                print 'nopa', nopa
                print 'len(Sposx)', len(Sposx)
                Spos=[]
                Sposx=[]
                Sposy=[]
                Sposz=[]
                totalno=0
                for i in range(nobin):
                        for j in range(nobin):
                                for k in range(nobin):
                                        if (H[i,j,k]>nopa):
                                                Sposx=np.append(Sposx,edges[0][i])
                                                Sposy=np.append(Sposy,edges[1][j])
                                                Sposz=np.append(Sposz,edges[2][k])
                                                totalno+=H[i,j,k]
                print 'len(Sposx)', len(Sposx)
                if ((len(Sposx)>3 and len(Sposx)<upno) or inell==1):
                        inell=1
                        #print 'fitting ell'
                        points = np.vstack(((Sposx+haloX),(Sposy+haloY),(Sposz+haloZ))).T
                        A, centroid = mvee(points)
                        print 'centroid', centroid
                        SXi=centroid[0]
                        SYi=centroid[1]
                        SZi=centroid[2]
                        U, D, V = la.svd(A)
                        rx, ry, rz = 1./np.sqrt(D)
                        u, v = np.mgrid[0:2*pi:20j, -pi/2:pi/2:10j]

                        def ellipse(u,v):
                            x = rx*cos(u)*cos(v)
                            y = ry*sin(u)*cos(v)
                            z = rz*sin(v)
                            return x,y,z

                        edgespoints=[]
                        E = np.dstack(ellipse(u,v))
                        E = np.dot(E,V) + centroid
                        x, y, z = np.rollaxis(E, axis = -1)
                        err=0
                        errlist=np.append(errlist,0)
                        inV=la.inv(V)
                        for i in range(1,len(edges[0])):
                                for j in range(1,len(edges[1])):
                                        for k in range(1, len(edges[2])):
                                                edgespoints=np.append(edgespoints,(((edges[0][i]+edges[0][i-1])/2-centroid[0]+haloX),((edges[1][j]+edges[1][j-1])/2-centroid[1]+haloY),((edges[2][k]+edges[2][k-1])/2-centroid[2]+haloZ)))
                        edgespoints=np.matrix(edgespoints.reshape((len(edgespoints)/3,3)))
                        rotback=np.dot(edgespoints,inV).T
                        sumHcrit=0
                        for i in range(0,len(edges[0])-1):
                                for j in range(0,len(edges[1])-1):
                                        for k in range(0, len(edges[2])-1):
                                                n=i*(len(edges[2])-1)*(len(edges[2])-1)+j*(len(edges[2])-1)+k
                                                if (np.power(rotback[0,n]/rx,2)+np.power(rotback[1,n]/ry,2)+np.power(rotback[2,n]/rz,2)<1):
                                                        sumHcrit=sumHcrit+H[i,j,k]
                        sumH=np.sum(H)
                        massrat=sumHcrit/sumH
                        print 'sum in ell/all', massrat
                        del x, y, z, E, points, A, centroid,U, D, V,u, v, inV, rotback
                        massratm=massrat
                        nopam=nopa
                        if massratm>upmr:
                                nopad=nopa*1.2
                        if massratm<dmr:
                                nopau=nopa*0.8
                        if (nit>3):
#                               print 'NominS', NominS[nit-1], NominS[nit-3]
                                if (np.abs(NominS[nit-1]-NominS[nit-3])<0.0001*NominS[nit-1]):
                                        err=5
                                        errlist=np.append(errlist,5)
                                        break
                        if (sumHcrit<1e-8):
                                err=6
                                errlist=np.append(errlist,6)
                                break
                elif (inell==0):
                        if(len(Sposx)<3):
                                nopa=nopa*0.9-1
                                nopau=nopa
                        if(len(Sposx)>upno):
                                nopa=nopa*1.1+1
                                nopad=nopa
        return SXi, SYi, SZi, rx, ry, rz, err, massrat, Sposx


def cosmichalo(runtodo):
        multifile='n'
        rundir=''
        subdir='hdf5/'
        beginno=400
        finalno=441
        halocolor='k'
        labelname='not set'
        xmax_of_box=40.0
        halostr='00'
        firever=1
        usepep=0
        maindir='scratch'
        highres=0
        fileno=4
        h0=0.702
        snumadd=0
        icolor=0
        snapsep=1
        suffixadd=''
        crlabel=''
        Sheaform=0
        if (runtodo=='m09'):
                rundir='m09_hr_Dec16_2013/'
                halostr='00'
        if (runtodo=='m10'):
                rundir='m10_hr_Dec9_2013/'
                halostr='00'
                subdir='/'
        if (runtodo=='m11'):
                rundir='m11_hhr_Jan9_2013/'
                halostr='01'
                halocolor='m'
                labelname='m11'
        if (runtodo=='m12v'):
                rundir='m12v_mr_Dec5_2013_3/'
                halo_to_do=[0]
        if (runtodo=='B1'):
                rundir='B1_hr_Dec5_2013_11/'
                halo_to_do=[0]
        if (runtodo=='m12qq'):
                rundir='m12qq_hr_Dec16_2013/'
                halo_to_do=[0]
                usepep=1
        if (runtodo=='383'):
                rundir='hydro_l8l_z3_383/'
                halostr='00'
                halocolor='y'
                labelname='m11h383'
        if (runtodo=='476'):
                rundir='hydro_l10l_z3_476/'
                halostr='00'
                halocolor='c'
                labelname='m10h476'
        if (runtodo=='553'):
                rundir='hydro_l10l_z3_553/'
                halostr='00'
                halocolor='r'
                labelname='m10h553'
        if (runtodo=='573'):
                rundir='hydro_l10l_z3_573/'
                halostr='00'
                halocolor='g'
                labelname='m10h573'
        if (runtodo=='f383'):
                rundir='FIRE_2_0_or_h383_criden1000_noaddm_sggs/'
                halostr='00'
                subdir='output'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                halocolor='y'
                icolor=0.85
                labelname='m11c'
                firever=2
                usepep=0
                maindir='oasis'
        if (runtodo=='f383sn199'):
                rundir='FIRE_2_0_h383_gasstrip_sn199/'
                halostr='00'
                subdir='output'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                halocolor='y'
                labelname='m11c'
                firever=2
                maindir='oasis'
        if (runtodo=='f383sn278'):
                rundir='FIRE_2_0_h383_gasstrip_sn278/'
                halostr='00'
                subdir='output'
                beginno=279
                finalno=600
                xmax_of_box=40.0
                halocolor='y'
                labelname='m11c'
                firever=2
                maindir='oasis'
                usepep=1
        if (runtodo=='f383sn472'):
                rundir='FIRE_2_0_h383_gasstrip_sn472/'
                halostr='00'
                subdir='output'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                halocolor='y'
                labelname='m11c'
                firever=2
                maindir='oasis'
        if (runtodo=='f476z0025'):
                rundir='FIRE_2_0_or_h476_gasstrip_z0_025/'
                halostr='00'
                subdir='output'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                halocolor='y'
                labelname='m11b'
                firever=2
                maindir='oasis'
        if (runtodo=='f476z0025'):
                rundir='FIRE_2_0_or_h476_gasstrip_z0_025/'
                halostr='00'
                subdir='output'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                halocolor='y'
                labelname='m11b'
                firever=2
                maindir='oasis'
        if (runtodo=='f476z005'):
                rundir='FIRE_2_0_or_h476_gasstrip_z0_05/'
                halostr='00'
                subdir='output'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                halocolor='y'
                labelname='m11b'
                firever=2
                maindir='oasis'
        if (runtodo=='f476z01'):
                rundir='FIRE_2_0_or_h476_gasstrip_z0_1/'
                halostr='00'
                subdir='output'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                halocolor='y'
                labelname='m11b'
                firever=2
                maindir='oasis'
        if (runtodo=='f476z015'):
                rundir='FIRE_2_0_or_h476_gasstrip_z0_15/'
                halostr='00'
                subdir='output'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                halocolor='y'
                labelname='m11b'
                firever=2
                maindir='oasis'
        if (runtodo=='f476z02'):
                rundir='FIRE_2_0_or_h476_gasstrip_z0_2/'
                halostr='00'
                subdir='output'
                beginno=487
                finalno=600
                xmax_of_box=40.0
                halocolor='y'
                labelname='m11b'
                firever=2
                maindir='oasis'
        if (runtodo=='f383_hv'):
                rundir='FIRE_2_0_or_h383_hv/'
                halostr='00'
                subdir='output'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                halocolor='y'
                labelname='m11c_hv'
                firever=2
                maindir='oasis'
                usepep=0
                highres=1
        if (runtodo=='fm11'):
                rundir='FIRE_2_0_m11/'
                halostr='01'
                subdir='output'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                multifile='y'
                halocolor='m'
                labelname='m11'
                icolor=0.7
                maindir='oasis'
                firever=2
        if (runtodo=='f573'):
                rundir='FIRE_2_0_or_h573_criden1000_noaddm_sggs/'
        #       halostr='00'
                halostr='00'
                subdir='output'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                halocolor='g'
                labelname='m10z'
                icolor=0.01
                firever=2
                maindir='oasis'

        if (runtodo=='f573z0'):
                rundir='FIRE_2_0_h573_2_21/'
                halostr='00'
                subdir='output'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                halocolor='g'
                labelname='m10z'
                icolor=0.01
                firever=2
                maindir='oasis'

        if (runtodo=='f573cr'):
                rundir='FIRE_2_0_h573_cr_2_21_M1_z0/'
                halostr='00'
                #halostr='03'
                subdir='output'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                halocolor='g'
                labelname='m10z'
                icolor=0.01
                firever=2
                maindir='/oasis/'
        if (runtodo=='f573crdc27'):
                rundir='FIRE_2_0_h573_cr_2_21_M1_z0_dc27/'
                halostr='00'
                #halostr='03'
                subdir='output'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                halocolor='g'
                labelname='m10z'
                icolor=0.01
                firever=2
                maindir='/oasis/'
        if (runtodo=='f573_hv'):
                rundir='FIRE_2_0_or_h573_hv/'
                halostr='00'
                subdir='output'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                halocolor='g'
                labelname='m10z_hv'
                firever=2
                maindir='oasis'
                usepep=1
                highres=1
                multifile='y'

        if (runtodo=='f476'):
                rundir='FIRE_2_0_or_h476_noaddm/'
                subdir='output'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                halocolor='c'
                icolor=0.5
                labelname='m11b'
                firever=2
                halostr='00'
                usepep=0
                maindir='/oasis/'

        if (runtodo=='fm11c'):
                rundir='m11c_res2100/'
                subdir='output'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                halocolor='c'
                icolor=0.5
                labelname='m11c'
                firever=2
                halostr='00'
                usepep=0
                maindir='/oasis/extra'
                snumadd=1
                
        if (runtodo=='m11bhv'):
                rundir='m11b_res260_FIRE/'
                halostr='02'
                subdir='output'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                halocolor='c'
                icolor=0.5
                labelname='m11bhv'
                firever=2
                maindir='oasis'

        if (runtodo=='f553'):
                rundir='FIRE_2_0_or_h553_criden1000_noaddm_sggs/'
                halostr='00'
                subdir='output'
                xmax_of_box=40.0
                beginno=100
                finalno=600
                halocolor='r'
                labelname='m11a'
                icolor=0.34
                firever=2
                maindir='oasis'
                usepep=0
                

        if (runtodo=='1146'):
                rundir='hydro_l10l_z3_1146_ssl/'
                halostr='01'
                halocolor='b'
                labelname='m10h1146'
        if (runtodo=='f1146m'):
                rundir='FIRE_2_0_or_h1146_mar/'
                halostr='01'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m10f1146_mar'
                firever=2
                maindir='oasis'
        if (runtodo=='f1146'):
                rundir='FIRE_2_0_or_h1146_criden1000_noaddm_sggs/'
                halostr='01'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m10y'
                firever=2
                maindir='oasis'
        if (runtodo=='f1146_hv'):
                rundir='FIRE_2_0_or_h1146_hv/'
                halostr='03'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m10yhr'
                firever=2
                maindir='oasis'
                highres=1
                multifile='y'

        if (runtodo=='fm12z'):
                maindir='oasis'
                rundir='m12z_res33000/'
                halostr='00'
                subdir='output'
                beginno=100
                finalno=440
                xmax_of_box=40.0
                multifile='y'
                halocolor='m'
                labelname='m12z'
                firever=2
        if (runtodo=='fm12c'):
                maindir='oasis'
                rundir='m12c_res56000/'
                halostr='02'
                subdir='output'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                multifile='n'
                halocolor='m'
                labelname='m12c'
                firever=2
        if (runtodo=='fm12i'):
                maindir='oasis'
                rundir='m12i_res7000/'
                halostr='00'
                subdir='output'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                multifile='n'
                halocolor='k'
                labelname='m12i'
                firever=2
                snumadd=1

        if (runtodo=='fm12imd'):
                rundir='m12i_res7100/'
                halostr='00'
                subdir='output'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                multifile='y'
                halocolor='k'
                labelname='m12i'
                firever=2
                maindir='oasis/metal_diffusion'
                highres=4
                crlabel='md'

        if (runtodo=='fm12fmd'):
                rundir='m12f_res7100/'
                halostr='00'
                beginno=600
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m12f'
                firever=2
                multifile='y'
                fileno=4
                maindir='oasis/metal_diffusion'
                usepep=1
                cosmo=1
                Sheaform=1
                highres=4
                crlabel='md'

        if (runtodo=='fm12_tmp'):
                maindir='scratch'
                rundir='m12_tmp/'
                halostr='00'
                subdir='output'
                beginno=100
                finalno=224
                xmax_of_box=40.0
                multifile='n'
                halocolor='k'
                labelname='m12i_M1CR'
                firever=2
                snumadd=0
        if (runtodo=='f1297'):
                rundir='FIRE_2_0_or_h1297_criden1000_noaddm_sggs/'
                halostr='03'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                labelname='m10f1297'
                subdir='output'
                firever=2
        if (runtodo=='1297'):
                rundir='hydro_l10l_z3_1297/'
                halostr='02'
        if (runtodo=='d573'):
                rundir='bare_l10l_z3_573/'
                halostr='00'
                beginno=100
                finalno=440
                xmax_of_box=40.0
                labelname='d573'
                subdir='output'
                firever=1
        if (runtodo=='dm11a_hv'):
                rundir='m11a_res260_dmo/'
                finalno=600
                xmax_of_box=40.0
                labelname='dm11a'
                subdir='output'
                firever=2
                maindir='oasis'
                usepep=0
                halostr='01'
        if (runtodo=='m11bhv'):
                rundir='m11b_res260_FIRE/'
                halostr='00'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                labelname='m11bhr'
                subdir='output'
                firever=2
                maindir='oasis'
                usepep=0
        if (runtodo=='dm11b_hv'):
                rundir='m11b_res260_dmo/'
                halostr='00'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                labelname='dm11b'
                subdir='output'
                firever=2
                maindir='oasis'
                usepep=0
        if (runtodo=='dm11bhv6rv'):
                rundir='m11b_res260_dmo_6rv/'
                halostr='00'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                labelname='dm11b'
                subdir='output'
                firever=2
                maindir='scratch'
                usepep=1
        if (runtodo=='dm11bhvf6rv'):
                rundir='m11b_res260_dmo_f6rv/'
                halostr='00'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                labelname='dm11b'
                subdir='output'
                firever=2
                maindir='scratch'
                usepep=1
        if (runtodo=='dm11bhvf8rvc'):
                rundir='m11b_res260_dmo_f8rv_convex/'
                halostr='00'
                beginno=100
                finalno=600
                xmax_of_box=40.0
                labelname='dm11b'
                subdir='output'
                firever=2
                maindir='oasis'
                usepep=1
        if (runtodo=='d553'):
                rundir='bare_l10l_z3_553/'
                halostr='00'
                beginno=100
                finalno=440
                xmax_of_box=40.0
                labelname='d553'
                subdir='output'
                firever=1
        if (runtodo=='d476'):
                rundir='bare_l10l_z3_476/'
                halostr='00'
                beginno=100
                finalno=440
                xmax_of_box=40.0
                labelname='d476'
                subdir='output'
                firever=1
                firever=1
        if (runtodo=='d383'):
                rundir='bare_l8l_z3_383/'
                halostr='00'
                beginno=100
                finalno=440
                xmax_of_box=40.0
                labelname='d383'
                subdir='output'
                firever=1
        if (runtodo=='df61'):
                rundir='FIRE_2_0_or_h61_dmo/'
                halostr='00'
                beginno=100
                finalno=600
                xmax_of_box=60.0
                labelname='dm11f61'
                subdir='output'
                firever=2
                multifile='y'
        if (runtodo=='f61'):
                rundir='FIRE_2_0_or_h61/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=60.0
                halocolor='b'
                labelname='m11f'
                firever=2
                maindir='oasis'
                multifile='y'
                icolor=0.93
        if (runtodo=='f46'):
                rundir='FIRE_2_0_or_h46/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=60.0
                halocolor='y'
                labelname='m11g'
                firever=2
                maindir='oasis'
                multifile='y'
        if (runtodo=='fm11e'):
                rundir='m11e_res7000/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m11e'
                firever=2
                h0=0.68
                maindir='oasis'

        if (runtodo=='fm11emd'):
                rundir='m11e_res7100/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m11e'
                firever=2
                h0=0.68
                maindir='oasis/metal_diffusion'
                highres=4
                crlabel='md'


        if (runtodo=='fm11h'):
                rundir='m11h_res7000/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m11h'
                firever=2
                h0=0.68
                maindir='oasis'

        if (runtodo=='fm11hmd'):
                rundir='m11h_res7100/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m11h'
                firever=2
                h0=0.68
                maindir='oasis/metal_diffusion'
                highres=4
                crlabel='md'

        if (runtodo=='fm11i'):
                rundir='m11i_res7000/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m11i'
                firever=2
                h0=0.68
                maindir='oasis'

        if (runtodo=='fm11imd'):
                rundir='m11i_res7100/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m11i'
                firever=2
                h0=0.68
                maindir='oasis/metal_diffusion'
                highres=4
                crlabel='md'


        if (runtodo=='fm11d'):
                rundir='m11d_res7000/'
                halostr='01'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m11d'
                firever=2
                h0=0.68
                maindir='oasis'

        if (runtodo=='fm11dmd'):
                rundir='m11d_res7100/'
                halostr='01'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m11d'
                firever=2
                h0=0.68
                maindir='oasis/metal_diffusion'
                highres=4
                crlabel='md'


        if (runtodo=='fm09'):
                rundir='m09/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='fm09'
                firever=2
                maindir='oasis/extra'
                multifile='y'
                fileno=8
                usepep=1
        if (runtodo=='fm10q'):
                rundir='m10q/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m10q'
                firever=2
                maindir='oasis/extra'
                multifile='y'
                fileno=8
                usepep=1

        if (runtodo=='fm10qmd'):
                rundir='m10q_res250/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m10q'
                firever=2
                maindir='oasis/metal_diffusion'
                multifile='y'
                fileno=8
                usepep=0
                highres=4
                crlabel='md'

        if (runtodo=='fm10v'):
                rundir='m10v/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='fm10v'
                firever=2
                maindir='oasis/extra'
                multifile='y'
                fileno=8
                usepep=1


        if (runtodo=='fm10vmd'):
                rundir='m10v_res250/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m10v'
                firever=2
                maindir='oasis/metal_diffusion'
                multifile='y'
                fileno=8
                usepep=1
                highres=4
                crlabel='md'


        if (runtodo=='fm10qCRK100'):
                rundir='m10q_mass250_MHDCR_K100_M1_tkFIX/'
                halostr='01'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='fm10qCRK100'
                firever=2
                maindir='scratch/'
                multifile='n'
                fileno=4
                usepep=0

        if (runtodo=='fm10qCRK1000'):
                rundir='m10q_mass250_MHDCR_K1000_M1_tkFIX/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='fm10qCRK1000'
                firever=2
                maindir='scratch/'
                multifile='n'
                fileno=4
                usepep=0


        if (runtodo=='fm11q'):
                rundir='m11q/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m11q'
                firever=2
                maindir='oasis/extra'
                multifile='y'
                fileno=4
                usepep=1



        if (runtodo=='fm11qmd'):
                rundir='m11q_res880/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m11q'
                firever=2
                maindir='oasis/metal_diffusion'
                multifile='y'
                fileno=4
                usepep=1
                Sheaform=1
                highres=4
                crlabel='md'

        if (runtodo=='fm11qCRK100'):
                rundir='m11q_mass7000_MHDCR_K100_M1_tkFIX/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='fm11qCRK100'
                firever=2
                maindir='scratch/'
                multifile='n'
                fileno=4
                usepep=1
                highres=4


        if (runtodo=='fm11qCRK1000'):
                rundir='m11q_mass7000_MHDCR_K1000_M1_tkFIX/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='fm11qCRK1000'
                firever=2
                maindir='scratch/'
                multifile='n'
                fileno=4
                usepep=1

        if (runtodo=='fm11qstd'):
                rundir='std_m11qFIRE2/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='k'
                labelname='fm11qstd'
                firever=2
                maindir='scratch/'
                multifile='n'
                fileno=4
                usepep=1




        if (runtodo=='fm11v1'):
                rundir='m11v/'
                halostr='01'
                beginno=600
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m11v1'
                firever=2
                maindir='oasis/extra'
                usepep=1


        if (runtodo=='fm11v'):
                rundir='m11v/'
                halostr='00'
                beginno=600
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m11v'
                firever=2
                maindir='oasis/extra'
                usepep=1
        if (runtodo=='fm11v2'):
                rundir='m11v/'
                halostr='02'
                beginno=600
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m11v2'
                firever=2
                maindir='oasis/extra'
                usepep=1
        if (runtodo=='fm11v3'):
                rundir='m11v/'
                halostr='03'
                beginno=600
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m11v3'
                firever=2
                maindir='oasis/extra'
                usepep=1
        if (runtodo=='fm11v4'):
                rundir='m11v/'
                halostr='04'
                beginno=600
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m11v4'
                firever=2
                maindir='oasis/extra'
                usepep=1

        if (runtodo=='fm12f'):
                rundir='m12f_ref13/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m12fhr'
                firever=2
                multifile='y'
                fileno=4
                maindir='oasis/extra'
                usepep=1

#        if (runtodo=='fm12i'):
 #               rundir='m12i_ref13/'
  #              halostr='00'
   #             beginno=100
    #            finalno=600
     #           subdir='output'
      #          xmax_of_box=40.0
       #         halocolor='b'
         #       labelname='fm12i_ref13'
        #        firever=2
          #      multifile='y'
           #     fileno=4
            #    maindir='oasis/extra'
             #   usepep=1

        if (runtodo=='fm12m'):
                rundir='m12m_ref13/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='fm12m'
                firever=2
                multifile='y'
                fileno=4
                maindir='oasis/extra'
                usepep=1

        if (runtodo=='fm12b'):
                rundir='m12b_ref12/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='m12b'
                firever=2
                maindir='oasis/extra'
                usepep=1


#        if (runtodo=='fm12c'):
 #               rundir='m12c_ref12/'
  #              halostr='00'
   #             beginno=100
    #            finalno=600
     #           subdir='output'
      #          xmax_of_box=40.0
       #         halocolor='b'
         #       labelname='fm12c_ref12'
        #        firever=2
          #      maindir='oasis/extra'
           #     usepep=1

        if (runtodo=='fm12q'):
                rundir='m12q_ref12/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='fm12q'
                firever=2
                maindir='oasis/extra'
                usepep=1

        if (runtodo=='b10lv'):
                rundir='bare_l8l_10lv/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='b10lv'
                firever=1
                maindir='scratch'
                usepep=0


        if (runtodo=='b11lv'):
                rundir='bare_l8l_11lv/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='b11lv'
                firever=2
                maindir='scratch'
                usepep=1

        if (runtodo=='m12ihrcr_b_70'):
                rundir='/m12i_mass7000_MHDCR_tkFIX/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                highres=2
                beginno=100
                finalno=601
                snapsep=10
                multifile=0
                labelname='m12ihr'
                halostr='00'
                halocolor='k'
                firever=2
                crlabel='cr70'


        if (runtodo=='m12ihrcr_700'):
                rundir='/m12i_mass7000_MHDCR_tkFIX/cr_700/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                highres=3
                usepep=0
                beginno=100
                finalno=601
                snapsep=10
                multifile=0
                labelname='m12ihr'
                halostr='00'
                halocolor='k'
                firever=2
                crlabel='cr700'


        if (runtodo=='m12ihrmhdcv'):
                rundir='/m12i_mass7000_MHDCR_tkFIX/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                beginno=0
                finalno=601
                snapsep=10
                highres=1
                multifile=0
                labelname='m12ihr'
                halostr='00'
                halocolor='k'
                firever=2
                crlabel='mhdcv'
                newlabel = 'MHD          no CR'


        if (runtodo=='m12icr_b_70'):
                rundir='/m12i_mass56000/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                highres=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=0
                labelname='m12i'
                halostr='00'
                halocolor='k'
                firever=2
                crlabel='cr70'

        if (runtodo=='m12icr_700'):
                rundir='/m12i_mass56000/cr_700/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                highres=3
                usepep=0
                beginno=520
                finalno=601
                snapsep=10
                multifile=0
                labelname='m12i'
                halostr='00'
                halocolor='k'
                firever=2
                crlabel='cr700'
                highres=3


        if (runtodo=='m12imhdcv'):
                rundir='/m12i_mass56000/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                beginno=0
                finalno=601
                snapsep=10
                highres=1
                multifile=0
                labelname='m12i'
                halostr='00'
                halocolor='k'
                firever=2
                crlabel='mhdcv'
                highres=1
                newlabel = 'MHD          no CR'


        if (runtodo=='m09mhdcv'):
                rundir='/m09/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=8
                firever=2
                beginno=0
                finalno=601
                snapsep=10
                multifile=1
                labelname='m09'
                halostr='00'
                halocolor='k'
                highres=1
                crlabel='mhdcv'
                newlabel = 'MHD          no CR'

        if (runtodo=='m09cr_b_70'):
                rundir='/m09/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=8
                firever=2
                beginno=0
                finalno=601
                snapsep=10
                multifile=1
                labelname='m09'
                halostr='00'
                halocolor='k'
                highres=2
                crlabel='cr70'

        if (runtodo=='m10vmhdcv'):
                rundir='/m10v/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=8
                firever=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=1
                labelname='m10v'
                halostr='04'
                halocolor='k'
                highres=1
                crlabel='mhdcv'
                newlabel = 'MHD          no CR'


        if (runtodo=='m10vcr_700'):
                rundir='/m10v/cr_700/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=8
                firever=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=1
                labelname='m10v'
                halostr='01'
                halocolor='k'
                highres=3
                crlabel='cr700'
                newlabel=r'MHD $\kappa$=3e29'

        if (runtodo=='m10vcr_b_70'):
                rundir='/m10v/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=8
                firever=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=1
                labelname='m10v'
                halostr='04'
                halocolor='k'
                highres=2
                crlabel='cr70'

        if (runtodo=='m11bcr_b_70'):
                rundir='/m11b/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                firever=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11b'
                halostr='00'
                halocolor='k'
                highres=2
                crlabel='cr70'



        if (runtodo=='m11bcr_700'):
                rundir='/m11b/cr_700/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                firever=2
                usepep=0
                beginno=520
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11b'
                halostr='00'
                halocolor='k'
                highres=3
                crlabel='cr700'
                newlabel=r'MHD $\kappa$=3e29'

        if (runtodo=='m11dcr_b_70'):
                rundir='/m11d/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                firever=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11d'
                halostr='00'
                halocolor='k'
                highres=2
                crlabel='cr70'

        if (runtodo=='m11dcr_700'):
                rundir='/m11d/cr_700/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                firever=2
                beginno=520
                finalno=601
                snapsep=10
                usepep=0
                multifile=0
                labelname='m11d'
                halostr='00'
                halocolor='k'
                highres=3
                crlabel='cr700'
                newlabel=r'MHD $\kappa$=3e29'


        if (runtodo=='m11hcr_700'):
                rundir='/m11h/cr_700/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                firever=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11h'
                halostr='00'
                halocolor='k'
                highres=3
                crlabel='cr700'
                newlabel=r'MHD $\kappa$=3e29'

        if (runtodo=='m11hcr_b_70'):
                rundir='/m11h/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                firever=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11h'
                halostr='00'
                halocolor='k'
                highres=2
                crlabel='cr70'

        if (runtodo=='m11fcr_700'):
                rundir='/m11f/cr_700/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=4
                firever=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=1
                labelname='m11f'
                halostr='00'
                halocolor='k'
                highres=3
                crlabel='cr700'
                newlabel=r'MHD $\kappa$=3e29'

        if (runtodo=='m11fcr_b_70'):
                rundir='/m11f/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=4
                firever=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=1
                labelname='m11f'
                halostr='00'
                halocolor='k'
                highres=2
                crlabel='cr70'

        if (runtodo=='m11gcr_700'):
                rundir='/m11g/cr_700/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=4
                firever=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=1
                labelname='m11g'
                halostr='00'
                halocolor='k'
                highres=3
                crlabel='cr700'
                newlabel=r'MHD $\kappa$=3e29'

        if (runtodo=='m11gcr_b_70'):
                rundir='/m11g/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=4
                firever=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=1
                labelname='m11g'
                halostr='00'
                halocolor='k'
                highres=2
                crlabel='cr70'

        if (runtodo=='f46'):
                rundir='/FIRE_2_0_or_h46/'
                subdir='/output/'
                maindir='/oasis/'
                nomultihdf5=4
                firever=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=1
                labelname='m11g'
                halostr='00'
                halocolor='k'
                highres=0
                crlabel='noCR'
                dclabel='m11g hydro'

        if (runtodo=='m11vcr_b_70'):
                rundir='/m11v/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                firever=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11v'
                halostr='00'
                halocolor='k'
                highres=2
                usepep=1
                crlabel='cr70'


        if (runtodo=='m11v1cr_b_70'):
                rundir='/m11v/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                firever=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11v1'
                halostr='01'
                halocolor='k'
                highres=2
                usepep=1
                crlabel='cr70'
        if (runtodo=='m11v2cr_b_70'):
                rundir='/m11v/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                firever=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11v2'
                halostr='02'
                halocolor='k'
                highres=2
                usepep=1
                crlabel='cr70'
        if (runtodo=='m11v3cr_b_70'):
                rundir='/m11v/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                firever=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11v'
                halostr='03'
                halocolor='k'
                highres=2
                usepep=1
                crlabel='cr70'
        if (runtodo=='m11v4cr_b_70'):
                rundir='/m11v/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                firever=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11v4'
                halostr='04'
                halocolor='k'
                highres=2
                usepep=1
                crlabel='cr70'

        if (runtodo=='m11bmhdcv'):
                rundir='/m11b/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                beginno=500
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11b'
                halostr='00'
                halocolor='r'
                firever=2
                highres=1
                crlabel='mhdcv'
                newlabel = 'MHD          no CR'

        if (runtodo=='m11cmhdcv'):
                rundir='/m11c/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                beginno=500
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11c'
                halostr='00'
                halocolor='r'
                firever=2
                highres=1
                crlabel='mhdcv'
                newlabel = 'MHD          no CR'




        if (runtodo=='m11dmhdcv'):
                rundir='/m11d/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                startno=0
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11d'
                halostr='00'
                halocolor='r'
                firever=2
                highres=1
                h0=0.68
                crlabel='mhdcv'
                newlabel = 'MHD          no CR'


        if (runtodo=='m11hmhdcv'):
                rundir='/m11h/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                startno=0
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11h'
                halostr='00'
                halocolor='r'
                firever=2
                highres=1
                h0=0.68
                crlabel='mhdcv'
                newlabel = 'MHD          no CR'


        if (runtodo=='m11fmhdcv'):
                rundir='/m11f/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=4
                beginno=500
                finalno=601
                snapsep=10
                multifile=1
                labelname='m11f'
                halostr='00'
                halocolor='r'
                firever=2
                highres=1
                crlabel='mhdcv'
                newlabel = 'MHD          no CR'


        if (runtodo=='m11gmhdcv'):
                rundir='/m11g/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=4
                startno=0
                finalno=601
                snapsep=10
                multifile=1
                labelname='m11g'
                halostr='00'
                halocolor='r'
                firever=2
                highres=1
                crlabel='mhdcv'
                newlabel = 'MHD          no CR'

        if (runtodo=='m11vmhdcv'):
                rundir='/m11v/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                startno=0
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11v'
                halostr='00'
                halocolor='r'
                firever=2
                highres=1
                crlabel='mhdcv'
                usepep=1
                newlabel = 'MHD          no CR'

        if (runtodo=='m11v1mhdcv'):
                rundir='/m11v/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                startno=0
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11v1'
                halostr='01'
                halocolor='r'
                firever=2
                highres=1
                crlabel='mhdcv'
                usepep=1
                newlabel = 'MHD          no CR'

        if (runtodo=='m11v2mhdcv'):
                rundir='/m11v/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                startno=0
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11v2'
                halostr='02'
                halocolor='r'
                firever=2
                highres=1
                crlabel='mhdcv'
                usepep=1
                newlabel = 'MHD          no CR'

        if (runtodo=='m11v3mhdcv'):
                rundir='/m11v/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                startno=0
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11v3'
                halostr='03'
                halocolor='r'
                firever=2
                highres=1
                crlabel='mhdcv'
                usepep=1
                newlabel = 'MHD          no CR'

        if (runtodo=='m11v4mhdcv'):
                rundir='/m11v/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                startno=0
                finalno=601
                snapsep=10
                multifile=0
                labelname='m11v4'
                halostr='04'
                halocolor='r'
                firever=2
                highres=1
                crlabel='mhdcv'
                usepep=1
                newlabel = 'MHD          no CR'


        if (runtodo=='m10qcr_b_70'):
                rundir='/m10q_mass250_MHDCR_M1_tkFIX/cr_b_70/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                firever=2
                beginno=500
                finalno=601
                snapsep=10
                multifile=0
                labelname='m10q'
                halostr='03'
                halocolor='k'
                highres=2
                crlabel='cr70'

        if (runtodo=='m10qmhdcv'):
                rundir='/m10q_mass250_MHDCR_M1_tkFIX/mhdcv/'
                subdir='/output/'
                maindir='/oasis/philruns/'
                nomultihdf5=1
                startno=0
                finalno=601
                snapsep=10
                multifile=0
                labelname='m10q'
                halostr='00'
                halocolor='r'
                firever=2
                highres=1
                crlabel='mhdcv'
                newlabel = 'MHD          no CR'



        return {'Sheaform':Sheaform,'crlabel':crlabel, 'suffixadd':suffixadd,'snapsep':snapsep, 'icolor':icolor, 'snumadd':snumadd,'h0':h0,'fileno':fileno,'rundir':rundir,'subdir':subdir,'halostr':halostr,'beginno':beginno,'finalno':finalno, 'multifile':multifile, 'halocolor':halocolor, 'labelname':labelname, 'xmax_of_box':xmax_of_box, 'firever':firever, 'usepep':usepep, 'maindir':maindir, 'highres':highres, 'newlabel':newlabel}

def muprofile(Slight, Sr, murl,bandneeded,UNITS_CGS=1,UNITS_SOLAR_BAND=0):
        Slight=np.array(Slight)
        Sr=np.array(Sr)
        murl=np.array(murl)
        Slin =[]
        for i in range(len(murl)):
                Slin.append(np.sum(Slight[Sr<murl[i]]))
        Slin=np.array(Slin)
        #print 'Slin', Slin
        Slbet = Slin[1:]-Slin[:-1]
        Areabet = np.pi*(np.square(murl[1:])-np.square(murl[:-1]))
        vmag = luminosity_to_magnitude(Slbet, \
        UNITS_CGS=UNITS_CGS, UNITS_SOLAR_BAND=UNITS_SOLAR_BAND,\
        BAND_NUMBER=bandneeded, \
        VEGA=0, AB=1 , \
        MAGNITUDE_TO_LUMINOSITY=0 )
        mul = np.log10(Areabet*1.425*1.425)*2.5+vmag+34.97
        magl = luminosity_to_magnitude(Slin, \
        UNITS_CGS=UNITS_CGS, UNITS_SOLAR_BAND=UNITS_SOLAR_BAND,\
        BAND_NUMBER=bandneeded, \
        VEGA=0, AB=1 , \
        MAGNITUDE_TO_LUMINOSITY=0 )
        return mul, magl



def reff_ell2d(H, edges, haloX, haloY, nobin, den):
        NominS=[]
        ellX=[]
        ellY=[]
        Spos=[]
        Sposx=[]
        Sposy=[]
        errlist=[]
        SXi=haloX
        SYi=haloY
        nit=0
        massrat=0
        upno=1000
        upmr=0.3
        dmr=0.7
        nopa=float(den)
        NominS=np.append(NominS,nopa)
        nopau=float(den*2)
        nopad=float(den*0.2)
        inell=0
        rx=0.
        ry=0.
        SXi=0.
        SYi=0.
        while (nit<100 and (massrat<dmr or massrat>upmr)):
                nit+=1
                if inell==1:
                        nopa=(nopau+nopad)/2.
                NominS=np.append(NominS,nopa)

                print 'nopa', nopa
                print 'len(Sposx)', len(Sposx)
                Spos=[]
                Sposx=[]
                Sposy=[]
                totalno=0
                for i in range(nobin):
                        for j in range(nobin):
                                if (H[i,j]>nopa):
                                        Sposx=np.append(Sposx,edges[0,i])
                                        Sposy=np.append(Sposy,edges[1,j])
                                        totalno+=H[i,j]
                print 'len(Sposx)', len(Sposx)
                print 'totalno', totalno
                if ((len(Sposx)>3 and len(Sposx)<upno) or inell==1):
                        inell=1
                        #print 'Sposx', Sposx
                        #print 'fitting ell'
                        points = np.vstack(((Sposx+haloX),(Sposy+haloY))).T
                        A, centroid = mvee(points)
                        print 'centroid', centroid
                        SXi=centroid[0]
                        SYi=centroid[1]
                        U, D, V = la.svd(A)
                        rx, ry = 1./np.sqrt(D)
                        u = np.mgrid[0:2*pi:20j]

                        def ellipse2d(u):
                            x = rx*cos(u)
                            y = ry*sin(u)
                            return x,y

                        edgespoints=[]
                        E = np.dstack(ellipse2d(u))
                        E = np.dot(E,V) + centroid
                        x, y= np.rollaxis(E, axis = -1)
                        err=0
                        errlist=np.append(errlist,0)
                        inV=la.inv(V)
                        for i in range(1,nobin):
                                for j in range(1,nobin):
                                        edgespoints=np.append(edgespoints,(((edges[0,i]+edges[0,i-1])/2-centroid[0]+haloX),((edges[1,j]+edges[1,j-1])/2-centroid[1]+haloY)))
                        edgespoints=np.matrix(edgespoints.reshape((len(edgespoints)/2,2)))
                        rotback=np.dot(edgespoints,inV).T
                        sumHcrit=0
                        for i in range(0,len(edges[0])-1):
                                for j in range(0,len(edges[1])-1):
                                        n=i*(len(edges[1])-1)+j
                                        if (np.power(rotback[0,n]/rx,2)+np.power(rotback[1,n]/ry,2)<1):
                                                sumHcrit=sumHcrit+H[i,j]
                        sumH=np.sum(H)
                        massrat=sumHcrit/sumH
                 #       print 'sum in ell/all', massrat
                        del x, y, E, points, A, centroid,U, D, V,u, inV, rotback
                        massratm=massrat
                        nopam=nopa
                        if massratm>upmr:
                                nopad=nopa*1.1
                        if massratm<dmr:
                                nopau=nopa*0.9
                        if (nit>3):
                #                print 'NominS', NominS[nit-1], NominS[nit-3]
                                if (np.abs(NominS[nit-1]-NominS[nit-3])<0.0001*NominS[nit-1]):
                                        err=5
                                        errlist=np.append(errlist,5)
                                        break
                        if (sumHcrit<1e-8):
                                err=6
                                errlist=np.append(errlist,6)
                                break
                elif (inell==0):
                        if(len(Sposx)<3):
                                nopa=nopa*0.97-1
                                nopau=nopa
                        if(len(Sposx)>upno):
                                nopa=nopa*1.03+1
                                nopad=nopa
        return SXi, SYi, rx, ry,  err, massrat, Sposx

def readtime(firever=2):
        file=open(programdir+'data/snapshot_times.txt','r')
        file.readline()
        file.readline()
        file.readline()
        dars = file.readlines()
        file.close()
        snap2list=[]
        a2list=[]
        time2list=[]
        red2list=[]
        for line in dars:
                xsd = line.split()
                snap2list.append(int(xsd[0]))
                a2list.append(float(xsd[1]))
                red2list.append(float(xsd[2]))
                time2list.append(float(xsd[3]))
        snap2list=np.array(snap2list)
        time2list=np.array(time2list)
        if firever==1:
                file=open(programdir+'/data/output_times.txt','r')
                dars = file.readlines()
                file.close()
                snaplist=[]
                alist=[]
                timelist=[]
                redlist=[]
                ncount=0
                for line in dars:
                        xsd = line.split()
                        alist.append(float(xsd[0]))
                        snaplist.append(ncount)
                        ncount+=1
                alist=np.array(alist)
                snaplist=np.array(snaplist)
                timelist=np.array(np.interp(alist,a2list,time2list))
                return {'snaplist':snaplist, 'timelist':timelist, 'alist':alist, 'redlist':redlist}
        if firever==2:
                return {'snaplist':snap2list, 'timelist':time2list, 'alist':a2list, 'redlist':red2list}

def outLgamma_nism(Grho,Neb,cregy): #Gadget unit (but physical) (no h; no comoving)
        Gn_in_cm_3 = (0.78)/proton_mass_in_g*Grho*1e10*solar_mass_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm
        tpi_in_sp = 1/Gn_in_cm_3/hadronicdecayrate
        Lgamma = cregy*1e10*solar_mass_in_g*km_in_cm*km_in_cm/tpi_in_sp*betapi/nopi_per_gamma #in erg/s
        return {'Lgamma':Lgamma, 'nism':Gn_in_cm_3}


def diffusionsol(Ecr0,kappa,time,rneed,offset=0.,dim=3.):
        ecr = Ecr0/np.power(2.*np.pi*(2.*kappa*time+offset*offset),dim/2.)*np.exp(-rneed*rneed/2./(2.*kappa*time+offset*offset))
        return ecr

def diffCRin(Ecr0,kappa,time,rneed,offset=0):
        Ecrin = Ecr0*(special.erf(rneed/np.sqrt(2.*(2.*kappa*time+offset*offset)))-rneed/np.power(np.pi*(kappa*time+0.5*offset*offset),1./2.)*np.exp(-rneed*rneed/2./(2.*kappa*time+offset*offset)))
        return Ecrin

def calNFW(rlist,reffp,meffp,zr=0): 
        rhocrit=127*(0.27*(1+zr)*(1+zr)*(1+zr)+0.73)
        xLam=-0.73/(0.27*(1+zr)*(1+zr)*(1+zr)+0.73)
        delta=18.0*np.pi*np.pi+82.0*(xLam)-39*(xLam)*(xLam)
        samlogMv=np.arange(7.0,14.0,0.05)
        logc=0.537+(1.025-0.537)*np.exp(-0.718*np.power(zr,1.08))+(0.024*zr-0.097)*(samlogMv-np.log10(0.7)-np.log10(1e12))
        cont=np.power(10,logc)
        Rvirn=np.power(np.power(10,samlogMv)*3.0/4.0/np.pi/rhocrit/delta,1.0/3.0)
#       finding the NFW profile (at z=0) that fits at reff and plot it
        #x=reffp/Rvirn
        #rho = rhocrit*delta/cont/cont/cont/x/(1/cont+x)/(1/cont+x)
        Rs = Rvirn/cont
        cdelta = cont*cont*cont*delta/3./(np.log(1+cont)-cont/(1+cont))
        #enclosed NFW halo mass within reffp with different total halo masses
        menc = 4*np.pi*rhocrit*cdelta*Rs*Rs*Rs*(np.log((Rs+reffp)/Rs)-reffp/(Rs+reffp))
        #infer the total halo mass with meffp
        logMveff = np.interp(meffp, menc, samlogMv)
        #calculate the enclosed mass with the total halo mass above
        logceff=0.537+(1.025-0.537)*np.exp(-0.718*np.power(zr,1.08))+(0.024*zr-0.097)*(logMveff-np.log10(0.7)-np.log10(1e12))
        conteff=np.power(10,logceff)
        Rvirneff=np.power(np.power(10,logMveff)*3.0/4.0/np.pi/rhocrit/delta,1.0/3.0)
        Rseff = Rvirneff/conteff
        cdeltaeff = conteff*conteff*conteff*delta/3./(np.log(1+conteff)-conteff/(1+conteff))
        #output the enclosed NFW halo mass in rlist that matches [reffp,meffp]
        mnfw = np.array(4*np.pi*rhocrit*cdeltaeff*Rseff*Rseff*Rseff*(np.log((Rseff+rlist)/Rseff)-rlist/(Rseff+rlist)))
        return {'mnfw':mnfw}

    
def findradwnism(Gextra, findradiusatnism):
    radlist = np.linspace(0.01,10,num=20)
    gasdenlist = radlist*0.0
    maxlength = 0.5
    for irad in range(len(radlist)-1):
        Gx = Gextra['p'][:,0]; Gy = Gextra['p'][:,1]; Gz = Gextra['p'][:,2];
        Gvx = Gextra['v'][:,0]; Gvy = Gextra['v'][:,1]; Gvz = Gextra['v'][:,2];
        Grho = Gextra['rho']; Gu = Gextra['u']; Gm = Gextra['m'];
        cutxy = ((Gx)*(Gx)+(Gy)*(Gy) > radlist[irad]*radlist[irad]) & ((Gx)*(Gx)+(Gy)*(Gy) < radlist[irad+1]*radlist[irad+1])
        cutz = (Gz)*(Gz) < maxlength*maxlength/4.
        cut = cutxy*cutz
        shellvol_in_cm3 = np.pi*(-np.power(radlist[irad],2)+np.power(radlist[irad+1],2))*kpc_in_cm*kpc_in_cm*kpc_in_cm*maxlength
        Gm_in_g=Gm[cut]*1e10*Msun_in_g
        Grho_in_g_cm_3 = np.sum(Gm_in_g)/shellvol_in_cm3
        Gnism_in_cm_3 = 1.0/proton_mass_in_g*Grho_in_g_cm_3
        gasdenlist[irad]=Gnism_in_cm_3
    radlist = np.flip(radlist)
    gasdenlist = np.flip(np.minimum.accumulate(gasdenlist))
    print 'radlist', radlist
    print 'gasdenlist', gasdenlist
    radneed = np.interp(findradiusatnism, gasdenlist, radlist)
    print 'radius that has nism = 1cm^-3', radneed
    return radneed
