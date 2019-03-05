import os
import numpy as np
import matplotlib
matplotlib.use('agg')
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
import matplotlib.colors
import matplotlib.cm
import matplotlib.pyplot as plt
import math
import h5py
import re
import sys
import glob
from numpy.linalg import inv
rcParams['figure.figsize'] = 8, 6
rcParams.update({'figure.autolayout': True})
rcParams.update({'font.size': 24})
rcParams['axes.unicode_minus'] = False
import matplotlib.patches as patches
rcParams['axes.linewidth'] = 2
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True
rcParams['xtick.major.width'] = 2
rcParams['ytick.major.width'] = 2
rcParams['xtick.minor.width'] = 2
rcParams['ytick.minor.width'] = 2
rcParams['text.usetex'] = True
import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='inferno', cmap=cmaps.inferno)
plt.register_cmap(name='plasma', cmap=cmaps.plasma)
from readsnap_samson import *
from Sasha_functions import *
from gadget_lib.cosmo import *
from samson_functions import *
from crtestfunction import *
#dirneed=['m12icr_b_70', 'm12imhdcv', 'm12icr_700', 'f553']
dirneed=[
'bwsmclrdc28mhd','bwmwmrdc28mhd','bwsbclrdc28mhd',
#'bwsmclrdc29','bwmwmrdc29','bwsbclrdc29',
#'bwsmclrdc28','bwmwmrdc28','bwsbclrdc28',
#'bwsmclr','bwmwmr','bwsbclrmhd',
#'bwsmclrmhd', 'bwmwmrmhd', 'bwsbclrmhd',
#'m11dmhdcv',\
#'m11bcr_b_70','m11bmhdcv','m11bcr_700',\
#'476',\
#'bwsmclr','bwsmclrdc0','bwsmclrdc27','bwsmclrdc28','bwsmclrdc29','bwsmclrstr','bwsmclrdc28mhd'
#'m09','m10','m11','m12v','B1','m12qq','383','476',\
#'fm10qmd','fm10vmd','fm11dmd','fm11emd', 'fm11hmd','fm11imd', 'fm11qmd','fm12fmd','fm12imd',\
#'m10qcr_b_70','m10vcr_b_70', 'm11bcr_b_70','m11dcr_b_70','m11vcr_b_70','m11v1cr_b_70',\
#'m12icr_b_70','m11fcr_b_70','m11hcr_b_70','m11gcr_b_70',\
#'m12fmhdcv','m12mmhdcv',\
#'m10qmhdcv','m10vmhdcv', 'm11bmhdcv','m11dmhdcv','m12imhdcv','m11fmhdcv','m11hmhdcv','m11gmhdcv',\
#'m12fcr_700','m12mcr_700',\
#'m10vcr_700', 'm11bcr_700','m11dcr_700','m11fcr_700','m11gcr_700','m11hcr_700','m12icr_700',
#'fm12m','fm12b','fm10q','fm11d','f553','f573','f476','f383','fm11q','fm11v','fm11v1','fm11v2','f1146','f46','f61',
]
galcen=0
hubble=0.702
#Nsnaplist=range(300,441,10)
#Nsnaplist=range(590,601,10)
#Nsnaplist=range(500,601,10)
Nsnaplist=range(450,500,1)
#Nsnaplist=[570,580,590,600]
#fmeat='test'
#fmeat='md'
#fmeat='all'
#fmeat='FIRE1'
#fmeat='FIRE2'
fmeat='dc28mhd'
#fmeat='hydro'
#fmeat='mhdcv'
#fmeat='cr_b_70'
#fmeat='cr_700'
#fmeat='crhydroave'
#fmeat='hydrosn600'
#Nsnaplist=[580,590,600]
tsep=10.0
needmsmv=1
needSFRMs=1
needmlms=1
avesnap=1

sfrl=np.array([])
mll =np.array([]) #massloading factor
mvl =np.array([])
mwl =np.array([]) #outflow rate
msl =np.array([])
runtodol = np.array([])
dcl = np.array([])
nsml = np.array([])
for runtodo in dirneed:
        for Nsnap in Nsnaplist:
                #haloinfo=cosmichalo(runtodo)
                haloinfo=outdirname(runtodo)
                rundir=haloinfo['rundir']
                maindir=haloinfo['maindir']
                print 'maindir', maindir
                subdir=haloinfo['subdir']
                halostr=haloinfo['halostr']
                snumadd=haloinfo['snumadd']
                usepep=haloinfo['usepep']
                hubble=haloinfo['h0']
                firever=haloinfo['firever']
                highres=haloinfo['highres']
                Rvir=haloinfo['Rvir']
                cosmo=haloinfo['cosmo']
                Ms = haloinfo['Msini']
                Mv = haloinfo['Mhini']
                if cosmo==1:
                        try:
                                if usepep==0:
                                        halosA = read_halo_history(rundir, halonostr=halostr,comoving=0,hubble=hubble,maindir=maindir,snumadd=snumadd)
                                else:
                                        halosA = read_halo_history_pep(rundir, Nsnap, singlesnap=1, firever=firever, halonostr=halostr, comoving=0, maindir=maindir)
                                redlist = halosA['redshift']
                                Rvirlist = halosA['R']
                                Mslist = halosA['Ms']
                                Mvlist = halosA['M']
                                #print 'Mslist', Mslist
                                #print 'Mvlist', Mvlist
                        except (IOError,KeyError):
                                continue
                        try:
                                Rvir = Rvirlist[-1]
                                Ms = Mslist[-1]
                                Mv = Mvlist[-1]
                        except TypeError:
                                Rvir=Rvirlist; Ms = Mslist; Mv = Mvlist;
                        
                Rup = Rvir*0.3
                Rdown = Rvir*0.2
                try:
                        emdata=gaswindphase(runtodo,Nsnap,rotface=0,withinr=100.0,zup=Rup,zdown=Rdown,userad=1)
                        sfrdata=outsfr(runtodo, Nsnap,tsep=tsep)
                        vr =  emdata['vr'] #in km/s
                        Gmass = emdata['Gmass'] #in 1e10Msun
                        print 'np.median(vr)', np.median(vr)
                        print 'np.sum(Gmass)', np.sum(Gmass)
                        SFR = sfrdata['SFR']# in Msun/yr
                        Nsm = sfrdata['Nsm']# in solar mass
                        print 'runtodo', runtodo
                        print 'dc', highres
                        print 'Gmass', np.sum(Gmass)
                        print 'SFR', SFR
                        print 'Rvir', Rvir
                        mwind = np.sum(Gmass*1e10*vr*km_in_cm)/(Rup-Rdown)/kpc_in_cm*yr_in_sec #Msun/yr
                        ml = mwind/SFR
                        msl=np.append(msl,Ms)
                        mvl=np.append(mvl,Mv)
                        nsml = np.append(nsml,Nsm)
                        mwl = np.append(mwl,mwind)
                        mll=np.append(mll,ml)
                        sfrl=np.append(sfrl,SFR)
                        runtodol= np.append(runtodol,runtodo)
                        dcl=np.append(dcl,highres)
                except (IOError,KeyError):
                        continue
                
print 'dcl', dcl
print 'dcl==0', dcl==0

if avesnap==1:
        print 'runtodol', runtodol
        print 'mvl', mvl
        print 'sfrl', sfrl
        cutinf = np.isfinite(mvl)
        mvlold=mvl[cutinf];mslold=msl[cutinf];mllold=mll[cutinf];sfrlold=sfrl[cutinf];
        runtodolold=runtodol[cutinf];mwlold=mwl[cutinf];dclold=dcl[cutinf];
        nsmoldl=nsml[cutinf]
        mvl=msl=mll=runtodol=mwl=dcl=sfrl=avensml=np.array([])
        for runtodo in dirneed:
                print 'runtodo', runtodo 
                print 'mvlold', mvlold[runtodolold==runtodo], 
                print 'ave', np.average(mvlold[runtodolold==runtodo])
                numoftimes = len(mvlold[runtodolold==runtodo])
                mvl=np.append(mvl,np.average(mvlold[runtodolold==runtodo]))
                msl=np.append(msl,np.average(mslold[runtodolold==runtodo]))
                mll=np.append(mll,np.average(mllold[runtodolold==runtodo]))
                mwl=np.append(mwl,np.average(mwlold[runtodolold==runtodo]))
                dcl=np.append(dcl,np.median(dclold[runtodolold==runtodo]))
                sfrl=np.append(sfrl,np.average(sfrlold[runtodolold==runtodo]))
                runtodol=np.append(runtodol,np.unique(runtodolold[runtodolold==runtodo]))
                avensml = np.append(avensml,np.average(nsmoldl[runtodolold==runtodo]))
        mlla = mwl/(avensml/numoftimes/tsep/1e6)
        print 'runtodol', runtodol
        print 'mvl', mvl
        print 'sfrl', sfrl
        print 'mlla', mlla

plt.plot(mvl[dcl==0],mll[dcl==0],marker='o',ls='none',color='y',label='Hydro')
plt.plot(mvl[dcl==1],mll[dcl==1],marker='s',ls='none',color='b',label='MHD')
plt.plot(mvl[dcl==2],mll[dcl==2],marker='^',ls='none',color='g',label=r'$\kappa$=3e28')
plt.plot(mvl[dcl==3],mll[dcl==3],marker='d',ls='none',color='k',label=r'$\kappa$=3e29')
plt.plot(mvl[dcl==4],mll[dcl==4],marker='*',ls='none',color='brown',label='metal diffusion')

plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$M_{\rm vir}\,[{\rm M}_{\odot}]$')
plt.ylabel(r'$\eta=\dot{M}_{\rm w}/\dot{M}_{*}$')
plt.legend(loc='best',fontsize=12)
filename='figures/Massloading_'+fmeat+'.pdf'
print 'filename', filename
plt.savefig(filename)
plt.clf()


if needmlms==1:
        plt.plot(msl[dcl==0],mll[dcl==0],marker='o',ls='none',color='y',label='Hydro')
        plt.plot(msl[dcl==1],mll[dcl==1],marker='s',ls='none',color='b',label='MHD')
        plt.plot(msl[dcl==2],mll[dcl==2],marker='^',ls='none',color='g',label=r'$\kappa$=3e28')
        plt.plot(msl[dcl==3],mll[dcl==3],marker='d',ls='none',color='k',label=r'$\kappa$=3e29')
        plt.plot(msl[dcl==4],mll[dcl==4],marker='*',ls='none',color='brown',label='metal diffusion')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'$M_{*}\,[{\rm M}_{\odot}]$')
        plt.ylabel(r'$\eta=\dot{M}_{\rm w}/\dot{M}_{*}$')
        plt.legend(loc='best',fontsize=12)
        filename='figures/MsMv_'+fmeat+'.pdf'
        print 'mvl', mvl
        print 'msl', msl
        print 'filename', filename
        plt.savefig(filename)
        plt.clf()

if needmsmv==1:
        plt.plot(mvl[dcl==0],msl[dcl==0],marker='o',ls='none',color='y',label='Hydro')
        plt.plot(mvl[dcl==1],msl[dcl==1],marker='s',ls='none',color='b',label='MHD')
        plt.plot(mvl[dcl==2],msl[dcl==2],marker='^',ls='none',color='g',label=r'$\kappa$=3e28')
        plt.plot(mvl[dcl==3],msl[dcl==3],marker='d',ls='none',color='k',label=r'$\kappa$=3e29')
        plt.plot(mvl[dcl==4],msl[dcl==4],marker='*',ls='none',color='brown',label='metal diffusion')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'$M_{\rm vir}\,[{\rm M}_{\odot}]$')
        plt.ylabel(r'$M_{*}\,[{\rm M}_{\odot}]$')
        plt.legend(loc='best',fontsize=12)
        filename='figures/MsMv_'+fmeat+'.pdf'
        print 'mvl', mvl
        print 'msl', msl
        print 'filename', filename
        plt.savefig(filename)
        plt.clf()



if needSFRMs==1:
        plt.plot(msl[dcl==0],sfrl[dcl==0],marker='o',ls='none',color='y',label='Hydro')
        plt.plot(msl[dcl==1],sfrl[dcl==1],marker='s',ls='none',color='b',label='MHD')
        plt.plot(msl[dcl==2],sfrl[dcl==2],marker='^',ls='none',color='g',label=r'$\kappa$=3e28')
        plt.plot(msl[dcl==3],sfrl[dcl==3],marker='d',ls='none',color='k',label=r'$\kappa$=3e29')
        plt.plot(msl[dcl==4],sfrl[dcl==4],marker='*',ls='none',color='brown',label='metal diffusion')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'$M_{*}\,[{\rm M}_{\odot}]$')
        plt.ylabel(r'SFR $[{\rm M}_{\odot}]$')
        plt.legend(loc='best',fontsize=12)
        filename='figures/SFRMs_'+fmeat+'.pdf'
        #print 'msl', msl
        #print 'sfrl', sfrl
        print 'filename', filename
        plt.savefig(filename)
        plt.clf()
