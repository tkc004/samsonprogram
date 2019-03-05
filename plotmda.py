import os
import pyfits
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
from readsnap_samson import *
from Sasha_functions import *
from gadget_lib.cosmo import *
from samson_functions import *
from samson_const import *

NewG = NewG_in_cgs/1e3
#NewG = 6.674e-11 #SI unit
Msun_in_kg = Msun_in_g/1e3
#Msun_in_kg = 2e30
kpc_in_m = kpc_in_cm/1e2
#kpc_in_m = 3.085e19

#dirneed=['f383']
#dirneed=['f476']
#dirneed=['f61']
#dirneed=['fm11']
#dirneed=['f573']
dirneed=['m11bcr_b_70','m11bmhdcv','m10qcr_b_70', 'm10qmhdcv']
#dirneed=['f573','f553','f476','fm11','f383','f61']
#dirneed=['1146']
#dirneed=['m11']
#dirneed=['573']
#dirneed=['573','553','476']
#dirneed=['553','476']
maxr=3.0
minr=1.0
Nsnap=600
#Nsnap=239 #t=5.0Gyr in FIRE2
#Nsnap=473 #t=11.0Gyr in FIRE2
#Nsnap=140
#Nsnap=239
#Nsnap=393
#Nsnap=600
#Nsnap=130 #t=2.0Gyr
#Nsnap=140 #t=2.2Gyr
#Nsnap=150 #t=2.4Gyr
#Nsnap=158 #2.5Gyr
#Nsnap=170 #t=2.8Gyr
#Nsnap=180 #t=3.0Gyr 
#Nsnap=190 #t=3.3Gyr
#Nsnap=200 #t=3.5Gyr
#Nsnap=225 #t=4.0Gyr
#Nsnap=230
#Nsnap=235
#Nsnap=250
#Nsnap=290 #t=5.9Gyr
#Nsnap=300 #6.4Gyr
#Nsnap=350 #t=9.0Gyr
#Nsnap=370 #t=9.8Gyr
#Nsnap=385 #t=10.5Gyr
#Nsnap=390 #t=10.8Gyr
#Nsnap=400 #t=11.3Gyr
#Nsnap=420 #t=12.4 Gyr
#Nsnap=430 #t=13.0Gyr
#Nsnap=440 #t=13.7 Gyr
galcen = 0
#reffp=2.33
#reffp=1.35
#reffp=2.59
reffp=4.0
#reffp=1.48




zr=0
rhocrit=127*(0.27*(1+zr)*(1+zr)*(1+zr)+0.73)
xLam=-0.73/(0.27*(1+zr)*(1+zr)*(1+zr)+0.73)
delta=18.0*np.pi*np.pi+82.0*(xLam)-39*(xLam)*(xLam)
samlogMv=np.arange(8.0,12.0,0.05)
logc=0.537+(1.025-0.537)*np.exp(-0.718*np.power(zr,1.08))+(0.024*zr-0.097)*(samlogMv-np.log10(0.7)-np.log10(1e12))
cont=np.power(10,logc)
Rvirn=np.power(np.power(10,samlogMv)*3.0/4.0/np.pi/rhocrit/delta,1.0/3.0)


for runtodo in dirneed:
	haloinfo=cosmichalo(runtodo)
	beginno=haloinfo['beginno']
	finalno=haloinfo['finalno']
	rundir=haloinfo['rundir']
	subdir=haloinfo['subdir']
	maindir=haloinfo['maindir']
	multifile=haloinfo['multifile']
	halocolor=haloinfo['halocolor']
	halostr=haloinfo['halostr']
	labelname=haloinfo['labelname']
	snapno=[]
	mwithin=[]
	halosA = read_halo_history(rundir, halonostr=halostr, maindir=maindir)
	redlist = halosA['redshift']
	xlist = halosA['x']
	ylist = halosA['y']
	zlist = halosA['z']
	halolist =halosA['ID']
	mvirlist = halosA['M']
	alist = np.array(1./(1.+np.array(redlist)))
	xlist = np.array(xlist)*np.array(alist)/0.702
        ylist = np.array(ylist)*np.array(alist)/0.702
        zlist = np.array(zlist)*np.array(alist)/0.702
#	print 'redlist', redlist
#	print 'xlist', xlist
#	print 'halolist', halolist


	if (int(Nsnap) < 10):
		Nsnapstring = '00'+str(Nsnap)
	elif (int(Nsnap) < 100):
		Nsnapstring = '0'+str(Nsnap)
	else:
		Nsnapstring = str(Nsnap)

	the_snapdir = '/home/tkc004/'+maindir+'/'+rundir+'/'+subdir+'/'
	the_prefix ='snapshot'
	if (multifile == 'y'):
		the_snapdir = '/home/tkc004/'+maindir+'/'+rundir+'/'+subdir+'/snapdir_'+Nsnapstring+'/'
	the_suffix = '.hdf5'
        if galcen == 1:
		finname = '/home/tkc004/'+maindir+'/'+rundir+'/center/galcen.txt'
		f = open(finname)
		f.readline()
		dars = f.readlines()
		snapD = []
		galXl = []
		galYl = []
		galZl = []
                for line in dars:
                        xsd = line.split()
                        snapD.append(int(float(xsd[0])))
                        galXl.append(float(xsd[2]))
                        galYl.append(float(xsd[3]))
                        galZl.append(float(xsd[4]))

	print the_snapdir, Nsnapstring
	header = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, header_only=1, h0=1,cosmological=1)
	G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, cosmological=1, h0=1)
	DM = readsnap(the_snapdir, Nsnapstring, 1, snapshot_name=the_prefix, extension=the_suffix, cosmological=1, h0=1)
	S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, cosmological=1, h0=1)
	ascale = header['time']
	thisz = 1./ascale-1.
	#print 'thisz', thisz
	hubble = header['hubble']
	#print 'hubble', hubble
	if galcen == 1:
                xcen=np.interp(Nsnap,snapD,galXl) #physical distance (kpc) (new program for FIRE 2.0 uses physical distance)
                ycen=np.interp(Nsnap,snapD,galYl)
                zcen=np.interp(Nsnap,snapD,galZl)
	else:
		xcen=np.interp(np.log(ascale),np.log(alist),xlist) #physical distance (kpc)
		ycen=np.interp(np.log(ascale),np.log(alist),ylist)
		zcen=np.interp(np.log(ascale),np.log(alist),zlist)

	halono=np.interp(np.log(ascale),np.log(alist),halolist)
	mvir = np.interp(np.log(ascale),np.log(alist),mvirlist)
	print 'cen', xcen, ycen, zcen
	Gpos = G['p'][:,:]
	Gmass = G['m'][:]
	Gx = Gpos[:,0]
	Gy = Gpos[:,1]
	Gz = Gpos[:,2]
	Gr = np.sqrt(np.square(Gx-xcen)+np.square(Gy-ycen)+np.square(Gz-zcen))

	Spos = S['p'][:,:]
	Smass = S['m'][:]
	Sx = Spos[:,0]
	Sy = Spos[:,1]
	Sz = Spos[:,2]
	Sr = np.sqrt(np.square(Sx-xcen)+np.square(Sy-ycen)+np.square(Sz-zcen))

	DMpos = DM['p'][:,:]
	DMmass = DM['m'][:]
	DMx = DMpos[:,0]
	DMy = DMpos[:,1]
	DMz = DMpos[:,2]
	DMr = np.sqrt(np.square(DMx-xcen)+np.square(DMy-ycen)+np.square(DMz-zcen))
	t_now=quick_lookback_time(thisz)

	Seff = np.sum(Smass[Sr<reffp])*1e10
	DMeff = np.sum(DMmass[DMr<reffp])*1e10

        print 't = '+str(t_now)+' Gyr'	
	print 'Mhalf, log', Seff+DMeff, np.log10(Seff+DMeff)

	logrlist = linspace(np.log10(minr),np.log10(maxr),num=50)
	Gmlist=[]
	Smlist=[]
	DMmlist=[]
	for logr in logrlist:
		withgr =  Gr<np.power(10,logr)
                withsr =  Sr<np.power(10,logr)
                withdr =  DMr<np.power(10,logr)
		Gwithinm = np.sum(Gmass[withgr]*1e10)
		Swithinm = np.sum(Smass[withsr]*1e10)
		DMwithinm = np.sum(DMmass[withdr]*1e10)
		Gmlist.append(Gwithinm)
		Smlist.append(Swithinm)
		DMmlist.append(DMwithinm)
	rlist=np.power(10,logrlist)
	vollist=np.power(rlist,3)*4.0*np.pi/3.0
	Gmlist=np.array(Gmlist)
	Smlist=np.array(Smlist)
	DMmlist=np.array(DMmlist)


#       finding the NFW profile (at z=0) that fits rho at reff and plot it
        #x=reffp/Rvirn
        #rho = rhocrit*delta/cont/cont/cont/x/(1/cont+x)/(1/cont+x)
	Rs = Rvirn/cont
	cdelta = cont*cont*cont*delta/3./(np.log(1+cont)-cont/(1+cont))
	menc = 4*np.pi*rhocrit*cdelta*Rs*Rs*Rs*(np.log((Rs+reffp)/Rs)-reffp/(Rs+reffp))
        Grho = (Gmlist[1:]-Gmlist[:-1])/(vollist[1:]-vollist[:-1])
        Srho = (Smlist[1:]-Smlist[:-1])/(vollist[1:]-vollist[:-1])
        DMrho = (DMmlist[1:]-DMmlist[:-1])/(vollist[1:]-vollist[:-1])
	Trho = np.array(Grho+Srho+DMrho)
        meff = np.interp(reffp, rlist, DMmlist+Smlist+Gmlist)
        logMveff = np.interp(meff, menc, samlogMv)
#        print 'meff', meff
#	print 'menc', menc
#	print 'logMveff', logMveff
#	print 'samlogMv', samlogMv
        logceff=0.537+(1.025-0.537)*np.exp(-0.718*np.power(zr,1.08))+(0.024*zr-0.097)*(logMveff-np.log10(0.7)-np.log10(1e12))
        conteff=np.power(10,logceff)
        Rvirneff=np.power(np.power(10,logMveff)*3.0/4.0/np.pi/rhocrit/delta,1.0/3.0)
        #rhonfw = rhocrit*delta/conteff/conteff/conteff/xeff/(1/conteff+xeff)/(1/conteff+xeff)
	Rseff = Rvirneff/conteff
        cdeltaeff = conteff*conteff*conteff*delta/3./(np.log(1+conteff)-conteff/(1+conteff))
        mnfw = np.array(4*np.pi*rhocrit*cdeltaeff*Rseff*Rseff*Rseff*(np.log((Rseff+rlist)/Rseff)-rlist/(Rseff+rlist)))
	rcrit = np.interp(200.*rhocrit, Trho[::-1], rlist[1:][::-1])
#	print '200.*rhocrit', 200.*rhocrit
#	print 'Grho+Srho+DMrho', Grho+Srho+DMrho
#	print 'rcrit', rcrit
#	print 'reffp', reffp
#	print 'mnfweff',   np.interp(reffp, rlist, mnfw)
#	plt.axvline(x=rcrit, ls='--',color='b')
#	plt.axvline(x=Rvirneff, ls='--', color='r')
#        plt.axvline(x=reffp, ls='--', color='k')
#	plt.plot(rlist,Gmlist,label='Gas',lw=2)
#	plt.plot(rlist,Smlist,label='Star',lw=2)
#        plt.plot(rlist,DMmlist,label='DM',lw=2)
#	plt.xlim([rlist[0],rlist[-1]])
	atot = (Gmlist+Smlist+DMmlist)/rlist/rlist*NewG*Msun_in_kg/kpc_in_m/kpc_in_m
	abar = (Gmlist+Smlist)/rlist/rlist*NewG*Msun_in_kg/kpc_in_m/kpc_in_m
	plt.plot(np.log10(abar),np.log10(atot),label=labelname,lw=2)
	plt.xlabel(r'${\rm log}[a({\rm bar})]\; [{\rm m/s^2}]$')
	plt.ylabel(r'${\rm log}[a({\rm tot})]\; [{\rm m/s^2}]$')
	plt.legend(loc='best',fontsize=14)
#	plt.title('t = '+str(t_now)+' Gyr')
#	print 't = '+str(t_now)+' Gyr'
x0 = np.arange(-14.,-8.,0.1)
plt.plot(x0,x0,linestyle=':',color='black',linewidth=2.)
gcross = 1.2e-10
obs = x0 - np.log10(1.-np.exp(-np.sqrt((10.**x0)/gcross)))
dobs = 0.13 + (-10.-x0) * 0.03
plt.plot(x0,obs,linestyle='--',color='black',linewidth=2.)
plt.plot(x0,obs+2.*dobs,linestyle='-.',color='black',linewidth=2.)
plt.plot(x0,obs-2.*dobs,linestyle='-.',color='black',linewidth=2.)
plt.savefig('figures/mda_'+str(Nsnap)+'.pdf')
plt.clf()

