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
from readsnap_samson import *
from Sasha_functions import *
from gadget_lib.cosmo import *
from samson_functions import *
rcParams['text.usetex'] = True
import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='inferno', cmap=cmaps.inferno)
plt.register_cmap(name='plasma', cmap=cmaps.plasma)
#dirneed=['f573z0']
#dirneed=['f573']
#dirneed=['m12icr_b_70']
#dirneed=['m11bmhdcv']
#dirneed=['f573crdc27']
#dirneed=['dm11bhvf8rvc']
#dirneed=['fm10qCRK100']
#dirneed=['fm11qstd']
#dirneed=['fm11qCRK100']
#dirneed=['fm11qCRK1000']
#dirneed=['m11fmhdcv']
#dirneed=['m12icr_b_70', 'm12icr_700']
#dirneed=['m12icr_700']
#dirneed=['fm12b','f476','f553']
#dirneed=['fm11imd']
#dirneed=['fm11q','fm11v']
#dirneed=['f1146']
#dirneed=['f476','f553','f573']
#dirneed=['m11hcr_b_70','m11hcr_700','m11hmhdcv']
#dirneed=['m11vcr_b_70','m11v1cr_b_70','m11v2cr_b_70']
#dirneed=['fm11v','fm11v1','fm11v2']
#dirneed=['m11vmhdcv','m11v1mhdcv','m11v2mhdcv']
#dirneed=['f46']
#dirneed=['f553']
#dirneed=['f553','f476','f383']
#dirneed=['fm12i','m12icr_700']
#dirneed=['m11fcr_700','m11gcr_700','m11hcr_700']
#dirneed=['m11bcr_700']
#dirneed=['m11bcr_b_70','m11bcr_700','m11dcr_b_70','m11dcr_700','m12icr_b_70', 'm12icr_700']
#dirneed=['m11dcr_b_70','m11dcr_700']
#dirneed=['m11dmhdcv']
#dirneed=['m11fcr_b_70']
#dirneed=['f383']
#dirneed=['m09cr_b_70']
#dirneed=['m10qcr_b_70']
#dirneed=['m10qmhdcv']
#dirneed=['m11bcr_b_70']
dirneed=['m11bcr_700']
#dirneed=['m11dcr_700']
#dirneed=['m11bmhdcv']
#dirneed=['f476']
#dirneed=['f61']
#dirneed=['fm11']
#dirneed=['f573']
#dirneed=['f573','fm11']
#dirneed=['1146']
#dirneed=['m11']
#dirneed=['573']
#dirneed=['573','553','476']
#dirneed=['553','476']#
maxr=100
#maxr=20.0
minr=1.0
#Nsnap=239 #t=5.0Gyr in FIRE2
#Nsnap=473 #t=11.0Gyr in FIRE2
#Nsnap=140
#Nsnap=239
#Nsnap=393
#Nsnap=550

Nsnap=600
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

usesingle=1

zr=0
rhocrit=127*(0.27*(1+zr)*(1+zr)*(1+zr)+0.73)
xLam=-0.73/(0.27*(1+zr)*(1+zr)*(1+zr)+0.73)
delta=18.0*np.pi*np.pi+82.0*(xLam)-39*(xLam)*(xLam)
samlogMv=np.arange(7.0,14.0,0.05)
logc=0.537+(1.025-0.537)*np.exp(-0.718*np.power(zr,1.08))+(0.024*zr-0.097)*(samlogMv-np.log10(0.7)-np.log10(1e12))
cont=np.power(10,logc)
Rvirn=np.power(np.power(10,samlogMv)*3.0/4.0/np.pi/rhocrit/delta,1.0/3.0)


for runtodo in dirneed:
	haloinfo=cosmichalo(runtodo)
	beginno=haloinfo['beginno']
	finalno=haloinfo['finalno']
	usepep= haloinfo['usepep']
	snapsep=haloinfo['snapsep']
	firever=haloinfo['firever']
	rundir=haloinfo['rundir']
	subdir=haloinfo['subdir']
	maindir=haloinfo['maindir']
	multifile=haloinfo['multifile']
	halocolor=haloinfo['halocolor']
	halostr=haloinfo['halostr']
	snumadd=haloinfo['snumadd']
	suffixadd=haloinfo['suffixadd']
	snapno=[]
	mwithin=[]

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
	print 'ascale', ascale
	print 'thisz', thisz
	hubble = header['hubble']
	print 'hubble', hubble
	print 'usepep', usepep
	print 'usesingle', usesingle
	print 'hubble', hubble
	if usepep==1:
		if usesingle==1:
                        halosA = read_halo_history_pep(rundir, Nsnap, singlesnap=1, firever=firever,\
			halonostr=halostr, hubble=hubble, comoving=0, maindir=maindir)
		else:
			halosA = read_halo_history_pep(rundir, Nsnap, snapsep=snapsep, singlesnap=0,\
			beginno=beginno, firever=firever,halonostr=halostr, hubble=hubble, comoving=0, maindir=maindir)
	else:
		halosA = read_halo_history(rundir, snumadd=snumadd,suffixadd=suffixadd, halonostr=halostr, maindir=maindir, hubble=hubble, comoving=0)
        redlist = halosA['redshift']
        xlist = halosA['x']
        ylist = halosA['y']
        zlist = halosA['z']
        halolist =halosA['ID']
        mvirlist = halosA['M']
	print 'xlist', xlist
        print 'mvirlist', mvirlist
        alist = np.array(1./(1.+np.array(redlist)))
	print 'alist', alist
        xlist = np.array(xlist)
        ylist = np.array(ylist)
        zlist = np.array(zlist)


	if usesingle==1 and usepep==1:
		xcen = xlist
		ycen = ylist
		zcen = zlist
		halono=halolist
		mvir = mvirlist
	else:
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
	Sm12 = np.amax(Smlist[-1])/2.
	rs12 = np.interp(Sm12, Smlist, rlist)
	print 'half stellar mass radius', rs12
	plt.axvline(x=rcrit, ls='--',color='b')
	plt.axvline(x=Rvirneff, ls='--', color='r')
        plt.axvline(x=reffp, ls='--', color='k')
	plt.plot(rlist,Gmlist,label='Gas',lw=4,ls='-.')
	plt.plot(rlist,Smlist,label='Star',lw=4)
        plt.plot(rlist,DMmlist,label='DM',lw=4,ls=':')
	plt.plot(rlist,mnfw, label='NFW',color='k',ls='--')
	plt.plot(rlist,Gmlist+Smlist+DMmlist,label='Tot',color='k',lw=2)
	plt.xlim([rlist[0],rlist[-1]])
	#plt.ylim([1e5,1e11])
	#plt.ylim([1e8,1e12])
	plt.ylim([1e5,1e12])
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel(r'$r\; [{\rm kpc}]$')
	plt.ylabel(r'${\rm M}_{\rm enc}\; [{\rm M}_{\odot}]$')
	plt.legend(loc='best',fontsize=14)
#	plt.title('t = '+str(t_now)+' Gyr')
#	print 't = '+str(t_now)+' Gyr'
	plt.savefig('figures/Mencgsdm_'+runtodo+'_'+str(Nsnap)+'.pdf')
	plt.clf()


	plt.axvline(x=reffp, ls='--', color='k')
        plt.plot(rlist[1:],Grho,label='Gas')
        plt.plot(rlist[1:],Srho,label='Star')
        plt.plot(rlist[1:],DMrho,label='DM')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'$r\;( {\rm kpc})$')
        plt.ylabel(r'$\rho ({\rm M}_{\odot}/{\rm kpc}^3)$')
        plt.legend(loc='best',fontsize=12)
        print 't = '+str(t_now)+' Gyr'
        plt.savefig('figures/rhogsdm_'+runtodo+'_'+str(Nsnap)+'.pdf')
        plt.clf()
	
