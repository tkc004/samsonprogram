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

NewG = 6.674e-11 #SI unit
Msun_in_kg = 2e30
kpc_in_m = 3.085e19

#dirneed=['f383']
#dirneed=['f476']
#dirneed=['f61']
#dirneed=['fm11']
#dirneed=['f573']
#dirneed=['f573','f553','f476','fm11','f383','f61']
#dirneed=['1146']
#dirneed=['m11']
#dirneed=['573']
#dirneed=['573','553','476']
#dirneed=['553','476']
dirneed=['dm11a_hv']
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
#reffp=1.48




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
	snapno=[]
	mwithin=[]
	halosA = read_halo_history(rundir, halonostr=halostr, maindir=maindir)
	redlist = halosA['redshift']
	xlist = halosA['x']
	ylist = halosA['y']
	zlist = halosA['z']
	halolist =halosA['ID']
	mvirlist = halosA['M']
	rvirlist =halosA['R']
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
	header = readsnap(the_snapdir, Nsnapstring, 1, snapshot_name=the_prefix, extension=the_suffix, header_only=1, h0=1,cosmological=1)
	print 'header', header
	DM = readsnap(the_snapdir, Nsnapstring, 1, snapshot_name=the_prefix, extension=the_suffix, cosmological=1, h0=1)
	contam = readsnap(the_snapdir, Nsnapstring, 2, snapshot_name=the_prefix, extension=the_suffix, cosmological=1, h0=1)
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
	Rvir = np.interp(np.log(ascale),np.log(alist),rvirlist)
	print 'cen', xcen, ycen, zcen

	DMpos = DM['p'][:,:]
	DMmass = DM['m'][:]
	DMx = DMpos[:,0]
	DMy = DMpos[:,1]
	DMz = DMpos[:,2]
	DMr = np.sqrt(np.square(DMx-xcen)+np.square(DMy-ycen)+np.square(DMz-zcen))
	cutDMr = DMr<Rvir
	DMcut = np.sum(DMmass[cutDMr])
	
        ctpos = contam['p'][:,:]
        ctmass = contam['m'][:]
        ctx = ctpos[:,0]
        cty = ctpos[:,1]
        ctz = ctpos[:,2]
        ctr = np.sqrt(np.square(ctx-xcen)+np.square(cty-ycen)+np.square(ctz-zcen))
	cutctr = ctr<Rvir
	ctmcut = np.sum(ctmass[cutctr])

	print 'mvir', mvir
	print 'Rvir', Rvir
	print 'DMcut', DMcut
	print 'ctmcut', ctmcut
	
	

	
