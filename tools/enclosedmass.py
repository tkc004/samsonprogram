import os
#import pyfits
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


def enclosedmass(runtodo, time, minr, maxr, usesnapno=0, galcen=0, unitRv=0,needDM=1,needS=1,needG=1,needHI=0):
		Gmlist=[]
		Smlist=[]
		DMmlist=[]
		HImlist=[]
		haloinfo=cosmichalo(runtodo)
		beginno=haloinfo['beginno']
		finalno=haloinfo['finalno']
		rundir=haloinfo['rundir']
		subdir=haloinfo['subdir']
		multifile=haloinfo['multifile']
		halocolor=haloinfo['halocolor']
		halostr=haloinfo['halostr']
		maindir=haloinfo['maindir']
		snumadd=haloinfo['snumadd']
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
                snapno=[]
                mwithin=[]
                halosA = read_halo_history(rundir, halonostr=halostr, maindir=maindir,snumadd=snumadd)
                redlist = halosA['redshift']
                xlist = halosA['x']
                ylist = halosA['y']
                zlist = halosA['z']
                halolist =halosA['ID']
                mvirlist = halosA['M']
		Rvirlist = halosA['R']
                alist = np.array(1./(1.+np.array(redlist)))
                xlist = np.array(xlist)*np.array(alist)/0.702
                ylist = np.array(ylist)*np.array(alist)/0.702
                zlist = np.array(zlist)*np.array(alist)/0.702
		Rvirlist = np.array(Rvirlist)*np.array(alist)/0.702
		mvirlist = np.array(mvirlist)/0.702
                #print 'redlist', redlist
                #print 'xlist', xlist
                #print 'halolist', halolist
		readtimelist = readtime()
                snaplist = readtimelist['snaplist']
		timelist = readtimelist['timelist']
		
		if usesnapno==0:
			Nsnap = int(np.interp(time, timelist, snaplist))
		else:
			Nsnap=time
			time = np.interp(Nsnap, snaplist, timelist)

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
                print 'the_snapdir' , the_snapdir
                header = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, header_only=1, h0=1,cosmological=1)
                ascale = header['time']
                thisz = 1./ascale-1.
                #print 'thisz', thisz
                hubble = header['hubble']
                #print 'hubble', hubble
                if galcen == 1:
                        xcen=np.interp(Nsnap,snapD,galXl) #physical distance (kpc)
                        ycen=np.interp(Nsnap,snapD,galYl)
                        zcen=np.interp(Nsnap,snapD,galZl)
                else:
                        xcen=np.interp(np.log(ascale),np.log(alist),xlist) #physical distance (kpc)
                        ycen=np.interp(np.log(ascale),np.log(alist),ylist)
                        zcen=np.interp(np.log(ascale),np.log(alist),zlist)
		Rvir=np.interp(np.log(ascale),np.log(alist),Rvirlist)
                halono=np.interp(np.log(ascale),np.log(alist),halolist)
                mvir = np.interp(np.log(ascale),np.log(alist),mvirlist)
		if unitRv==1:
			minr*=Rvir
			maxr*=Rvir
                #print 'cen', xcen, ycen, zcen
		if needG==1:
			G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, cosmological=1, h0=1)
			Gpos = G['p'][:,:]
			Gmass = G['m'][:]
			Gx = Gpos[:,0]
			Gy = Gpos[:,1]
			Gz = Gpos[:,2]
			Gr = np.sqrt(np.square(Gx-xcen)+np.square(Gy-ycen)+np.square(Gz-zcen))
			if needHI==1:
				Gmetal = G['z'][:,:]
				Gnh = G['nh'][:]
				H_fraction = 1.0-Gmetal[:,1]-Gmetal[:,0]
				print 'H_fraction', H_fraction
				H_mass = Gmass*H_fraction
				print 'neutral abundance', Gnh
				H_p0_mass = H_mass*Gnh

		if needS==1:
			S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, cosmological=1, h0=1)
			Spos = S['p'][:,:]
			Smass = S['m'][:]
			Sx = Spos[:,0]
			Sy = Spos[:,1]
			Sz = Spos[:,2]
			Sr = np.sqrt(np.square(Sx-xcen)+np.square(Sy-ycen)+np.square(Sz-zcen))

		if needDM==1:
			DM = readsnap(the_snapdir, Nsnapstring, 1, snapshot_name=the_prefix, extension=the_suffix, cosmological=1, h0=1)
			DMpos = DM['p'][:,:]
			DMmass = DM['m'][:]
			DMx = DMpos[:,0]
			DMy = DMpos[:,1]
			DMz = DMpos[:,2]
			DMr = np.sqrt(np.square(DMx-xcen)+np.square(DMy-ycen)+np.square(DMz-zcen))
                t_now=quick_lookback_time(thisz)

                logrlist = linspace(np.log10(minr),np.log10(maxr),num=20)
		if needG==1:
			Gmlist=[]
			for logr in logrlist:
				withgr =  Gr<np.power(10,logr)
				Gwithinm = np.sum(Gmass[withgr]*1e10)
				Gmlist.append(Gwithinm)
		        Gmlist=np.array(Gmlist)
			if needHI==1:
				mHIlist=[]
				for logr in logrlist:
					withgr =  Gr<np.power(10,logr)
					HIwithinm = np.sum(H_p0_mass[withgr]*1e10)
					HImlist.append(HIwithinm)
				HImlist=np.array(HImlist)
		if needS==1:
			Smlist=[]
			for logr in logrlist:
				withsr =  Sr<np.power(10,logr)
				Swithinm = np.sum(Smass[withsr]*1e10)
				Smlist.append(Swithinm)
			Smlist=np.array(Smlist)
		if needDM==1:
			DMmlist=[]
			for logr in logrlist:
				withdr =  DMr<np.power(10,logr)
				DMwithinm = np.sum(DMmass[withdr]*1e10)
				DMmlist.append(DMwithinm)
			DMmlist=np.array(DMmlist)
                rlist=np.power(10,logrlist)
                vollist=np.power(rlist,3)*4.0*np.pi/3.0
		return {'rlist':rlist,'vollist':vollist,'Gmlist':Gmlist, 'HImlist':HImlist, 'Smlist':Smlist,'DMmlist':DMmlist,'halocolor':halocolor}






