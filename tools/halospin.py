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

def gethalospin(runtodo, time, rfactor=1, galcen=0,useAHF=0,hubble=0.702):
		print 'runtodo', runtodo
                haloinfo=cosmichalo(runtodo)
                beginno=haloinfo['beginno']
                finalno=haloinfo['finalno']
                rundir=haloinfo['rundir']
                subdir=haloinfo['subdir']
                multifile=haloinfo['multifile']
                halocolor=haloinfo['halocolor']
                halostr=haloinfo['halostr']
                maindir=haloinfo['maindir']
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
                halosA = read_halo_history(rundir, halonostr=halostr, maindir=maindir)
                redlist = halosA['redshift']
                xlist = halosA['x']
                ylist = halosA['y']
                zlist = halosA['z']
		xvlist = halosA['xv']
                yvlist = halosA['yv']
                zvlist = halosA['zv']
                halolist =halosA['ID']
                mvirlist = halosA['M']
                Rvirlist = halosA['R']
		lambdalist = halosA['lambdalist']
                alist = np.array(1./(1.+np.array(redlist)))
                xlist = np.array(xlist)*np.array(alist)/0.702
                ylist = np.array(ylist)*np.array(alist)/0.702
                zlist = np.array(zlist)*np.array(alist)/0.702
		xvlist = np.array(xvlist)*np.array(alist)
                yvlist = np.array(yvlist)*np.array(alist)
                zvlist = np.array(zvlist)*np.array(alist)
                Rvirlist = np.array(Rvirlist)*np.array(alist)/0.702
		#Rvirlist = np.array(Rvirlist)/0.702
		mvirlist = np.array(mvirlist)/0.702
		lambdalist = np.array(lambdalist)	
                #print 'redlist', redlist
                #print 'xlist', xlist
                #print 'halolist', halolist
                readtimelist = readtime()
                snaplist = readtimelist['snaplist']
                timelist = readtimelist['timelist']
		atlist = readtimelist['alist']
                Nsnap = int(np.interp(time, timelist, snaplist))

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
		if useAHF==0:
			header = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, header_only=1, h0=1,cosmological=1)
			ascale = header['time']
			thisz = 1./ascale-1.
			#print 'thisz', thisz
			hubble = header['hubble']
			print 'hubble', hubble
		else:
			ascale = np.interp(time, timelist, atlist)
			thisz = 1./ascale-1.
		
                if galcen == 1:
                        xcen=np.interp(Nsnap,snapD,galXl) #physical distance (kpc)
                        ycen=np.interp(Nsnap,snapD,galYl)
                        zcen=np.interp(Nsnap,snapD,galZl)
                else:
                        xcen=np.interp(np.log(ascale),np.log(alist),xlist) #physical distance (kpc)
                        ycen=np.interp(np.log(ascale),np.log(alist),ylist)
                        zcen=np.interp(np.log(ascale),np.log(alist),zlist)
		xvcen=np.interp(np.log(ascale),np.log(alist),xvlist)
                yvcen=np.interp(np.log(ascale),np.log(alist),yvlist)
                zvcen=np.interp(np.log(ascale),np.log(alist),zvlist)
                Rvir=np.interp(np.log(ascale),np.log(alist),Rvirlist)
                halono=np.interp(np.log(ascale),np.log(alist),halolist)
                mvir = np.interp(np.log(ascale),np.log(alist),mvirlist)
		spinlambda = np.interp(np.log(ascale),np.log(alist),lambdalist)
		#print 'lambdalist', lambdalist
		#print 'Rvirlist', Rvirlist
		if useAHF==0:
			DM = readsnap(the_snapdir, Nsnapstring, 1, snapshot_name=the_prefix, extension=the_suffix, cosmological=1, h0=1)
			DMpos = np.array(DM['p'][:,:]) #in kpc
			DMmass = np.array(DM['m'][:]*1e10) #in Msun
			DMvel = np.array(DM['v'][:,:]) # in km/s
			DMrelposraw = DMpos - np.array([xcen,ycen,zcen])
			cut = np.sqrt(DMrelposraw[:,0]*DMrelposraw[:,0]+DMrelposraw[:,1]*DMrelposraw[:,1]+DMrelposraw[:,2]*DMrelposraw[:,2])<Rvir*rfactor
			DMrelpos = DMrelposraw[cut]
			DMrelvel = DMvel[cut] - np.array([xvcen,yvcen,zvcen])
			DMmasscut = DMmass[cut]
			DMmassinR = np.sum(DMmasscut)
			DMmom = DMmasscut[:,np.newaxis]*DMrelvel # in Msun km/s
			#print 'DMpos', DMpos
			#print 'DMrelpos', DMrelpos
			#print 'DMrelvel', DMrelvel
			#print 'DMmass', DMmass
			print 'DMmassvir', DMmassinR
			print 'mvir', mvir 
			#print 'DMmassnew', DMmass[:,np.newaxis]
			#print 'DMmom', DMmom
			print 'DM cm vel', np.sum(DMmom, axis=0)/DMmassinR
			print 'vcen', xvcen, yvcen, zvcen
			DMangmom = np.cross(DMrelpos,DMmom) # in Msun kpc km/s
			sumDMangmom = np.linalg.norm(np.sum(DMangmom, axis=0))
			vcir = np.sqrt(NewtonG_cgs*DMmassinR*solar_mass_in_g/Rvir/rfactor/kpc_in_cm) # in cm/s
			spinlambda = sumDMangmom*solar_mass_in_g*kpc_in_cm*km_in_cm/np.sqrt(2)/DMmassinR/solar_mass_in_g/vcir/Rvir/rfactor/kpc_in_cm
		return spinlambda


                     
