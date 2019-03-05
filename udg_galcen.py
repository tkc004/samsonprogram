from struct import *
import sys
import os
import gc
import matplotlib.pyplot as plt
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
import errno
import numpy.linalg as la
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
mpl.use('Agg')
from colors_table import *
from samson_functions import *
from lum_mag_conversions import luminosity_to_magnitude
from cosmo import *
from readsnap_origin import *
pi  = np.pi
sin = np.sin
cos = np.cos

#runtodo='m11' #'m12v' 'm12qq' 'B1' 'm09' 'm10'
runtodo=str(sys.argv[1])


#Rfrac=0.2
Rfrac=0.3 #TK test

haloinfo=cosmichalo(runtodo)
beginno=haloinfo['beginno']
finalno=haloinfo['finalno']
rundir=haloinfo['rundir']
subdir=haloinfo['subdir']
multifile=haloinfo['multifile']
halocolor=haloinfo['halocolor']
halostr=haloinfo['halostr']
maindir = haloinfo['maindir']
snapsep = haloinfo['snapsep']


path = '/home/tkc004/'+maindir+'/'+rundir+'/center'
try:
  os.mkdir(path)
except OSError as exc:
    if exc.errno == errno.EEXIST and os.path.isdir(path):
	pass

spmname='/home/tkc004/'+maindir+'/'+rundir+'/center/galcen.txt'
spmfile=open(spmname,"w")

dirname='/home/tkc004/'+maindir+'/'+rundir
hubble=0.702
halosA = read_halo_history(rundir, maindir=maindir,halonostr=halostr, hubble=hubble, comoving=0)
halored = halosA['redshift']
#print 'Mstarl', np.log10(Mstarl)
ascale = 1./(halored+1.)
for Nsnap in range(beginno,finalno,snapsep):

	if (int(Nsnap) < 10):
		Nsnapstring = '00'+str(Nsnap)
	elif (int(Nsnap) < 100):
		Nsnapstring = '0'+str(Nsnap)
	else:
		Nsnapstring = str(Nsnap)

	snapdir='/home/tkc004/'+maindir+'/'+rundir+'/'+subdir
        if (multifile == 'y'):
                snapdir = '/home/tkc004/'+maindir+'/'+rundir+'/'+subdir+'/snapdir_'+Nsnapstring+'/'
	print snapdir
        ppp_head=readsnap(snapdir,Nsnapstring,0,h0=1,cosmological=1,header_only=1)
        time=ppp_head['time']	
	hubble=ppp_head['hubble'] 
	if Nsnap==beginno:
		halosA = read_halo_history(rundir, maindir=maindir,halonostr=halostr, hubble=hubble, comoving=0)
		halored = halosA['redshift']
		haloXl = halosA['x']
		haloYl = halosA['y']
		haloZl = halosA['z']
		haloid =halosA['ID']
		mvirlist = halosA['M']
		Mstarl = halosA['Ms']
		Rvirl = halosA['R']
        #print 'hubble',ppp_head['hubble'] 
	S = readsnap(snapdir, Nsnapstring, 4, h0=1,cosmological=1)
	posS = S['p']
	xS = posS[:,0]
	yS = posS[:,1]
	zS = posS[:,2]
	Smass = S['m']
        Sage=get_stellar_ages(S,ppp_head,cosmological=1); #should be in gyr  
#	print 'Sage', Sage
	#print 'log(Smass)', np.log10(np.sum(Smass))+10
        Smet=S['z'];
        if(checklen(Smet.shape)>1): Smet=Smet[:,0]
        ztime=1./time-1.; t_now=quick_lookback_time(ztime); # should be in gyr  
#	print 'ztime', ztime
	haloX = np.interp(np.log10(time),np.log10(ascale),haloXl)
        haloY = np.interp(np.log10(time),np.log10(ascale),haloYl)
        haloZ = np.interp(np.log10(time),np.log10(ascale),haloZl)
        Rvir = np.interp(np.log10(time),np.log10(ascale),Rvirl)
	Mstar = np.interp(np.log10(time),np.log10(ascale),Mstarl)
	subhaloN = np.interp(np.log10(time),np.log10(ascale),haloid)
	#print 'haloXYZ', haloX, haloY, haloZ
	xSrel=xS-haloX
	ySrel=yS-haloY
	zSrel=zS-haloZ
	print 'xSrel', xSrel
	ScloseHx = abs(xSrel) < Rfrac*Rvir
	ScloseHy = abs(ySrel) < Rfrac*Rvir
	ScloseHz = abs(zSrel) < Rfrac*Rvir
	ScloseHxy = ScloseHx * ScloseHy * ScloseHz
	Sdists = calcdist2(haloX, haloY, haloZ, xS[ScloseHxy], yS[ScloseHxy], zS[ScloseHxy])
	Sclosest = Sdists < 0.25*Rvir
	SrelposX = xSrel[:][ScloseHxy][Sclosest]
	SrelposY = ySrel[:][ScloseHxy][Sclosest]
	SrelposZ = zSrel[:][ScloseHxy][Sclosest]
	Sraddists=Sdists[Sclosest]
	Smassclosest = Smass[ScloseHxy][Sclosest]*1e10
	Smetclosest = Smet[ScloseHxy][Sclosest]
	SageGyrcloset = Sage[ScloseHxy][Sclosest]
	Sco=np.array([SrelposX,SrelposY,SrelposZ])
	Scor=Sco.T
	nobin=30
	H, edges = np.histogramdd(Scor, bins = (nobin, nobin, nobin), weights=Smassclosest)
	#den=np.sum(H)/(pi*4./3.)*8e-6
	den=np.sum(H)/(pi*4./3.)*8e-5 #TK test
	# do an ellipsoidal fit to locate the galaxy center
	SXi, SYi, SZi, rx, ry, rz, err, massrat, Sposx = reff_ell(H, edges, haloX, haloY, haloZ, nobin, den)
        spmfile.write(str(Nsnap)+'     '+str(subhaloN)+'    '+str(SXi)+'   '+str(SYi)+ '    '+str(SZi)+'        '+str(err)+'      '+str(time)+'\n')


spmfile.close()


