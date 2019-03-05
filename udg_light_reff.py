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
runtodo = 'm11dcr_700'
#runtodo = 'm11bcr_700'
#runtodo='m10vcr_700'
#runtodo='m10qcr_b_70'
#runtodo='m11qcr_b_70'
#runtodo='m11hcr_b_70'
#runtodo='m11fcr_700'
#runtodo='m11dcr_700'
#runtodo='m11bcr_700'
#runtodo='m11hcr_700'
#runtodo='m11bcr_b_70'
#runtodo='m11hmhdcv'
#runtodo='m11dmhdcv'
#runtodo='m11dcr_b_70'
#runtodo='m11fcr_b_70'
#runtodo='m11fmhdcv'
#runtodo='m11bmhdcv'
#runtodo='m11bcr_b_70'
#runtodo='m11' #'m12v' 'm12qq' 'B1' 'm09' 'm10'
#runtodo=str(sys.argv[1])
#bandneeded = int(sys.argv[2])
bandneeded = 3 #Johnson V band
#bandneeded = 2 #Johnson B band
surmagcut = 26
#bandneeded = 10 #SDSS g band
#bandneeded = 11 #SDSS r band
snapsep=10
haloinfo=cosmichalo(runtodo)
beginno=haloinfo['beginno']
finalno=haloinfo['finalno']
rundir=haloinfo['rundir']
subdir=haloinfo['subdir']
multifile=haloinfo['multifile']
halocolor=haloinfo['halocolor']
halostr=haloinfo['halostr']
maindir=haloinfo['maindir']
if bandneeded == 2:
	spmname='data/'+runtodo+'_gallightinfo_Bband.txt'
elif bandneeded == 3:
        spmname='data/'+runtodo+'_gallightinfo_Vband.txt'
spmfile=open(spmname,"w")
spmfile.write('Snapshot no   ' + 'Halo number   ' + 'X     '+'Y      '+'Z        '+'rx       '+'ry        '+'rz       '+'errlist   '+'\n')

dirname='/home/tkc004/'+maindir+'/'+rundir
hubble=haloinfo['h0']
halosA = readhalos(dirname,halostr,hubble=hubble)
halored = halosA['z']
Rvirl = halosA['rv']
haloXl = halosA['haloX']
haloYl = halosA['haloY']
haloZl = halosA['haloZ']
haloid = halosA['id']
Mstarl = halosA['ms']
print 'haloid', haloid
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
	if multifile=='y':
		snapdir+='/snapdir_'+Nsnapstring
	print 'snapdir', snapdir
        ppp_head=readsnap(snapdir,Nsnapstring,0,h0=1,cosmological=1,header_only=1)
        time=ppp_head['time']
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
	#print 'haloXYZ', haloX, haloY, haloZ
	xSrel=xS-haloX
	ySrel=yS-haloY
	zSrel=zS-haloZ
	print 'xSrel', xSrel
	ScloseHx = abs(xSrel) < 0.1*Rvir
	ScloseHy = abs(ySrel) < 0.1*Rvir
	ScloseHz = abs(zSrel) < 0.1*Rvir
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
	Slight = Smassclosest*colors_table(SageGyrcloset, Smetclosest/0.019, \
            BAND_ID=bandneeded, KROUPA_IMF=1)
	Sco=np.array([SrelposX,SrelposY,SrelposZ])
	Scor=Sco.T
	nobin=30
	H, edges = np.histogramdd(Scor, bins = (nobin, nobin, nobin), weights=Slight)
	den=np.sum(H)/(pi*4./3.)*8e-6
	# do an ellipsoidal fit to locate the galaxy center
	SXi, SYi, SZi, rx, ry, rz, err, massrat, Sposx = reff_ell(H, edges, haloX, haloY, haloZ, nobin, den)

	#find half light radius
        xSrels=xS-SXi
        ySrels=yS-SYi
        zSrels=zS-SZi
	Sr = np.sqrt(np.square(xSrels)+np.square(ySrels)+np.square(zSrels))
        Sclosest = Sr < 0.25*Rvir
        SrelposXs = xSrels[:][Sclosest]
        SrelposYs = ySrels[:][Sclosest]
        SrelposZs = zSrels[:][Sclosest]
        Smassclosest = Smass[Sclosest]*1e10 
	Srclosest = Sr[Sclosest]
       	Smassneed = np.sum(Smassclosest)*0.5; ri=rx; si=0.1; ncountmax=100; relerr=0.1
        rhalf,shalf = find_enclosed_radius(Srclosest,Smassclosest,Smassneed,ri,si,ncountmax,relerr)
	Smetclosest = Smet[Sclosest]
        SageGyrcloset = Sage[Sclosest]
        Slight = Smassclosest*colors_table(SageGyrcloset, Smetclosest/0.019, \
            BAND_ID=bandneeded, KROUPA_IMF=1)
	Slightneed=np.sum(Slight)*0.5; ri=rx; si=0.1; ncountmax=100; relerr=0.1
	rlhalf,slhalf = find_enclosed_radius(Srclosest,Slight,Slightneed,ri,si,ncountmax,relerr)
	Sxy = np.sqrt(np.square(SrelposXs)+np.square(SrelposYs))
        zhalf,szhalf = find_enclosed_radius(Sxy,Slight,Slightneed,ri,si,ncountmax,relerr)
        Sxz = np.sqrt(np.square(SrelposXs)+np.square(SrelposZs))
        yhalf,syhalf = find_enclosed_radius(Sxz,Slight,Slightneed,ri,si,ncountmax,relerr)
        Syz = np.sqrt(np.square(SrelposZs)+np.square(SrelposYs))
        xhalf,sxhalf = find_enclosed_radius(Syz,Slight,Slightneed,ri,si,ncountmax,relerr)
	print 'Smassneed, Mstar/2', Smassneed, Mstar*0.5
	print 'rhalf,shalf', rhalf, shalf
	print 'rlhalf,slhalf',rlhalf, slhalf
	print 'xhalf,sxhalf',xhalf, sxhalf
	print 'yhalf,syhalf',yhalf, syhalf
	print 'zhalf,szhalf',zhalf, szhalf
	vmag = luminosity_to_magnitude( Slightneed*2., \
        UNITS_CGS=1, \
        BAND_NUMBER=bandneeded, \
        VEGA=0, AB=1 , \
        MAGNITUDE_TO_LUMINOSITY=0 )
	print 'vmag', vmag
	reff = (xhalf+yhalf+zhalf)/3.0
        effsur = np.log10(reff*reff*2.0*np.pi*1.425*1.425)*2.5+vmag+34.97
	print 'effective surface brightness', effsur
	ryz =  np.sqrt(SrelposYs*SrelposYs+SrelposZs*SrelposZs)
	withinreff = ryz < xhalf
	slightwithinreff = np.sum(Slight[withinreff])
        vmagwithinreff = luminosity_to_magnitude( slightwithinreff, \
        UNITS_CGS=1, \
        BAND_NUMBER=bandneeded, \
        VEGA=0, AB=1 , \
        MAGNITUDE_TO_LUMINOSITY=0 )
	meansurx = np.log10(xhalf*xhalf*np.pi*1.425*1.425)*2.5+vmagwithinreff+34.97
	murl = np.linspace(1.0,Rvir*0.2,num=50)
	mul,magl = muprofile(Slight, Srclosest, murl, bandneeded)
        muyzl,magyzl = muprofile(Slight, ryz, murl, bandneeded)
#	print 'murl', murl
#	print 'mul', mul
#	print 'magl', magl
	magat26 = np.interp(surmagcut,mul,magl[:-1])
	rat26 = np.interp(surmagcut, mul, murl[:-1])
	muateff = np.interp(reff, murl[:-1],mul)
        muateffyz = np.interp(xhalf, murl[:-1],muyzl)
	Slightneed26 = np.sum(Slight[Srclosest<rat26])*0.5
	r26lhalf,s26lhalf = find_enclosed_radius(Srclosest,Slight,Slightneed26,ri,si,ncountmax,relerr)
        z26half,s26zhalf = find_enclosed_radius(Sxy,Slight,Slightneed26,ri,si,ncountmax,relerr)
        y26half,s26yhalf = find_enclosed_radius(Sxz,Slight,Slightneed26,ri,si,ncountmax,relerr)
        x26half,s26xhalf = find_enclosed_radius(Syz,Slight,Slightneed26,ri,si,ncountmax,relerr)
	print 'murl', murl
	print 'mul', mul
	print 'muateff', muateff
	print 'magat26', magat26
	print 'rat26', rat26
	print 'r26lhalf', r26lhalf
	print 'muateffyz', muateffyz
	print 'meansurx', meansurx
	spmfile.write(str(Nsnap)+'     '+halostr+'    '+str(SXi)+'   '+str(SYi) \
+ '    '+str(SZi)+'        '+str(rx)+'         '+str(ry)+'           '+str(rz) \
+'               '+str(err)+'           '+str(massrat) +'          '+str(len(Sposx)) \
+'        '+str(rhalf)+'       '+str(rlhalf)+'          '+str(xhalf)+'         ' \
+str(yhalf)+'             '+str(zhalf)+'             '+str(vmag)+'             ' \
+str(effsur)+'        '+str(muateff)+'               '+str(magat26)+'        '+str(rat26) \
+'            '+str(r26lhalf)+'           '+str(x26half)+'          '+str(y26half)+'            '+str(z26half)+'      '+str(meansurx)+'     '+str(muateffyz)+'\n')



spmfile.close()


