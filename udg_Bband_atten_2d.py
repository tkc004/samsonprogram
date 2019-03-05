import os
import pyfits
import numpy as np
import matplotlib
matplotlib.use('agg')
from pylab import *
import matplotlib.colors
import matplotlib.cm
import matplotlib.pyplot as plt
import math
import h5py
from samson_functions import *
from lum_mag_conversions import *
beginno=615
finalno=678
haloname=str(sys.argv[1])
nametail=str(sys.argv[2])
dir=str(sys.argv[3])
#haloname='f61'
#nametail='aBxr'
halonamel=haloname+nametail+dir
icount = 0
surmagcut = 26
bandneeded = 2 # B band


spmname=halonamel+'r_galatten_Bband.txt'
spmfile=open(spmname,"w")
for i in range(beginno,finalno,1):
	haloinfo=cosmichalo(haloname)
	xmax_of_box=haloinfo['xmax_of_box']
	movdir=haloinfo['subdir']
	labelname=haloinfo['labelname']
	halocolor=haloinfo['halocolor']
	if dir=='z':
		filemeat = movdir+'_s0'+str(i)+'_t000_star'
	else:
		filemeat = movdir+'_s0'+str(i)+'_t090_star'
	
        filename11='/home/tkc004/scratch/test_movie_maker/'+halonamel+'r/'+movdir+'/'+filemeat+'.dat'
	print 'filename', filename11
	infi11=h5py.File(filename11,'r')
	image24_11 = np.array(infi11["image24"])
	massmap_11 = np.array(infi11["massmap"])
	infi11.close()
	oimage24_11=np.power(10,image24_11)
	print 'image24_11.shape', image24_11.shape
	xno, yno = image24_11.shape
	print 'total luminosity (Lsun)', np.log10(np.sum(oimage24_11))
	absmag = 4.77-2.5*np.log10(np.sum(oimage24_11))
	apmag = absmag-5.0*(1.0-np.log10(145e6))
	print 'absolute magnitude', absmag
	xmax = np.mod(np.argmax(oimage24_11),xno)
	ymax = np.argmax(oimage24_11)/xno
	nobin=200
	xarr = np.linspace(-xmax_of_box/2.,xmax_of_box/2.,num=xno)
	haloX=0.0
	haloY=0.0
#	print 'xmax, ymax', xmax, ymax
	l_bol_sun = 3.9e33
	oimg_cgs = oimage24_11*l_bol_sun
	H = oimg_cgs[xmax-100:xmax+100,ymax-100:ymax+100]
	edges = np.array([xarr[xmax-100:xmax+100],xarr[ymax-100:ymax+100]],np.dtype(float))
        den=np.sum(H)/(pi*4./3.)*8e-4
#	print 'edges', edges
#	print 'H.shape', H.shape
	#print 'len(H[H>den])', len(H[H>den])
	SXi,SYi,rx,ry,err,massrat,Sposx=reff_ell2d(H,edges,haloX,haloY,nobin,den)
	print 'finish! center is ', SXi, SYi
	Sr=[]
	Sl=[]
	for i in range(xno):
		for j in range(xno):
			if (oimg_cgs[i,j]>1e-8):
				Sr=np.append(Sr, np.sqrt(np.square(xarr[i]-SXi)+np.square(xarr[j]-SYi)))
				Sl=np.append(Sl, oimg_cgs[i,j])
	Srclosest = np.array(Sr); Slight = np.array(Sl)
	Slightneed = np.sum(Slight)/2.
	ri=2.0; si=0.1; ncountmax=100; relerr=0.1
        rlhalf,slhalf = find_enclosed_radius(Srclosest,Slight,Slightneed,ri,si,ncountmax,relerr)
	print 'r12', rlhalf, slhalf
        murl = np.linspace(1.0,10.0,num=50)
	#print 'Srclosest', Srclosest
        mul,magl = muprofile(Slight, Srclosest, murl, bandneeded, UNITS_SOLAR_BAND=0, UNITS_CGS=1)
	print 'mul', mul
	print 'magl', magl
        magat26 = np.interp(surmagcut,mul,magl[:-1])
        rat26 = np.interp(surmagcut, mul, murl[:-1])
	muatreff = np.interp(rlhalf,murl[:-1],mul)
        Slightneed26 = np.sum(Slight[Srclosest<rat26])*0.5
        r26lhalf,s26lhalf = find_enclosed_radius(Srclosest,Slight,Slightneed26,ri,si,ncountmax,relerr)
	print 'r26lhalf', r26lhalf, s26lhalf
	print 'muatreff', muatreff
        vmag = luminosity_to_magnitude( Slightneed*2., \
        UNITS_CGS=1, \
        BAND_NUMBER=bandneeded, \
        VEGA=0, AB=1 , \
        MAGNITUDE_TO_LUMINOSITY=0 )
        print 'vmag', vmag
        reff = rlhalf
        effsur = np.log10(reff*reff*2.0*np.pi*1.425*1.425)*2.5+vmag+34.97
        spmfile.write(str(vmag)+'       '+str(rlhalf)+'           '+str(effsur)+'       '+str(muatreff) \
+'               '+str(magat26)+'        '+str(rat26) +'            '+str(r26lhalf)+'\n')

spmfile.close()
	
