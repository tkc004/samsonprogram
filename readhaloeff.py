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

#dirneed=['1146']
#dirneed=['383']
#dirneed=['476']
#dirneed=['573']
#dirneed=['573','553','476']
#dirneed=['553','476']
maxr=20.0
minr=1.0
#Nsnap=130 #t=2.0Gyr
#Nsnap=140 #t=2.2Gyr
#Nsnap=150 #t=2.4Gyr
#Nsnap=158 #2.5Gyr
#Nsnap=170 #t=2.8Gyr
#Nsnap=180 #t=3.0Gyr 
#Nsnap=190 #t=3.3Gyr
#Nsnap=200 #t=3.5Gyr
#Nsnap=225 #t=4.0Gyr
#Nsnap=230 #t=4.1Gyr
#Nsnap=235 #t=4.2Gyr
#Nsnap=250 #t=4.6Gyr
#Nsnap=290 #t=5.9Gyr
#Nsnap=300 #6.4Gyr
#Nsnap=350 #t=9.0Gyr
#Nsnap=385 #t=10.5Gyr
#Nsnap=390 #t=10.8Gyr
#Nsnap=400 #t=11.3Gyr
#Nsnap=420 #t=12.4 Gyr
#Nsnap=430 #t=13.0Gyr
#Nsnap=440 #t=13.7 Gyr
galcen = 1
#reffp=2.33
#reffp=1.35
#reffp=2.59
#reffp=1.48
#datalist=[[,],[,],[,]]
#datalist=[[199,1.41],[350,2.74],[430,3.21]]
#datalist=[[225,2.63],[291,2.55],[430,1.99]]
#datalist=[[180,1.16],[225,2.09],[385,3.13]]
#datalist=[[158,1.48],[225,1.88]]
#datalist=[[140,2.12],[235,1.98]]
#dirneed=['f573']
dirneed=['f573','f553','f476','fm11','f383','f61']
#tlist=[[13.4,8.3,13.5],[10.4,5.8,13.5],[10.7,2.3,13.5],[2.7,2.0,6.0],[2.0,2.0,6.5]]
#tlist=[[13.4,8.3,13.5],[10.4,5.8,13.5],[10.7,2.3,13.5],[2.7,2.0,6.0],[2.0,2.0,6.5],[2.0,2.0,2.4]]
tlist=[[13.4,8.2,13.6],[10.4,5.7,13.6],[10.7,2.6,13.6],[2.7,2.3,6.0],[2.0,2.0,6.5],[2.2,2.2,2.4]]


refflist=[[2.9,1.4,3.2],[1.4,1.3,1.4],[2.0,1.3,2.0],[1.7,1.3,4.0],[1.3,1.3,1.5],[1.6,1.6,2.5]] #x
axislist=[[0.35,0.61,0.32],[0.66,0.49,0.72],[0.63,0.52,0.91],[0.58,0.56,0.4],[0.35,0.35,0.75],[0.71,0.71,0.34]] #x
#refflist=[[1.9,1.4,2.0],[1.3,1.3,1.3],[1.3,1.3,2.2],[2.6,1.4,3.2],[1.3,1.3,1.4],[1.5,1.5,1.5]] #y
#axislist=[[0.64,0.47,0.68],[0.57,0.42,0.64],[0.54,0.54,0.28],[0.44,0.64,0.61],[0.48,0.48,0.52],[0.63,0.63,0.63]] #y
#refflist=[[1.9,1.8,1.7],[1.3,1.3,1.3],[1.4,1.3,1.6],[2.2,1.4,2.3],[1.3,1.3,1.3],[1.3,1.3,1.5]] #z
#axislist=[[0.35,0.3,0.37],[0.64,0.66,0.59],[0.53,0.5,0.54],[0.36,0.16,0.75],[0.58,0.58,0.55],[0.73,0.73,0.64]] #z
zr=0

rhocrit=127*(0.27*(1+zr)*(1+zr)*(1+zr)+0.73)
xLam=-0.73/(0.27*(1+zr)*(1+zr)*(1+zr)+0.73)
delta=18.0*np.pi*np.pi+82.0*(xLam)-39*(xLam)*(xLam)
samlogMv=np.arange(8.0,12.0,0.05)
logc=0.537+(1.025-0.537)*np.exp(-0.718*np.power(zr,1.08))+(0.024*zr-0.097)*(samlogMv-np.log10(0.7)-np.log10(1e12))
cont=np.power(10,logc)
Rvirn=np.power(np.power(10,samlogMv)*3.0/4.0/np.pi/rhocrit/delta,1.0/3.0)



def meff(runtodo,time,reffp,unitRv=0):
	haloinfo=cosmichalo(runtodo)
	rundir=haloinfo['rundir']
	subdir=haloinfo['subdir']
	maindir=haloinfo['maindir']
	multifile=haloinfo['multifile']
	halocolor=haloinfo['halocolor']
	halostr=haloinfo['halostr']
        snapno=[]
        mwithin=[]
        halosA = read_halo_history(rundir, halonostr=halostr,maindir=maindir)
        redlist = halosA['redshift']
        xlist = halosA['x']
        ylist = halosA['y']
        zlist = halosA['z']
        halolist =halosA['ID']
        mvirlist = halosA['M']
	mstarlist = halosA['Ms']
	rvirlist = halosA['R']
        alist = np.array(1./(1.+np.array(redlist)))
        xlist = np.array(xlist)*np.array(alist)/0.702
        ylist = np.array(ylist)*np.array(alist)/0.702
        zlist = np.array(zlist)*np.array(alist)/0.702
	rvirlist = np.array(rvirlist)*np.array(alist)/0.702
        print 'halolist', halolist
	snaplist, timelist = readtime()
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
		xcen=np.interp(Nsnap,snapD,galXl) #physical distance (kpc)
		ycen=np.interp(Nsnap,snapD,galYl)
		zcen=np.interp(Nsnap,snapD,galZl)
	else:
		xcen=np.interp(np.log(ascale),np.log(alist),xlist) #physical distance (kpc)
		ycen=np.interp(np.log(ascale),np.log(alist),ylist)
		zcen=np.interp(np.log(ascale),np.log(alist),zlist)
	rvir = np.interp(np.log(ascale),np.log(alist),rvirlist)

#	print 'alist', alist
	halono=np.interp(np.log(ascale),np.log(alist),halolist)
	#print 'mvirlist', mvirlist
	mvir = np.interp(np.log(ascale),np.log(alist),mvirlist/0.702)
	#print 'mstarlist', mstarlist
        mstar = np.interp(np.log(ascale),np.log(alist),mstarlist/0.702)
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

	if unitRv==1:
		reffp*=rvir
	mlist=[]
	rlist=np.linspace(0.1,20,num=200)
        for i in range(len(rlist)):
                withinrlist = Sr<rlist[i]
                withinrm = np.sum(Smass[withinrlist]*1e10)
                mlist = np.append(mlist,withinrm)
        r12 = np.interp(mlist[-1]*0.5,mlist,rlist)

	Seff = np.sum(Smass[Sr<reffp])*1e10
	DMeff = np.sum(DMmass[DMr<reffp])*1e10
	Geff = np.sum(Gmass[Gr<reffp])*1e10

	print 'Mhalf, log', Seff+DMeff, np.log10(Seff+DMeff)
	return {'Seff':Seff,'DMeff':DMeff,'Geff':Geff,'mvir':mvir,'mstar':mstar, 'r12':r12}

Mvirstr=' '
MsMvstr=' '
Meffstr=' '
M12str=' '
Meffwgstr=' '
MsM12str=' '
Ms02str=' '
Ms02Mvstr=' '
for ncount in range(len(dirneed)):
	runtodo = dirneed[ncount]
	Selist=[]
	DMelist=[]
	Gelist=[]
	Mvirlist=[]
	Mstarlist=[]
	ms02list=[]
	m12dmslist=[]
	for icount in range(len(tlist[ncount])):
		time=tlist[ncount][icount]
		reffp=refflist[ncount][icount]*np.sqrt(axislist[ncount][icount])*4.0/3.0
		eff=meff(runtodo,time,reffp)
		r12=eff['r12']
		m02 = meff(runtodo,time,0.2,unitRv=1)
		m12eff = meff(runtodo,time, r12)
		m12dmslist.append(m12eff['Seff']+m12eff['DMeff'])
		ms02list.append(m02['Seff'])
		Selist.append(eff['Seff'])
		DMelist.append(eff['DMeff'])
		Gelist.append(eff['Geff'])
		Mvirlist.append(eff['mvir'])
		Mstarlist.append(eff['mstar'])
	m12dmslist=np.array(m12dmslist)
	ms02list=np.array(ms02list)
	Selist=np.array(Selist)
	DMelist=np.array(DMelist)
	Gelist=np.array(Gelist)
	Mvirlist=np.array(Mvirlist)
	Mstarlist=np.array(Mstarlist)
	Mvirround = np.around(Mvirlist/1e10, decimals=1)
	Mvirstr = Mvirstr+'&      $'+str(Mvirround[0])+'$& $\,_{'+str(Mvirround[1])+'}^{'+str(Mvirround[2])+'} $'
	MsMvround = np.around(Mstarlist*1e3/Mvirlist, decimals=2)
	MsMvstr = MsMvstr+'&      $'+str(MsMvround[0])+'$& $\,_{'+str(MsMvround[1])+'}^{'+str(MsMvround[2])+'} $'
        Ms02Mvround = np.around(ms02list*1e3/Mvirlist, decimals=2)
        Ms02Mvstr = Ms02Mvstr+'&      $'+str(Ms02Mvround[0])+'$& $\,_{'+str(Ms02Mvround[1])+'}^{'+str(Ms02Mvround[2])+'} $'
        Ms02round = np.around(ms02list/1e8, decimals=2)
        Ms02str = Ms02str+'&      $'+str(Ms02round[0])+'$& $\,_{'+str(Ms02round[1])+'}^{'+str(Ms02round[2])+'} $'
	Meffround = np.around((Selist+DMelist)/1e8, decimals=2)
        Meffstr = Meffstr+'&      $'+str(Meffround[0])+'$& $\,_{'+str(Meffround[1])+'}^{'+str(Meffround[2])+'} $'
	M12round = np.around((m12dmslist)/1e8, decimals=2)
	M12str = M12str+'&      $'+str(M12round[0])+'$& $\,_{'+str(M12round[1])+'}^{'+str(M12round[2])+'} $'
        Meffwground = np.around((Selist+DMelist+Gelist)/1e8, decimals=2)
        Meffwgstr = Meffwgstr+'&      $'+str(Meffwground[0])+'$& $\,_{'+str(Meffwground[1])+'}^{'+str(Meffwground[2])+'} $'
	MsM12round = np.around(10*Selist/(Selist+DMelist), decimals=2)
        MsM12str = MsM12str+'&      $'+str(MsM12round[0])+'$& $\,_{'+str(MsM12round[1])+'}^{'+str(MsM12round[2])+'} $'
	
print 'Ms02', Ms02str+'\\\\ [1.5mm]'
print 'Mvir', Mvirstr+'\\\\ [1.5mm]'
print 'MsMv', MsMvstr+'\\\\ [1.5mm]'
print 'Ms02Mv', Ms02Mvstr+'\\\\ [1.5mm]'
print 'Meff', Meffstr+'\\\\ [1.5mm]'
print 'M12', M12str+'\\\\ [1.5mm]'
print 'Meffwg', Meffwgstr+'\\\\ [1.5mm]'
print 'MsM12', MsM12str+'\\\\ [1.5mm]'
