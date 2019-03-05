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
rcParams['text.usetex'] = True
import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='inferno', cmap=cmaps.inferno)
plt.register_cmap(name='plasma', cmap=cmaps.plasma)
from readfittxtdat import *
from samson_functions import cosmichalo
filesuffix='_gallightinfo_Vband.txt'
#halolist=['1146', '573','m11', '383']
#halolist=['m10','f1297','1297','f1146','1146','f573','573','f553','553','f476','476','m11','f383','383']
#halolist=['573','m11', '383']
#halolist=['m11bcr_b_70', 'm11bmhdcv', 'm11fcr_b_70','m11fmhdcv','m11hmhdcv']
#halolist=['f573','f553','f476']
#halolist=['f573','f553','f476','f383']
halolist=['f476']
#halolist=['f573','f553','f476','fm11','f383']
#halolist=['f573','f553','f476','fm11','f383','f61','fm11h','fm11i']
#halolist=['f573','f553','f476','fm11','f383']
snapd=277
snapsep=5


for haloname in halolist:
	print 'haloname', haloname
	ms='h'
        haloinfo=cosmichalo(haloname)
        rundir=haloinfo['rundir']
        subdir=haloinfo['subdir']
        maindir=haloinfo['maindir']
        multifile=haloinfo['multifile']
        hc=haloinfo['halocolor']
        halostr=haloinfo['halostr']
	labelname=haloinfo['labelname']
	icolor=haloinfo['icolor']
	snapl=[]
	vmagl=[]
#	effsur=[]
	xhalf=[]
	yhalf=[]
	zhalf=[]
	flog = open('/home/tkc004/samsonprogram/data/'+haloname+filesuffix, 'r')
	flog.readline()
	dars = flog.readlines()
	flog.close()
	for line in dars:
		xsd = line.split()
#		print 'xsd', xsd
		snapl=np.append(snapl, float(xsd[0]))
		xhalf=np.append(xhalf, float(xsd[13]))
                yhalf=np.append(yhalf, float(xsd[14]))
                zhalf=np.append(zhalf, float(xsd[15]))
		vmagl=np.append(vmagl, float(xsd[16]))
#		effsur=np.append(effsur, float(xsd[17]))
	snapl=np.array(snapl)
	xhalf=np.array(xhalf[snapl>snapd])
	yhalf=np.array(yhalf[snapl>snapd])
	zhalf=np.array(zhalf[snapl>snapd])
	vmagl=np.array(vmagl[snapl>snapd])
	xhalf=xhalf[::snapsep]
	yhalf=yhalf[::snapsep]
	zhalf=zhalf[::snapsep]
	vmagl=vmagl[::snapsep]
#	reff = np.sqrt(np.array(yhalf)*np.array(zhalf))
        reff = np.sqrt((xhalf*xhalf+yhalf*yhalf+zhalf*zhalf)/3.)
	maxhalf = np.fmax(xhalf,yhalf)
	maxhalf = np.fmax(maxhalf,zhalf)
        minhalf = np.fmin(xhalf,yhalf)
        minhalf = np.fmin(minhalf,zhalf)
	axisratio = minhalf/maxhalf
	ell = np.sqrt(maxhalf*maxhalf-minhalf*minhalf)/maxhalf 
	#print 'xhalf', xhalf
	#print 'yhalf', yhalf
	#print 'zhalf', zhalf
	#print 'maxhalf', maxhalf
	#print 'minhalf', minhalf
#        effsur = np.log10(reff*reff*2.0*np.pi*1.425*1.425)*2.5+vmagl+34.97
#	effsur = np.array(effsur)
#	print 'vmagl, effsur', vmagl, effsur
#	if haloname[0]=='f':
#		plt.plot(vmagl, effsur, label=labelname,marker=ms,ls='none', mfc=cmaps.plasma(icolor),ms=12)
#	else:
#		plt.plot(vmagl, effsur, label=labelname,marker=ms,ls='none',mfc=cmaps.plasma(icolor),ms=12)
#plt.plot(vmagLG, effsurLG, label='LG', marker='^', ls='none',ms=8)
#plt.plot(vmagnb, effsurnb, label='Nearby', marker='*',ls='none',ms=8)

#	print 'axisratio, reff', axisratio, reff

	plt.plot(reff, ell, label=labelname,marker=ms,ls='none',ms=12)

plt.xlabel(r'$r_{\rm mean}$')
#plt.ylabel(r'$r_{\rm min}/r_{\rm max}$')
plt.ylabel('Ellipticity')
plt.legend(loc='best', fontsize=14, ncol=4,numpoints=1)
plt.savefig('figures/reffellMs_m11b.pdf')
plt.clf()
