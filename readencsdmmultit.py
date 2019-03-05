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
rcParams['figure.figsize'] = 5, 5
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
rcParams['xtick.minor.width'] = 2
rcParams['ytick.minor.width'] = 2
rcParams['text.usetex'] = True
import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='inferno', cmap=cmaps.inferno)
plt.register_cmap(name='plasma', cmap=cmaps.plasma)
from readsnap_samson import *
from Sasha_functions import *
from gadget_lib.cosmo import *
from samson_functions import *
from enclosedmass import enclosedmass
lsl = ['solid','dashed','dashdot','dotted']
dirneed=[['f573',2.0,6.5,10.0,13.8],['f553',2.0,6.5,10.0,13.8],['f476',2.0,6.5,10.0,13.8],['fm11',2.0,6.5,10.0,13.8],['f383',2.0,6.5,10.0,13.8],['f61',2.0,6.5,10.0,13.8]]
#dirneed=[['f573',8.2,13.5],['f553',5.7,13.5],['f476',2.6,13.5],['fm11',2.3,6.0],['f383',2.0,6.5],['f61',2.2,2.4]]
#dirneed=[['fm12z',1.8,2.0]]
#dirneed=[['fm12i',1.8,2.0]]
#dirneed=[['fm12i',1.8,2.0],['fm12c',1.8,2.0]]
maxr=120.0
minr=0.5
galcen=0
#Nsnap=440
DMonly=1

for runtodol in dirneed:
	runtodo=runtodol[0]
	for icount in range(len(runtodol)-1):
		time=runtodol[icount+1]
		print 'runtodo, time', runtodo, time
		emdata = enclosedmass(runtodo, time,minr,maxr, galcen=galcen)
		rlist =  emdata['rlist']
		vollist = emdata['vollist']
		Gmlist = emdata['Gmlist']
		Smlist = emdata['Smlist']
		DMmlist = emdata['DMmlist']
		haloinfo=cosmichalo(runtodo)
		labelname=haloinfo['labelname']
		halocolor=haloinfo['halocolor']
		icolor=haloinfo['icolor']
		print 'halocolor', halocolor
		print 'labelname', labelname
		if DMonly==1:
			plt.plot(rlist,DMmlist,ls=lsl[icount], color=cmaps.plasma(icolor),lw=2,label='t = '+str(time)+'Gyr')
		else:
			plt.plot(rlist,Smlist+DMmlist,ls=ilsl[icount], color=cmaps.plasma(icolor), lw=2)
		f = open('/home/tkc004/samsonprogram/data/DMSMenc'+runtodol[0]+'_info.txt', 'w')
		for i in range(len(rlist)):
			f.write(str(rlist[i])+'      '+str(Smlist[i])+'      '+str(DMmlist[i])+'  \n')
		f.close()
	plt.errorbar([3.2,8.1], [2.6e9,4.5e9], yerr=[[1.3e9,3.4e9],[2.8e9,2.8e9]],label='VCC1287', fmt='o')
	plt.errorbar(4.7, 7.0e9, xerr=[[0.2],[0.2]], yerr=[[2.0e9],[3.0e9]],label='Dragonfly44',  fmt='s')
	reffp = 4.7; meffp = 7.0e9;
	dataNFW = calNFW(rlist,reffp,meffp,zr=0)
	mnfw = dataNFW['mnfw']
	#coreect for baryon:
	plt.plot(rlist,mnfw,ls='dotted', color='k',lw=2,label='NFW (t = 13.6 Gyr)')
        dataNFW = calNFW(rlist,reffp,meffp,zr=3.24)
        mnfw = dataNFW['mnfw']
	plt.plot(rlist,mnfw,ls='solid', color='k',lw=2,label='NFW (t = 2.0 Gyr)')
	plt.xscale('log')
	plt.yscale('log')
	plt.xlim([minr,maxr])
	plt.ylim([2e7,6e11])
	plt.xlabel(r'$r\; ({\rm kpc})$')
	if DMonly==1:
		plt.ylabel(r'${\rm M}_{\rm enc,DM} ({\rm M}_{\odot})$')
	else:
		plt.ylabel(r'${\rm M}_{\rm enc} ({\rm M}_{\odot})$')
	plt.legend(loc='best',fontsize=12)
	plt.title(labelname)
	if DMonly==1:
		plt.savefig('figures/Mencdmmultitime_udg_dmonly'+runtodo+'.pdf')
	else:
		plt.savefig('figures/Mencdmmultitime_udg_sdmonly'+runtodo+'.pdf')
	plt.clf()
