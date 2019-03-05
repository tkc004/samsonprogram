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
runtodo='fm11v'
galcen=0
#Nsnap=440

gl = []
sl = []
dml = []
snapl = []
snaplist=range(300,601,10)


haloinfo=cosmichalo(runtodo)
rundir=haloinfo['rundir']
maindir=haloinfo['maindir']
subdir=haloinfo['subdir']
halostr=haloinfo['halostr']
snumadd=haloinfo['snumadd']
halosA = read_halo_history(rundir, halonostr=halostr, maindir=maindir,snumadd=snumadd)
redlist = halosA['redshift']
Rvirlist = halosA['R']
alist = np.array(1./(1.+np.array(redlist)))
Rvirlist = np.array(Rvirlist)*np.array(alist)/0.702
Rvir = Rvirlist[-1]
print 'Rvir', Rvir


for snapno in snaplist:
	try:
		emdata = enclosedmass(runtodo, snapno, 0.001, 0.02*Rvir, usesnapno=1,needDM=1,needS=1,needG=1)
		rlist =  emdata['rlist']
		vollist = emdata['vollist']
		Gmlist = emdata['Gmlist']
		Smlist = emdata['Smlist']
		DMmlist = emdata['DMmlist']
		haloinfo=cosmichalo(runtodo)
		labelname=haloinfo['labelname']
		halocolor=haloinfo['halocolor']
		icolor=haloinfo['icolor']
		gl.append(Gmlist[-1])
		sl.append(Smlist[-1])
		dml.append(DMmlist[-1])
		snapl.append(snapno)
	except (IOError,KeyError):
		continue
	
plt.plot(snapl,gl,marker='o',label='Gas')
plt.plot(snapl,sl,marker='*',label='Star')
plt.plot(snapl,dml,marker='s',label='DM')
	
plt.yscale('log')
plt.xlabel('Snapshot number')
plt.ylabel(r'${\rm M}_{\rm enc} ({\rm M}_{\odot})$')
plt.legend(loc='best',fontsize=12)
filename='figures/Mwithin_'+runtodo+'.pdf'
print 'filename', filename
plt.savefig(filename)
plt.clf()
