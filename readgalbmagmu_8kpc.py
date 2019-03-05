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
from samson_functions import cosmichalo
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
filesuffix='nBxr_galatten_Bband.txt'
#halolist=['1146', '573','m11', '383']
#halolist=['f573','f553','f476','f383','f61']
#halolist=['f573','f553','f476','fm11','f383','f61']
halolist=['f61']
#attenlist=['f383','f61']
attenlist=['f61']
dirlist=['x','y','z']
#halolist=['f573','f553','f476','fm11','f383','1146','573','553','476','m11','383']
#halolist=['f1146','f573','f476','f383','1146','573','553','476','m11','383']
#halolist=['573','m11', '383']
snapsep=5

ftxt = open('/home/tkc004/scratch/galfit-example/EXAMPLE/JANSEN_2000_data.txt', 'r')
ftxt.readline()
dars = ftxt.readlines()
ftxt.close()
Mbl=[]
mube=[]

for line in dars:
	xsd = line.split()
	Mbl=np.append(Mbl, float(xsd[10]))
	mube=np.append(mube, float(xsd[18]))

plt.plot(Mbl, mube, label='Nearby dwarfs',marker='.',ls='none',color='0.5')
for dir in dirlist:
	filesuffix='nB'+dir+'r_galatten_Bband.txt'
	for haloname in halolist:
		print haloname
		ms='s'
		snapl=[]
		bmagl=[]
		muateff=[]
		xhalf=[]
		yhalf=[]
		zhalf=[]
		flog = open('/home/tkc004/samsonprogram/data/'+haloname+filesuffix, 'r')
		flog.readline()
		dars = flog.readlines()
		flog.close()
		haloinfo=cosmichalo(haloname)
		beginno=haloinfo['beginno']
		finalno=haloinfo['finalno']
		rundir=haloinfo['rundir']
		subdir=haloinfo['subdir']
		multifile=haloinfo['multifile']
		hc=haloinfo['halocolor']
		labelname=haloinfo['labelname']
		firever=haloinfo['firever']
		if dir=='x':
			ms='s'
		elif dir=='y':
			ms='<'
			labelname=''
		elif dir=='z':
			ms='o'
			labelname=''
		for line in dars:
			xsd = line.split()
	#		print 'xsd', xsd
			reff=np.append(xhalf, float(xsd[1]))
			bmagl=np.append(bmagl, float(xsd[0]))
			muateff=np.append(muateff, float(xsd[3]))
	#	reff = np.sqrt(np.array(yhalf)*np.array(zhalf))
	#        reff = np.array((xhalf+yhalf+zhalf)/3.0)
	 #       effsur = np.log10(reff*reff*2.0*np.pi*1.425*1.425)*2.5+vmagl+34.97
		muateff = np.array(muateff)
		bmagl = np.array(bmagl)
	#	print 'vmagl, effsur', vmagl, effsur
		plt.plot(bmagl[::snapsep], muateff[::snapsep], label=labelname,marker=ms,ls='none', mfc='none',markeredgecolor=hc,markersize=8)

	for haloname in attenlist:
		ms='s'
		snapl=[]
		bmagl=[]
		muateff=[]
		xhalf=[]
		yhalf=[]
		zhalf=[]
		flog = open('/home/tkc004/samsonprogram/data/'+haloname+'aBin'+dir+'r_galatten_Bband.txt', 'r')
		flog.readline()
		dars = flog.readlines()
		flog.close()
		haloinfo=cosmichalo(haloname)
		beginno=haloinfo['beginno']
		finalno=haloinfo['finalno']
		rundir=haloinfo['rundir']
		subdir=haloinfo['subdir']
		multifile=haloinfo['multifile']
		hc=haloinfo['halocolor']
		labelname=haloinfo['labelname']
		firever=haloinfo['firever']
		if dir=='x':
			ms='s'
		elif dir=='y':
			ms='<'
			labelname=''
		elif dir=='z':
			ms='o'
			labelname=''
		for line in dars:
			xsd = line.split()
	#               print 'xsd', xsd
			reff=np.append(xhalf, float(xsd[1]))
			bmagl=np.append(bmagl, float(xsd[0]))
			muateff=np.append(muateff, float(xsd[3]))
		muateff = np.array(muateff)
		plt.plot(bmagl[::snapsep], muateff[::snapsep],marker=ms,ls='none',mfc=hc,markersize=8)

plt.xlabel(r'$M_{\rm B}$')
plt.ylabel(r'$\mu({\rm B},r_{1/2})$')
plt.ylim([29,16])
plt.xlim([-13.5,-22.5])
plt.legend(loc='best',ncol=3, fontsize=14,numpoints=1)
plt.savefig('figures/bmag_effsur_8kpc.pdf')
plt.clf()
