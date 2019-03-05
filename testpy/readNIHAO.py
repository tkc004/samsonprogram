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
rcParams['text.usetex'] = True
import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='inferno', cmap=cmaps.inferno)
plt.register_cmap(name='plasma', cmap=cmaps.plasma)
ftxt = open('/home/tkc004/samsonprogram/data/NIHAOXII.txt', 'r')
ftxt.readline()
dars = ftxt.readlines()
ftxt.close()
Mbl=[]
Msl=[]
Mngl=[]

for line in dars:
	xsd = line.split()
	Mbl=np.append(Mbl, float(xsd[9]))
	Msl=np.append(Msl, float(xsd[6]))
	Mngl=np.append(Mngl,float(xsd[7]))

Mgl = np.log10(np.power(10,Mbl)-np.power(10,Msl))
gfrac = np.power(10,Mgl-Mbl)
Msngl = np.log10(np.power(10,Msl)+np.power(10,Mngl))
ngfrac = np.power(10,Mngl-Msngl)
plt.scatter(Msl, gfrac, label='NIHAOII',marker='o',color='b')
plt.scatter(Msl, ngfrac, label='neutral only', marker='^',color='r')
plt.xlabel(r'$\log(M_{\rm s}/M_\odot)$')
plt.ylabel(r'$\log(M_{\rm gas}/M_{\rm bar})$')
plt.ylim([0,1])
plt.xlim([6,11])
plt.legend(loc='best',ncol=3, fontsize=14,numpoints=1)
plt.savefig('figures/msgfrac_NIHAO.pdf')
plt.clf()
