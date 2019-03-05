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

NewG = 6.674e-11 #SI unit
Msun_in_kg = 2e30
kpc_in_m = 3.085e19

#dirneed=['f383']
#dirneed=['f476']
#dirneed=['f61']
#dirneed=['fm11']
#dirneed=['f573']
#dirneed=['f573','f553','f476','fm11','f383','f61']
#dirneed=['1146']
#dirneed=['m11']
#dirneed=['573']
#dirneed=['573','553','476']
#dirneed=['553','476']
dirneed=['m11bhv']
Nsnap=346




for runtodo in dirneed:
	haloinfo=cosmichalo(runtodo)
	beginno=haloinfo['beginno']
	finalno=haloinfo['finalno']
	rundir=haloinfo['rundir']
	subdir=haloinfo['subdir']
	maindir=haloinfo['maindir']
	multifile=haloinfo['multifile']
	halocolor=haloinfo['halocolor']
	halostr=haloinfo['halostr']
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
	header = readsnap(the_snapdir, Nsnapstring, 1, snapshot_name=the_prefix, extension=the_suffix, header_only=1, h0=1,cosmological=1)
	G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, cosmological=1, h0=1)
	ascale = header['time']
	thisz = 1./ascale-1.
	#print 'thisz', thisz
	hubble = header['hubble']
	#print 'hubble', hubble
	Gmass = G['m'][:]
	print 'np.amin(Gmass)', np.amin(Gmass)
	print 'np.amax(Gmass)', np.amax(Gmass)
