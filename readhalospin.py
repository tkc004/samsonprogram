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
rcParams['xtick.major.width'] = 2
rcParams['ytick.major.width'] = 2
rcParams['xtick.minor.width'] = 2
rcParams['ytick.minor.width'] = 2
from readsnap_samson import *
from Sasha_functions import *
from gadget_lib.cosmo import *
from samson_functions import *
from halospin import gethalospin
#dirneed=['m10h1146','m10h573','m10h553','m10h476','m11h383','m11']
#dirneed=['f1146','f573','f553','f476','f383','fm11']
#dirneed=[['f573',13.8],['f553',13.8],['f476',13.8],['f383',13.8],['f61',13.8]]
#dirneed=[['d573',13.8],['d553',13.8],['d476',13.8],['d383',13.8],['df61',13.8]]
#dirneed=[['d573',6],['d553',6],['d476',6],['d383',6],['df61',6]]
#dirneed=[['f61',13.8]]
dirneed=[['f573',6],['f553',6],['f476',6],['f383',6],['f61',6]]
#dirneed=[['f573',2],['f553',2],['f476',2],['f383',2],['f61',2]]
#dirneed=['m11']
#dirneed=['553','476']
useAHF=0
rfactor=1


for runtodol in dirneed:
	runtodo=runtodol[0]
	time = runtodol[1]
	spinlambda = gethalospin(runtodo, time, rfactor=rfactor,useAHF=useAHF)
	print 'spinlambda', spinlambda 
