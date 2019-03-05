import os
#import pyfits
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
from enclosedmass import enclosedmass
#dirneed=['m10h1146','m10h573','m10h553','m10h476','m11h383','m11']
#dirneed=['f1146','f573','f553','f476','f383','fm11']
dirneed=[['f476',92.0],['m11bmhdcv',92.0],['m11bcr_b_70',92.0],['m11bcr_700',92.0]]
#dirneed=[['f573',85.0],['f553',90.0],['f476',92.0],['fm11',140.0],['f383',140.0],['f61',210.0],['m11bmhdcv',92.0],['m11gmhdcv',210]]
#dirneed=[['f573',8.2,13.5],['f553',5.7,13.5],['f476',2.6,13.5],['fm11',2.3,6.0],['f383',2.0,6.5],['f61',2.2,2.4]]
#dirneed=['m11']
#dirneed=['553','476']
minr=1.0
time=13.8
galcen=0
#Nsnap=440
runlist=[]
Gmelist=[]
Smelist=[]
HImllist=[]
fHIlist=[]
for runtodol in dirneed:
	runtodo=runtodol[0]
	maxr=runtodol[1]*0.5
	emdata = enclosedmass(runtodo, time,minr,maxr, galcen=galcen,needHI=1)
	rlist =  emdata['rlist']
	vollist = emdata['vollist']
	Gmlist = emdata['Gmlist']
	Smlist = emdata['Smlist']
	DMmlist = emdata['DMmlist']
	HImlist = emdata['HImlist']
	haloinfo=cosmichalo(runtodo)
	labelname=haloinfo['labelname']
	halocolor=haloinfo['halocolor']
	Gmrv25 = np.interp(runtodol[1]*0.25,rlist,Gmlist)
	Smrv25 = np.interp(runtodol[1]*0.25,rlist,Smlist)
        HImrv25 = np.interp(runtodol[1]*0.25,rlist,HImlist)
	fHI = HImrv25/(HImrv25+Smrv25)
	runlist = np.append(runlist,runtodo)
	Gmelist = np.append(Gmelist,Gmrv25)
	Smelist = np.append(Smelist,Smrv25)
	HImllist = np.append(HImllist,HImrv25)
	fHIlist = np.append(fHIlist,fHI)
	#print 'runtodo, Gm<r25, mHI<r25', runtodo, Gmrv25, HImrv25
	#print 'Sm<r25', Smrv25
	#print 'fHI', HImrv25/(HImrv25+Smrv25)
print 'runlist', runlist
print 'Smelist', Smelist
print 'HImllist', HImllist
print 'fHIlist', fHIlist

