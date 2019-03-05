from struct import *
import sys
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pylab
import numpy as np
import math
from distcalcs import *
from zip_em_all import *
import time
from scipy.interpolate import interp1d
from datetime import date
import scipy.stats as stats
import scipy.optimize as optimize
import errno
import collections
import binascii
from Sasha_functions import read_halo_history_pep
from samson_functions import cosmichalo
from readsnap_origin import readsnap
#from collections import Counter
# this is a a code that finds the Lagrangian region in the initial condition by particle tracking

today = date.today()
print today

haloname='dm11bhv6rv'
Nsnapstring='600'
the_prefix='snapshot'
the_suffix='.hdf5'
halostr='00'

finname = 'ICs/ics_halo_t_476_hv_3digit.dat'
#we only need non zero shiftback if we use a zoom in box not at the center and shift_back = yes in output of the MUSIC config
#search for shift_x, shift_y, shift_z in  conf_log.txt for the value
#shiftbackx=-32./256.
#shiftbacky=-64./256.
#shiftbackz=0.

shiftbackx=0.
shiftbacky=-55./256.
shiftbackz=-37./256.


single_file = 'y'
m10_type = 'n'
nfiles = 1














#much is copied from snap2tip.c
pc = 3.08567758e18
kpc = pc * 1e3
year = 3.15569e7

H_MASSFRAC = 0.76
PROTONMASS  = 1.6726e-24
GAMMA = 5.0/3.0
BOLTZMANN  = 1.3806e-16

fbar = 0.044 #from paramvalue

if (single_file != 'y'):
	finname_base = finname


Gpos_ari = []
DMpos_ari = []
Wpos_ari = []
Spos_ari = []

Gvel_ari = []
DMvel_ari = []
Wvel_ari = []
Svel_ari = []

Gpid_ari = []
DMpid_ari = []
Wpid_ari = []
Spid_ar = []

Gmass_ari = []
DMmass_ari = []
Wmass_ari = []
Smass_ari = []

Tb_ari = []
rho_ari = []
Neb_ari = []
Nh_ari = []
SML_ari = []
SFRb_ari = []
age_ari = []

Gmet_ari = []
Smet_ari = []

if (m10_type == 'y'):
	P3pos_ari = []
	P5pos_ari = []
	P3vel_ari = []
	P5vel_ari = []
	P3pid_ari = []
	P5pid_ari = []
	P3mass_ari = []
	P5mass_ari = []


#making up some bogus crap
a = 1.0
thetime = 0
theredshift = 1.0 / a - 1.0
redshiftstring = ''
comoving_volume = 1
physical_volume = 1
gadget_volume = 1
boxsize = 60000 #in kpc/h comoving
h = 0.7
omega_m = 0.27
omega_L = 0.73

filecount = 0

while (filecount < nfiles):
	#read the header
	if  (single_file != 'y'):
		finname = finname_base+'.'+str(filecount)
	print 'doing ',finname
	f = open(finname, 'rb')
	headertypei = 6*'i' + 6*'d' + 'd' + 'd' + 'i' + 'i' + 6*'i' + 2*'i' + 4*'d'
	headersizei = 6*4 +   6*8 +    8 +   8 +   4 +   4 +   6*4 + 2*4  +  4*8
	remnanti = 256 - headersizei
	junk = f.read(4)
	headeri = f.read(headersizei)
	foo3i = unpack(headertypei, headeri)
	print 'header: ', foo3i
	Npartsi = foo3i[0:6]
	Mpartsi = foo3i[6:12]
	thetimei = foo3i[12]
	redshifti = foo3i[13]
	flagSFRi = foo3i[14]
	flagfeedbacki = foo3i[15]
	Nparts2i = foo3i[16:22]
	moreflagsi = foo3i[22:24]
	cosmosi = foo3i[24:28]
	boxsizei = cosmosi[0]
	omega_mi = cosmosi[1]
	omega_Li = cosmosi[2]
	hi = cosmosi[3]
	junk = f.read(remnanti)
	junk = f.read(4)
	junk = f.read(4)
	NSPHi = Npartsi[0]
	Ndarki = Npartsi[1]
	Nweirdi = Npartsi[2]
	Nstari = Npartsi[4]
	NP3i = Npartsi[3]
	NP5i = Npartsi[5]

	posdi = '(3,)f4'
	posdtypei = np.dtype(posdi)

	DMpos_ari.append(  np.fromfile(f, dtype=posdtypei, count = Ndarki))
	Wpos_ari.append(  np.fromfile(f, dtype=posdtypei, count = Nweirdi))
	junk = unpack('i',f.read(4))
	junk = unpack('i',f.read(4))

	DMvel_ari.append( np.fromfile(f, dtype=posdtypei, count = Ndarki))
        Wvel_ari.append(  np.fromfile(f, dtype=posdtypei, count = Nweirdi))
	junk = unpack('i',f.read(4))
	junk = unpack('i',f.read(4))

#	piddi = 'i4'
#	piddtypei = np.dtype(piddi)

#	piddtypei = np.dtype(np.int64)

	piddtypei = np.dtype(np.int32)

	DMpid_ari.append(   np.fromfile(f, dtype=piddtypei, count = Ndarki))
        Wpid_ari.append(  np.fromfile(f, dtype=piddtypei, count = Nweirdi))
	junk = unpack('i',f.read(4))
	massdi = 'f4'
	massdtypei = np.dtype(massdi)

	DMmass_ari.append( np.zeros(Ndarki))
	DMmass_ari[filecount][0:Ndarki]=Mpartsi[1]
        Wmass_ari.append(  np.fromfile(f, dtype=massdtypei, count = Nweirdi))

	f.close()
	filecount+=1


print np.count_nonzero(Gpos_ari)

DMposi = zip_em_all(DMpos_ari, nfiles) 
del(DMpos_ari)

DMveli = zip_em_all(DMvel_ari, nfiles)
del(DMvel_ari)

DMpidi = zip_em_all(DMpid_ari, nfiles)


del(DMpid_ari)


DMmassi = zip_em_all(DMmass_ari, nfiles)

del(DMmass_ari)


Wposi = zip_em_all(Wpos_ari, nfiles)
del(Wpos_ari)

Wveli = zip_em_all(Wvel_ari, nfiles)
del(Wvel_ari)

Wpidi = zip_em_all(Wpid_ari, nfiles)
del(Wpid_ari)

Wmassi = zip_em_all(Wmass_ari, nfiles)
del(Wmass_ari)


totalDmassi = np.sum(DMmassi)


xDi = DMposi[:,0]
yDi = DMposi[:,1]
zDi = DMposi[:,2]

xWi = Wposi[:,0]
yWi = Wposi[:,1]
zWi = Wposi[:,2]


#print 'Counter(DMpidi).most_common(3)',Counter(DMpidi).most_common(3)
#print 'np.amax(DMpidi)',np.amax(DMpidi)
#print 'np.amin(DMpidi)',np.amin(DMpidi)
#print 'len(DMpidi)', len(DMpidi)
#print 'NSPHi', NSPHi
#print 'Ndarki', Ndarki
#print 'Nweirdi', Nweirdi
#print 'np.amax(DMmassi)', np.amax(DMmassi)
#print 'np.amin(DMmassi)', np.amin(DMmassi)
#print 'DMmassi', DMmassi

ai = float(thetimei)
print 'using a = ',ai
theredshifti = np.absolute(1.0/ai - 1)
redshiftstringi = "{0:.3f}".format(theredshifti)
print 'z = ',"{0:.3f}".format(theredshifti)

haloinfo=cosmichalo(haloname)
beginno=haloinfo['beginno']
finalno=haloinfo['finalno']
rundir=haloinfo['rundir']
subdir=haloinfo['subdir']
halostr=haloinfo['halostr']
labelname=haloinfo['labelname']
maindir=haloinfo['maindir']
firever=haloinfo['firever']
multifile=haloinfo['multifile']
usepep=haloinfo['usepep']
h0=haloinfo['h0']
the_snapdir = '/home/tkc004/'+maindir+'/'+rundir+'/'+subdir

print 'the_snapdir', the_snapdir
print 'Nsnapstring', Nsnapstring
DM = readsnap(the_snapdir, Nsnapstring, 1, snapshot_name=the_prefix, extension=the_suffix,h0=1,cosmological=1)
W = readsnap(the_snapdir, Nsnapstring, 2, snapshot_name=the_prefix, extension=the_suffix,h0=1,cosmological=1)
header=readsnap(the_snapdir, Nsnapstring, 1, snapshot_name=the_prefix, extension=the_suffix,h0=1,cosmological=1,header_only=1)
hubble = header['hubble']
atime = header['time']
halosA = read_halo_history_pep(rundir, finalno, beginno=beginno, singlesnap=1, firever=firever,halonostr=halostr, hubble=hubble, comoving=0, maindir=maindir)
redlist = halosA['redshift']
haloid = halosA['ID']
a_scale = 1.0/(1.0+redlist)
xcen = halosA['x']
ycen = halosA['y']
zcen = halosA['z']
Rvir = halosA['R']
Mvir = halosA['M']
DMpos = DM['p'][:,:]
DMpid = DM['id']
DMm = DM['m']
Wpos = W['p'][:,:]
Wpid = W['id']
Wm = W['m']
DMx = DMpos[:,0]
DMy = DMpos[:,1]
DMz = DMpos[:,2]
Wx = Wpos[:,0]
Wy = Wpos[:,1]
Wz = Wpos[:,2]
DMxrel = DMx-xcen
DMyrel = DMy-ycen
DMzrel = DMz-zcen
DMr = np.sqrt(DMxrel*DMxrel+DMyrel*DMyrel+DMzrel*DMzrel)
cutrv = DMr<Rvir
DMmcut = DMm[cutrv]
Wxrel = Wx-xcen
Wyrel = Wy-ycen
Wzrel = Wz-zcen
Wr = np.sqrt(Wxrel*Wxrel+Wyrel*Wyrel+Wzrel*Wzrel)
cutrvW = Wr<Rvir
Wmcut = Wm[cutrvW]

Wmsum =  np.sum(Wmcut)
DMmsum = np.sum(DMmcut)

print 'Wmass', Wmsum
print 'DMmass', DMmsum
print 'Wmass/DMmass', Wmsum/DMmsum


DMpid=np.array(DMpid[cutrv])
idint0=np.in1d(DMpidi,DMpid)
print 'DMpidi[idint0]', np.sort(DMpidi[idint0])
if len(DMpidi[idint0])==0:
	print 'no id'
idneeded=np.in1d(DMpidi,DMpid)
print 'xmax', np.amax(xDi[idneeded])
print 'ymax', np.amax(yDi[idneeded])
print 'zmax', np.amax(zDi[idneeded])

print 'xmin', np.amin(xDi[idneeded])
print 'ymin', np.amin(yDi[idneeded])
print 'zmin', np.amin(zDi[idneeded])

Wpid=np.array(Wpid[cutrvW])
idint0W=np.in1d(Wpidi,Wpid)
if len(Wpidi[idint0W])==0:
        print 'no id'
idneededW=np.in1d(Wpidi,Wpid)
print 'xmaxW', np.amax(xWi[idneededW])
print 'ymaxW', np.amax(yWi[idneededW])
print 'zmaxW', np.amax(zWi[idneededW])

print 'xminW', np.amin(xWi[idneededW])
print 'yminW', np.amin(yWi[idneededW])
print 'zminW', np.amin(zWi[idneededW])



#random = np.vstack((xDi[idneeded]/40000.-shiftbackx,yDi[idneeded]/40000.-shiftbacky,zDi[idneeded]/40000.-shiftbackz)).T

xcor = xDi[idneeded]/40000.-shiftbackx
ycor = yDi[idneeded]/40000.-shiftbacky
zcor = zDi[idneeded]/40000.-shiftbackz

print 'xcor', xcor
print 'ycor', ycor
print 'zcor', zcor


xcorW = xWi[idneededW]/40000.-shiftbackx
ycorW = yWi[idneededW]/40000.-shiftbacky
zcorW = zWi[idneededW]/40000.-shiftbackz


#xcor = random[0]
#ycor = random[1]
#zcor = random[2]


#filename='lagregion/parlist_halo_dm.dat'
#infoname='lagregion/info.txt'
#text_file = open(filename, "w")
#ncount=0
#for x in random:
#	for z in x:
#		y='%.6f' % float(z)
#		text_file.write(y)
#		text_file.write(' ')
#		if (ncount==2):
#			text_file.write('\n')
#			ncount=-1
#		ncount=ncount+1
#text_file.close()

#info_file = open(infoname,"w")
#info_file.write('Mvir   '+'Rvir   '+'\n')
#info_file.write(str(Mvir)+'  '+str(Rvir)+'\n')
#info_file.close()


xedges=np.linspace(0.,1.,num=1000)
yedges=np.linspace(0.,1.,num=1000)
#print 'xedges', xedges
#print 'yedges', yedges
#histogram reverse x y coordinate
Hxy, xedges, yedges = np.histogram2d(ycor, xcor, bins=(xedges, yedges))
plt.imshow(Hxy, interpolation='nearest', origin='low',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('figures/zooms_xy_his.png')
plt.clf()

Hxz, xedges, yedges = np.histogram2d(zcor, xcor, bins=(xedges, yedges))
plt.imshow(Hxz, interpolation='nearest', origin='low',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.xlabel('x')
plt.ylabel('z')
plt.savefig('figures/zooms_xz_his.png')
plt.clf()

Hyz, xedges, yedges = np.histogram2d(zcor, ycor, bins=(xedges, yedges))
plt.imshow(Hyz, interpolation='nearest', origin='low',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.xlabel('y')
plt.ylabel('z')
plt.savefig('figures/zooms_yz_his.png')
plt.clf()


Hxy, xedges, yedges = np.histogram2d(ycorW, xcorW, bins=(xedges, yedges))
plt.imshow(Hxy, interpolation='nearest', origin='low',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('figures/zooms_xy_hisW.png')
plt.clf()

Hxz, xedges, yedges = np.histogram2d(zcorW, xcorW, bins=(xedges, yedges))
plt.imshow(Hxz, interpolation='nearest', origin='low',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.xlabel('x')
plt.ylabel('z')
plt.savefig('figures/zooms_xz_hisW.png')
plt.clf()

Hyz, xedges, yedges = np.histogram2d(zcorW, ycorW, bins=(xedges, yedges))
plt.imshow(Hyz, interpolation='nearest', origin='low',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.xlabel('y')
plt.ylabel('z')
plt.savefig('figures/zooms_yz_hisW.png')
plt.clf()
