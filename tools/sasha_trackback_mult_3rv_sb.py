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
#from collections import Counter


today = date.today()
print today





doplots = 'n' #do any plots at all
w='n' #'w' for write to file
showplots = 'n' #'y' to display plots
weighter ='s' #'m' or 's' apply histogram weights?]
write_outflow_file = 'y' #write the outlfows to txt file
use_halo_tree = 'y' #AHF Merger Trace info instead of the actual halo number
Rvir_must_grow = 'n'
use_mbp = 'n' # use most bound particle position and velocity for halo center.
single_file = 'y'
m10_type = 'n'
debug = 'n'
fun_with_metals = 'y' #what should i do with metals?  OK i made a note about it.
nfiles = 1

galmode = 'n'  #'g' for galaxy mode - use SKID centers
starmode = 'S' #  'S' to use info about young stars
print_all_shells = 'y' #write outflow in each shell?
use_tempcut = 'n'
use_metalcut = 'n'
use_rhocut = 'n' 
velocity_diagnostics = 'y'
use_inner_velocity = 'n'
find_exact_SFR = 'y'
write_exact_SFR = 'n'
use_just_outflows ='n'

use_galactic_gas = 'y' # 'y' to use galaxy above the overdensity creterion, 'n' for halo gas instead 
cooltemp = 'n'  #reverses T cut
use_metal_poor = 'n' #reverses metalcut

snapno = 440
halo_to_do = range(1200,1600) #halos to do
#halo_to_do = [141]



#THRESHOLDS - adjust as needed
Inner_SF_thresh = 0.2 #distinguish between "inner" and "outer" SF
Inner_vel_thresh = 0.1
MetRichThresh = 1e-3 #for the purpose of cutting "metal rich gas"
rho_cut_value = 100000 #separates 'galactic' gas. 0 for using all gas
Temp_cut_value = 100000 #seperates 'cool' gas 
SFR_time_range = 3e7 #timescale in years used for SFR calculations 
fixed_outflow_cut = -9999 #fixed cut for determining what an 'outflow'. If less than -1, uses 1-D velocity dispersion of halo
fixed_inflow_cut = -9999
R_bins = 10 # this is the number of bins you'll get from 0 to Rvir. IN reality the outflows are computed starting with bin #2 




shiftbackx=-32./256.
shiftbacky=-64./256.
shiftbackz=0.












finname = 'P-Gadget3/ics_l8l_10lv.dat'

#finname = 'snaps/snapshot_440'


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


totalDmassi = np.sum(DMmassi)


xDi = DMposi[:,0]
yDi = DMposi[:,1]
zDi = DMposi[:,2]


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





















if (int(snapno) < 10):
	Nsnapstring = '00'+str(snapno)
elif (int(snapno) < 100):
	Nsnapstring = '0'+str(snapno)
else:
	Nsnapstring = str(snapno)
finname = 'snaps/snapshot_' + Nsnapstring
if (single_file != 'y'):
	finname = 'snapdir_'+Nsnapstring+'/snap_convertedshot_' + Nsnapstring


directory_name='figures/snapshot'+ Nsnapstring
path='figures'

try:
  os.mkdir(directory_name)
except OSError as exc:
    if exc.errno == errno.EEXIST and os.path.isdir(path):
	pass




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


Gpos_ar = []
DMpos_ar = []
Wpos_ar = []
Spos_ar = []

Gvel_ar = []
DMvel_ar = []
Wvel_ar = []
Svel_ar = []

Gpid_ar = []
DMpid_ar = []
Wpid_ar = []
Spid_ar = []

Gmass_ar = []
DMmass_ar = []
Wmass_ar = []
Smass_ar = []

Tb_ar = []
rho_ar = []
Neb_ar = []
Nh_ar = []
SML_ar = []
SFRb_ar = []
age_ar = []

Gmet_ar = []
Smet_ar = []

if (m10_type == 'y'):
	P3pos_ar = []
	P5pos_ar = []
	P3vel_ar = []
	P5vel_ar = []
	P3pid_ar = []
	P5pid_ar = []
	P3mass_ar = []
	P5mass_ar = []


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
	headertype = 6*'i' + 6*'d' + 'd' + 'd' + 'i' + 'i' + 6*'i' + 2*'i' + 4*'d'
	headersize = 6*4 +   6*8 +    8 +   8 +   4 +   4 +   6*4 + 2*4  +  4*8
	remnant = 256 - headersize
	junk = f.read(4)
	header = f.read(headersize)
	foo3 = unpack(headertype, header)
	print 'header: ', foo3
	Nparts = foo3[0:6]
	Mparts = foo3[6:12]
	thetime = foo3[12]
	redshift = foo3[13]
	flagSFR = foo3[14]
	flagfeedback = foo3[15]
	Nparts2 = foo3[16:22]
	moreflags = foo3[22:24]
	cosmos = foo3[24:28]
	boxsize = cosmos[0]
	omega_m = cosmos[1]
	omega_L = cosmos[2]
	h = cosmos[3]
	junk = f.read(remnant)
	junk = f.read(4)
	junk = f.read(4)
	NSPH = Nparts[0]
	Ndark = Nparts[1]
	Nweird = Nparts[2]
	Nstar = Nparts[4]
	NP3 = Nparts[3]
	NP5 = Nparts[5]

	posd = '(3,)f4'
	posdtype = np.dtype(posd)

	DMpos_ar.append(  np.fromfile(f, dtype=posdtype, count = Ndark))
        Wpos_ar.append(  np.fromfile(f, dtype=posdtype, count = Nweird))

	junk = unpack('i',f.read(4))
	junk = unpack('i',f.read(4))

	DMvel_ar.append( np.fromfile(f, dtype=posdtype, count = Ndark))
        Wvel_ar.append(  np.fromfile(f, dtype=posdtype, count = Nweird))

	junk = unpack('i',f.read(4))
	junk = unpack('i',f.read(4))

#	pidd = 'i4'
#	piddtype = np.dtype(np.int64)
	piddtype = np.dtype(np.int32)
#	piddtype = np.dtype(pidd)

	DMpid_ar.append(   np.fromfile(f, dtype=piddtype, count = Ndark))
        Wpid_ar.append(  np.fromfile(f, dtype=piddtype, count = Nweird))
	

	junk = unpack('i',f.read(4))
	junk = unpack('i',f.read(4))

	massd = 'f4'
	massdtype = np.dtype(massd)

	DMmass_ar.append( np.zeros(Ndark))
	DMmass_ar[filecount][0:Ndark]=Mparts[1]
        Wmass_ar.append(  np.fromfile(f, dtype=massdtype, count = Nweird))
	

	
	
	f.close()
	filecount+=1


print np.count_nonzero(Gpos_ar)

DMpos = zip_em_all(DMpos_ar, nfiles) 
del(DMpos_ar)

DMvel = zip_em_all(DMvel_ar, nfiles)
del(DMvel_ar)

DMpid = zip_em_all(DMpid_ar, nfiles)





del(DMpid_ar)


DMmass = zip_em_all(DMmass_ar, nfiles)

del(DMmass_ar)


totalDmass = np.sum(DMmass)


#print 'Counter(DMpid).most_common(3)',Counter(DMpid).most_common(3)
#print 'np.amax(DMpid)',np.amax(DMpid)
#print 'np.amin(DMpid)',np.amin(DMpid)
#print 'len(DMpid)', len(DMpid)
#print 'NSPH', NSPH
#print 'Ndark', Ndark
#print 'Nweird', Nweird
#print 'np.amax(DMmass)', np.amax(DMmass)
#print 'np.amin(DMmass)', np.amin(DMmass)
#print 'DMmass', DMmass



a = float(thetime)
print 'using a = ',a
theredshift = np.absolute(1.0/a - 1)
redshiftstring = "{0:.3f}".format(theredshift)
print 'z = ',"{0:.3f}".format(theredshift)



#print MeanWeights
#print MeanWeights_con  
#print (1.0/MeanWeights_con)
#converted_rho = rho * UnitMass_in_g * 6.02e23 / MeanWeights_con
#converted_rho /= (UnitLength_in_cm * UnitLength_in_cm *UnitLength_in_cm)
#print Neb
#print converted_rho
#print rho
#print Nh 



#now for halo catalog

LU = 1
MU = 1
LU2 = 1
finname = 'Pep/snap'+ Nsnapstring+ 'RPep.z'+redshiftstring+'.AHF_halos'
f = open(finname)
f.readline()
dars = f.readlines()
 
halos  =[]
halID = []

#most massive galaxy in each halo from skid sorted catalog
#as of right now this is bogus and obsolete for what i'm doing
halogalmatch = [1, 0, 2, 6, 4, 10, 7, 6, 5, 17]

#i wish i didn't have to use obtuse variable names but i do.
for line in dars:
	xsd = line.split()
	halID.append(int(xsd[0]))
	centX = float(xsd[5]) 
	centY = float(xsd[6]) 
	centZ = float(xsd[7]) 
	velX = float(xsd[8]) 
	velY = float(xsd[9]) 
	velZ = float(xsd[10]) 
	M = float(xsd[3]) 
	Vmax = float(xsd[16]) 
	Vsig = float(xsd[18]) 
	Rvir = float(xsd[11]) 
	#12/09/13 using Pep positions for Mgas and Mstar 
	Mgas = float(xsd[53])
	Mstar = float(xsd[73])
	mbpcentX = float(xsd[46])
	mbpcentY = float(xsd[47])
	mbpcentZ = float(xsd[48])
	mbpvelX = float(xsd[43])
	mbpvelY = float(xsd[44])
	mbpvelZ = float(xsd[45])
	fMhires = float(xsd[37])
	if (use_mbp == 'y'):
		centX = mbpcentX 
		centY = mbpcentY
		centZ = mbpcentZ
		velX = mbpvelX 
		velY = mbpvelY
		velZ = mbpvelZ
	halos.append([centX, centY, centZ, velX, velY, velZ, M, Rvir, Vmax, Vsig, Mgas, Mstar,fMhires])
f.close()
#do i have to do this?
halosA = np.array(halos)
halomax = len(halos)

xD = DMpos[:,0]
yD = DMpos[:,1]
zD = DMpos[:,2]



subhaloN=0

#print 'halos center = ', halosA[subhaloN][0], halosA[subhaloN][1], halosA[subhaloN][2]
#print 'halos vel    = ', halosA[subhaloN][3], halosA[subhaloN][4], halosA[subhaloN][5]
#print 'Rvir         = ', halosA[subhaloN][7]


halodata=[]

alistno=[]
for subhaloN in halo_to_do:
	Nhalostring=str(subhaloN)

	haloX=halosA[subhaloN][0]
	haloY=halosA[subhaloN][1]
	haloZ=halosA[subhaloN][2]
	haloVX = halosA[subhaloN][3]
	haloVY = halosA[subhaloN][4]
	haloVZ = halosA[subhaloN][5]


	print 'haloXYZ', haloX, haloY, haloZ


	xDMrel=xD-haloX
	yDMrel=yD-haloY
	zDMrel=zD-haloZ






	Rvir=halosA[subhaloN][7]
	Mvir=halosA[subhaloN][6]
	print 'Rvir', Rvir
	print 'Mvir', Mvir
	if ((Mvir>5e9) and (Mvir<1e12)):
		inthaloN=int(subhaloN)
		alistno=np.append(alistno,inthaloN)
		fM=halosA[subhaloN][12]
		DcloseHx = abs(xDMrel) < Rvir*3.0
		DcloseHy = abs(yDMrel) < Rvir*3.0
		DcloseHz = abs(zDMrel) < Rvir*3.0
		DcloseHxy = DcloseHx * DcloseHy * DcloseHz
		DMdists = calcdist2(haloX, haloY, haloZ, xD[DcloseHxy], yD[DcloseHxy], zD[DcloseHxy])

	#	print 'DcloseHx', DcloseHx


		if (halosA[subhaloN][12]<0.9999): 
			continue

		DMclosest = DMdists < Rvir*3.0



		#print 'Rvir', Rvir
		#print 'DMdists', DMdists
		#print 'xD[DcloseHxy]', xD[DcloseHxy]

		DMrelposX = xDMrel[:][DcloseHxy][DMclosest]
		DMrelposY = yDMrel[:][DcloseHxy][DMclosest]
		DMrelposZ = zDMrel[:][DcloseHxy][DMclosest]


		DMvelX = DMvel[:,0][DcloseHxy][DMclosest] * np.sqrt(a)
		DMvelY = DMvel[:,1][DcloseHxy][DMclosest] * np.sqrt(a)
		DMvelZ = DMvel[:,2][DcloseHxy][DMclosest] * np.sqrt(a)


		DMrelvelX = DMvelX-haloVX
		DMrelvelY = DMvelY-haloVY
		DMrelvelZ = DMvelZ-haloVZ



		DMveldists=calcdist2(haloVX,haloVY,haloVZ,DMvelX,DMvelY,DMvelZ)

		DMraddists=DMdists[DMclosest]/h

		DMparno=len(DMraddists)



		DMmassclosest = DMmass[DcloseHxy][DMclosest]*1e10/h
		
		DMwid=np.array(DMpid[DcloseHxy][DMclosest])
		print 'sumxD', np.count_nonzero(xD)
		print 'sumid', np.count_nonzero(DMpid)
		print 'sumDMrelposX', np.count_nonzero(DMrelposX)
		print 'sumDMmassclosest', np.sum(DMmassclosest)
		print 'sum_wid', np.count_nonzero(DMwid)




#		plt.plot(xD[DcloseHxy][DMclosest],yD[DcloseHxy][DMclosest],ls='none',marker='o')
#		plt.show()


#		print 'xmax', np.amax(xD[DcloseHxy][DMclosest])
#		print 'ymax', np.amax(yD[DcloseHxy][DMclosest])
#		print 'zmax', np.amax(zD[DcloseHxy][DMclosest])

#		print 'xmin', np.amin(xD[DcloseHxy][DMclosest])
#		print 'ymin', np.amin(yD[DcloseHxy][DMclosest])
#		print 'zmin', np.amin(zD[DcloseHxy][DMclosest])




#		DMpidi=np.array(DMpidi)
#		DMwid=np.array(DMwid)

		

#		counter = Counter(DMpidi)
#		counterw = Counter(DMwid)

#		print 'counter.most_common DMpidi', counter.most_common(3)
#		print 'counter.most_common DMwid', counterw.most_common(3)

		idint0=np.in1d(DMpidi,DMwid)
		print 'DMpidi[idint0]', np.sort(DMpidi[idint0])
		if len(DMpidi[idint0])==0:
			continue
		idneeded=np.in1d(DMpidi,DMwid)
		print 'sumidint0', np.count_nonzero(DMpidi[idint0])
		print 'idneeded', np.sort(DMpidi[idneeded])
		print 'sum_id', np.count_nonzero(DMpidi[idneeded])
		print 'DMwid', np.sort(DMwid)
		print 'sum_wid', np.count_nonzero(DMwid)
		print 'pid', np.sort(DMpidi)
		print 'sum_pid', np.count_nonzero(DMpidi)
		print 'DMmass[-1]', DMmassi[-1]
		print 'sumDMmass', np.sum(DMmassi[idneeded])*1e10/h

	#	plt.plot(xDi[idint0],yDi[idint0],ls='none',marker='o')
	#	plt.show()

	#	plt.plot(xDi[idint0],zDi[idint0],ls='none',marker='o')
	#	plt.show()

		print 'xmax', np.amax(xDi[idneeded])
		print 'ymax', np.amax(yDi[idneeded])
		print 'zmax', np.amax(zDi[idneeded])

		print 'xmin', np.amin(xDi[idneeded])
		print 'ymin', np.amin(yDi[idneeded])
		print 'zmin', np.amin(zDi[idneeded])

		#print 'original'
		#print 'xmax', np.amax(xDi)
		#print 'ymax', np.amax(yDi)
		#print 'zmax', np.amax(zDi)

		#print 'xmin', np.amin(xDi)
		#print 'ymin', np.amin(yDi)
		#print 'zmin', np.amin(zDi)

		#print [x for x, y in collections.Counter(DMpidi).items() if y > 1]

		random = np.vstack((xDi[idneeded]/40000.-shiftbackx,yDi[idneeded]/40000.-shiftbacky,zDi[idneeded]/40000.-shiftbackz)).T
		#asran=binascii.b2a_uu(random)
		#strran.savetxt("ellipsoid.dat")


		filename='zooms_3rv/parlist_halo_'+str(subhaloN)+'.dat'
		infoname='zooms_3rv/info_'+str(subhaloN)+'.txt'
		text_file = open(filename, "w")
		ncount=0
		for x in random:
			for z in x:
				y='%.3f' % float(z)
				text_file.write(y)
				text_file.write(' ')
				if (ncount==2):
					text_file.write('\n')
					ncount=-1
				ncount=ncount+1
		text_file.close()

		info_file = open(infoname,"w")
		info_file.write('Mvir   '+'Rvir   '+'\n')
		info_file.write(str(Mvir)+'  '+str(Rvir)+'\n')
		info_file.close()

alistf = open('iclist.txt', "w")
for number in alistno:
	alistf.write(str(int(number))+'\n')
alistf.close()




