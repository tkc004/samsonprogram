from struct import *
import sys
import os
import matplotlib as mpl
#mpl.use('Agg')
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
Nsnapstring = str(snapno)
halo_to_do = range(1600) #halos to do
#halo_to_do = [0]



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







#finname = 'iclist.txt'
#f = open(finname)
#dars = f.readlines()
#f.close()
#alist = []
#Nlist = []

#Nhalo = 0
#for line in dars:
#        xsd = line.split()
#        Nlist.append(str(xsd[0]))
#        Nhalo = Nhalo + 1



nloop=0
nNlist=[]
Mlist=[]
Rlist=[]
parnolist=[]
pmratiolist=[]
totalenclist=[]

for haloN in halo_to_do:
	strhaloN=str(haloN)
	finname = 'zooms_3rv/ics_halo_'+strhaloN+'.dat'
	infoname = 'zooms_3rv/info_'+strhaloN+'.txt'
	if os.path.isfile(finname):	
		try:	
			f = open(infoname)
			dars = f.readlines()
			f.close()

			for line in dars:
				xsd = line.split()
				Mhalo=str(xsd[0])
				Rhalo=str(xsd[1])

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

				piddtypei = np.dtype(np.int32)

				DMpid_ari.append(   np.fromfile(f, dtype=piddtypei, count = Ndarki))
				Wpid_ari.append(  np.fromfile(f, dtype=piddtypei, count = Nweirdi))
				junk = unpack('i',f.read(4))
				massdi = 'f4'
				massdtypei = np.dtype(massdi)

				DMmass_ari.append( np.zeros(Ndarki))
				DMmass_ari[filecount][0:Ndarki]=Mpartsi[1]
                                Wmass_ari.append( np.zeros(Nweirdi))
				Wmass_ari[filecount][0:Nweirdi]=Mpartsi[2]
				#Wmass_ari.append(  np.fromfile(f, dtype=massdtypei, count = Nweirdi))
				f.close()
				filecount+=1






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




			ai = float(1.0)
			print 'using a = ',ai
			theredshifti = np.absolute(1.0/ai - 1)
			redshiftstringi = "{0:.3f}".format(theredshifti)
			print 'z = ',"{0:.3f}".format(theredshifti)

			LU = 1
			MU = 1
			LU2 = 1
			finname = 'Pep/snap'+ Nsnapstring+ 'RPep.z'+redshiftstringi+'.AHF_halos'
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


			
			fMhires=halos[haloN][12]
			haloX=halos[haloN][0]
			haloY=halos[haloN][1]
			haloZ=halos[haloN][2]
			if (fMhires<0.999999):
				continue
					
			xDMrel=xDi-haloX
			yDMrel=yDi-haloY
			zDMrel=zDi-haloZ
                        xWrel=xWi-haloX
                        yWrel=yWi-haloY
                        zWrel=zWi-haloZ

			xDc=xDMrel<1000
			yDc=yDMrel<1000
			zDc=zDMrel<1000
			xWc=xWrel<1000
			yWc=yWrel<1000
			zWc=zWrel<1000
			critD=xDc*yDc*zDc
			critW=xWc*yWc*zWc
			DMenc=np.sum(DMmassi[critD])
			Wenc=np.sum(Wmassi[critW])
			totalenc=DMenc+Wenc

			parno=len(xDi)
			nNlist=np.append(nNlist,strhaloN)
			Mlist=np.append(Mlist,str(np.log10(float(Mhalo))))
			Rlist=np.append(Rlist,Rhalo)
			parnolist=np.append(parnolist,parno)
			pmratio=parno/float(Mhalo)*np.amin(DMmassi)*1e10
			pmratiolist=np.append(pmratiolist,pmratio)
			totalenclist=np.append(totalenclist,totalenc)
			Mhalo=str(np.log10(float(Mhalo)))
			print 'strhaloN, Mhalo,     Rhalo,   parno,  ell mass/Mvir, enc(1Mpc^3)'
			print strhaloN, Mhalo,     Rhalo,   parno,  pmratio, totalenc
			nloop=nloop+1
		except error:
			print 'StructError'































pmname='zooms_3rv/pmratio.txt'
pmfile=open(pmname,"w")
pmfile.write('halo no.   ' + 'Mvir   ' + 'Rvir   '+'Part no.   '+'ell mass/Mvir'+'            enc(1Mpc)'+'\n')
for ncount in range(nloop):
	pmfile.write(nNlist[ncount]+'    '+Mlist[ncount]+'  '+Rlist[ncount]+'   '+str(parnolist[ncount])+ '    '+str(pmratiolist[ncount])+'      '+str(totalenclist[ncount])+'\n')
pmfile.close()


orderp=pmratiolist.argsort()
snNlist=nNlist[orderp]
sMlist=Mlist[orderp]
sRlist=Rlist[orderp]
sparnolist=parnolist[orderp]
spmratiolist=pmratiolist[orderp]
stotalenclist=totalenclist[orderp]


spmname='zooms_3rv/sort_pmratio.txt'
spmfile=open(spmname,"w")
spmfile.write('halo no.   ' + 'Mvir   ' + 'Rvir   '+'Part no.   '+'ell mass/Mvir'+'        enc(1Mpc)'+'\n')
for ncount in range(nloop):
        spmfile.write(snNlist[ncount]+'    '+sMlist[ncount]+'  '+sRlist[ncount]+'   '+str(sparnolist[ncount])+ '    '+str(spmratiolist[ncount])+'        '+str(stotalenclist)+'\n')
spmfile.close()








