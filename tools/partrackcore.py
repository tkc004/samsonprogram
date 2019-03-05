from general_particle_tracking_function import *
from Sasha_functions import *
from readsnap import readsnap
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys

multifile = 'n'

Nsnap = 300


if (int(Nsnap) < 10):
	Nsnapstring = '00'+str(Nsnap)
elif (int(Nsnap) < 100):
	Nsnapstring = '0'+str(Nsnap)
else:
	Nsnapstring = str(Nsnap)

the_snapdir = './hdf5/'
the_prefix ='snapshot'
if (multifile == 'y'):
	the_snapdir = './snapdir_'+Nsnapstring+'/'
the_suffix = '.hdf5'

#get the first 20 DM particles and see where they came from
the_ptype = 1
gthe_ptype = 0
sthe_ptype = 4


S = readsnap(the_snapdir, Nsnapstring, the_ptype, snapshot_name=the_prefix, extension=the_suffix)
gS = readsnap(the_snapdir, Nsnapstring, gthe_ptype, snapshot_name=the_prefix, extension=the_suffix)
sS = readsnap(the_snapdir, Nsnapstring, sthe_ptype, snapshot_name=the_prefix, extension=the_suffix)

PID_list = S['id'][:]
DMpos = S['p'][:,:]
partX = DMpos[:,0]
partY = DMpos[:,1]
partZ = DMpos[:,2]
gPID_list = gS['id'][:]
Gpos = gS['p'][:,:]
gpartX = Gpos[:,0]
gpartY = Gpos[:,1]
gpartZ = Gpos[:,2]
sPID_list = sS['id'][:]
Spos = sS['p'][:,:]
spartX = Spos[:,0]
spartY = Spos[:,1]
spartZ = Spos[:,2]

header = S['header'][:]
redshift = header[3]
boxsize = header[9]
a = 1/(1+redshift)

haloN=0

Nhalostring=str(haloN)

halostats = find_halo_now(haloN, a)
#print 'halo stats ',halostats
haloX = halostats[2]
haloY = halostats[3]
haloZ = halostats[4]

print 'haloXYZ', haloX, haloY, haloZ
dists = calcdist2(haloX, haloY, haloZ, partX, partY, partZ, boxsize)
gdists = calcdist2(haloX, haloY, haloZ, gpartX, gpartY, gpartZ, boxsize)
sdists = calcdist2(haloX, haloY, haloZ, spartX, spartY, spartZ, boxsize)


close2 = dists<2
close5 = dists<5
close10 = dists<10

gclose2 = gdists<2
gclose5 = gdists<5
gclose10 = gdists<10

sclose2 = sdists<2
sclose5 = sdists<5
sclose10 = sdists<10

idneeded2 = PID_list[close2]
idneeded5 = PID_list[close5]
idneeded10 = PID_list[close10]

gidneeded2 = gPID_list[gclose2]
gidneeded5 = gPID_list[gclose5]
gidneeded10 = gPID_list[gclose10]

sidneeded2 = sPID_list[sclose2]
sidneeded5 = sPID_list[sclose5]
sidneeded10 = sPID_list[sclose10]

totpar2=len(idneeded2)
totpar5=len(idneeded5)
totpar10=len(idneeded10)

gtotpar2=len(gidneeded2)
gtotpar5=len(gidneeded5)
gtotpar10=len(gidneeded10)

stotpar2=len(sidneeded2)
stotpar5=len(sidneeded5)
stotpar10=len(sidneeded10)


totno=len(DMpos[close2])


print 'Total no of DM within 2kpc', len(DMpos[close2])

#print 'dists[close2]', dists[close2]


Nsnapneeded=Nsnap-10
tableadist=[]
gtableadist=[]
stableadist=[]
tabledc=[]
tablegc=[]
tablesc=[]
while (Nsnapneeded<440):
	Nsnapneeded = Nsnapneeded+1

	TTbl2, na, boxsize = match_particles(idneeded2, Nsnapneeded, ptype=the_ptype)
        TTbl5, na, boxsize = match_particles(idneeded5, Nsnapneeded, ptype=the_ptype)
        TTbl10, na, boxsize = match_particles(idneeded10, Nsnapneeded, ptype=the_ptype)

        gTTbl2, na, boxsize = match_particles(gidneeded2, Nsnapneeded, ptype=gthe_ptype)
        gTTbl5, na, boxsize = match_particles(gidneeded5, Nsnapneeded, ptype=gthe_ptype)
        gTTbl10, na, boxsize = match_particles(gidneeded10, Nsnapneeded, ptype=gthe_ptype)

        sTTbl2, na, boxsize = match_particles(sidneeded2, Nsnapneeded, ptype=sthe_ptype)
        sTTbl5, na, boxsize = match_particles(sidneeded5, Nsnapneeded, ptype=sthe_ptype)
        sTTbl10, na, boxsize = match_particles(sidneeded10, Nsnapneeded, ptype=sthe_ptype)




	npartX2 = TTbl2[:,1]
	npartY2 = TTbl2[:,2]
	npartZ2 = TTbl2[:,3]
        npartX5 = TTbl5[:,1]
        npartY5 = TTbl5[:,2]
        npartZ5 = TTbl5[:,3]
        npartX10 = TTbl10[:,1]
        npartY10 = TTbl10[:,2]
        npartZ10 = TTbl10[:,3]


        gnpartX2 = gTTbl2[:,1]
        gnpartY2 = gTTbl2[:,2]
        gnpartZ2 = gTTbl2[:,3]
        gnpartX5 = gTTbl5[:,1]
        gnpartY5 = gTTbl5[:,2]
        gnpartZ5 = gTTbl5[:,3]
        gnpartX10 = gTTbl10[:,1]
        gnpartY10 = gTTbl10[:,2]
        gnpartZ10 = gTTbl10[:,3]


        snpartX2 = sTTbl2[:,1]
        snpartY2 = sTTbl2[:,2]
        snpartZ2 = sTTbl2[:,3]
        snpartX5 = sTTbl5[:,1]
        snpartY5 = sTTbl5[:,2]
        snpartZ5 = sTTbl5[:,3]
        snpartX10 = sTTbl10[:,1]
        snpartY10 = sTTbl10[:,2]
        snpartZ10 = sTTbl10[:,3]




	halostats = find_halo_now(haloN, na)
	#print 'halo stats ',halostats
	nhaloX = halostats[2]
	nhaloY = halostats[3]
	nhaloZ = halostats[4]

	print 'Nsnap, na, haloXYZ', Nsnapneeded, na,  haloX, haloY, haloZ
	ndists2 = calcdist2(nhaloX, nhaloY, nhaloZ, npartX2, npartY2, npartZ2, boxsize)
        ndists5 = calcdist2(nhaloX, nhaloY, nhaloZ, npartX5, npartY5, npartZ5, boxsize)
        ndists10 = calcdist2(nhaloX, nhaloY, nhaloZ, npartX10, npartY10, npartZ10, boxsize)

        gndists2 = calcdist2(nhaloX, nhaloY, nhaloZ, gnpartX2, gnpartY2, gnpartZ2, boxsize)
        gndists5 = calcdist2(nhaloX, nhaloY, nhaloZ, gnpartX5, gnpartY5, gnpartZ5, boxsize)
        gndists10 = calcdist2(nhaloX, nhaloY, nhaloZ, gnpartX10, gnpartY10, gnpartZ10, boxsize)

        sndists2 = calcdist2(nhaloX, nhaloY, nhaloZ, snpartX2, snpartY2, snpartZ2, boxsize)
        sndists5 = calcdist2(nhaloX, nhaloY, nhaloZ, snpartX5, snpartY5, snpartZ5, boxsize)
        sndists10 = calcdist2(nhaloX, nhaloY, nhaloZ, snpartX10, snpartY10, snpartZ10, boxsize)


#	histmass=[]
#	for jn in range(100):
#		 nclose=ndists2<jn
		 
		 
	nclose2 = ndists2<2
	nclose5 = ndists2<5
	nclose10 = ndists2<10
	nclose20 = ndists2<20
	nclose100 = ndists2<100

        gnclose2 = gndists2<2
        gnclose5 = gndists2<5
        gnclose10 = gndists2<10
        gnclose20 = gndists2<20
        gnclose100 = gndists2<100

        snclose2 = sndists2<2
        snclose5 = sndists2<5
        snclose10 = sndists2<10
        snclose20 = sndists2<20
        snclose100 = sndists2<100


	avedist2 = np.sum(ndists2)/totpar2
        avedist5 = np.sum(ndists5)/totpar5
        avedist10 = np.sum(ndists10)/totpar10

        gavedist2 = np.sum(gndists2)/gtotpar2
        gavedist5 = np.sum(gndists5)/gtotpar5
        gavedist10 = np.sum(gndists10)/gtotpar10

        savedist2 = np.sum(sndists2)/stotpar2
        savedist5 = np.sum(sndists5)/stotpar5
        savedist10 = np.sum(sndists10)/stotpar10

	nredshift = 1/na-1


	dc2 = len(npartX2[nclose2])
	dc5 = len(npartX2[nclose5])
	dc10= len(npartX2[nclose10])
	dc20= len(npartX2[nclose20])
	dc100= len(npartX2[nclose100])
        gc2 = len(gnpartX2[gnclose2])
        gc5 = len(gnpartX2[gnclose5])
        gc10= len(gnpartX2[gnclose10])
        gc20= len(gnpartX2[gnclose20])
        gc100= len(gnpartX2[gnclose100])
        sc2 = len(snpartX2[snclose2])
        sc5 = len(snpartX2[snclose5])
        sc10= len(snpartX2[snclose10])
        sc20= len(snpartX2[snclose20])
        sc100= len(snpartX2[snclose100])


	print 'average dists', avedist2
	tableadist.append([nredshift,avedist2,avedist5,avedist10])
	gtableadist.append([nredshift,gavedist2,gavedist5,gavedist10])
	stableadist.append([nredshift,savedist2,savedist5,savedist10])
	tabledc.append([nredshift,dc2,dc5,dc10,dc20,dc100])
	tablegc.append([nredshift,gc2,gc5,gc10,gc20,gc100])
	tablesc.append([nredshift,sc2,sc5,sc10,sc20,sc100])

arradist=np.array(tableadist)
garradist=np.array(gtableadist)
sarradist=np.array(stableadist)
tdc=np.array(tabledc)
tgc=np.array(tablegc)
tsc=np.array(tablesc)
print arradist[:,0]


finname = 'parwi'+Nsnapstring+'_'+Nhalostring+'.eps'
plt.plot(tdc[:,0],tdc[:,1],label='2kpc', color='r')
plt.plot(tdc[:,0],tdc[:,2],label='5kpc', color='g')
plt.plot(tdc[:,0],tdc[:,3],label='10kpc', color='y')
plt.plot(tdc[:,0],tdc[:,2],label='20kpc', color='b')
plt.plot(tdc[:,0],tdc[:,3],label='100kpc', color='k')
plt.xlabel(r'z')
plt.ylabel(r'par no')
plt.legend(loc=3)
plt.savefig(finname)
plt.clf()


finname = 'gparwi'+Nsnapstring+'_'+Nhalostring+'.eps'
plt.plot(tgc[:,0],tgc[:,1],label='2kpc', color='r')
plt.plot(tgc[:,0],tgc[:,2],label='5kpc', color='g')
plt.plot(tgc[:,0],tgc[:,3],label='10kpc', color='y')
plt.plot(tgc[:,0],tgc[:,2],label='20kpc', color='b')
plt.plot(tgc[:,0],tgc[:,3],label='100kpc', color='k')
plt.xlabel(r'z')
plt.ylabel(r'par no')
plt.legend(loc=3)
plt.savefig(finname)
plt.clf()


finname = 'sparwi'+Nsnapstring+'_'+Nhalostring+'.eps'
plt.plot(tsc[:,0],tsc[:,1],label='2kpc', color='r')
plt.plot(tsc[:,0],tsc[:,2],label='5kpc', color='g')
plt.plot(tsc[:,0],tsc[:,3],label='10kpc', color='y')
plt.plot(tsc[:,0],tsc[:,2],label='20kpc', color='b')
plt.plot(tsc[:,0],tsc[:,3],label='100kpc', color='k')
plt.xlabel(r'z')
plt.ylabel(r'par no')
plt.legend(loc=3)
plt.savefig(finname)
plt.clf()



finname = 'ave_dist_from_2kpc'+Nsnapstring+'_'+Nhalostring+'.eps'
plt.plot(arradist[:,0],arradist[:,1], label='2kpc', color='r')
plt.xlabel(r'z')
plt.ylabel(r'Ave dist')
plt.savefig(finname)
plt.clf()

finname = 'gas_ave_dist_from_2kpc'+Nsnapstring+'_'+Nhalostring+'.eps'
plt.plot(garradist[:,0],garradist[:,1], label='2kpc', color='r')
plt.xlabel(r'z')
plt.ylabel(r'Ave dist')
plt.savefig(finname)
plt.clf()

finname = 'star_ave_dist_from_2kpc'+Nsnapstring+'_'+Nhalostring+'.eps'
plt.plot(sarradist[:,0],sarradist[:,1], label='2kpc', color='r')
plt.xlabel(r'z')
plt.ylabel(r'Ave dist')
plt.savefig(finname)
plt.clf()

finname = 'dm_g_s_ave_dist_from_2kpc'+Nsnapstring+'_'+Nhalostring+'.eps'
plt.plot(arradist[:,0],arradist[:,1], label='DM', color='r')
plt.plot(garradist[:,0],garradist[:,1], label='Gas', color='b')
plt.plot(sarradist[:,0],sarradist[:,1], label='Star', color='g')
plt.xlabel(r'z')
plt.ylabel(r'Ave dist')
plt.legend(ncol=3)
plt.savefig(finname)
plt.clf()



finname = 'ave_dist'+Nsnapstring+'_'+Nhalostring+'.eps'
plt.plot(arradist[:,0],arradist[:,1],label='2kpc', color='r')
plt.plot(arradist[:,0],arradist[:,2],label='5kpc', color='g')
plt.plot(arradist[:,0],arradist[:,3],label='10kpc', color='y')
plt.xlabel(r'z')
plt.ylabel(r'Ave dist')
plt.legend(loc=3)
plt.savefig(finname)
plt.clf()


finname = 'gas_ave_dist'+Nsnapstring+'_'+Nhalostring+'.eps'
plt.plot(garradist[:,0],garradist[:,1],label='2kpc', color='r')
plt.plot(garradist[:,0],garradist[:,2],label='5kpc', color='g')
plt.plot(garradist[:,0],garradist[:,3],label='10kpc', color='y')
plt.xlabel(r'z')
plt.ylabel(r'Ave dist')
plt.legend(loc=3)
plt.savefig(finname)
plt.clf()


finname = 'star_ave_dist'+Nsnapstring+'_'+Nhalostring+'.eps'
plt.plot(sarradist[:,0],sarradist[:,1],label='2kpc', color='r')
plt.plot(sarradist[:,0],sarradist[:,2],label='5kpc', color='g')
plt.plot(sarradist[:,0],sarradist[:,3],label='10kpc', color='y')
plt.xlabel(r'z')
plt.ylabel(r'Ave dist')
plt.legend(loc=3)
plt.savefig(finname)
plt.clf()

