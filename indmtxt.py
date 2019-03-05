from readsnap_samson import *
from Sasha_functions import *
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from gadget_lib.cosmo import *
from samson_functions import *
from pylab import *
rcParams['xtick.major.width'] = 2
rcParams['ytick.major.width'] = 2

#dirneed=['383','g383_nag','g0320_383','g0301_383']

#wistr='8_1'
wistr=str(sys.argv[2])
if wistr=='8_1':
	withinkpc=8.1
if wistr=='4_6':
	withinkpc=4.6
galcen=1
runtodo=str(sys.argv[1])
#dirneed=['f383']
snapno=[]
mwithin=[]
haloinfo=cosmichalo(runtodo)
beginno=haloinfo['beginno']
finalno=haloinfo['finalno']
rundir=haloinfo['rundir']
subdir=haloinfo['subdir']
maindir=haloinfo['maindir']
multifile=haloinfo['multifile']
halocolor=haloinfo['halocolor']
halostr=haloinfo['halostr']

if galcen == 1:
	spmname='/home/tkc004/'+maindir+'/'+rundir+'/dmwithin'+wistr+'kpc_gc.txt'
else:
	spmname='/home/tkc004/'+maindir+'/'+rundir+'/dmwithin'+wistr+'kpc.txt'
spmfile=open(spmname,"w")
spmfile.write('snapshot number'+'     '+'time (Gyr)'+'    '+'DM mass within '+wistr+' kpc (Msun)'+'    '+'halo number'+'       '+'virial mass'+'        '+'gas mass within '+wistr+' kpc (Msun)'+'\n')
halosA = read_halo_history(rundir, halonostr=halostr, maindir=maindir)
redlist = halosA['redshift']
xlist = halosA['x']
ylist = halosA['y']
zlist = halosA['z']
halolist =halosA['ID']
mvirlist = halosA['M']
alist = np.array(1./(1.+np.array(redlist)))
xlist = np.array(xlist)*np.array(alist)/0.702
ylist = np.array(ylist)*np.array(alist)/0.702
zlist = np.array(zlist)*np.array(alist)/0.702
print 'redlist', redlist
print 'xlist', xlist
print 'halolist', halolist

if galcen == 1:
	finname = '/home/tkc004/'+maindir+'/'+rundir+'/center/galcen.txt'
	f = open(finname)
	f.readline()
	dars = f.readlines()
	snapD = []
	galXl = []
	galYl = []
	galZl = []
	for line in dars:
		xsd = line.split()
		snapD.append(int(float(xsd[0])))
		galXl.append(float(xsd[2]))
		galYl.append(float(xsd[3]))
		galZl.append(float(xsd[4]))

for Nsnap in range(beginno,finalno):
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
	header = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, header_only=1, h0=1,cosmological=1)
	G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, cosmological=1, h0=1)
	DM = readsnap(the_snapdir, Nsnapstring, 1, snapshot_name=the_prefix, extension=the_suffix, cosmological=1, h0=1)
	S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, cosmological=1, h0=1)
	ascale = header['time']
	thisz = 1./ascale-1.
	print 'thisz', thisz
	hubble = header['hubble']
	print 'hubble', hubble
	if galcen == 1:
		xcen=np.interp(Nsnap,snapD,galXl) #physical distance (kpc) new fire 2.0 uses physical distance for center
		ycen=np.interp(Nsnap,snapD,galYl)
		zcen=np.interp(Nsnap,snapD,galZl)
	else:

		xcen=np.interp(np.log(ascale),np.log(alist),xlist) #physical distance (kpc)
		ycen=np.interp(np.log(ascale),np.log(alist),ylist)
		zcen=np.interp(np.log(ascale),np.log(alist),zlist)
	halono=np.interp(np.log(ascale),np.log(alist),halolist)
	mvir = np.interp(np.log(ascale),np.log(alist),mvirlist)
	print 'cen', xcen, ycen, zcen
	Gpos = G['p'][:,:]
	Gmass = G['m'][:]
	Gx = Gpos[:,0]
	Gy = Gpos[:,1]
	Gz = Gpos[:,2]
	Gr = np.sqrt(np.square(Gx-xcen)+np.square(Gy-ycen)+np.square(Gz-zcen))
	withinr = Gr<withinkpc
	Gwithinm = np.sum(Gmass[withinr]*1e10)

	Spos = S['p'][:,:]
	Smass = S['m'][:]
	Sx = Spos[:,0]
	Sy = Spos[:,1]
	Sz = Spos[:,2]
	Sr = np.sqrt(np.square(Sx-xcen)+np.square(Sy-ycen)+np.square(Sz-zcen))
	Stot = np.sum(Smass)
	withinr = Sr<withinkpc
	Swithinm = np.sum(Smass[withinr]*1e10)


	DMpos = DM['p'][:,:]
	DMmass = DM['m'][:]
	DMx = DMpos[:,0]
	DMy = DMpos[:,1]
	DMz = DMpos[:,2]
	DMr = np.sqrt(np.square(DMx-xcen)+np.square(DMy-ycen)+np.square(DMz-zcen))
	withinr = DMr<withinkpc
	DMtot = np.sum(DMmass)
	DMwithinm = np.sum(DMmass[withinr]*1e10) #in solar mass
	masswithin = Gwithinm+Swithinm+DMwithinm
	print 'Nsnap, DM within, DMtot, Stot', Nsnap, DMwithinm, DMtot, Stot
	print 'Gm within, Sm within, tot m within', Gwithinm, Swithinm, np.log10(masswithin)
	t_now=quick_lookback_time(thisz)
	print Nsnap,t_now,masswithin
	spmfile.write(str(Nsnap)+'       '+str(t_now)+'       '+str(masswithin)+'        '+str(halono)+'     '+ str(mvir)+'     '+str(Gwithinm)+'       '+str(DMwithinm)+'       '+str(Swithinm)+'\n')
spmfile.close()	
