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
galcen=1
#runtodo=str(sys.argv[1])
Nsnap=546
dirneed=['f573','f553','f476']
snapno=[]
mwithin=[]
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
	labelname=haloinfo['labelname']

	if galcen == 1:
		spmname='/home/tkc004/samsonprogram/data/starinfo_'+labelname+'_'+str(Nsnap)+'_gc.txt'
	else:
		spmname='/home/tkc004/samsonprogram/data/starinfo_'+labelname+'_'+str(Nsnap)+'.txt'
	spmfile=open(spmname,"w")
	spmfile.write('[0] star forming time [Gyr]'+'    '+'[1] stellar mass [Msun]'+'      '+'[2] x [kpc]'+'    '+'[3] y [kpc]'+'       '+'[4] z [kpc]'+'\n')
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
	S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, cosmological=1, h0=1)
	ascale = header['time']
	thisz = 1./ascale-1.
	print 'thisz', thisz
	t_now=quick_lookback_time(thisz)

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

	Spos = S['p'][:,:]
	Sm = S['m'][:]*1e10
	Sage = S['age'][:]
	Sx = Spos[:,0]-xcen
	Sy = Spos[:,1]-ycen
	Sz = Spos[:,2]-zcen
	Sr = np.sqrt(np.square(Sx)+np.square(Sy)+np.square(Sz))
	withinr = Sr<10. #kpc
	Sxcut = Sx[withinr]; Sycut = Sy[withinr]; Szcut = Sz[withinr];
	Smcut = Sm[withinr];
	Sredshift = 1./Sage[withinr]-1.
	SagecutGyr = quick_lookback_time(Sredshift)
	for i in range(len(Smcut)):
		#'[0] stellar age [Gyr]'+'    '+'[1] stellar mass [Msun]'+'      '+'[2] x [kpc]'+'    '+'[3] y [kpc]'+'       '+'[4] z [kpc]'
		spmfile.write(str(SagecutGyr[i])+'       '+str(Smcut[i])+'       '+str(Sxcut[i])+'        '+str(Sycut[i])+'     '+ str(Szcut[i])+'\n')
	spmfile.close()	
	print 'filename', spmname
