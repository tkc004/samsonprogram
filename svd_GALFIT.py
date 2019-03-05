from readsnap_samson import *
from Sasha_functions import *
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from gadget_lib.cosmo import *
from samson_functions import *
#dirneed=['383','g383_nag','g0320_383','g0301_383']
withinkpc=20.0
jtime=0
#reff=3.0
#reff=1.98
#r12 = reff*4.0/3.0
#dirneed=['f573']
#dirneed=['f383sn472']
#dirneed=['f383sn199','f383sn278','f338sn472']
dirneed=['f573','f553','f476','fm11','f383','f61']
tlist=[[13.4,8.2,13.5],[10.4,5.7,13.5],[10.7,2.6,13.5],[2.7,2.3,6.0],[2.0,2.0,6.5],[2.2,2.2,2.4]]
#tlist=[[13.4,8.3,13.5],[10.4,5.8,13.5],[10.7,2.3,13.5],[2.7,2.0,6.0],[2.0,2.0,6.5],[2.0,2.0,2.4]]
#tlist=[[13.6,13.6,13.6],[13.6,13.6,13.6],[13.6,13.6,13.6],[13.6,13.6,13.6],[13.6,13.6,13.6],[13.6,13.6,13.6]]
#tlist=[[13.3,8.3,13.6],[10.4,4.8,13.6],[11.4,11.4,13.6],[2.6,2.0,6.0],[2.0,2.0,6.5],[2.2,2.2,2.4]]
reffxlist=[[2.9,1.4,3.2],[1.4,1.3,1.4],[2.0,1.3,2.0],[1.7,1.3,4.0],[1.3,1.3,1.5],[1.6,1.6,2.5]]
axisxlist=[[0.35,0.61,0.32],[0.66,0.49,0.72],[0.63,0.52,0.91],[0.58,0.56,0.4],[0.35,0.35,0.75],[0.71,0.71,0.34]]
reffylist=[[1.9,1.4,2.0],[1.3,1.3,1.3],[1.3,1.3,2.2],[2.6,1.4,3.2],[1.3,1.3,1.4],[1.5,1.5,1.5]]
axisylist=[[0.64,0.47,0.68],[0.57,0.42,0.64],[0.54,0.54,0.28],[0.44,0.64,0.61],[0.48,0.48,0.52],[0.63,0.63,0.63]]
reffzlist=[[1.9,1.8,1.7],[1.3,1.3,1.3],[1.4,1.3,1.6],[2.2,1.4,2.3],[1.3,1.3,1.3],[1.3,1.3,1.5]]
axiszlist=[[0.35,0.3,0.37],[0.64,0.66,0.59],[0.53,0.5,0.54],[0.36,0.16,0.75],[0.58,0.58,0.55],[0.73,0.73,0.64]]

#tlist=[[13.6,13.6,13.6]]
#reffxlist=[[1.85,1.85,1.85]]
#axisxlist=[[0.56,0.56,0.56]]
#reffylist=[[2.22,2.22,2.22]]
#axisylist=[[0.45,0.45,0.45]]
#reffzlist=[[1.8,1.8,1.8]]
#axiszlist=[[0.6,0.6,0.6]]



#Nsnap=234

#Nsnap=157
#dirneed=['553','476']
def stellarveldis(runtodo,time,rxeff,ryeff,rzeff):
	gasveldisneed=0
	rlist=np.linspace(0.1,20,num=200)
	mlist=[]
        snaplist, timelist = readtime()
        Nsnap = int(np.interp(time, timelist, snaplist))
	haloinfo=cosmichalo(runtodo)
	rundir=haloinfo['rundir']
	subdir=haloinfo['subdir']
	maindir=haloinfo['maindir']
	multifile=haloinfo['multifile']
	halocolor=haloinfo['halocolor']
	halostr=haloinfo['halostr']
        halosA = read_halo_history(rundir, halonostr=halostr,maindir=maindir)
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
	print 'halolist', halolist


	if (int(Nsnap) < 10):
		Nsnapstring = '00'+str(Nsnap)
	elif (int(Nsnap) < 100):
		Nsnapstring = '0'+str(Nsnap)
	else:
		Nsnapstring = str(Nsnap)

	the_snapdir = '/home/tkc004/'+maindir+'/'+rundir+'/'+subdir
	the_prefix ='snapshot'
	the_suffix = '.hdf5'
	header = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, header_only=1, h0=1,cosmological=1)
	S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, h0=1,cosmological=1)
	ascale = header['time']
	thisz = 1./ascale-1.
	print 'thisz', thisz
	hubble = header['hubble']
	print 'hubble', hubble
	xcen=np.interp(np.log(ascale),np.log(alist),xlist) #physical distance (kpc)
	ycen=np.interp(np.log(ascale),np.log(alist),ylist)
	zcen=np.interp(np.log(ascale),np.log(alist),zlist)
	halono=np.interp(np.log(ascale),np.log(alist),halolist)
	mvir = np.interp(np.log(ascale),np.log(alist),mvirlist)
	print 'cen', xcen, ycen, zcen

	Spos = S['p'][:,:]
	Smass = S['m'][:]
	Svel = S['v'][:,:]
	Sx = Spos[:,0]
	Sy = Spos[:,1]
	Sz = Spos[:,2]
	Svx=Svel[:,0]
	Svy=Svel[:,1]
	Svz=Svel[:,2]
	Sr = np.sqrt(np.square(Sx-xcen)+np.square(Sy-ycen)+np.square(Sz-zcen))
	Srx = np.sqrt(np.square(Sy-ycen)+np.square(Sz-zcen))
        Sry = np.sqrt(np.square(Sx-xcen)+np.square(Sz-zcen))
        Srz = np.sqrt(np.square(Sx-xcen)+np.square(Sy-ycen))
	for i in range(len(rlist)):
		withinrlist = Sr<rlist[i]
		withinrm = np.sum(Smass[withinrlist]*1e10)
		mlist = np.append(mlist,withinrm)
	r12 = np.interp(mlist[-1]*0.5,mlist,rlist)
	reff = r12*3.0/4.0
	withinr = Sr<r12
	withinx = Srx<withinkpc
        withiny = Sry<withinkpc
        withinz = Srz<withinkpc
	withinrx = Sr<rxeff*4./3.
	withinry = Sr<ryeff*4./3.
	withinrz = Sr<rzeff*4./3.
	Swithinm = np.sum(Smass[withinr]*1e10)
        Swithinmx = np.sum(Smass[withinrx]*1e10)
        Swithinmy = np.sum(Smass[withinry]*1e10)
        Swithinmz = np.sum(Smass[withinrz]*1e10)
	print 'rlist', rlist[::10]
	print 'mlist', mlist[::10]
	print 'reff', reff
	xdis=np.std(Svx[withinx])
	ydis=np.std(Svy[withiny])
	zdis=np.std(Svz[withinz])
	print 'x dispersion', xdis
        print 'y dispersion', ydis
        print 'z dispersion', zdis
	print 'average dispersion', (np.std(Svx[withinx])+ np.std(Svy[withinx])+ np.std(Svz[withinx]))/3.0
	print 'm1/2 from sigmax', np.power(np.std(Svx[withinx]),2)*9.3e5*reff
        print 'm1/2 from sigmay', np.power(np.std(Svy[withiny]),2)*9.3e5*reff
        print 'm1/2 from sigmaz', np.power(np.std(Svz[withinz]),2)*9.3e5*reff
	print 'average m1/2', np.power((np.std(Svx[withinx])+ np.std(Svy[withinx])+ np.std(Svz[withinx]))/3.0,2)*9.3e5*reff
	del S, Spos, Smass, Sx, Sy, Sz, Svx, Svy, Svz, Sr, Srx, Sry, Srz
	G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, cosmological=1, h0=1)
	Gpos = G['p'][:,:]
	Gmass = G['m'][:]
	Gx = Gpos[:,0]
	Gy = Gpos[:,1]
	Gz = Gpos[:,2]
	Gr = np.sqrt(np.square(Gx-xcen)+np.square(Gy-ycen)+np.square(Gz-zcen))
	withinr = Gr<r12
        withinrx = Gr<rxeff*4./3.
        withinry = Gr<ryeff*4./3.
        withinrz = Gr<rzeff*4./3.
	Gwithinm = np.sum(Gmass[withinr]*1e10)
        Gwithinmx = np.sum(Gmass[withinrx]*1e10)
        Gwithinmy = np.sum(Gmass[withinry]*1e10)
        Gwithinmz = np.sum(Gmass[withinrz]*1e10)
	if gasveldisneed==1:
		Gv=G['v'][:,:]
		Gvx=Gv[:,0]
		Gvy=Gv[:,1]
		Gvz=Gv[:,2]
		withingx = np.abs(Gx-xcen)<withinkpc
		withingy = np.abs(Gy-ycen)<withinkpc
		withingz = np.abs(Gz-zcen)<withinkpc
		within= withingx*withingy*withingz
		xgdis=np.std(Gvx[within])
		ygdis=np.std(Gvy[within])
		zgdis=np.std(Gvz[within])
		print 'xgdis', xgdis
		print 'ygdis', ygdis
		print 'zgdis', zgdis
		
	del G, Gpos, Gmass, Gx, Gy, Gz, Gr
        DM = readsnap(the_snapdir, Nsnapstring, 1, snapshot_name=the_prefix, extension=the_suffix, cosmological=1, h0=1)
	DMpos = DM['p'][:,:]
	DMmass = DM['m'][:]
	DMx = DMpos[:,0]
	DMy = DMpos[:,1]
	DMz = DMpos[:,2]
	DMr = np.sqrt(np.square(DMx-xcen)+np.square(DMy-ycen)+np.square(DMz-zcen))
	withinr = DMr<r12
        withinrx = DMr<rxeff*4./3.
        withinry = DMr<ryeff*4./3.
        withinrz = DMr<rzeff*4./3.
	DMwithinm = np.sum(DMmass[withinr]*1e10) #in solar mass
        DMwithinmx = np.sum(DMmass[withinrx]*1e10)
        DMwithinmy = np.sum(DMmass[withinry]*1e10)
        DMwithinmz = np.sum(DMmass[withinrz]*1e10)
	masswithin = Gwithinm+Swithinm+DMwithinm
        masswithinx = Gwithinmx+Swithinmx+DMwithinmx
        masswithiny = Gwithinmy+Swithinmy+DMwithinmy
        masswithinz = Gwithinmz+Swithinmz+DMwithinmz
	print 'masswithin', masswithin
	print 'log(masswithin)', np.log10(masswithin)
	print 'expected dispersion', np.sqrt(masswithin/9.3e5/reff)
	return reff, xdis, ydis, zdis, masswithin, masswithinx, masswithiny, masswithinz, Gwithinm

adstr=' '
Mexpstr=' '
Mexpxstr=' '
Mexpystr=' '
Mexpzstr=' '
Mexpavestr=' '
M12str=' '
M12xstr=' '
M12ystr=' '
M12zstr=' '
M12avestr=' '
r12str=' '
G12str=' '
for ncount in range(len(dirneed)):
        runtodo = dirneed[ncount]
        refflist=[]
        xdislist=[]
        ydislist=[]
        zdislist=[]
	minlist=[]
	time=tlist[ncount][jtime]
	reffx=reffxlist[ncount][jtime]*np.sqrt(axisxlist[ncount][jtime])
	reffy=reffylist[ncount][jtime]*np.sqrt(axisylist[ncount][jtime])
	reffz=reffzlist[ncount][jtime]*np.sqrt(axiszlist[ncount][jtime])
	reff, xdis, ydis, zdis, masswithin, masswithinx, masswithiny, masswithinz, Gwithinm=stellarveldis(runtodo,time,reffx,reffy,reffz)
	r12 = reff*4.0/3.0
	avedis=(xdis+ydis+zdis)/3.0
        adround = np.around(avedis, decimals=1)
	dmin = np.around(np.amin([xdis,ydis,zdis]), decimals=1)
        dmax = np.around(np.amax([xdis,ydis,zdis]), decimals=1)
        adstr = adstr+'&      $'+str(adround)+'$& $\,_{'+str(dmin)+'}^{'+str(dmax)+'} $'
	Mexp=np.power(avedis,2)*9.3e5*reff
	Mexpround = np.around(Mexp/1e9, decimals=2)
	Mexpstr = Mexpstr+'&\multicolumn{2}{l}{'+str(Mexpround)+'}'
        Mexpx=np.power(xdis,2)*9.3e5*reffx
        Mexpxround = np.around(Mexpx/1e9, decimals=2)
        Mexpxstr = Mexpxstr+'&\multicolumn{2}{l}{'+str(Mexpxround)+'}'
        Mexpy=np.power(ydis,2)*9.3e5*reffy
        Mexpyround = np.around(Mexpy/1e9, decimals=2)
        Mexpystr = Mexpystr+'&\multicolumn{2}{l}{'+str(Mexpyround)+'}'
        Mexpz=np.power(zdis,2)*9.3e5*reffz
        Mexpzround = np.around(Mexpz/1e9, decimals=2)
        Mexpzstr = Mexpzstr+'&\multicolumn{2}{l}{'+str(Mexpzround)+'}'
        Mexpaveround = np.around((Mexpx+Mexpy+Mexpz)/3e9, decimals=2)
        Mexpavestr = Mexpavestr+'&\multicolumn{2}{l}{'+str(Mexpaveround)+'}'
        M12round = np.around(masswithin/1e9, decimals=2)
        M12str = M12str+'&\multicolumn{2}{l}{'+str(M12round)+'}'
        M12xround = np.around(masswithinx/1e9, decimals=2)
        M12xstr = M12xstr+'&\multicolumn{2}{l}{'+str(M12xround)+'}'
        M12yround = np.around(masswithiny/1e9, decimals=2)
        M12ystr = M12ystr+'&\multicolumn{2}{l}{'+str(M12yround)+'}'
        M12zround = np.around(masswithinz/1e9, decimals=2)
        M12zstr = M12zstr+'&\multicolumn{2}{l}{'+str(M12zround)+'}'
	M12averound = np.around((masswithinx+masswithiny+masswithinz)/3e9, decimals=2)
        M12avestr = M12avestr+'&\multicolumn{2}{l}{'+str(M12averound)+'}'
        G12round = np.around(Gwithinm/1e9, decimals=2)
        G12str = G12str+'&\multicolumn{2}{l}{'+str(G12round)+'}'
        r12round = np.around(r12, decimals=2)
        r12str = r12str+'&\multicolumn{2}{l}{'+str(r12round)+'}'



print 'ave dis', adstr+'\\\\'
print 'Mexp', Mexpstr+'\\\\'
print 'Mexpx', Mexpxstr+'\\\\'
print 'Mexpy', Mexpystr+'\\\\'
print 'Mexpz', Mexpzstr+'\\\\'
print 'Mexpave', Mexpavestr+'\\\\'
print 'M12', M12str+'\\\\'
print 'M12x', M12xstr+'\\\\'
print 'M12y', M12ystr+'\\\\'
print 'M12z', M12zstr+'\\\\'
print 'M12ave', M12avestr+'\\\\'
print 'r12', r12str+'\\\\'
print 'G12', G12str+'\\\\'
