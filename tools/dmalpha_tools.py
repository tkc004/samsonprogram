from readsnap_samson import *
from Sasha_functions import *
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from gadget_lib.cosmo import *
from samson_functions import *
mpl.use('Agg')
from matplotlib import rcParams
rcParams['axes.linewidth'] = 2
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True
#rcParams['text.usetex'] = True
rcParams['axes.unicode_minus']=False

rcParams.update({'figure.autolayout': True})
#dirneed=['fm11']
#tlist=[13.8]
#dirneed=['m09','m10','1297','1146','573','553','476','m11','383']
#dirneed=['f1146','f573','f553','f476','fm11','f383']
#tlist=[[13.4,8.3,13.5],[10.4,5.8,13.5],[10.7,2.3,13.5],[2.7,2.0,6.0],[2.0,2.0,6.5]]
#plotname='Oh15_fire2'

def get_text_positions(text, x_data, y_data, txt_width, txt_height):
    a = zip(y_data, x_data)
    text_positions = list(y_data)
    for index, (y, x) in enumerate(a):
        local_text_positions = [i for i in a if i[0] > (y - txt_height) 
                            and (abs(i[1] - x) < txt_width * 2) and i != (y,x)]
        if local_text_positions:
            sorted_ltp = sorted(local_text_positions)
            if abs(sorted_ltp[0][0] - y) < txt_height: #True == collision
                differ = np.diff(sorted_ltp, axis=0)
                a[index] = (sorted_ltp[-1][0] + txt_height, a[index][1])
                text_positions[index] = sorted_ltp[-1][0] + txt_height*1.01
                for k, (j, m) in enumerate(differ):
                    #j is the vertical distance between words
                    if j > txt_height * 2: #if True then room to fit a word in
                        a[index] = (sorted_ltp[k][0] + txt_height, a[index][1])
                        text_positions[index] = sorted_ltp[k][0] + txt_height
                        break
    return text_positions

def text_plotter(text, x_data, y_data, text_positions, txt_width,txt_height):
    for z,x,y,t in zip(text, x_data, y_data, text_positions):
        plt.annotate(str(z), xy=(x-txt_width/2, t), xytext=(x-txt_width/2,t),\
	size=8, bbox=dict(boxstyle='round,pad=0.1', fc='yellow', alpha=0.5))
        if y != t:
            plt.arrow(x, t,0,y-t, color='k',alpha=0.3, width=txt_width*0.01, 
                head_width=txt_width*0.02, head_length=txt_height*0.5, 
                zorder=0,length_includes_head=True)


def plotMvMsRvalpha(Mlist,alpha12,alpha51,zr=0.,xaxis='mv',labl=[],marker='o',needlabel=1,label='',alpha12only=0,color='y'): 
	print Mlist,alpha12,alpha51,zr,xaxis,marker,needlabel,label,alpha12only
	if len(labl)==0:
		needanno=0
	else:
		needanno=1
	redshiftstr=str(zr)
	if xaxis=='mv':
		xlabel=r'${\rm \log_{10}(M_{\rm h}/M_{\odot})}$'
		xrange=[8,13]
		redshiftstr=str(zr)
		rhocrit=127*(0.27*(1+zr)*(1+zr)*(1+zr)+0.73)
		xLam=-0.73/(0.27*(1+zr)*(1+zr)*(1+zr)+0.73)
		delta=18.0*np.pi*np.pi+82.0*(xLam)-39*(xLam)*(xLam)
		samlogMv=np.arange(8.0,13.0,0.2)
		logc=0.537+(1.025-0.537)*np.exp(-0.718*np.power(zr,1.08))+(0.024*zr-0.097)*(samlogMv-np.log10(0.7)-np.log10(1e12))
		cont=np.power(10,logc)
		Rvirn=np.power(np.power(10,samlogMv)*3.0/4.0/np.pi/rhocrit/delta,1.0/3.0)
		xr12=0.015
		xr51=0.0075
		alphanfw12=-3.0+2.0/(1.0+cont*xr12)
		alphanfw51=-3.0+2.0/(1.0+cont*xr51)
		xr12=0.015
		xr51=0.0075
		alphanfw12=-3.0+2.0/(1.0+cont*xr12)
		alphanfw51=-3.0+2.0/(1.0+cont*xr51)
	if xaxis=='ms':
		xlabel=r'${\rm \log_{10}(M_{*}/M_{\odot})}$'
		xrange=[4,12]
	if xaxis=='ms_mv':
                xlabel=r'${\rm \log_{10}(M_{*}/M_{\rm h})}$'
                xrange=[-4,-1]

	log10Mlist=np.log10(np.array(Mlist))
	if needlabel==1:
		label12 = '1-2%Rv0'
		label51 = '0.5-1%Rv0'
		if not label == '':
			label12=label
	if needlabel==1:
		plt.plot(log10Mlist,alpha12, markersize=10, label=label12, color=color, marker=marker,linestyle='none')
		if alpha12only==0:
			plt.plot(log10Mlist,alpha51, markersize=10, label=label51, color='r', marker=marker,linestyle='none')
	else:
                plt.plot(log10Mlist,alpha12, markersize=10, color=color, marker=marker,linestyle='none')
                if alpha12only==0:
			plt.plot(log10Mlist,alpha51, markersize=10, color='r', marker=marker,linestyle='none')
	if xaxis=='mv':
		plt.plot(samlogMv,alphanfw12, ls='dotted', color=color, lw=3)
		if alpha12only==0:
			plt.plot(samlogMv,alphanfw51, ls='--', color='r')
	if alpha12only==0:
		for i,j,k in zip(log10Mlist,alpha12,alpha51):
		       plt.plot([i,i],[j,k],color = 'g')
	if needanno==1:
		text = labl; x_data = log10Mlist; y_data = alpha12;	
		txt_width = 0.3; txt_height = 0.1;
		text_positions = get_text_positions(text, x_data, y_data, txt_width, txt_height)
		text_plotter(text, x_data, y_data, text_positions, txt_width, txt_height)
#		for label, x, y in zip(labl,log10Mlist,alpha12):
#			plt.annotate(label,xy=(x, y), xytext=(-2, 4),
#			textcoords='offset points', size=6, ha='right', va='bottom',
#			bbox=dict(boxstyle='round,pad=0.1', fc='yellow', alpha=0.5),
#			arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
	plt.legend(loc=1, ncol=3, shadow=True, numpoints=1)
        plt.xlim(xrange)
        plt.ylim([-2,0.5])
#        finname = 'figures/'+plotname+'.pdf'
        plt.xlabel(xlabel, fontsize=24)
        plt.ylabel(r'$\alpha$', fontsize=24)
        plt.text(8.5, 0.1, r'$z=$'+redshiftstr, fontsize=30)
        plt.tick_params(axis='both', labelsize=16)
#        plt.savefig(finname)
#        plt.clf()

                                                                  

def plotvsTHINGS():
 	txtdir='/home/tkc004/python/'
        filename=txtdir+'littlethingspts.txt'
        ff=open(filename,'r')
        dars = ff.readlines()
        ff.close()
        log10Ms = []
        alpha_pts = []

        for line in dars:
                xsd = line.split()
                log10Ms.append(float(xsd[0]))
                alpha_pts.append(float(xsd[1]))

        filename=txtdir+'littlethingsup.txt'
        ff=open(filename,'r')
        dars = ff.readlines()
        ff.close()
        log10Ms_up = []
        alpha_up = []

        for line in dars:
                xsd = line.split()
                log10Ms_up.append(float(xsd[0]))
                alpha_up.append(float(xsd[1]))

        filename=txtdir+'littlethingslow.txt'
        ff=open(filename,'r')
        dars = ff.readlines()
        ff.close()
        log10Ms_low = []
        alpha_low = []

        for line in dars:
                xsd = line.split()
                log10Ms_low.append(float(xsd[0]))
                alpha_low.append(float(xsd[1]))
        filename=txtdir+'thingspts.txt'
        ff=open(filename,'r')
        dars = ff.readlines()
        ff.close()
        log10Mt = []
        alpha_ptt = []

        for line in dars:
                xsd = line.split()
                log10Mt.append(float(xsd[0]))
                alpha_ptt.append(float(xsd[1]))

        filename=txtdir+'thingsptsup.txt'
        ff=open(filename,'r')
        dars = ff.readlines()
        ff.close()
        log10Mt_up = []
        alpha_upt = []

        for line in dars:
                xsd = line.split()
                log10Mt_up.append(float(xsd[0]))
                alpha_upt.append(float(xsd[1]))

        filename=txtdir+'thingsptslow.txt'
        ff=open(filename,'r')
        dars = ff.readlines()
        ff.close()
        log10Mt_low = []
        alpha_lowt = []
        for line in dars:
                xsd = line.split()
                log10Mt_low.append(float(xsd[0]))
                alpha_lowt.append(float(xsd[1]))


        yed=np.absolute(np.array(alpha_pts)-np.array(alpha_low))
        yeu=np.absolute(np.array(alpha_up)-np.array(alpha_pts))

        ytd=np.absolute(np.array(alpha_ptt)-np.array(alpha_lowt))
        ytu=np.absolute(np.array(alpha_upt)-np.array(alpha_ptt))

        plt.errorbar(log10Mt, alpha_ptt, yerr=[ytd,ytu], markersize=16, label='THINGS', color='r', marker='s',linestyle='none',lw=0.5,capthick=1, fillstyle='none',mew=1.5)
        plt.errorbar(log10Ms, alpha_pts, yerr=[yed,yeu], markersize=16, label='LITTLE THINGS', color='k', marker='<',linestyle='none',lw=0.5,capthick=1, fillstyle='none',mew=1.5)



def plotalpha37(Mstar37,alpha37,marker='o',needlabel=1):
	if needlabel==1:
		plt.plot(Mstar37,alpha37, markersize=16, label='simulations', color='b',marker=marker,linestyle='none')
	else:
		plt.plot(Mstar37,alpha37, markersize=16, color='b',marker=marker,linestyle='none')
        plt.legend(loc=4, ncol=1, shadow=True, numpoints=1)
        plt.xlim([4, 12])
        plt.ylim([-2,1])
        plt.axhline(y=-1,linewidth=1, color='r',ls='--')
        plt.xlabel(r'${\rm \log_{10}(M_{*}/M_{\odot})}$', fontsize=24)
        plt.ylabel(r'$\alpha$', fontsize=24)
        plt.tick_params(axis='both', labelsize=20)




def findalpha(DMx,DMy,DMz,DMmass,cen=[0.0,0.0,0.0],Rbound=200,fitrange=[1.0,2.0]):
        DMr = np.sqrt(np.square(DMx-cen[0])+np.square(DMy-cen[1])+np.square(DMz-cen[2]))
        closest=DMr<Rbound
        logDMraddists=np.log10(DMr[closest])
        DMmassclosest=DMmass[closest]
        hist, bins=np.histogram(logDMraddists, bins=40,weights=DMmassclosest)
        radbins=pow(10,bins)
        vol=radbins*radbins*radbins*np.pi*4./3.
        difvol=diffarray(vol)
        newhist=hist/difvol
        newbins=middlepointarray(radbins)
        minfitregion=newbins>fitrange[0]
        maxfitregion=newbins<fitrange[1]
        truefitreg=minfitregion*maxfitregion
        xdata=np.log(newbins[truefitreg])
        ydata=np.log(newhist[truefitreg])
        x0=np.array([0.0,0.0])
        def func(x, a, b):
                return a+b*x

        try:
                outfit=optimize.curve_fit(func, xdata, ydata, x0)
                Rviralpha=outfit[0][1]
        except (RuntimeError, TypeError, NameError):
                print "centralalpha_1_2Rvir not exist"
                Rviralpha=-3
        return Rviralpha, xdata, ydata

def dmalpha(runtodo,time,Rrange=[0.3,0.7],timesRv=0,usesnapno=0):
        haloinfo=cosmichalo(runtodo)
        beginno=haloinfo['beginno']
        finalno=haloinfo['finalno']
        rundir=haloinfo['rundir']
        subdir=haloinfo['subdir']
        multifile=haloinfo['multifile']
        halocolor=haloinfo['halocolor']
        halostr=haloinfo['halostr']
	firever=haloinfo['firever']
	maindir=haloinfo['maindir']
	highres=haloinfo['highres']
	usepep=haloinfo['usepep']
	suffixadd=haloinfo['suffixadd']
	snumadd=haloinfo['snumadd']
	Sheaform=haloinfo['Sheaform']
	rlist=np.linspace(0.1,20,num=200)
	mlist=[]
	if usesnapno==0:
		timedata = readtime(firever=firever)
		snaplist = timedata['snaplist']
		timelist = timedata['timelist']
		Nsnap = int(np.interp(time, timelist, snaplist))
	else:
		Nsnap = time
	if (int(Nsnap) < 10):
		Nsnapstring = '00'+str(Nsnap)
	elif (int(Nsnap) < 100):
		Nsnapstring = '0'+str(Nsnap)
	else:
		Nsnapstring = str(Nsnap)

	the_snapdir = '/home/tkc004/'+maindir+'/'+rundir+'/'+subdir
	print 'snapdir', the_snapdir
	print 'Nsnap', Nsnap
	the_prefix ='snapshot'
	the_suffix = '.hdf5'
	header = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, header_only=1, h0=1,cosmological=1)
	ascale = header['time']

	thisz = 1./ascale-1.
	print 'thisz', thisz
	if thisz<1e-7:
		ascale=1.
	hubble = header['hubble']
	print 'hubble', hubble
	if usepep==0:
		halosA = read_halo_history(rundir, snumadd=snumadd,suffixadd=suffixadd,halonostr=halostr,hubble=hubble,comoving=0,maindir=maindir) #to output comoving distance comoving = 1
		redlist = halosA['redshift']
		xlist = halosA['x'] #center of the halo physical kpc without h
		ylist = halosA['y']
		zlist = halosA['z']
		halolist =halosA['ID']
		mvirlist = halosA['M'] #Msun
		mslist = halosA['Ms'] # stellar mass
		rvirlist = halosA['R'] #physical kpc without h
		alist=1.0/(1.0+redlist)
		xcen=np.interp(np.log(ascale),np.log(alist),xlist) #physical distance (kpc)
		ycen=np.interp(np.log(ascale),np.log(alist),ylist)
		zcen=np.interp(np.log(ascale),np.log(alist),zlist)
		halono=np.interp(np.log(ascale),np.log(alist),halolist)
		mvir = np.interp(np.log(ascale),np.log(alist),mvirlist)
		Rvir = np.interp(np.log(ascale),np.log(alist),rvirlist)
		ms = np.interp(np.log(ascale),np.log(alist),mslist)

	else:
                halosA = read_halo_history_pep(rundir, Nsnap, Sheaform=Sheaform, halonostr=halostr,hubble=hubble,comoving=1,maindir=maindir,firever=firever)
                xcen = halosA['x'] #center of the halo physical kpc without h
                ycen = halosA['y']
                zcen = halosA['z']
                halono =halosA['ID']
                mvir = halosA['M'] #Msun
                ms = halosA['Ms'] # stellar mass
                Rvir = halosA['R'] #physical kpc without h
        DM = readsnap(the_snapdir, Nsnapstring, 1, snapshot_name=the_prefix, extension=the_suffix, cosmological=1, h0=1)
	DMpos = DM['p'][:,:]
	DMmass = DM['m'][:]
	DMx = DMpos[:,0]
	DMy = DMpos[:,1]
	DMz = DMpos[:,2]
#	Rviralpha, xdata, ydata = findalpha(DMx,DMy,DMz,DMmass,cen=[xcen,ycen,zcen],Rbound=Rvir,fitrange=[0.01*Rvir,0.02*Rvir])
	fitrange=np.array(Rrange)
	if timesRv==1:
		fitrange*=Rvir
        Rvir37alpha, xdata, ydata = findalpha(DMx,DMy,DMz,DMmass,cen=[xcen,ycen,zcen],Rbound=Rvir,fitrange=fitrange)
        return Rvir37alpha, xdata, ydata, mvir, ms
