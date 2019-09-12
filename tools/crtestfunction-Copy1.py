from samson_const import *
import matplotlib as mpl
from readsnap_cr import readsnapcr
import Sasha_functions as SF
import graphics_library as GL
import gas_temperature as GT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from samson_functions import *
from matplotlib import rcParams
from pylab import *




def weighted_quantile(values, quantiles, sample_weight=None, values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of initial array
    :param old_style: if True, will correct output to be consistent with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)/100.
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), 'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)


def dirout(runtodo,startno,Nsnap,snapsep, Rfrac,\
        physicalR=-1.0, shiftz=1): #if physicalR >0, physicalR is the radius to consider gamma rays 
        pretime=0
        presnap=startno
        snaplist=[]
        enclist=[]
        englist=[]
        enllist=[]
        sml=[]
        diskml=[]
        bulgeml=[]
        nsml=[]
        timel=[]
        Lgcalist=[]
        Lgamma=[]
        Lgamma_sfr=[]
        Lgamma_sfr_noIa=[]
        for i in range(startno,Nsnap, snapsep):
                info=outdirname(runtodo, i)
                rundir=info['rundir']
                maindir=info['maindir']
                halostr=info['halostr']
                runtitle=info['runtitle']
                slabel=info['slabel']
                snlabel=info['snlabel']
                dclabel=info['dclabel']
                resolabel=info['resolabel']
                the_snapdir=info['the_snapdir']
                Nsnapstring=info['Nsnapstring']
                havecr=info['havecr']
                Fcal=info['Fcal']
                iavesfr=info['iavesfr']
                timestep=info['timestep']
                cosmo=info['cosmo']
                color=info['color']
                withinRv=info['withinRv']
                usepep=info['usepep']
                beginno=info['beginno']
                finalno=info['finalno']
                snapsep=info['snapsep']
                firever=info['firever']         
                initsnap=info['initsnap']       
                haveB=info['haveB']
                M1speed=info['M1speed']
                Rvirguess=info['Rvir']
                the_prefix=info['the_prefix']
                the_suffix=info['the_suffix']
                correctIa=info['correctIa']
                halostr=info['halostr']
                galcen = info['galcen']
                if cosmo==0:
                        inittime=initsnap
                else:
                        inittime=0
                        print 'set initial time'        
                print 'the_snapdir', the_snapdir
                print 'Nsnapstring', Nsnapstring
                print 'havecr', havecr
                print 'withinRv', withinRv
                print 'cosmo', cosmo
                try:
                        if cosmo==1:
                                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                                S = readsnapcr(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                                if correctIa==1:
                                        Disk = readsnapcr(the_snapdir, Nsnapstring, 2, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                                        Bulge = readsnapcr(the_snapdir, Nsnapstring, 3, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                        else:
                                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                                S = readsnapcr(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                                if correctIa==1:
                                        Disk = readsnapcr(the_snapdir, Nsnapstring, 2, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                                        Bulge = readsnapcr(the_snapdir, Nsnapstring, 3, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        if havecr>4:
                                cregyl = G['cregyl']*1e10*solar_mass_in_g*km_in_cm*km_in_cm #in erg (originally in code unit: 1e10Msun*(km/s)^2)
                                cregyg = G['cregyg']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                        if havecr > 0:
                                cregy_codeunit = G['cregy']
                                cregy  = cregy_codeunit*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                        Grho = G['rho']
                        Neb = G['ne']
                except KeyError:
                        print 'Keyerror'
                        break
                try:
                        header=S['header']
                        timeneed=header[2]
                        print 'timeneed', timeneed
                        Smi=S['m']
                        Sage=S['age']
                        if withinRv ==1 and cosmo==1:
                                Sp = S['p']
                                Sx = Sp[:,0]
                                Sy = Sp[:,1]
                                Sz = Sp[:,2]
                                header=readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1,header_only=1)
                                h0 = header['hubble']
                                atime = header['time']
                                print 'halostr', halostr
                                if usepep==1:
                                        halosA = SF.read_halo_history_pep(rundir, finalno, snapsep=snapsep, beginno=beginno, singlesnap=0, firever=firever,halonostr=halostr, hubble=h0, comoving=0, maindir=maindir)
                                        afactor=atime
                                else:
                                        halosA = SF.read_halo_history(rundir, maindir=maindir, halonostr=halostr, hubble=h0, comoving=0)
                                        afactor=1.0
                                redlist = halosA['redshift']
                                haloid = halosA['ID']
                                a_scale = 1.0/(1.0+redlist)
                                xcenl = halosA['x']*afactor
                                ycenl = halosA['y']*afactor
                                zcenl = halosA['z']*afactor
                                Rvirl = halosA['R']*afactor
                                if galcen == 1:
                                        finname = '/home/tkc004/'+maindir+'/'+rundir+'/center/galcen.txt'
                                        f = open(finname)
                                        f.readline()
                                        dars = f.readlines()
                                        f.close()
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
                                        print 'xcenl', xcenl
                                        print 'galXl', galXl
                                        xcen=np.interp(i,snapD,galXl) #physical distance (kpc) (new program for FIRE 2.0 uses physical distance)
                                        ycen=np.interp(i,snapD,galYl)
                                        zcen=np.interp(i,snapD,galZl)
                                else:
                                        xcen = np.interp(atime,a_scale,xcenl)
                                        ycen = np.interp(atime,a_scale,ycenl)
                                        zcen = np.interp(atime,a_scale,zcenl)
                                Rvir = np.interp(atime,a_scale,Rvirl)
                                Sxrel = Sx-xcen
                                Syrel = Sy-ycen
                                Szrel = Sz-zcen
                                Sr = np.sqrt(Sxrel*Sxrel+Syrel*Syrel+Szrel*Szrel)
                                cutrvs = Sr<Rvir*0.1
                                Smi = Smi[cutrvs]
                                Sage = Sage[cutrvs]
                        Sm = np.sum(Smi)*1e10 #in solar mass
                        tcut=Sage>pretime
                        Nsm = np.sum(Smi[tcut])*1e10
                        if correctIa==1:
                                diskm = np.sum(Disk['m'])*1e10
                                bulgem = np.sum(Bulge['m'])*1e10
                except KeyError:
                        print 'key error'
                        Sm = 0.
                        Nsm = 0.
                        timeneed=0
                        if correctIa==1:
                                diskm = 0.
                                bulgem = 0.
                if withinRv ==1:
                        Gp = G['p']
                        Gx = Gp[:,0]
                        Gy = Gp[:,1]
                        Gz = Gp[:,2]
                        header=readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1,header_only=1)
                        h0 = header['hubble']
                        atime = header['time']
                        if cosmo==1:
                                if usepep==1:
                                        halosA = SF.read_halo_history_pep(rundir, finalno, beginno=beginno, snapsep=snapsep, singlesnap=0, firever=firever,halonostr=halostr, hubble=h0, comoving=0, maindir=maindir)
                                        afactor=atime
                                else:
                                        halosA = SF.read_halo_history(rundir, maindir=maindir, halonostr=halostr, hubble=h0, comoving=0)
                                        afactor=1.0
                                redlist = halosA['redshift']
                                haloid = halosA['ID']
                                a_scale = 1.0/(1.0+redlist)
                                xcenl = halosA['x']*afactor
                                ycenl = halosA['y']*afactor
                                zcenl = halosA['z']*afactor
                                Rvirl = halosA['R']*afactor
                                if galcen == 1:
                                        finname = '/home/tkc004/'+maindir+'/'+rundir+'/center/galcen.txt'
                                        f = open(finname)
                                        f.readline()
                                        dars = f.readlines()
                                        f.close()
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
                                        xcen=np.interp(i,snapD,galXl) #physical distance (kpc) (new program for FIRE 2.0 uses physical distance)
                                        ycen=np.interp(i,snapD,galYl)
                                        zcen=np.interp(i,snapD,galZl)
                                else:
                                        xcen = np.interp(atime,a_scale,xcenl)
                                        ycen = np.interp(atime,a_scale,ycenl)
                                        zcen = np.interp(atime,a_scale,zcenl)
                                Rvir = np.interp(atime,a_scale,Rvirl)
                        else:
                                Gm = G['m']
                                xcen=ycen=zcen=0;
                                if shiftz==1:
                                        datasup=1
                                        Sm = S['m']
                                        Sp = S['p']
                                        Sx = Sp[:,0]
                                        Sy = Sp[:,1]
                                        Sz = Sp[:,2]
                                        xcen, ycen, zcen = findcennew(runtodo,Nsnap,withinr=5.,datasup=datasup,Gx=Sx,Gy=Sy,Gz=Sz,Gm=Sm)
                                Rvir=Rvirguess
                        Gxrel = Gx-xcen
                        Gyrel = Gy-ycen
                        Gzrel = Gz-zcen
                        Gr = np.sqrt(Gxrel*Gxrel+Gyrel*Gyrel+Gzrel*Gzrel)
                        if physicalR>0.0:
                                cutrv = Gr<physicalR
                        else:
                                cutrv = Gr<Rvir*Rfrac
                        Grho = Grho[cutrv]
                        Neb = Neb[cutrv]
                        if havecr>0:
                                cregy = cregy[cutrv]
                                cregy_codeunit = cregy_codeunit[cutrv]
                        if havecr>4:
                                cregyl = cregyl[cutrv]
                                cregyg = cregyg[cutrv]
                if cosmo==1:
                        readtimelist=readtime(firever=2)
                        snap2list=readtimelist['snaplist']
                        time2list=readtimelist['timelist']
                        a2list=readtimelist['alist']
                        tnow = np.interp(timeneed,a2list,time2list)*1e9
                        pret = np.interp(pretime,a2list,time2list)*1e9
                if havecr>4:
                        cregygt=np.sum(cregyg)
                        cregylt=np.sum(cregyl)
                if havecr>0:
                        cregyt =np.sum(cregy)
                        Lout = outLgamma_nism(Grho,Neb,cregy_codeunit)
                        Lgcal = np.sum(Lout['Lgamma'])
                        print 'Lgcal', Lgcal
                snaplist.append(float(i))
                if cosmo==1:
                        timel.append(tnow)
                else:
                        timel.append(float(i)*0.98*1e6)
                if havecr>0:
                        enclist.append(cregyt)
                if havecr>4:
                        englist.append(cregygt)
                        enllist.append(cregylt)
                if havecr>0:
                        Lgcalist.append(Lgcal)
                sml.append(Sm)
                if correctIa==1:
                        diskml.append(diskm)
                        bulgeml.append(bulgem)
                nsml.append(Nsm)
                pretime=timeneed
                del G, S
        sml=np.array(sml)
        if correctIa==1:
                diskml=np.array(diskml)
                bulgeml=np.array(bulgeml)
        nsml=np.array(nsml)
        if havecr>0:
                enclist=np.array(enclist)
        if havecr>4:
                englist=np.array(englist)
                enllist=np.array(enllist)
        snaplist=np.array(snaplist)
        if havecr>0:
                Lgcalist=np.array(Lgcalist)
        timel=np.array(timel) #in yr
        #above is the coefficient for Salpeter only; for Kroupa, the coefficient is 50% larger:
        avesfrl=(nsml[1:])/(timel[1:]-timel[:-1]) #in Msun/yr
        print 'nsml',nsml
        print 'avesfrl', avesfrl
        Lsfr = Kroupa_Lsf*avesfrl*solar_mass_in_g/yr_in_sec*cspeed_in_cm_s*cspeed_in_cm_s #times 1.5? 6.2e-4 for Kroupa 3.8e-4 for Salpeter
        if correctIa==1:
                IaeffSFR=5.3e-8/3.0e-4*(sml[1:]+diskml[1:]+bulgeml[1:])/0.03753/1e9
                print 'IaeffSFR', IaeffSFR
                LsfrnoIa = Kroupa_Lsf*(avesfrl+IaeffSFR)*solar_mass_in_g/yr_in_sec*cspeed_in_cm_s*cspeed_in_cm_s
        print 'Lsfr', Lsfr
        print 'timel', timel
        if havecr>4:
                Lgamma = (enllist[1:]-enllist[:-1])/((timel[1:]-timel[:-1])*yr_in_sec)/(hadronicdecayrate+coulombdecayrate)*hadronicdecayrate*betapi/nopi_per_gamma
                Lgamma_sfr = np.absolute(Lgamma/Lsfr)
        if havecr>0:
                Lgcal_sfr = Lgcalist[1:]/Lsfr
        else:
                Lgcal_sfr = 0.0
        if correctIa==1:
                Lgamma_sfr_noIa = Lgamma/LsfrnoIa
        return {'snaplist':snaplist, 'timel':timel, 'avesfrl':avesfrl, 'sml':sml, 'Lgamma':np.absolute(Lgamma),\
 'Lgamma_sfr':Lgamma_sfr, 'Lgcal_sfr':Lgcal_sfr, 'Lgamma_sfr_noIa':Lgamma_sfr_noIa,\
'Lgcal':np.absolute(Lgcalist),'Lsfr':Lsfr}


def gassurout(runtodo,Nsnap,projection,shiftz=1): 
        withinr=10.0
        nogrid=100
        maxlength=10.0
        info=outdirname(runtodo, Nsnap)
        rundir=info['rundir']
        maindir=info['maindir']
        halostr=info['halostr']
        runtitle=info['runtitle']
        slabel=info['slabel']
        snlabel=info['snlabel']
        dclabel=info['dclabel']
        resolabel=info['resolabel']
        the_snapdir=info['the_snapdir']
        Nsnapstring=info['Nsnapstring']
        havecr=info['havecr']
        Fcal=info['Fcal']
        iavesfr=info['iavesfr']
        timestep=info['timestep']
        color=info['color']
        cosmo=info['cosmo']
        usepep=info['usepep']
        haveB=info['haveB']
        Rvir=info['Rvir']
        beginno = info['beginno']
        finalno = info['finalno']
        snapsep = info['snapsep']
        firever = info['firever']
        the_prefix=info['the_prefix']
        the_suffix=info['the_suffix']
        if cosmo==1:
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                header=readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1,header_only=1)
                h0 = header['hubble']
                atime = header['time']
                print 'snapsep', snapsep
                print 'Nsnap', Nsnap
                print 'surout eginno', beginno
                if usepep==1:
                        halosA = SF.read_halo_history_pep(rundir, Nsnap, snapsep=snapsep,\
                         singlesnap=1, beginno=beginno, firever=firever,halonostr=halostr, hubble=h0, comoving=0, maindir=maindir)
                        afactor=atime
                else:
                        halosA = SF.read_halo_history(rundir, maindir=maindir, halonostr=halostr, hubble=h0, comoving=0)
                        afactor=1.0
                redlist = halosA['redshift']
                haloid = halosA['ID']
                a_scale = 1.0/(1.0+redlist)
                xcenl = halosA['x']*afactor
                print 'xcenl', xcenl
                ycenl = halosA['y']*afactor
                zcenl = halosA['z']*afactor
                Rvirl = halosA['R']*afactor
                xcen = np.interp(atime,a_scale,xcenl)
                ycen = np.interp(atime,a_scale,ycenl)
                zcen = np.interp(atime,a_scale,zcenl)
                Rvir = np.interp(atime,a_scale,Rvirl)
                print 'xcen', xcen
        else:
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                Gpos = G['p']
                Gx = Gpos[:,0]
                Gy = Gpos[:,1]
                Gz = Gpos[:,2]
                Gm = G['m']
                xcen=ycen=zcen=0;
                if shiftz==1:
                        datasup=1
                        xcen, ycen, zcen = findcennew(runtodo,Nsnap,withinr=5.,datasup=1,Gx=Gx,Gy=Gy,Gz=Gz,Gm=Gm)
        Gpos = G['p']
        Gx = Gpos[:,0]-xcen
        Gy = Gpos[:,1]-ycen
        Gz = Gpos[:,2]-zcen
        Gr = np.sqrt(Gx*Gx+Gy*Gy+Gz*Gz)
        Gcut = Gr<0.5*Rvir
        Gx = Gx[Gcut]
        Gy = Gy[Gcut]
        Gz = Gz[Gcut]
        print 'np.amax(Gx),np.amin(Gx)', np.amax(Gx),np.amin(Gx)
        print 'np.amax(Gy),np.amin(Gy)', np.amax(Gy),np.amin(Gy)
        print 'np.amax(Gz),np.amin(Gz)', np.amax(Gz),np.amin(Gz)
        Gm =np.array(G['m'])
        Gm = Gm[Gcut]
        
        #Gr = np.sqrt(Gx*Gx+Gy*Gy)
        if projection=='y':
                Gr = np.sqrt(Gx*Gx+Gz*Gz)
        if projection=='x':
                Gr = np.sqrt(Gy*Gy+Gz*Gz)
        if projection=='z':
                Gr = np.sqrt(Gy*Gy+Gx*Gx)
        dr = withinr/nogrid
        gassurdenl=[]
        rneed=[]
        for i in range(nogrid):
                rcut = Gr<dr*(i+1)
                Gmwithin = np.sum(Gm[rcut])*1e10*Msun_in_g
                gassurden = Gmwithin/(np.pi*dr*(i+1)*dr*(i+1)*kpc_in_cm*kpc_in_cm)
                gassurdenl = np.append(gassurdenl, gassurden)
                rneed=np.append(rneed,dr*(i+1))
        print 'rneed', rneed
        print 'gassurdenl', gassurdenl
        return {'rneed':rneed, 'gassurdenl':gassurdenl}

def gassurlout(runtodo,startno,Nsnap,snapsep, rad, projection,shiftz=1):
        snapl=[]
        gassurl=[]
        for i in range(startno,Nsnap,snapsep):
                outdata = gassurout(runtodo,i,projection,shiftz=shiftz)
                rneed=outdata['rneed']
                gassurdenl=outdata['gassurdenl']
                gassur=np.interp(rad,rneed,gassurdenl)
                print 'gassur in gassurlout', gassur
                snapl=np.append(snapl,i)
                gassurl = np.append(gassurl,gassur)
        snapl=np.array(snapl)
        gassurl=np.array(gassurl)
        return {'snapl': snapl, 'gassurl': gassurl}



def outcrdata(runtodo, Nsnap, startno,snapsep,shiftz=1):
        withinr=10.0
        maxlength=0.5
        nbin=10
        Rfrac=0.1
        nooftimes=0
        Smnew=0.
        totdiskgm=0.
        diskcoldgm=0.
        Lgamma=0.
        z12gas=0.
        r12cr=0.
        z12dengas=0.
        gassurden=0.
        crden=0.
        for ino in range(startno,Nsnap,snapsep):
                info=outdirname(runtodo, ino)
                haveB=info['haveB']
                rundir=info['rundir']
                runtitle=info['runtitle']
                slabel=info['slabel']
                snlabel=info['snlabel']
                dclabel=info['dclabel']
                resolabel=info['resolabel']
                the_prefix=info['the_prefix']
                the_suffix=info['the_suffix']
                the_snapdir=info['the_snapdir']
                Nsnapstring=info['Nsnapstring']
                havecr=info['havecr']
                Fcal=info['Fcal']
                iavesfr=info['iavesfr']
                timestep=info['timestep']
                maindir=info['maindir']
                haveB=info['haveB']
                usepep=info['usepep']
                Rvir=info['Rvir']
                print 'the_snapdir', the_snapdir
                print 'Nsnapstring', Nsnapstring
                print 'havecr', havecr
                cosmo=info['cosmo']
                if cosmo==1:
                        h0=1
                else:
                        h0=0
                header = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo, header_only=1)
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo)
                if cosmo==1:
                        halosingle = SF.read_halo_history(rundir, halonostr='00',hubble=hubble, comoving=0, maindir=maindir, singlesnap=1, atime=ascale)
                        xcen = halosingle['x']
                        ycen = halosingle['y']
                        zcen = halosingle['z']
                        xvcen = halosingle['xv']
                        yvcen = halosingle['yv']
                        zvcen = halosingle['zv']
                        Rvirnow = halosingle['R']
                        MgAHF = halosingle['Mg']
                else:
                        xcen=0
                        ycen=0
                        if shiftz==1:
                                zcen=findcenz(runtodo,ino)
                        else:
                                zcen=0.
                        print 'zcen', zcen
                        xvcen=0
                        yvcen=0
                        zvcen=0
                Gpos = G['p']
                Gvel = G['v']
                Gu = G['u'] #internal energy per unit mass # in km^2/s^2
                Gx = Gpos[:,0]-xcen
                Gy = Gpos[:,1]-ycen
                Gz = Gpos[:,2]-zcen
                Gvx = Gvel[:,0]-xvcen   #km/s
                Gvy = Gvel[:,1]-yvcen
                Gvz = Gvel[:,2]-zvcen
                Grho = G['rho']*1e10*Msun_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm #g/cm^3
                Gm = G['m']*1e10 #Msun
                Neb = G['ne']
                if havecr>0:
                        cregy = G['cregy']*1e10*Msun_in_g*km_in_cm*km_in_cm #erg
                TrueTemp, converted_rho  = SF.convertTemp(Gu, Neb, Grho)
                cutxy = Gx*Gx+Gy*Gy < withinr*withinr
                cutz = Gz*Gz < maxlength*maxlength
                cut = cutxy*cutz
                Gmcut = Gm[cut]
                Tcut = TrueTemp[cut]
                coldgascut = Tcut < 1.0e4 #K
                totdiskgm += np.sum(Gmcut) #Msun
                diskcoldgm += np.sum(Gmcut[coldgascut]) #Msun

                cutr = Gx*Gx+Gy*Gy < withinr*withinr/4.
                Gmcutr = Gm[cutr]
                Gzcutr = Gz[cutr]
                mlist=[]
                zlist=np.linspace(0.1,50,num=1000)
                for i in range(len(zlist)):
                        withinzlist = np.absolute(Gzcutr)<zlist[i]
                        withinzm = np.sum(Gmcutr[withinzlist])
                        mlist = np.append(mlist,withinzm)
                z12gas += np.interp(mlist[-1]*0.5,mlist,zlist)


                mlist=[]
                zlist=np.linspace(0.1,50,num=1000)
                Gnism_in_cm_3 = (0.78+0.22*Neb*0.76)/proton_mass_in_g*Grho
                dencut = Gnism_in_cm_3[cutr]>10.0
                Gmden = Gmcutr[dencut]
                Gzden = Gzcutr[dencut]
                for i in range(len(zlist)):
                        withinzlist = np.absolute(Gzden)<zlist[i]
                        withinzm = np.sum(Gmden[withinzlist])
                        mlist = np.append(mlist,withinzm)
                z12dengas += np.interp(mlist[-1]*0.5,mlist,zlist)

                Grr = np.sqrt(Gx*Gx+Gy*Gy+Gz*Gz)
                rvcut = Grr<Rvir*0.25
                Gmrv = Gm[rvcut]
                Gxrv = Gx[rvcut]
                Gyrv = Gy[rvcut]
                Gzrv = Gz[rvcut]
                Gr = np.sqrt(Gxrv*Gxrv+Gzrv*Gzrv)
                rgsd = 1.0 # consider gas mass only within
                rcut = Gr<rgsd #kpc
                Gmwithin = np.sum(Gmrv[rcut])*Msun_in_g
                gassurden += Gmwithin/(np.pi*rgsd*rgsd*kpc_in_cm*kpc_in_cm)


                if havecr>0:
                        crden += np.sum(cregy[cut]*erg_in_eV)/(np.pi*withinr*withinr*maxlength*2.0*kpc_in_cm*kpc_in_cm*kpc_in_cm)
                        crlist=[]
                        rlist=np.linspace(0.1,50,num=1000)
                        for i in range(len(rlist)):
                                withinrlist = np.absolute(Gz)<rlist[i]
                                withinrcr = np.sum(cregy[withinrlist])
                                crlist = np.append(crlist,withinrcr)
                        r12cr+= np.interp(crlist[-1]*0.5,crlist,rlist)
                else:
                        r12cr+=0.0
                nooftimes+=1



        S = readsnapcr(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
        Sm = S['m']*1e10 #Msun
        Sp = S['p']
        Sage = S['age']
        Sx = Sp[:,0]
        Sy = Sp[:,1]
        Sz = Sp[:,2]
        header=readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1,header_only=1)
        h0 = header['hubble']
        atime = header['time']
        if cosmo==1:
                if usepep==1:
                        halosA = SF.read_halo_history_pep(rundir, finalno, beginno=beginno,\
                         singlesnap=0, firever=firever,halonostr=halostr, hubble=h0, comoving=1, maindir=maindir)
                        afactor=atime
                else:
                        halosA = SF.read_halo_history(rundir, maindir=maindir, halonostr=halostr, hubble=h0, comoving=0)
                        afactor=1.0
                redlist = halosA['redshift']
                haloid = halosA['ID']
                a_scale = 1.0/(1.0+redlist)
                xcenl = halosA['x']*afactor
                ycenl = halosA['y']*afactor
                zcenl = halosA['z']*afactor
                Rvirl = halosA['R']*afactor
                xcen = np.interp(atime,a_scale,xcenl)
                ycen = np.interp(atime,a_scale,ycenl)
                zcen = np.interp(atime,a_scale,zcenl)
                Rvir = np.interp(atime,a_scale,Rvirl)
                Sxrel = Sx-xcen
                Syrel = Sy-ycen
                Szrel = Sz-zcen
        print 'Sage', np.amin(Sage), np.amax(Sage)
        agecut = Sage>0
        Smnew = np.sum(Sm[agecut])


        outdata = dirout(runtodo,startno,Nsnap,snapsep, Rfrac)
        Lgamma=0.
        if havecr>4:
                Lgamma=np.average(outdata['Lgamma'])
        elif havecr>0:
                Lgamma=np.average(outdata['Lgcalist'])
        avesfr = np.average(outdata['avesfrl'])

        return {'avesfr':avesfr, 'Smnew':Smnew, 'totdiskgm':totdiskgm/nooftimes, 'diskcoldgm':diskcoldgm/nooftimes, 'Lgamma':Lgamma,\
'crden':crden/nooftimes, 'z12gas':z12gas/nooftimes, 'r12cr':r12cr/nooftimes, 'z12dengas': z12dengas/nooftimes, 'gassurden':gassurden/nooftimes}


def findcenz(runtodo,Nsnap,dir='z',withinr=15.,datasup=0,Gx=[],Gy=[],Gz=[],Gm=[]):
        maxlength = 10.
        nogrid = 50
        info=outdirname(runtodo, Nsnap)
        rundir=info['rundir']
        maindir=info['maindir']
        halostr=info['halostr']
        runtitle=info['runtitle']
        slabel=info['slabel']
        snlabel=info['snlabel']
        dclabel=info['dclabel']
        resolabel=info['resolabel']
        the_snapdir=info['the_snapdir']
        Nsnapstring=info['Nsnapstring']
        havecr=info['havecr']
        Fcal=info['Fcal']
        iavesfr=info['iavesfr']
        timestep=info['timestep']
        cosmo=info['cosmo']
        color=info['color']
        withinRv=info['withinRv']
        usepep=info['usepep']
        beginno=info['beginno']
        finalno=info['finalno']
        firever=info['firever']
        initsnap=info['initsnap']
        haveB=info['haveB']
        M1speed=info['M1speed']
        Rvirguess=info['Rvir']
        the_prefix=info['the_prefix']
        the_suffix=info['the_suffix']
        if datasup==0:
                G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                Gp = G['p']
                Gm = G['m']
                Gz = Gp[:,2]
                Gx = Gp[:,0]
                Gy = Gp[:,1]
        
        if dir=='z':
                px=Gx;py=Gy;pz=Gz
        if dir=='y':
                px=Gz;py=Gx;pz=Gy
        if dir=='x':
                px=Gy;py=Gz;pz=Gx
        dz = maxlength/nogrid
        Gnism_in_cm_3l=[]
        zl =[]
        for iz in range(nogrid):
                cutxy = px*px+py*py < withinr*withinr
                cutz = (pz> dz*iz-maxlength/2.) & (pz<dz*(iz+1)-maxlength/2.)
                cut = cutxy*cutz
                Gmcut = Gm[cut]
                Gm_in_g = Gmcut*1e10*2e33
                shellvol_in_cm3 = np.pi*dz*np.power(withinr,2)*3.086e21*3.086e21*3.086e21
                Grho_in_g_cm_3 = Gm_in_g/shellvol_in_cm3
                Gnism_in_cm_3 = np.sum(1./protonmass_in_g*Grho_in_g_cm_3)
                Gnism_in_cm_3l = np.append(Gnism_in_cm_3l, Gnism_in_cm_3)
                zl = np.append(zl, dz*(iz+0.5)-maxlength/2.)
        print 'Gnism_in_cm_3l', Gnism_in_cm_3l
        zcen = zl[np.argmax(Gnism_in_cm_3l)]
        print 'zcen',zcen
        print 'Gnism_in_cm_3l', Gnism_in_cm_3l[np.argmax(Gnism_in_cm_3l)]
        return zcen
    
    
def findcennew(runtodo,Nsnap,withinr=15.,datasup=0,Gx=[],Gy=[],Gz=[],Gm=[]):
        maxlength = 10.
        nogrid = 50
        info=outdirname(runtodo, Nsnap)
        rundir=info['rundir']
        maindir=info['maindir']
        halostr=info['halostr']
        runtitle=info['runtitle']
        slabel=info['slabel']
        snlabel=info['snlabel']
        dclabel=info['dclabel']
        resolabel=info['resolabel']
        the_snapdir=info['the_snapdir']
        Nsnapstring=info['Nsnapstring']
        havecr=info['havecr']
        Fcal=info['Fcal']
        iavesfr=info['iavesfr']
        timestep=info['timestep']
        cosmo=info['cosmo']
        color=info['color']
        withinRv=info['withinRv']
        usepep=info['usepep']
        beginno=info['beginno']
        finalno=info['finalno']
        firever=info['firever']
        initsnap=info['initsnap']
        haveB=info['haveB']
        M1speed=info['M1speed']
        Rvirguess=info['Rvir']
        the_prefix=info['the_prefix']
        the_suffix=info['the_suffix']
        if datasup==0:
                G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                Gp = G['p']
                Gm = G['m']
                Gz = Gp[:,2]
                Gx = Gp[:,0]
                Gy = Gp[:,1]
        Gr = np.sqrt(Gx*Gx+Gy*Gy+Gz*Gz)
        cutr = Gr<withinr
        Gx=Gx[cutr]; Gy=Gy[cutr]; Gz=Gz[cutr]; Gm=Gm[cutr];
        H, (xedges, yedges, zedges) = np.histogramdd((Gx, Gy, Gz), weights=Gm, bins = nogrid)
        maxindex = np.unravel_index(H.argmax(), H.shape)
        xcor = (xedges[:-1]+xedges[1:])/2.; ycor = (yedges[:-1]+yedges[1:])/2.; zcor = (zedges[:-1]+zedges[1:])/2.
        xmax = xcor[maxindex[0]]; ymax = ycor[maxindex[1]]; zmax = zcor[maxindex[2]]
        return xmax, ymax, zmax


def gaswindphase(runtodo,Nsnap,rotface=1,Tcut=1.0,highTcut=1e10,vcut=0.0,
                 vhcut=1e10,withinr=20.0,zup=25.0,zdown=20.0,userad=0):
        info=outdirname(runtodo, Nsnap)
        rundir=info['rundir']
        maindir=info['maindir']
        halostr=info['halostr']
        runtitle=info['runtitle']
        slabel=info['slabel']
        snlabel=info['snlabel']
        dclabel=info['dclabel']
        resolabel=info['resolabel']
        the_snapdir=info['the_snapdir']
        Nsnapstring=info['Nsnapstring']
        havecr=info['havecr']
        Fcal=info['Fcal']
        iavesfr=info['iavesfr']
        timestep=info['timestep']
        color=info['color']
        cosmo=info['cosmo']
        usepep=info['usepep']
        haveB=info['haveB']
        Rvir=info['Rvir']
        the_prefix=info['the_prefix']
        the_suffix=info['the_suffix']
        finalno=info['finalno']
        beginno=info['beginno']
        firever=info['firever']
        h0 = info['h0']
        snapsep=info['snapsep'] 
        if havecr==0:
                cregy=[]
        vmax=8.5
        vmin=0
        if runtitle=='MW':
                vmax = 7.0
                vmin = -1.0
        if runtitle=='SMC':
                vmax = 5.0
                vmin = 0.5
        if runtitle=='SBC':
                vmax = 9.0
                vmin=-1.0
        if cosmo==1:
                print 'the_snapdir', the_snapdir
                print 'Nsnapstring', Nsnapstring
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                header=readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1,header_only=1)
                h0 = header['hubble']
                atime = header['time']
                if usepep==1:
                        halosA = SF.read_halo_history_pep(rundir, finalno, snapsep=snapsep,\
                        beginno=beginno, singlesnap=0, firever=firever,halonostr=halostr, hubble=h0, comoving=1, maindir=maindir)
                        afactor=atime
                else:
                        halosA = SF.read_halo_history(rundir, maindir=maindir, halonostr=halostr, hubble=h0, comoving=0)
                        afactor=1.0
                redlist = halosA['redshift']
                haloid = halosA['ID']
                a_scale = 1.0/(1.0+redlist)
                xcenl = halosA['x']*afactor
                ycenl = halosA['y']*afactor
                zcenl = halosA['z']*afactor
                xvcenl = halosA['xv']
                yvcenl = halosA['yv']
                zvcenl = halosA['zv']
                Rvirl = halosA['R']*afactor
                Lxstarl = halosA['Lxstar']
                Lystarl = halosA['Lystar']
                Lzstarl = halosA['Lzstar']
                xcen = np.interp(atime,a_scale,xcenl)
                ycen = np.interp(atime,a_scale,ycenl)
                zcen = np.interp(atime,a_scale,zcenl)
                xvcen = np.interp(atime,a_scale,xvcenl)
                yvcen = np.interp(atime,a_scale,yvcenl)
                zvcen = np.interp(atime,a_scale,zvcenl)
                Rvir = np.interp(atime,a_scale,Rvirl)
                Lxstar = np.interp(atime,a_scale,Lxstarl)
                Lystar = np.interp(atime,a_scale,Lystarl)
                Lzstar = np.interp(atime,a_scale,Lzstarl)
        else:
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                xcen =0; ycen=0; zcen=0; xvcen=0; yvcen=0; zvcen=0;
        Gpos = G['p']
        Gmass = G['m'][:]
        Gvel = G['v'][:,:]
        Gu = G['u'][:]
        rho = G['rho'][:]
        Neb = G['ne'][:]
        Gid = G['id'][:]
        if havecr>0:
                cregy=G['cregy']
        partX = Gpos[:,0]-xcen
        partY = Gpos[:,1]-ycen
        partZ = Gpos[:,2]-zcen
        partR = np.sqrt(partX*partX+partY*partY+partZ*partZ)
        vx = Gvel[:,0]-xvcen
        vy = Gvel[:,1]-yvcen
        vz = Gvel[:,2]-zvcen
        if rotface==1:
                Gr =  np.sqrt(partX*partX+partY*partY+partZ*partZ)
                cutr = Gr < 6. #kpc
                Gmcutr = Gmass[cutr]
                Gxcutr = partX[cutr]; Gycutr = partY[cutr]; Gzcutr = partZ[cutr];
                Gvxcutr = vx[cutr]; Gvycutr = vy[cutr]; Gvzcutr = vz[cutr];
                Lang = [0.,0.,0.]
                for i in range(len(Gxcutr)):
                        Lang += Gmcutr[i]*np.cross([Gxcutr[i],Gycutr[i],Gzcutr[i]],[Gvxcutr[i],Gvycutr[i],Gvzcutr[i]])
                partX,partY,partZ = SF.rotateL_to_z(partX,partY,partZ,Lang[0],Lang[1],Lang[2])
                vx, vy, vz = SF.rotateL_to_z(vx,vy,vz,Lang[0],Lang[1],Lang[2])
        vr = (vx*partX+vy*partY+vz*partZ)/partR
        header = G['header'][:]
        redshift = header[3]
        boxsize = header[9]
        if userad==1:
                cutr = (np.absolute(partR)>zdown) & (np.absolute(partR)<zup)
                cutv = (vr>vcut)*(vr<vhcut)
                cut = cutr*cutv
        else:
                cutz = (np.absolute(partZ)>zdown) & (np.absolute(partZ)<zup) #kpc
                cutxy = partX*partX+partY*partY<withinr*withinr
                absvz = vz*partZ/np.absolute(partZ)
                cutv = (absvz>vcut)*(absvz<vhcut) #outflowing gas
                cut = cutz*cutxy*cutv
        if Tcut>0.0:
                TrueTemp, converted_rho  = SF.convertTemp(Gu, Neb, rho)
                cutT = TrueTemp>Tcut
                cuthT = TrueTemp<highTcut
                cut = cut*cutT*cuthT
        Gu = Gu[cut]
        rho = rho[cut]
        Neb = Neb[cut]
        Gmass = Gmass[cut]
        partX = partX[cut]
        partY = partY[cut]
        partZ = partZ[cut]
        partR = partR[cut]
        vx = vx[cut]
        vy = vy[cut]
        vz = vz[cut]
        vr = vr[cut]
        Gid = Gid[cut]
        if havecr>0:
                cregy = cregy[cut]
        TrueTemp, converted_rho  = SF.convertTemp(Gu, Neb, rho)
        del G, Gpos, Gvel
        return {'vmax':vmax,'vmin':vmin,'Gid':Gid,'TrueTemp':TrueTemp,'convertedrho':converted_rho,\
 'vx':vx, 'vy':vy, 'vz':vz, 'vr':vr, 'partX':partX,'partY':partY, 'partZ':partZ, 'partR':partR,\
'Gu':Gu, 'Gmass':Gmass,'rho':rho,'cregy':cregy}



def gastrack(runtodo,snaptrack,snapstart,Tcut_t=1.0,highTcut_t=1e10,\
        vcut_t=0.0,vhcut_t=1e10,withinr_t=20.0,zup_t=25.0,zdown_t=20.0,userad=0):
        data = gaswindphase(runtodo,snaptrack,Tcut=Tcut_t,highTcut=highTcut_t,\
        vcut=vcut_t,vhcut=vhcut_t,withinr=withinr_t,zup=zup_t,zdown=zdown_t)
        Gid_t = data['Gid']
        #print 'Gid_t', np.sort(Gid_t)
        vcut = -1e10 #if we track gas particles, we should consider all velocities, temp, etc
        withinr=1e10
        zup=1e10
        zdown=0.001
        Tcut = 1.0
        highTcut =1e15
        vhcut = 1e15
        data = gaswindphase(runtodo,snapstart,Tcut=Tcut,highTcut=highTcut,\
        vcut=vcut,vhcut=vhcut,withinr=withinr,zup=zup,zdown=zdown,userad=userad)
        partX = data['partX']
        partY = data['partY']
        partZ = data['partZ']
        partR = data['partR']
        vx = data['vx']
        vy = data['vy']
        vz = data['vz']
        vr = data['vr']
        TrueTemp = data['TrueTemp']
        converted_rho = data['convertedrho']
        Gmass = data['Gmass']
        vmax = data['vmax']
        vmin = data['vmin']
        Gid = data['Gid']
        idint0=np.in1d(Gid,Gid_t)
        partR=partR[idint0]
        vr=vr[idint0]
        partX=partX[idint0]
        partY=partY[idint0]
        partZ=partZ[idint0]
        vx=vx[idint0]
        vy=vy[idint0]
        vz=vz[idint0]
        TrueTemp=TrueTemp[idint0]
        converted_rho=converted_rho[idint0]
        Gmass = Gmass[idint0]
        Gid = Gid[idint0]
        #print 'Gid_after', np.sort(Gid)
        return {'vmax':vmax,'vmin':vmin,'Gid':Gid,'TrueTemp':TrueTemp,'convertedrho':converted_rho,\
 'vx':vx, 'vy':vy, 'vz':vz, 'vr':vr, 'partX':partX,'partY':partY, 'partZ':partZ, 'partR':partR, 'Gmass':Gmass}



def outsfr(runtodo, Nsnap, tsep=10.0):
        info=outdirname(runtodo, Nsnap)
        rundir=info['rundir']
        maindir=info['maindir']
        halostr=info['halostr']
        runtitle=info['runtitle']
        slabel=info['slabel']
        snlabel=info['snlabel']
        dclabel=info['dclabel']
        resolabel=info['resolabel']
        the_snapdir=info['the_snapdir']
        Nsnapstring=info['Nsnapstring']
        havecr=info['havecr']
        Fcal=info['Fcal']
        iavesfr=info['iavesfr']
        timestep=info['timestep']
        cosmo=info['cosmo']
        color=info['color']
        withinRv=info['withinRv']
        usepep=info['usepep']
        beginno=info['beginno']
        finalno=info['finalno']
        snapsep=info['snapsep']
        firever=info['firever']
        initsnap=info['initsnap']
        haveB=info['haveB']
        M1speed=info['M1speed']
        Rvirguess=info['Rvir']
        the_prefix=info['the_prefix']
        the_suffix=info['the_suffix']
        correctIa=info['correctIa']
        halostr=info['halostr']
        galcen = info['galcen']
        if cosmo==0:
                inittime=initsnap
        else:
                inittime=0
                print 'set initial time'
        print 'the_snapdir', the_snapdir
        print 'Nsnapstring', Nsnapstring
        print 'havecr', havecr
        print 'withinRv', withinRv
        print 'cosmo', cosmo
        try:
                if cosmo==1:
                        S = readsnapcr(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                else:
                        S = readsnapcr(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
        except KeyError:
                print 'Keyerror'
        try:
                header=S['header']
                timeneed=header[2]
                print 'timeneed', timeneed
                Smi=S['m']
                Sage=S['age']
                apre = timeneed-1.0e-3*tsep
                if cosmo==1:
                        readtimelist=readtime(firever=2)
                        snap2list=readtimelist['snaplist']
                        time2list=readtimelist['timelist']
                        a2list=readtimelist['alist']
                        tnow_in_yr = np.interp(timeneed,a2list,time2list)*1e9
                        tpre = tnow_in_yr-1.0e6*tsep
                        apre = np.interp(tpre/1.0e9,time2list,a2list)
                        Sp = S['p']
                        Sx = Sp[:,0]
                        Sy = Sp[:,1]
                        Sz = Sp[:,2]
                        header=readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1,header_only=1)
                        h0 = header['hubble']
                        atime = header['time']
                        print 'halostr', halostr
                        if usepep==1:
                                halosA = SF.read_halo_history_pep(rundir, finalno, snapsep=snapsep, beginno=beginno, singlesnap=0, firever=firever,halonostr=halostr, hubble=h0, comoving=0, maindir=maindir)
                                afactor=atime
                        else:
                                halosA = SF.read_halo_history(rundir, maindir=maindir, halonostr=halostr, hubble=h0, comoving=0)
                                afactor=1.0
                        redlist = halosA['redshift']
                        haloid = halosA['ID']
                        a_scale = 1.0/(1.0+redlist)
                        xcenl = halosA['x']*afactor
                        ycenl = halosA['y']*afactor
                        zcenl = halosA['z']*afactor
                        Rvirl = halosA['R']*afactor
                        if galcen == 1:
                                finname = '/home/tkc004/'+maindir+'/'+rundir+'/center/galcen.txt'
                                f = open(finname)
                                f.readline()
                                dars = f.readlines()
                                f.close()
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
                                xcen=np.interp(Nsnap,snapD,galXl) #physical distance (kpc) (new program for FIRE 2.0 uses physical distance)
                                ycen=np.interp(Nsnap,snapD,galYl)
                                zcen=np.interp(Nsnap,snapD,galZl)
                        else:
                                xcen = np.interp(atime,a_scale,xcenl)
                                ycen = np.interp(atime,a_scale,ycenl)
                                zcen = np.interp(atime,a_scale,zcenl)
                        Rvir = np.interp(atime,a_scale,Rvirl)
                        Sxrel = Sx-xcen
                        Syrel = Sy-ycen
                        Szrel = Sz-zcen
                        Sr = np.sqrt(Sxrel*Sxrel+Syrel*Syrel+Szrel*Szrel)
                        cutrvs = Sr<Rvir*0.1
                        Smi = Smi[cutrvs]
                        Sage = Sage[cutrvs]
                Sm = np.sum(Smi)*1e10 #in solar mass
                tcut=Sage>apre
                Nsm = np.sum(Smi[tcut])*1e10
                SFR = Nsm/1.0e6/tsep #Msun/yr
                #del Sx, Sy, Sz, Sr, S, Smi, Sage, Sm, Sp
        except KeyError:
                print 'Keyerror'
        return {'SFR':SFR, 'Nsm':Nsm}

    
    
def gasphase(runtodo,Nsnap,rotface=1,Tcut=1.0,highTcut=1e10,vcut=0.0,vhcut=1e10,withinr=20.0,
             withoutr=0.01,zup=25.0,zdown=20.0,userad=0):
        loccen=0
        info=outdirname(runtodo, Nsnap)
        rundir=info['rundir']
        maindir=info['maindir']
        halostr=info['halostr']
        runtitle=info['runtitle']
        slabel=info['slabel']
        snlabel=info['snlabel']
        dclabel=info['dclabel']
        resolabel=info['resolabel']
        the_snapdir=info['the_snapdir']
        Nsnapstring=info['Nsnapstring']
        havecr=info['havecr']
        Fcal=info['Fcal']
        iavesfr=info['iavesfr']
        timestep=info['timestep']
        color=info['color']
        cosmo=info['cosmo']
        usepep=info['usepep']
        haveB=info['haveB']
        Rvir=info['Rvir']
        the_prefix=info['the_prefix']
        the_suffix=info['the_suffix']
        finalno=info['finalno']
        beginno=info['beginno']
        firever=info['firever']
        h0 = info['h0']
        snapsep=info['snapsep'] 
        snumadd=info['snumadd']
        if havecr==0:
            cregy=[]
        if haveB==0:
            Brms=[]
        vmax=8.5
        vmin=0
        if runtitle=='MW':
                vmax = 7.0
                vmin = -1.0
        if runtitle=='SMC':
                vmax = 5.0
                vmin = 0.5
        if runtitle=='SBC':
                vmax = 9.0
                vmin=-1.0
        if cosmo==1:
                print 'the_snapdir', the_snapdir
                print 'Nsnapstring', Nsnapstring
        G = readsnapwcen(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix,\
         havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
         loccen=loccen,runtodo=runtodo,rundir=rundir,halostr=halostr)
        if havecr>0:
            cregy=G['cregy']
        if haveB>0:
            GB = G['B']
            Bx = GB[:,0]; By = GB[:,1]; Bz = GB[:,2];
            Brms = np.sqrt(Bx*Bx+By*By+Bz*Bz)
        Gpos = G['p']
        Gvel = G['v']
        Gu = G['u']
        rho = G['rho']
        Neb = G['ne']
        Gmass = G['m']
        Gid = G['id']
        partX = Gpos[:,0]
        partY = Gpos[:,1]
        partZ = Gpos[:,2]
        partR = np.sqrt(partX*partX+partY*partY+partZ*partZ)
        vx = Gvel[:,0]
        vy = Gvel[:,1]
        vz = Gvel[:,2]
        vr = (vx*partX+vy*partY+vz*partZ)/partR
        header = G['header'][:]
        redshift = header[3]
        boxsize = header[9]
        if userad==1:
                cutr = (np.absolute(partR)>zdown) & (np.absolute(partR)<zup)
                cutv = (vr>vcut)*(vr<vhcut)
                cut = cutr*cutv
        else:
                cutz = (np.absolute(partZ)>zdown) & (np.absolute(partZ)<zup) #kpc
                cutxyin = partX*partX+partY*partY<withinr*withinr
                cutxyout = partX*partX+partY*partY>withoutr*withoutr
                cutxy = cutxyin*cutxyout
                absvz = vz*partZ/np.absolute(partZ)
                cutv = (absvz>vcut)*(absvz<vhcut) #outflowing gas
                cut = cutz*cutxy*cutv
        if Tcut>0.0:
                TrueTemp, converted_rho  = SF.convertTemp(Gu, Neb, rho)
                cutT = TrueTemp>Tcut
                cuthT = TrueTemp<highTcut
                cut = cut*cutT*cuthT
        Gu = Gu[cut]
        rho = rho[cut]
        Neb = Neb[cut]
        Gmass = Gmass[cut]
        partX = partX[cut]
        partY = partY[cut]
        partZ = partZ[cut]
        partR = partR[cut]
        vx = vx[cut]
        vy = vy[cut]
        vz = vz[cut]
        vr = vr[cut]
        Gid = Gid[cut]
        if havecr>0:
            cregy = cregy[cut]
        if haveB>0:
            Brms = Brms[cut]
        TrueTemp, converted_rho  = SF.convertTemp(Gu, Neb, rho)
        del G, Gpos, Gvel
        return {'vmax':vmax,'vmin':vmin,'Gid':Gid,'TrueTemp':TrueTemp,'convertedrho':converted_rho,\
 'vx':vx, 'vy':vy, 'vz':vz, 'vr':vr, 'partX':partX,'partY':partY, 'partZ':partZ, 'partR':partR,\
'Gu':Gu, 'Gmass':Gmass,'rho':rho,'cregy':cregy, 'Brms':Brms}
    
    
    
def calstrloss(G, withinr = 20., nobin = 30):
    '''
    A rough calculation of streaming loss vA dot gradP:
    divide the volume radially; then calculate the dp/dr and vA; 
    streaming loss ~ vA * dp/dr
    
    G: data structure from readsnap after re-center (in physical units)
    i.e. coordinate in kpc, B in gauss, m in 1e10 msun, etc.
    '''
    
    Gp = G['p']; Gm = G['m']; GB = G['B']; cregy = G['cregy']; Grho = G['rho'];
    Gu = G['u'];
    vth = np.sqrt(Gu)*km_in_cm
    Gx = Gp[:,0]; Gy = Gp[:,1]; Gz = Gp[:,2];
    Gr = np.sqrt(Gx*Gx+Gy*Gy+Gz*Gz)
    GB2 = GB[:,0]*GB[:,0]+GB[:,1]*GB[:,1]+GB[:,2]*GB[:,2];
    GBmag = np.sqrt(GB2);
    cregy_in_cgs = cregy*cregy_codetocgs
    Grho_in_cgs = Grho*rho_codetocgs
    
    radl = np.linspace(0.01, withinr, num=nobin)
    Pcrl = np.array([])
    vA2l = np.array([])
    for i in range(nobin-1):
        cutrin = Gr<radl[i+1]
        cutrout = Gr>radl[i]
        cut = cutrin*cutrout
        vol = 4./3.*np.pi*(np.power(radl[i+1],3)-np.power(radl[i],3))
        vol_in_cgs = vol*vol_codetocgs
        creden_in_cgs = np.sum(cregy_in_cgs[cut])/vol_in_cgs
        Pcr_in_cgs = creden_in_cgs*(CRgamma-1.)
        #vA2_in_cgs = np.average(GB2[cut]/(4.*np.pi*Grho_in_cgs[cut])+vth[cut]*vth[cut],weights=Gm[cut])
        vA2_in_cgs = np.average(GB2[cut]/(4.*np.pi*Grho_in_cgs[cut]),weights=Gm[cut])
        Pcrl = np.append(Pcrl,Pcr_in_cgs); vA2l = np.append(vA2l,vA2_in_cgs);
    radlm = (radl[1:]+radl[:-1])/2.
    volm = 4./3.*np.pi*(np.power(radlm[1:],3)-np.power(radlm[:-1],3))*vol_codetocgs
    dPcr_dr = (Pcrl[1:]-Pcrl[:-1])/(radlm[1:]-radlm[:-1])/x_codetocgs
    vA2lm = (vA2l[1:]+vA2l[:-1])/2.
    vAlm = np.sqrt(vA2lm)
    strl =  -np.absolute(vAlm*dPcr_dr) #streaming loss in erg/s/cm^3
    totstrl = np.sum(strl*volm)
    return {'totstrl':totstrl, 'strl':strl, 'vAlm':vAlm, 'dPcr_dr':dPcr_dr, 'radlm':radlm}




def energyoutput(runtodo,Nsnap,shiftz=1,usesolarcircle=0,rotface=1,usecentral=0,cylr=9.0,cylrin=0.01,cylz=0.5,nbin = 5):
        #see Su 2016 for the definition of turbulent energy
        #cylr=10. #kpc
        #cylz=2. #kpc #thickness of the cylinder we consider
        
        info=outdirname(runtodo, Nsnap)
        rundir=info['rundir']
        runtitle=info['runtitle']
        slabel=info['slabel']
        snlabel=info['snlabel']
        dclabel=info['dclabel']
        resolabel=info['resolabel']
        the_snapdir=info['the_snapdir']
        Nsnapstring=info['Nsnapstring']
        havecr=info['havecr']
        Fcal=info['Fcal']
        iavesfr=info['iavesfr']
        timestep=info['timestep']
        maindir=info['maindir']
        haveB=info['haveB']
        the_prefix=info['the_prefix']
        the_suffix=info['the_suffix']
        usepep=info['usepep']
        snumadd=info['snumadd']
        halostr=info['halostr']
        print 'the_snapdir', the_snapdir
        print 'Nsnapstring', Nsnapstring
        print 'havecr', havecr
        cosmo=info['cosmo']
        if cosmo==1:
                h0=1
        else:
                h0=0
        header = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo, header_only=1)
        ascale = header['time']
        #print 'this time', ascale
        thisred = 1./ascale-1.
        hubble = header['hubble']
        #print 'hubble', hubble
        G = readsnapwcen(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix,\
        havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
        loccen=shiftz,runtodo=runtodo,rundir=rundir,halostr=halostr)
        egydata = calenergyoutput(G, usesolarcircle=usesolarcircle,usecentral=usecentral,cylr=cylr,cylrin=cylrin,cylz=cylz,nbin = nbin,havecr=havecr,haveB=haveB,cosmo=cosmo,Nsnap=Nsnap,ascale=ascale)
        return egydata
        
def calenergyoutput(G, usesolarcircle=0,usecentral=0,cylr=9.0,cylrin=0.01,cylz=0.5,nbin = 5,havecr=0, haveB=0,cosmo=0,Nsnap=600,ascale=1.0):        
        if usesolarcircle==1:
                cylrin=7.0; cylr=9.0
        if usecentral==1:
                cylrin=0.01; cylr=3.0
        fraccut = 68.0 #in percentage
        turl = []
        gml = []
        gvl = []
        gminl = []
        cregyl = []
        therml = []
        Begyl = []
        Begy_2Bzl = []
        Evzl = []
        xcell = []
        ycell = []
        zcell = []
        Omegal = [] #angular rotation frequency in rad/s
        Gpos = G['p']
        Gvel = G['v']*km_in_cm #now in cm/s
        Gx = Gpos[:,0]
        Gy = Gpos[:,1]
        Gz = Gpos[:,2]
        Gvx = Gvel[:,0]
        Gvy = Gvel[:,1]
        Gvz = Gvel[:,2]
        Grho = G['rho']*1e10*solar_mass_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm
        if havecr>0:
                cregy  = G['cregy']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
        Gm = G['m']*1e10*solar_mass_in_g #in g
        Gtotm = np.sum(Gm)
        GEint = G['u']*km_in_cm*km_in_cm*Gm #in erg
        Gvol = Gm/Grho #in cm3
        if haveB>0:
                GB = G['B']
                Bx = GB[:,0]
                By = GB[:,1]
                Bz = GB[:,2]
                B2 = Bx*Bx+By*By+Bz*Bz
                Begy = B2/np.pi/8.*Gm/Grho
                Begy_2Bz = Begy-Bz*Bz/np.pi/4.*Gm/Grho
        Evz = Gm*Gvz*Gvz        
        Grxy = np.sqrt(Gx*Gx+Gy*Gy)
        dr = (cylr-cylrin)/float(nbin)
        cutz = np.absolute(Gz)<cylz/2.
        radl=[]
        for nrad in range(nbin+1):
                radl = np.append(radl,cylrin+nrad*dr)
        for nrad in range(nbin):
                cutr = (Grxy<radl[nrad+1]) & (Grxy>radl[nrad])
                #print 'radl[nrad+1]', radl[nrad+1]
                #print 'radl[nrad]', radl[nrad]
                cut=cutr*cutz
                Gmc=Gm[cut]
                Grhoc=Grho[cut]
                Grc=Grxy[cut]
                Gxc=Gx[cut]
                Gyc=Gy[cut]
                Gzc=Gz[cut]
                Gvxc = Gvx[cut]
                Gvyc = Gvy[cut]
                Gvzc = Gvz[cut]
                #print 'len(Gmc)', len(Gmc)
                Gvzcave = np.average(Gvz[cut],weights=Gmc)
                #calculate average rotational velocity:
                #theta hat dot v = vtheta
                Grxyc = np.sqrt(Gxc*Gxc+Gyc*Gyc)
                vth = -Gyc/Grxyc*Gvxc+Gxc/Grxyc*Gvyc
                vthave = np.average(vth,weights=Gmc)
                Omega = vthave/(radl[nrad+1]+radl[nrad])*2.0/kpc_in_cm
                #print 'Gvxc', Gvxc
                Gvxc = Gvxc+Gyc/Grxyc*vthave
                Gvyc = Gvyc-Gxc/Grxyc*vthave
                vzlog = np.log10(np.absolute(Gvzc-Gvzcave))
                logvzsigma = weighted_quantile(vzlog,fraccut,sample_weight=Gmc)
                cutvz = np.log10(vzlog)<logvzsigma ###
                Grcz = Grc[cutvz]
                Gmcz = Gmc[cutvz]
                Grhoz = Grhoc[cutvz]
                Gzcz = Gzc[cutvz]
                Gxcz = Gxc[cutvz]
                Gycz = Gyc[cutvz]
                Gvzcz = Gvzc[cutvz]
                Gvxcz = Gvxc[cutvz]
                Gvycz = Gvyc[cutvz]
                totnum = len(Gmcz)
                print 'totnum', totnum
                vol = np.pi*(radl[nrad+1]*radl[nrad+1]-radl[nrad]*radl[nrad])*cylz
                numden = totnum/vol
                vol15 = 15.0/numden #TK test
                dl = np.power(vol15,1./3.)
                print 'dl', dl
                #print 'int(dr/dl)+1', int(cylz/dl)+1
                numcount=0
                rscl = []
                for nrs in range(int(dr/dl)+2):
                        prsc = 100.0/(int(dr/dl)+1)*nrs
                        if prsc>100.0:
                                prsc =100.0
                        rsc = np.percentile(Grcz,prsc)
                        rscl = np.append(rscl,rsc)
                for nrs in range(int(dr/dl)+1):
                        rcyls = (nrad)*dr+nrs*dl
                        cutrs = (Grcz<rscl[nrs+1]) & (Grcz>rscl[nrs])
                        Gmcrs = Gmcz[cutrs]
                        Gzcrs = Gzcz[cutrs]
                        Gxcrs = Gxcz[cutrs]
                        Gycrs = Gycz[cutrs]
                        Gvzcrs = Gvzcz[cutrs]
                        Gvxcrs = Gvxcz[cutrs]
                        Gvycrs = Gvycz[cutrs]
                        cutrsa = (Grc<rscl[nrs+1]) & (Grc>rscl[nrs])
                        Gmcsa = Gmc[cutrsa]
                        Grhocsa = Grhoc[cutrsa]
                        Gxcsa = Gxc[cutrsa]
                        Gycsa = Gyc[cutrsa]
                        Gzcsa = Gzc[cutrsa]
                        rzscl = []
                        for nzs in range(int(cylz/dl)+2):
                                przsc = 100.0/(int(cylz/dl)+1)*nzs
                                if przsc>100.0:
                                        przsc=100.0
                                rzsc = np.percentile(Gzcrs,przsc)
                                rzscl = np.append(rzscl,rzsc)
                        for nzs in range(int(cylz/dl)+1):
                                cutzs = (Gzcrs<rzscl[nzs+1]) & (Gzcrs>rzscl[nzs])
                                Gmcrzs = Gmcrs[cutzs]
                                Gxcrzs = Gxcrs[cutzs]
                                Gycrzs = Gycrs[cutzs]
                                Gzcrzs = Gzcrs[cutzs]
                                Gvxcrzs = Gvxcrs[cutzs]
                                Gvycrzs = Gvycrs[cutzs]
                                Gvzcrzs = Gvzcrs[cutzs]
                                cutzsa = (Gzcsa<rzscl[nzs+1]) & (Gzcsa>rzscl[nzs])
                                Gmcsa = Gmc[cutrsa]
                                Grhocsa = Grhoc[cutrsa]
                                Gxcsa = Gxc[cutrsa]
                                Gycsa = Gyc[cutrsa]
                                Gzcsa = Gzc[cutrsa]
                                Gvxcsa = Gvxc[cutrsa]
                                Gvycsa = Gvyc[cutrsa]
                                Gvzcsa = Gvzc[cutrsa]                                
                                totnphi=2.0*np.pi*rcyls/dl
                                phi = np.arctan2(Gycrzs,Gxcrzs)
                                phi[phi<0] += 2.0*np.pi
                                phil=[]
                                phia = np.arctan2(Gycsa,Gxcsa)
                                phia[phia<0] += 2.0*np.pi
                                phial=[]
                                for nphi in range(int(totnphi)+2):
                                        pphi = 100.0*nphi/(int(totnphi)+1)
                                        if pphi>100.0:
                                                pphi=100.0
                                        phineed = np.percentile(phi,pphi)
                                        phil = np.append(phil, phineed)
                                        phianeed = np.percentile(phia,pphi)
                                        phial = np.append(phial, phianeed)                                        
                                #print 'phil', phil
                                for nphi in range(int(totnphi)+1):
                                        phis = nphi*2.0*np.pi/totnphi
                                        cutphi = (phi>phil[nphi]) & (phi<phil[nphi+1])
                                        Gms = Gmcrzs[cutphi]
                                        Gvxs = Gvxcrzs[cutphi]
                                        Gvys = Gvycrzs[cutphi]
                                        Gvzs = Gvzcrzs[cutphi]
                                        cutphia = (phia>phial[nphi]) & (phia<phial[nphi+1])
                                        Gmsa = Gmcsa[cutphia]
                                        Grhosa = Grhocsa[cutphia]
                                        Gvxsa = Gvxcsa[cutphia]
                                        Gvysa = Gvycsa[cutphia]
                                        Gvzsa = Gvzcsa[cutphia]
                                        Gvxsrel = Gvxs-np.average(Gvxs,weights=Gms)
                                        Gvysrel = Gvys-np.average(Gvys,weights=Gms)
                                        Gvzsrel = Gvzs-np.average(Gvzs,weights=Gms)
                                        logGv2s = np.log10(Gvxsrel*Gvxsrel+Gvysrel*Gvysrel+Gvzsrel*Gvzsrel)
                                        if len(Gms)<2:
                                                Gmcell = Gms
                                                Gvxcell = Gvxsrel
                                                Gvycell = Gvysrel
                                                Gvzcell = Gvzsrel
                                        else:
                                                logv2sigma = weighted_quantile(logGv2s,fraccut,sample_weight=Gms)
                                                cutv2 = logGv2s<logv2sigma ###
                                                Gmcell = Gms[cutv2]
                                                Gvxcell = Gvxsrel[cutv2]
                                                Gvycell = Gvysrel[cutv2]
                                                Gvzcell = Gvzsrel[cutv2]
                                        turenergy = np.sum(0.5*Gmcell*(Gvxcell*Gvxcell+Gvycell*Gvycell+Gvzcell*Gvzcell)) #in erg
                                        turl = np.append(turl,turenergy)
                                        gml = np.append(gml,np.sum(Gmcell))
                                        Omegal = np.append(Omegal,Omega)
                                        gvl = np.append(gvl,np.sum(Gmsa/Grhosa))
                                        xcell = np.append(xcell,0.5*(rscl[nrs+1]+rscl[nrs])*np.cos(0.5*(phil[nphi]+phil[nphi+1])))
                                        ycell = np.append(ycell,0.5*(rscl[nrs+1]+rscl[nrs])*np.sin(0.5*(phil[nphi]+phil[nphi+1])))
                                        zcell = np.append(zcell,0.5*(rzscl[nzs+1]+rzscl[nzs]))
                                        
        Gxcutz = Gx[cutz]
        Gycutz = Gy[cutz]
        Gzcutz = Gz[cutz]
        if havecr>0:
                cregyl = cregy[cutz]
        therml = GEint[cutz]
        voll = Gvol[cutz]
        if haveB>0:
                Begyl = Begy[cutz]
                Begy_2Bzl = Begy_2Bz[cutz]
        Evzl = Evz[cutz]
        Gmcutz = Gm[cutz]
        Grxycz = np.sqrt(Gxcutz*Gxcutz+Gycutz*Gycutz)
        rxycell = np.sqrt(xcell*xcell+ycell*ycell)
        print 'np.amax(radl)', np.amax(radl)
        print 'np.amax(rscl)', np.amax(rscl)
        print 'np.amax(xcell)', np.amax(xcell)
        print 'np.amax(rxycell)', np.amax(rxycell)
        cutr = (Grxycz<cylr) & (Grxycz>cylrin)
        Gmcutzr = Gmcutz[cutr]
        print 'Gmcutzr', np.sum(Gmcutzr)
        if haveB>0:
                Begyl=Begyl[cutr]
                Begy_2Bzl = Begy_2Bzl[cutr]
        Evzl=Evzl[cutr]
        therml=therml[cutr]
        voll=voll[cutr]
        if havecr>0:
                cregyl=cregyl[cutr]
        Gxpl=Gxcutz[cutr]
        Gypl=Gycutz[cutr]
        Gzpl=Gzcutz[cutr]
        if cosmo==1:
                readtimelist=readtime(firever=2)
                snap2list=readtimelist['snaplist']
                time2list=readtimelist['timelist']
                a2list=readtimelist['alist']
                tnow = np.interp(ascale,a2list,time2list)*1e9
        if cosmo==1:
                timen = tnow
        else:
                timen = float(Nsnap)*0.98*1e6
        return {'turl':turl,'therml':therml,'Begyl':Begyl,'cregyl':cregyl,\
                'Evzl':Evzl,'Begy_2Bzl':Begy_2Bzl, 'voll':voll,\
        'cylrin':cylrin,'cylr':cylr,'cylz':cylz,'nbin':nbin,\
        'Gxpl':Gxpl, 'Gypl':Gypl, 'Gzpl':Gzpl, 'xcell':xcell,\
         'ycell':ycell, 'zcell':zcell, 'Omegal':Omegal,'gminl':Gmcutzr,
                'gvcell':gvl,'gml':gml, 'timen':timen, 'G':G}

    
    
    
    
def caltur_Su():
        radl=[]
        for nrad in range(nbin+1):
                radl = np.append(radl,cylrin+nrad*dr)
        for nrad in range(nbin):
                cutr = (Grxy<radl[nrad+1]) & (Grxy>radl[nrad])
                #print 'radl[nrad+1]', radl[nrad+1]
                #print 'radl[nrad]', radl[nrad]
                cut=cutr*cutz
                Gmc=Gm[cut]
                Grhoc=Grho[cut]
                Grc=Grxy[cut]
                Gxc=Gx[cut]
                Gyc=Gy[cut]
                Gzc=Gz[cut]
                Gvxc = Gvx[cut]
                Gvyc = Gvy[cut]
                Gvzc = Gvz[cut]
                #print 'len(Gmc)', len(Gmc)
                Gvzcave = np.average(Gvz[cut],weights=Gmc)
                #calculate average rotational velocity:
                #theta hat dot v = vtheta
                Grxyc = np.sqrt(Gxc*Gxc+Gyc*Gyc)
                vth = -Gyc/Grxyc*Gvxc+Gxc/Grxyc*Gvyc
                vthave = np.average(vth,weights=Gmc)
                Omega = vthave/(radl[nrad+1]+radl[nrad])*2.0/kpc_in_cm
                #print 'Gvxc', Gvxc
                Gvxc = Gvxc+Gyc/Grxyc*vthave
                Gvyc = Gvyc-Gxc/Grxyc*vthave
                vzlog = np.log10(np.absolute(Gvzc-Gvzcave))
                logvzsigma = weighted_quantile(vzlog,fraccut,sample_weight=Gmc)
                cutvz = np.log10(vzlog)<logvzsigma ###
                Grcz = Grc[cutvz]
                Gmcz = Gmc[cutvz]
                Grhoz = Grhoc[cutvz]
                Gzcz = Gzc[cutvz]
                Gxcz = Gxc[cutvz]
                Gycz = Gyc[cutvz]
                Gvzcz = Gvzc[cutvz]
                Gvxcz = Gvxc[cutvz]
                Gvycz = Gvyc[cutvz]
                totnum = len(Gmcz)
                print 'totnum', totnum
                vol = np.pi*(radl[nrad+1]*radl[nrad+1]-radl[nrad]*radl[nrad])*cylz
                numden = totnum/vol
                vol15 = 15.0/numden #TK test
                dl = np.power(vol15,1./3.)
                print 'dl', dl
                #print 'int(dr/dl)+1', int(cylz/dl)+1
                numcount=0
                rscl = []
                for nrs in range(int(dr/dl)+2):
                        prsc = 100.0/(int(dr/dl)+1)*nrs
                        if prsc>100.0:
                                prsc =100.0
                        rsc = np.percentile(Grcz,prsc)
                        rscl = np.append(rscl,rsc)
                for nrs in range(int(dr/dl)+1):
                        rcyls = (nrad)*dr+nrs*dl
                        cutrs = (Grcz<rscl[nrs+1]) & (Grcz>rscl[nrs])
                        Gmcrs = Gmcz[cutrs]
                        Gzcrs = Gzcz[cutrs]
                        Gxcrs = Gxcz[cutrs]
                        Gycrs = Gycz[cutrs]
                        Gvzcrs = Gvzcz[cutrs]
                        Gvxcrs = Gvxcz[cutrs]
                        Gvycrs = Gvycz[cutrs]
                        cutrsa = (Grc<rscl[nrs+1]) & (Grc>rscl[nrs])
                        Gmcsa = Gmc[cutrsa]
                        Grhocsa = Grhoc[cutrsa]
                        Gxcsa = Gxc[cutrsa]
                        Gycsa = Gyc[cutrsa]
                        Gzcsa = Gzc[cutrsa]
                        rzscl = []
                        for nzs in range(int(cylz/dl)+2):
                                przsc = 100.0/(int(cylz/dl)+1)*nzs
                                if przsc>100.0:
                                        przsc=100.0
                                rzsc = np.percentile(Gzcrs,przsc)
                                rzscl = np.append(rzscl,rzsc)
                        for nzs in range(int(cylz/dl)+1):
                                cutzs = (Gzcrs<rzscl[nzs+1]) & (Gzcrs>rzscl[nzs])
                                Gmcrzs = Gmcrs[cutzs]
                                Gxcrzs = Gxcrs[cutzs]
                                Gycrzs = Gycrs[cutzs]
                                Gzcrzs = Gzcrs[cutzs]
                                Gvxcrzs = Gvxcrs[cutzs]
                                Gvycrzs = Gvycrs[cutzs]
                                Gvzcrzs = Gvzcrs[cutzs]
                                cutzsa = (Gzcsa<rzscl[nzs+1]) & (Gzcsa>rzscl[nzs])
                                Gmcsa = Gmc[cutrsa]
                                Grhocsa = Grhoc[cutrsa]
                                Gxcsa = Gxc[cutrsa]
                                Gycsa = Gyc[cutrsa]
                                Gzcsa = Gzc[cutrsa]
                                Gvxcsa = Gvxc[cutrsa]
                                Gvycsa = Gvyc[cutrsa]
                                Gvzcsa = Gvzc[cutrsa]                                
                                totnphi=2.0*np.pi*rcyls/dl
                                phi = np.arctan2(Gycrzs,Gxcrzs)
                                phi[phi<0] += 2.0*np.pi
                                phil=[]
                                phia = np.arctan2(Gycsa,Gxcsa)
                                phia[phia<0] += 2.0*np.pi
                                phial=[]
                                for nphi in range(int(totnphi)+2):
                                        pphi = 100.0*nphi/(int(totnphi)+1)
                                        if pphi>100.0:
                                                pphi=100.0
                                        phineed = np.percentile(phi,pphi)
                                        phil = np.append(phil, phineed)
                                        phianeed = np.percentile(phia,pphi)
                                        phial = np.append(phial, phianeed)                                        
                                #print 'phil', phil
                                for nphi in range(int(totnphi)+1):
                                        phis = nphi*2.0*np.pi/totnphi
                                        cutphi = (phi>phil[nphi]) & (phi<phil[nphi+1])
                                        Gms = Gmcrzs[cutphi]
                                        Gvxs = Gvxcrzs[cutphi]
                                        Gvys = Gvycrzs[cutphi]
                                        Gvzs = Gvzcrzs[cutphi]
                                        cutphia = (phia>phial[nphi]) & (phia<phial[nphi+1])
                                        Gmsa = Gmcsa[cutphia]
                                        Grhosa = Grhocsa[cutphia]
                                        Gvxsa = Gvxcsa[cutphia]
                                        Gvysa = Gvycsa[cutphia]
                                        Gvzsa = Gvzcsa[cutphia]
                                        Gvxsrel = Gvxs-np.average(Gvxs,weights=Gms)
                                        Gvysrel = Gvys-np.average(Gvys,weights=Gms)
                                        Gvzsrel = Gvzs-np.average(Gvzs,weights=Gms)
                                        logGv2s = np.log10(Gvxsrel*Gvxsrel+Gvysrel*Gvysrel+Gvzsrel*Gvzsrel)
                                        if len(Gms)<2:
                                                Gmcell = Gms
                                                Gvxcell = Gvxsrel
                                                Gvycell = Gvysrel
                                                Gvzcell = Gvzsrel
                                        else:
                                                logv2sigma = weighted_quantile(logGv2s,fraccut,sample_weight=Gms)
                                                cutv2 = logGv2s<logv2sigma ###
                                                Gmcell = Gms[cutv2]
                                                Gvxcell = Gvxsrel[cutv2]
                                                Gvycell = Gvysrel[cutv2]
                                                Gvzcell = Gvzsrel[cutv2]
                                        turenergy = np.sum(0.5*Gmcell*(Gvxcell*Gvxcell+Gvycell*Gvycell+Gvzcell*Gvzcell)) #in erg
                                        turl = np.append(turl,turenergy)
                                        gml = np.append(gml,np.sum(Gmcell))
                                        Omegal = np.append(Omegal,Omega)
                                        gvl = np.append(gvl,np.sum(Gmsa/Grhosa))
                                        xcell = np.append(xcell,0.5*(rscl[nrs+1]+rscl[nrs])*np.cos(0.5*(phil[nphi]+phil[nphi+1])))
                                        ycell = np.append(ycell,0.5*(rscl[nrs+1]+rscl[nrs])*np.sin(0.5*(phil[nphi]+phil[nphi+1])))
                                        zcell = np.append(zcell,0.5*(rzscl[nzs+1]+rzscl[nzs]))
                                        
    
    
    


def calasyn(p):
    from scipy.special import gamma
    asyn= np.sqrt(np.pi)/2.*gamma(p/4.+19./12.)*gamma(p/4.-1./12.)*gamma(p/4.+5./4.)/(p+1.)/gamma(p/4.+7./4.)
    return asyn
    

def calsynke(ecre,p,Emax,Emin):
    # ecre = cosmic ray electron density in GeV/m^3
    # Emin #Pfrommer+04 used 0.1GeV; Longair used 
    # We should adapt 1 GeV such that the ratio of CR proton to electron energy is aorund 100 (see Beck & Krause 2005)
    # p is the cosmic ray electron spectral index (assume to be 2.5)
    # N(E) = ke*E^(-p)
    # or ke = Int( E N(E) dE )/Int( E^(-p+1) dE) = (2-p)*ecre/(Emax^(2-p)-Emin^(2-p))
    ke = (2.-p)*ecre/(np.power(Emax,2.-p)-np.power(Emin,2.-p)) #ke in the above unit (GeV based), depending on p
    return ke
    
def calsynchrotron(alpha,B,ke,nu,V):
    #alpha = synchrotron spectral index (assume to be 0.75)
    #B = B field strength in Tesla
    #p = electron spectral index
    #nu = synchrotron frequency in Hz 
    #V = m^3
    #energy of electron measured in GeV; i.e. N(E) in GeV^-1 m^-1
    # N(E) = ke*E^(-p), or ke = Int( E N(E) dE )/Int( E^(-p+1) dE) 
    p=2.*alpha+1.
    asyn = calasyn(p)
    bigA = 2.344e-25*np.power(3.217e17,alpha)*asyn
    radiolum_in_W_Hz = bigA*V*ke*np.power(B,1+alpha)*np.power(nu,-alpha)
    return radiolum_in_W_Hz


def calBmin(eta,radiolum_in_W_Hz,V,nu):
    #in SI unit
    # only for alpha (synchrotron spectral index) = 0.75 
    # and assuming numin = nu and neglecting numax in the calculation
    # see Longair High energy astrophysics 3rd edition
    #fnu
    Bmin_in_T = 1.8*np.power(eta*radiolum_in_W_Hz/V,2./7.)*np.power(nu,1./7.)
    Bmin_in_G = Bmin_in_T*Tesla_in_Gauss
    return Bmin_in_G


def calsyndefault(B,ecre,V,nu=1e6,alpha=0.75,Emax=100,Emin=1.0):
    #in SI unit, except ecre in GeV/m^3
    p = 2.*alpha+1.
    asyn = calasyn(p)
    ke = calsynke(ecre,p,Emax,Emin)
    radiolum_in_W_Hz = calsynchrotron(alpha,B,ke,nu,V)
    return radiolum_in_W_Hz

def convertsyncetoHz(apara,B,reverse=0):
    # calculate the maximum intensity in synchrotron radiation (see Longair High energy astrophysics 3rd edition)
    #if reverse=1; input nu in Hz and output Eelectron in GeV
    #input eelectron in GeV
    #B in Tesla
    Cprop = 0.29*3./2.*echarge_in_C/2./np.pi/electronmass_in_kg\
    /electronmass_in_kg/electronmass_in_kg/cspeed_in_m_s/cspeed_in_m_s/cspeed_in_m_s/cspeed_in_m_s
    if reverse==0:
        Eelectron=apara
        Eelectron=Eelectron/Joule_in_GeV #convert to SI unit
        nu = Eelectron*Eelectron*Cprop*B #nu is the frequency of synchrotron in Hz (with the maximum intensity)
        outpara = nu
    if reverse==1:
        nu=apara
        Eelectron = np.sqrt(nu/Cprop/B) #Eelectron in SI
        Eelectron=Eelectron*Joule_in_GeV #convert to GeV
        outpara = Eelectron       
    return outpara



def pressureXYZ(G, xlist, ylist, zlist, dx, dy, dz,havecr=0,haveB=0): #calculate density with r and z (cylindrical) coordinate 
    #xlist, ylist and zlist must have the same size; (x,y,z) is the coordinate
    Gp = G['p']; Gm = G['m']*1e10*Msun_in_g; #now cgs
    Gx = Gp[:,0]; Gy = Gp[:,1]; Gz = Gp[:,2];  #now kpc
    Gvel = G['v']*km_in_cm #now in cm/s
    Gvx = Gvel[:,0]; Gvy = Gvel[:,1]; Gvz = Gvel[:,2]    
    Grho = G['rho']*1e10*solar_mass_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm
    GEint = G['u']*km_in_cm*km_in_cm*Gm #in erg #thermal energy
    Gvol = Gm/Grho
    if havecr>0: cregy = G['cregy']*1e10*solar_mass_in_g*km_in_cm*km_in_cm #in erg
    if haveB>0:
            GB = G['B']
            Bx = GB[:,0]
            By = GB[:,1]
            Bz = GB[:,2]
            B2 = Bx*Bx+By*By+Bz*Bz
            Begy = B2/np.pi/8.*Gm/Grho
            Begy_2Bz = Begy-Bz*Bz/np.pi/4.*Gm/Grho    
    rhol = np.array([])
    pthl = np.array([])
    pturl = np.array([])
    pcrl = np.array([])
    pBl = np.array([])
    for (x,y,z) in zip(xlist,ylist,zlist):
        cutx = np.logical_and(Gx>x-dx/2.,Gx<x+dx/2.)
        cuty = np.logical_and(Gy>y-dy/2.,Gy<y+dy/2.)
        cutz = np.logical_and(Gz>z-dz/2.,Gz<z+dz/2.)
        volume = dx*dy*dz*kpc_in_cm*kpc_in_cm*kpc_in_cm
        cut = cutx*cuty*cutz
        Gmcut = np.sum(Gm[cut])
        Gvolcut = np.sum(Gvol[cut])
        Gthcut = np.sum(GEint[cut])
        pth = Gthcut/Gvolcut*(GAMMA-1.0)
        vzave = np.average(Gvz[cut])
        Gztur = np.sum(Gm[cut]*np.square(Gvz[cut]-vzave))
        ptur = Gztur/Gvolcut
        if havecr>0: 
            Gcregycut = np.sum(cregy[cut])
            pcr = Gcregycut/Gvolcut*(CRgamma-1.0)
            pcrl = np.append(pcrl,pcr)
        if haveB>0: 
            Begy_2Bzcut = np.sum(Begy_2Bz[cut])
            pB = np.sum(Begy_2Bz[cut])/Gvolcut
            pBl = np.append(pBl,pB)
        rho = Gmcut/volume
        rhol = np.append(rhol,rho)
        pthl = np.append(pthl,pth)
        pturl = np.append(pturl,ptur)
    return {'rhol':rhol, 'pthl':pthl, 'pturl':pturl, 'pcrl':pcrl, 'pBl':pBl}