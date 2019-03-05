from samson_const import *
from pathloc import *
import matplotlib as mpl
mpl.use('Agg')
from readsnap_cr import readsnapcr
import Sasha_functions as SF
import graphics_library as GL
import gas_temperature as GT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from samson_functions import *
from crtestfunction import *
from matplotlib import rcParams
from pylab import *
from textwrap import wrap
from scipy.optimize import curve_fit
#rcParams['figure.figsize'] = 5, 5
rcParams['figure.figsize'] = 10, 5
rcParams['font.size']=12
rcParams['font.family']='serif'
rcParams['text.usetex']=True
#rcParams.update({'figure.autolayout': True})
import matplotlib.patches as patches
rcParams['axes.linewidth'] = 2
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True
rcParams['axes.unicode_minus']=False
colortable = [ 'b', 'g', 'r']
#dirneed=['bwmwmr','bwmwmrdc27','bwmwmrdc28', 'bwmwmrdc29']
dirneed=['bwmwmrdc29']
#dirneed=['bwsmclrdc28']
#dirneed=['bwmwmrdc29']
#dirneed=['m11gcr_700']
#dirneed=['bwsbclrdc27']
#dirneed=['m11dcr_b_70']
#wanted='rhog'
wanted='nismcre'
startno=400
Nsnap=500
snapsep=10
the_prefix='snapshot'
the_suffix='.hdf5'
fmeat=dirneed[-1]


if wanted=='nismgamma' or wanted=='nismcre' or wanted=='nismcrad' or wanted=='nismcrg' or wanted=='nismcrl' or wanted == 'rcrad':
    rcParams['figure.figsize'] = 5,4
    nogrid=15
    Rvcut = 0
    atstarburst=0
    trackneed=0
    newlabelneed=1
    for runtodo in dirneed:
        timel=[]
        #enllist=[]
        #precount=0
        info=outdirname(runtodo, Nsnap)
        havecr=info['havecr']
        if havecr==0:
            continue
        if atstarburst==1:
            if runtodo=='bwsbclr':
                    Nsnap=523
            if runtodo=='bwsbclrmhd':
                    Nsnap=589
            if runtodo=='bwsbclrdc0':
                    Nsnap=590
            if runtodo=='bwsbclrdc27':
                    Nsnap=138
            if runtodo=='bwsbclrdc28':
                    Nsnap=520
            if runtodo=='bwsbclrdc29':
                    Nsnap=433
            if runtodo=='bwsbclrstr':
                    Nsnap=245
            if runtodo=='bwsbclrdc28mhd':
                    Nsnap=558
            if runtodo=='bwsbclrdc28str':
                    Nsnap=370
        for i in [Nsnap]:
            info=outdirname(runtodo, i)
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
            cosmo=info['cosmo']
            color=info['color']
            haveB=info['haveB']
            Rvir = info['Rvir']
            newlabel=info['newlabel']
            if runtitle=='SMC':
                labelneed=dclabel
            if newlabelneed==1:
                labelneed="\n".join(wrap(newlabel,17))
            else:
                labelneed=''
            ptitle=title
            if runtitle=='SMC':
                    ptitle='Dwarf'
            elif runtitle=='SBC':
                    ptitle='Starburst'
            elif runtitle=='MW':
                    ptitle=r'$L\star$ Galaxy'
            print 'the_snapdir', the_snapdir
            print 'Nsnapstring', Nsnapstring
            print 'havecr', havecr
            G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
            header = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix,\
             extension=the_suffix, havecr=havecr, header_only=1)
            print 'time', header['time']
            timeneed=i*0.98*1e6 #in yr
            Grho_t = G['rho'] #1e10Msun per kpc^3
            Gid_t = G['id']
            cregy_codeunit = G['cregy']
            cregy = cregy_codeunit*1e10*solar_mass_in_g*km_in_cm*km_in_cm #original cosmic ray energy in 1e10Msun km^2/sec^2
            Neb = G['ne']
            cregyl = G['cregyl']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
            if wanted == 'nismcrad' or wanted == 'rcrad':
                cregyad_t = G['cregya']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
            if wanted == 'nismcrg':
                cregyg_t = G['cregyg']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
            Gp_t = G['p']
            Gx_t = Gp_t[:,0]
            Gy_t = Gp_t[:,1]
            Gz_t = Gp_t[:,2]
            Gr_t = np.sqrt(Gx_t*Gx_t+Gy_t*Gy_t+Gz_t*Gz_t)
            if Rvcut==1:
                Grtcut = Gr_t<1.0*Rvir
                Grho_t=Grho_t[Grtcut]
                cregy_codeunit=cregy_codeunit[Grtcut]
                cregy=cregy[Grtcut]
                Neb = Neb[Grtcut]
                cregyl = cregyl[Grtcut]
                Gid_t = Gid_t[Grtcut]
                Gr_t=Gr_t[Grtcut]
                if wanted == 'nismcrad' or wanted == 'rcrad':
                    cregyad_t = cregyad_t[Grcut]
                if wanted == 'nismcrg':
                    cregyg_t = cregyg_t[Grcut]
            if wanted == 'nismcrad' or wanted == 'nismcrg' or wanted == 'nismcrl' or wanted == 'rcrad':
                snaptrack=i-snapsep
                info=outdirname(runtodo, snaptrack)
                Nsnapstring = info['Nsnapstring']
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                timetrack=snaptrack*0.98*1e6
                if trackneed==1:
                    Gp= G['p']
                    Gx= Gp[:,0]
                    Gy= Gp[:,1]
                    Gz= Gp[:,2]
                    Gr= np.sqrt(Gx*Gx+Gy*Gy+Gz*Gz)
                    Grho = G['rho']
                    Gid = G['id']
                    Gid, unindex = np.unique(Gid, return_index=True)
                    Gid_t, unindex_t = np.unique(Gid_t, return_index=True)
                    idint0=np.in1d(Gid,Gid_t)
                    # there are some duplicate IDs, from particle splitting? I drop those repetitions.
                    Gid = Gid[idint0]
                    Gidinds = Gid.argsort()
                    Gid_tinds = Gid_t.argsort()
                    Grho = Grho[unindex]
                    Grho = Grho[idint0]
                    Grho = Grho[Gidinds]
                    Grho_t = Grho_t[unindex_t]
                    Grho_t = Grho_t[Gid_tinds]
                    Gr = Gr[unindex]
                    Gr = Gr[idint0]
                    Gr = Gr[Gidinds]
                    Gr_t = Gr_t[unindex_t]
                    Gr_t = Gr_t[Gid_tinds]
                if wanted == 'nismcrad' or wanted == 'rcrad':
                    cregyad = G['cregya']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                    if trackneed==1:
                        cregyad = cregyad[unindex]
                        cregyad = cregyad[idint0]
                        cregyad_t = cregyad_t[unindex_t]
                        dcrad = (cregyad_t[Gid_tinds]-cregyad[Gidinds])/(timeneed-timetrack)/yr_in_sec
                if wanted == 'nismcrg':
                    cregyg = G['cregyg']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                    cregyg = cregyg[idint0]
                    dcrg = (cregyg_t[Gid_tinds]-cregyg[Gidinds])/(timeneed-timetrack)/yr_in_sec
                if wanted == 'nismcrl':
                    cregyl_p = G['cregyl']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                    cregyl_p = cregyl_p[idint0]
                    dcrl = (cregyl[Gid_tinds]-cregyl_p[Gidinds])/(timeneed-timetrack)/yr_in_sec
            Gnism_in_cm_3 = 1.0/proton_mass_in_g*Grho_t*1e10*Msun_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm
            if wanted=='nismgamma':
                Lout = outLgamma_nism(Grho,Neb,cregy_codeunit)
                Lgamma = Lout['Lgamma']
            LogGnism = np.log10(Gnism_in_cm_3)
            LogGr = np.log10(Gr_t)
            if trackneed==1:
                LogGnxaxis = np.linspace(-6,2.5,num=nogrid)
            else:
                LogGnxaxis = np.linspace(-4,4.0,num=nogrid)
            if wanted == 'rcrad':
                LogGnxaxis = np.linspace(-1,2.0,num=nogrid)
            dx =LogGnxaxis[1]-LogGnxaxis[0]
            Lgammal = []
            crel=[]
            dcradl=[]
            dcrgl=[]
            dcrll=[]
            for inism in range(nogrid-1):
                if wanted == 'rcrad':
                    cutg = (LogGr > LogGnxaxis[inism]) & (LogGr < LogGnxaxis[inism+1])
                else:
                    cutg = (LogGnism > LogGnxaxis[inism]) & (LogGnism < LogGnxaxis[inism+1])
                if wanted=='nismgamma':
                    Lgammacut = Lgamma[cutg]
                    Lgammal = np.append(Lgammal, np.sum(Lgammacut))
                if wanted=='nismcre':
                    crecut = cregy[cutg]
                    crel = np.append(crel, np.sum(crecut))  
                if wanted=='nismcrad' or wanted == 'rcrad':
                    dcradcut = dcrad[cutg]
                    dcradl = np.append(dcradl, np.sum(dcradcut)) 
                if wanted=='nismcrg':  
                    dcrgcut = dcrg[cutg] 
                    dcrgl = np.append(dcrgl, np.sum(dcrgcut))
                if wanted=='nismcrl':
                    dcrlcut = dcrl[cutg]
                    dcrll = np.append(dcrll, np.sum(dcrlcut))
            if haveB>0:
                ls='dashed'
            else:
                ls='solid'
            if wanted == 'rcrad':
                plt.plot((LogGnxaxis[1:]+LogGnxaxis[:-1])/2.,dcradl/dx, label=labelneed, color=color,lw=2,ls=ls)
            if wanted=='nismgamma':
                plt.plot((LogGnxaxis[1:]+LogGnxaxis[:-1])/2.,Lgammal/dx, label=labelneed, color=color,lw=2,ls=ls)
            if wanted=='nismcre':
                plt.plot((LogGnxaxis[1:]+LogGnxaxis[:-1])/2.,crel/dx, label=labelneed, color=color,lw=2,ls=ls)
            if wanted=='nismcrad':
                    plt.plot((LogGnxaxis[1:]+LogGnxaxis[:-1])/2.,dcradl/dx, label=labelneed, color=color,lw=2,ls=ls)
            if wanted=='nismcrg':
                    plt.plot((LogGnxaxis[1:]+LogGnxaxis[:-1])/2.,dcrgl/dx, label=labelneed, color=color,lw=2,ls=ls)
            if wanted=='nismcrl':
                    plt.plot((LogGnxaxis[1:]+LogGnxaxis[:-1])/2.,dcrll/dx, label=labelneed, color=color,lw=2,ls=ls)
    if runtitle=='SMC':
        plt.legend(loc='best',fontsize=9,ncol=1)
    if trackneed==0:
        plt.yscale('log')
    if wanted=='rcrad':
            plt.ylabel(r'$\mathrm{d} \dot{E}_{\rm Ad}/\mathrm{d} \log r\;[{\rm erg/s}]$',fontsize=16)
            filename = plotloc+'CRplot/rcrad/rcrad_'+fmeat+'.pdf'
    if wanted=='nismgamma':
        plt.ylabel(r'$\mathrm{d} L_{\gamma}/\mathrm{d} \log (n)[{\rm erg/s}]$',fontsize=16)
        filename = plotloc+'CRplot/nismgamma/nismgamma_'+fmeat+'.pdf'
    if wanted=='nismcre':
        plt.ylim(ymin=1e48)
        plt.ylabel(r'$\mathrm{d} E_{\rm cr}/\mathrm{d} \log (n)[{\rm erg}]$',fontsize=16)
        filename = plotloc+'CRplot/nismcre/nismcre_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
    if wanted=='nismcrad':
            plt.ylabel(r'$\mathrm{d} \dot{E}_{\rm Ad}/\mathrm{d} \log (n)[{\rm erg/s}]$',fontsize=16)
            filename = plotloc+'CRplot/nismcrad/nismcrad_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
    if wanted=='nismcrg':
            plt.ylabel(r'$\mathrm{d} \dot{E}_{\rm SNe}/\mathrm{d} \log (n)[{\rm erg/s}]$',fontsize=16)
            filename = plotloc+'CRplot/nismcrg/nismcrg_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
    if wanted=='nismcrl':
            plt.ylabel(r'$\mathrm{d} \dot{E}_{\rm Loss}/\mathrm{d} \log (n)[{\rm erg/s}]$',fontsize=16)
            filename = plotloc+'CRplot/nismcrl/nismcrl_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
    print 'filename', filename
    if wanted=='rcrad':
        plt.xlabel(r'$\log (r\;[{\rm kpc}])$',fontsize=16)
    else:
        plt.xlabel(r'$\log (n[{\rm cm^{-3}}])$',fontsize=16)
    plt.title(ptitle,fontsize=21)
    plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
    plt.savefig(filename,bbox_inches='tight')
    plt.clf()
