from samson_const import *
import matplotlib as mpl
from pathloc import *
mpl.use('Agg')
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
from textwrap import wrap
#rcParams['figure.figsize'] = 5, 5
rcParams['figure.figsize'] = 10, 5
rcParams['font.size']=18
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
#dirneed=['bwsbclrdc28']
#dirneed=['m12imhdcv','m12icr_b_70','m12icr_700','m12icr_70va']
#dirneed=['m12mmhdcv','m12mcr_b_70','m12mcr_700','m12mcr_70va']
#dirneed=['m11dmhdcv','m11dcr_b_70','m11dcr_700','m11dcr_70va']
#dirneed = ['bwmwlrstr','bwmwlrstrc1000']
#dirneed=['bwmwlrdc27ds','bwmwlrdc27']
#dirneed=['mw_mr_2_21','bwmwmrdc0','bwmwmrdc28']
#dirneed=['m09mhdcv','m09cr_b_70','m10qmhdcv','m10qcr_b_70','m11bmhdcv','m11bcr_b_70','m11fmhdcv','m11fcr_b_70','m12imhdcv','m12icr_b_70']
#dirneed=['m10vmhdcv','m10vcr_b_70', 'm10vcr_700']
#dirneed=['m11gcr_b_70','m11gcr_700','m11gmhdcv','f46']
#dirneed=['m12icr_b_70','m12icr_700','m12imhdcv','fm12i']
#dirneed=['m12icr_b_70']
#dirneed=['m12mcr_b_70']
#dirneed=['m11bcr_700']
#dirneed=['bwsmclrdc28']
dirneed=['bwmwmrdc28']
#dirneed=['m12mcr_700']
#startno=10
#startno=30
#startno=400
#startno=480
#startno=300
#startno=350
#startno=100
#startno=200
#startno=490
#startno=400
#startno=570
#startno=490
#startno=500
#startno=11
#startno=300
#startno=440
#startno=600
#startno=10
startno=490
#startno=0
#startno=480
#startno=450
#startno=506
#startno=200
#startno=350
#startno=499
#startno=600
#startno=550
#startno=590
#Nsnap=10
#Nsnap=50
#Nsnap=15
#Nsnap=20
#fmeat='mwmr'
#fmeat='mwdc27ds'
#fmeat='smc10kpc'
#fmeat='smc'
#fmeat='sbc'
#fmeat='M1str_'+str(Nsnap)
#fmeat='dc27g_sn'+str(Nsnap)
#fmeat='M1dc29g_2_21_smallts_sn'+str(Nsnap)
#fmeat='M1dc27g_2_21_smallts_sn'+str(Nsnap)
#fmeat='M1dc27_2_21_smallts_sn'+str(Nsnap)
#Nsnap=100
#Nsnap=601
#Nsnap=999
#Nsnap=200
#Nsnap=0
#Nsnap=210
#Nsnap=250
#Nsnap=300
#Nsnap=460
Nsnap=500
#Nsnap=651
#Nsnap=10
#Nsnap=801
#Nsnap=430
#Nsnap=600
#Nsnap=650
#Nsnap=590
#Nsnap=10
#Nsnap=560
#Nsnap=506
#Nsnap=200
#Nsnap=15
#Nsnap=10
#Nsnap=30
#Nsnap=440
#Nsnap=1000
#Nsnap=250
#Nsnap=400
#Nsnap=601
#Nsnap=300
#Nsnap=700
#Nsnap=500
#Nsnap=506
#Nsnap=518
#Nsnap=8
#Nsnap=555
#Nsnap=600
#Nsnap=1
#Nsnap=161
#Nsnap=160
#Nsnap=580
#Nsnap=30
#Nsnap=10
snapsep=10
#snapsep=1
#snapsep=5
#snapsep=490
#snapsep=20
fmeat='m12m'
#wanted='phaseTwindz'
#wanted='Twindbar'
#wanted='Tbar'
#wanted='phaseTz'
#wanted='rhoT'
#wanted='rhoTwind'
#wanted='Tzmwed'
#wanted='vzz'
#wanted='vzztrack'
#wanted='Tvz'
#wanted='Tvztrack'
#wanted='testCRcum'
#wanted='testCRdis'
#wanted='Tz'
#wanted='rhoztrack'
#wanted='Tztrack'
#wanted='rhoz'
#wanted = 'nismB'
#wanted ='dpcrz_rhog'
#wanted = 'dpcrz'
#wanted = 'massloadingSasha'
#wanted = 'outflowSasha'
#wanted = 'massloadingBooth'
#wanted = 'outflowBooth'
#wanted = 'gassurden'
#wanted = 'gassurdentime'
#wanted = 'outflowall'
#wanted = 'gasdenmidplane'
#wanted = 'outflow'
#wanted = 'sfrv'
#wanted = 'cumnsmrad'
#wanted = 'nsrad'
#wanted = 'sfrrad'
#wanted = 'sfrarad'
#wanted='sncretime'
#wanted='crdenplanes'
#wanted='crdenv'
#wanted='crdenmidplane'
#wanted='gammadensph'
#wanted='credensph'
#wanted='decayratiosph'
#wanted='decaytimesph'
#wanted='gammasph'
#wanted='Gmencsph'
#wanted='gasdensph'
#wanted='gasden'
#wanted='avecrdenr'
#wanted='avekedenr'
#wanted='avethdenr'
#wanted='aveedenr'
#wanted='denrtime'
#wanted='avedenr'
#wanted = 'FpFgr'
#wanted='gmr'
#wanted = 'rhovr'
#wanted='denr'
#wanted = 'rhovr_denr'
#wanted='crdenr'
#wanted='crer'
#wanted='crdpr'
#wanted='crecumz'
wanted='nismcrad'
#wanted='crecumr'
#wanted='Begyz'
#wanted='kez'
#wanted='ethz'
#wanted='vz'
#wanted='pz'
#wanted='dpz'
#wanted='gammar'
#wanted='gammaz'
#wanted='cresph'
#wanted='gmz'
#wanted='gasdenz'
#wanted='gasdeny'
#wanted='gasdenx'
#wanted='nismcrg'
#wanted = 'rcrad'
#wanted='nismcrl'
#wanted='nismcrad'
#wanted='nismcre'
#wanted='nismgamma'
#wanted='nismcumgamma'
#wanted='Begytime'
#wanted='Brmstime'
#wanted='Begy_mtime'
#wanted='Bcumgm'
#wanted='nismcumgm'
#wanted='gammadecaytime'
#wanted='gammacumdecaytime'
#wanted='credecaytime'
#wanted='crecumdecaytime'
#wanted='nismcumcre'
#wanted='pvolcumcre'
#wanted='pdxcumcre'
#wanted='pdxcumN'
#wanted='pvolcumN'
#wanted = 'crasnap'
#wanted='testCRcum'
#wanted='rhoT'
#wanted='cratime'
#wanted='gasdenv'
#wanted='crdrad'
#wanted='cramap'
#wanted='cramapv'
#wanted='Vquiverxy'
#wanted='Vquiverxz'
#wanted='Bquiverxy'
#wanted='Bquiverxz'
#wanted='Smassxy'
#wanted='Smassxz'
#wanted='testquiver'
#wanted='Twindm'
#wanted='vzwindm'
#wanted='Tz'
#wanted='pcrpth'
#wanted='Alfven_sound'
#wanted='cramapv'
#wanted='crturmap'
#wanted='outflowwrite'
#wanted = 'gassurden'
#wanted='gasfrac'
#wanted='printparmass'
#wanted='enchange'
#wanted='dirssfr10_200Myr'
#wanted='dirage10_200Myr'
#wanted='diragehis'
#wanted='dirage200Myr'
#wanted='diragestd'
#wanted='diragestd10_200Myr'
#wanted='diragemax'
#wanted='dirsm'
#wanted='dirsfr'
#wanted='testcr'
#wanted='dirheader'
#wanted='dirgammasfr'
#wanted='dirgamma'
#wanted='cosmic_raywrite'
#wanted='crdensol' #cosmic ray energy density at solar circle (observation ~ 1 eV/cm^3)
#wanted='smchange'
#wanted='cumgamma'
#wanted='smass'
#wanted='crdtime'
#wanted='gaskechange'
#wanted='gaske' #evolution of total gas kinetic energy 
#wanted='crchange'
#wanted='crtime'
#wanted='crcooling'
#wanted='lgratio'
#wanted='dg_sfr'
#wanted='gamma_partit'
#wanted='gamma_aveCR'
#wanted='dlgratio'
#wanted='gsmratio'
#title='LSG'
title='MW'
#title='Cosmo'
#title='Dwarf'
#titleneed='Dwarf'
#title='SBG'
titleneed=title
needlog=0
dclabelneed=1
correctIa=0
useM1=1
the_prefix='snapshot'
the_suffix='.hdf5'
withinr=15.0
nogrid=40
maxlength=10.0
med = -0.1 
wholeboxgas=1
diskgas=1
Rfrac = 0.5
nosum=0
xaxis_snapno=0

print 'wanted', wanted
print 'fmeat', fmeat
print 'runtodo', dirneed

if wanted=='nismgamma' or wanted=='nismcre' or wanted=='nismcrad' or wanted=='nismcrg' or wanted=='nismcrl' or wanted == 'rcrad':
    rcParams['figure.figsize'] = 5,4
    nogrid=15
    Rvcut = 0
    atstarburst=0
    trackneed=1
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
