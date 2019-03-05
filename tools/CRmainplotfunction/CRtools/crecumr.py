## from samson_const import *
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
wanted='crecumr'
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
#betapi = 0.7 ##should be between 0.5-0.9 from Lacki 2011 #fraction of pi has enough energy (>1GeV)
#yr_in_sec = 3.2e7
#nopi_per_gamma = 3.0
#solar_mass_in_g = 2e33
#km_in_cm = 1e5
#kpc_in_cm = 3.086e21
#erg_in_eV = 6.242e11
#cspeed_in_cm_s = 3e10
#proton_mass_in_g = 1.67e-24
#Kroupa_Lsf=6.2e-4
#Salpeter_Lsf=3.8e-4
#pidecay_fac = 2.0e5*250.0*yr_in_sec #over neff in cm_3 for pi decay time in s
#hadronicdecayrate = 5.86e-16 #from Guo 2008. The hadronic decay rate of CR
#coulombdecayrate = 1.65e-16 #Coulomb loss of CR


print 'wanted', wanted
print 'fmeat', fmeat
print 'runtodo', dirneed


if wanted=='gmz' or wanted=='crecumz' or wanted=='crecumr' or wanted=='Begyz' or wanted=='ethz' or wanted=='kez':
    rcParams['figure.figsize'] = 5,4
    withinr=100.0
    #withinr=20.
    #maxlength=10.0
    minlength=0.5
    maxlength=100.
    nogrid=500
    normalized =1
    newlabelneed=0
    strlabelneed=0
    legendneed=0
    massloss=0.3
    rotface=0
    for runtodo in dirneed:
        info=outdirname(runtodo, Nsnap)
        havecr = info['havecr']
        haveB = info['haveB']
        Nsnapstring = info['Nsnapstring']
        the_prefix = info['the_prefix']
        the_suffix = info['the_suffix']
        the_snapdir = info['the_snapdir']
        cosmo = info['cosmo']
        if wanted=='crecumz':
            if havecr==0:
                continue
        if wanted=='Begyz':
            if haveB==0:
                continue
        zlist=np.linspace(minlength,maxlength,num=nogrid)
        if wanted=='crecumz' or wanted=='crecumr':
            crel=zlist*0.
        if wanted=='Begyz':
            Begyl=zlist*0.
        if wanted=='ethz':
            ethl = zlist*0.
        if wanted=='kez':
            kel = zlist*0.
        Gm_in_sunl=zlist*0.
        nooftimes=0
        Esncr=1.
        if (wanted=='crecumz' or wanted=='crecumr') and normalized == 1:
            Nsnapstring0 = '000'
            S0 = readsnapcr(the_snapdir, Nsnapstring0, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
            S0mass = np.sum(S0['m']*1e10)
            print 'Nsnapstring', Nsnapstring
            S = readsnapcr(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
            Smass = np.sum(S['m']*1e10)
            Esncr = 0.1*(Smass-S0mass)*1e51*0.01*(1.+massloss) #total CR energy from SNII in erg
            print 'Esncr', Esncr
        for i in range(startno,Nsnap,snapsep):
                info=outdirname(runtodo, i)
                M1speed=info['M1speed']
                rundir=info['rundir']
                runtitle=info['runtitle']
                slabel=info['slabel']
                snlabel=info['snlabel']
                dclabel=info['dclabel']
                resolabel=info['resolabel']
                the_snapdir=info['the_snapdir']
                Nsnapstring=info['Nsnapstring']
                Fcal=info['Fcal']
                avesfr=info['iavesfr']
                timestep=info['timestep']
                haveB=info['haveB']
                havecr = info['havecr']
                color=info['color']
                newlabel=info['newlabel']
                strlabel=info['strlabel']
                ptitle=title
                if runtitle=='SMC':
                    ptitle='Dwarf'
                elif runtitle=='SBC':
                    ptitle='Starburst'
                elif runtitle=='MW':
                    ptitle=r'$L\star$ Galaxy'
                    labelneed=dclabel
                if newlabelneed==1:
                    labelneed="\n".join(wrap(newlabel,17)) 
                if strlabelneed==1:
                    labelneed="\n".join(wrap(strlabel,17)) 
                print 'labelneed,newlabel', labelneed,newlabel
                if cosmo==1:
                    datasup=0;
                else:
                    datasup=1;
                G = readsnapfromrun(runtodo,Nsnap,1,rotface=rotface,datasup=datasup)
                Gp = G['p']
                Gv = G['v']
                Gu = G['u']
                Gm = G['m']
                Grho = G['rho']
                if wanted=='crecumz' or wanted=='crecumr':
                    cregy = G['cregy']
                    print 'CR energy in erg', np.sum(cregy)*1e10*Msun_in_g*km_in_cm*km_in_cm
                if wanted=='Begyz':
                    Bfield = G['B']
                    Bx = Bfield[:,0]
                    By = Bfield[:,1]
                    Bz = Bfield[:,2]
                    B2 = Bx*Bx+By*By+Bz*Bz
                    Begy = B2/8./np.pi*Gm/Grho*kpc_in_cm*kpc_in_cm*kpc_in_cm
                if wanted=='kez':
                    Gvx = Gv[:,0]
                    Gvy = Gv[:,1]
                    Gvz = Gv[:,2]
                    ke = 0.5*Gm*(Gvx*Gvx+Gvy*Gvy+Gvz*Gvz)*1e10*Msun_in_g*km_in_cm*km_in_cm
                Gz = Gp[:,2]
                Gx = Gp[:,0]
                Gy = Gp[:,1]
                dz = maxlength/nogrid
                if normalized==1:
                    cutxy = Gx*Gx+Gy*Gy < withinr
                    cutz = np.absolute(Gz)<dz*nogrid
                    cut= cutxy*cutz
                    Gmtot = np.sum(Gm[cut]*1e10)
                for iz in range(nogrid):
                    if wanted=='crecumr':
                        cut = np.sqrt(Gx*Gx+Gy*Gy+Gz*Gz)<zlist[iz]
                    else:
                        cutxy = Gx*Gx+Gy*Gy < withinr
                        cutz = np.absolute(Gz)<zlist[iz]
                        cut = cutxy*cutz
                        Gmcut = Gm[cut]
                        Gm_in_sun = Gmcut*1e10
                    if wanted=='crecumz' or wanted=='crecumr':
                        cregycut = cregy[cut]*1e10*Msun_in_g*km_in_cm*km_in_cm
                        crel[iz] += np.sum(cregycut)
                    if wanted=='Begyz':
                        Begycut = Begy[cut]
                        Begyl[iz] += np.sum(Begycut)
                    if wanted=='ethz':
                        ethcut = Gu[cut]*km_in_cm*km_in_cm*Gm_in_sun*Msun_in_g
                        ethl[iz] += np.sum(ethcut)
                    if wanted=='kez':
                        kecut = ke[cut]
                        kel[iz] += np.sum(kecut)
                        Gm_in_sunl[iz] += np.sum(Gm_in_sun)
                nooftimes+=1
        if wanted=='crecumz' or wanted=='crecumr':
            crel=crel/nooftimes
        if wanted=='Begyz':
            Begyl=Begyl/nooftimes
        if wanted=='ethz':
            ethl=ethl/nooftimes
        if wanted=='kez':
            kel=kel/nooftimes
        Gm_in_sunl=Gm_in_sunl/nooftimes
        if haveB==0:
            lsn='solid'
        else:
            lsn='dashed'
        if wanted == 'gmz':
            if normalized==0:
                plt.plot(zlist, Gm_in_sunl, label=labelneed, ls=lsn,lw=2,color=color)
            else:
                plt.plot(zlist, Gm_in_sunl/Gmtot/(1.+massloss), label=labelneed, ls=lsn,lw=2,color=color)
        if wanted == 'crecumz' or wanted=='crecumr':
            plt.plot(zlist, crel/Esncr, label=labelneed, ls=lsn,lw=2,color=color)
        if wanted == 'Begyz':
            plt.plot(zlist, Begyl, label=labelneed, ls=lsn,lw=2,color=color)
            if wanted == 'ethz':
                        plt.plot(zlist, ethl, label=labelneed, ls=lsn,lw=2,color=color)
        if wanted == 'kez':
            plt.plot(zlist, kel, label=labelneed, ls=lsn,lw=2,color=color)
    plt.yscale('log')
    plt.xscale('log')
    if wanted=='crecumr':
                plt.xlabel(r'$r\; [{\rm kpc}]$',fontsize=18)
    else:
        plt.xlabel(r'$z\; [{\rm kpc}]$',fontsize=18)
    if wanted == 'gmz':
        if normalized==1:
            plt.ylabel(r'$M_{\rm gas}(<z)/M_0$',fontsize=18)
        else:
            plt.ylabel(r'$M_{\rm gas}(<z) [ M_\odot]$',fontsize=18)
        if runtitle=='SMC' and legendneed==1:
            #plt.legend(loc=2,fontsize=8,ncol=3,frameon=False)
            plt.legend(loc=4,fontsize=8,ncol=3)
    if wanted == 'crecumz' or wanted=='crecumr':
        if normalized==1:
            plt.ylabel(r'$E_{\rm cr}(<z)/E_{\rm cr,SN}(tot)$',fontsize=18)
        else:
            plt.ylabel(r'$E_{\rm cr}(<z)\;[{\rm erg}]$',fontsize=18)
        if wanted=='crecumr':
            if normalized==1:
                plt.ylabel(r'$E_{\rm cr}(<r)/E_{\rm cr,SN}(tot)$',fontsize=18)
            else:   
                plt.ylabel(r'$E_{\rm cr}(<r)\;[{\rm erg}]$',fontsize=18)
                if runtitle=='SMC' and legendneed==1:
                        #plt.legend(loc=2,fontsize=8,ncol=3,frameon=False)
                        plt.legend(loc=4,fontsize=8,ncol=3)
    if wanted == 'Begyz':
        plt.ylabel(r'$E_{\rm B}(<z)\;[{\rm erg}]$',fontsize=18)
        if wanted == 'ethz':
                plt.ylabel(r'$E_{\rm th}(<z)\;[{\rm erg}]$',fontsize=18)
        if wanted == 'kez':
                plt.ylabel(r'$E_{\rm KE}(<z)\;[{\rm erg}]$',fontsize=18)
    plt.tight_layout()
    if wanted == 'gmz':
        if normalized==0:
            filename = 'CRplot/gmz/gasmass_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        else:
            filename = 'CRplot/gmz/gasmass_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'_normal.pdf'
    if wanted == 'crecumz':
        filename = 'CRplot/crecumz/crecumz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    if wanted == 'crecumr':
                filename = 'CRplot/crecumr/crecumr_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    if wanted == 'Begyz':
                filename = 'CRplot/Begyz/Begyz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    if wanted == 'ethz':
                filename = 'CRplot/ethz/ethz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    if wanted == 'kez':
                filename = 'CRplot/kez/kez_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    print 'filename', filename
    plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
    plt.title(ptitle)
    plt.savefig(filename)
    plt.clf()
