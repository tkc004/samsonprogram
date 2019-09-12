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
import collections


def avedentestinput(subdict):
        dirneed=subdict['dirneed']
        wanted=subdict['wanted']
        print 'wanted', wanted
        startno=subdict['startno']
        Nsnap=subdict['Nsnap']
        snapsep=subdict['snapsep']
        the_prefix=subdict['the_prefix']
        the_suffix=subdict['the_suffix']
        fmeat=subdict['fmeat']


        print 'wanted', wanted
        print 'fmeat', fmeat
        print 'runtodo', dirneed



        if wanted=='gmr' or wanted=='avedenr' or wanted=='avecrdenr'\
        or wanted=='avethdenr' or wanted=='aveedenr' or wanted=='FpFgr'\
        or wanted=='avekedenr':
            plt.figure(figsize=(5,4))
            def funcab(x, a, b):
                    return a+b*x
            def func(x, a, b, c):
                    return a+b*x+c*x*x
            def d10func(r, a, b, c):
                    return np.power(10.,func(np.log10(r), a, b, c))*(b+2*c*np.log(r)/np.log(10.))/r
            nogrid=20
            withinr=2.0
            atstarburst=1
            linelabelneed=0
            if wanted=='aveedenr':
                linelabelneed=1
            labcount=0
            newlabelneed=subdict['newlabelneed']
            legendneed=subdict['legendneed']
            for runtodo in dirneed:
                if atstarburst==1:
                    if runtodo=='bwsbclr':
                            Nsnap=284
                            #Nsnap=523
                    if runtodo=='bwsbclrmhd':
                            Nsnap=600
                            #Nsnap=589
                    if runtodo=='bwsbclrdc0':
                            Nsnap=286
                            #Nsnap=590
                    if runtodo=='bwsbclrdc27':
                            Nsnap=369
                            #Nsnap=138
                    if runtodo=='bwsbclrdc28':
                            Nsnap=504
                            #Nsnap=520
                    if runtodo=='bwsbclrdc29':
                            Nsnap=430
                            #Nsnap=433
                    if runtodo=='bwsbclrstr':
                            Nsnap=358
                            #Nsnap=245
                    if runtodo=='bwsbclrdc28mhd':
                            Nsnap=433
                            #Nsnap=558
                    if runtodo=='bwsbclrdc28str':
                            Nsnap=532
                            #Nsnap=370
                for i in [Nsnap]:
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
                    firever=info['firever']
                    initsnap=info['initsnap']
                    haveB=info['haveB']
                    M1speed=info['M1speed']
                    Rvirguess=info['Rvir']
                    newlabel=info['newlabel']
                    snumadd=info['snumadd']
                    ptitle=title
                    labelneed=dclabel
                    if newlabelneed==1:
                            labelneed="\n".join(wrap(newlabel,17))
                if runtitle=='SMC':
                    ptitle='Dwarf'
                elif runtitle=='SBC':
                    ptitle='Starburst'
                elif runtitle=='MW':
                    ptitle=r'$L\star$ Galaxy'
                lsn='solid'
                if haveB>0:
                    lsn='dashed'
                if havecr==0 and wanted=='avecrdenr':
                    continue
                if cosmo==1:
                    h0=1
                    rotface=1
                    datasup=0;
                else:
                    h0=0
                    rotface=0
                    datasup=1;
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,cosmological=cosmo,h0=cosmo)
                Gp = G['p']
                Gv = G['v']
                Gu = G['u']
                Gm = G['m']
                Gvx = Gv[:,0]
                Gvy = Gv[:,1]
                Gvz = Gv[:,2]
                GEint = Gu*km_in_cm*km_in_cm*Gm*1e10*Msun_in_g
                cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                ke = 0.5*(Gvx*Gvx+Gvy*Gvy+Gvz*Gvz)*km_in_cm*km_in_cm*Gm*1e10*Msun_in_g  
                #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                Gz = Gp[:,2]
                Gx = Gp[:,0]
                Gy = Gp[:,1]
                dr = withinr/nogrid
                from crtestfunction import findcenz
                datasup=1
                xcen,ycen,zcen = findcennew(runtodo,Nsnap,withinr=2.5,datasup=datasup,Gx=Gx,Gy=Gy,Gz=Gz,Gm=Gm)
                Gz = Gz-zcen; Gy = Gy-ycen; Gx = Gx-xcen;
                GEint = Gu*km_in_cm*km_in_cm*Gm*1e10*Msun_in_g
                cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                ke = 0.5*(Gvx*Gvx+Gvy*Gvy+Gvz*Gvz)*km_in_cm*km_in_cm*Gm*1e10*Msun_in_g   
                Gr = np.sqrt(Gx*Gx+Gy*Gy+Gz*Gz)
                dr = withinr/nogrid
                Gm_in_sunl=[]
                GEl = []
                crel = []
                radl =[]
                kel = []
                for irad in range(nogrid):
                    cut = (Gr < dr*(irad+1)) 
                    crcut = cregy[cut]*1e10*Msun_in_g*km_in_cm*km_in_cm #erg
                    crel = np.append(crel, np.sum(crcut))
                    Gmcut = Gm[cut]
                    Gm_in_sun = Gmcut*1e10
                    Gm_in_sunl = np.append(Gm_in_sunl, np.sum(Gm_in_sun))
                    GEcut = GEint[cut]
                    kecut = ke[cut]
                    GEl = np.append(GEl,np.sum(GEcut)) #erg
                    kel = np.append(kel,np.sum(kecut))
                    radl = np.append(radl, dr*(irad+1))
                if wanted == 'gmr':
                    plt.plot(radl, Gm_in_sunl, label=runtodo)
                if wanted == 'avedenr':
                    vol = np.power(radl*kpc_in_cm,3.)*4./3.*np.pi
                    den = Gm_in_sunl*Msun_in_g/vol/protonmass_in_g
                    plt.plot(radl,den,label=labelneed,lw=2,ls=lsn,color=color)
                    print 'radl, den', radl, den
                    x0 = [1.0,-0.1]
                    xdata = np.log10(radl[1:]); ydata = np.log10(den[1:])
                    outfit=optimize.curve_fit(funcab, xdata, ydata,x0)
                    afit=outfit[0][0]
                    bfit=outfit[0][1]
                    print 'afit, bfit', afit, bfit
                if wanted == 'avecrdenr' or wanted == 'aveedenr':
                    if havecr>0:
                        vol = np.power(radl*kpc_in_cm,3.)*4./3.*np.pi
                        crden = crel*erg_in_eV/vol
                        if linelabelneed==1:
                            if labcount==2:
                                plt.plot(radl,crden,label='CR',lw=2,ls=lsn,color=color)
                            else:
                                plt.plot(radl,crden,lw=2,ls=lsn,color=color)
                        else:
                            plt.plot(radl,crden,label=dclabel,lw=2,ls=lsn,color=color)
                if wanted == 'avethdenr' or wanted == 'aveedenr':
                    vol = np.power(radl*kpc_in_cm,3.)*4./3.*np.pi
                    thden = GEl*erg_in_eV/vol
                    if linelabelneed==1:
                        if labcount==2:
                            plt.plot(radl,thden,label='Thermal',lw=1,ls=lsn,color=color)
                        else:
                            plt.plot(radl,thden,lw=1,ls=lsn,color=color)
                    else:
                        plt.plot(radl,thden,label=dclabel,lw=2,ls=lsn,color=color)
                if wanted == 'avekedenr':
                    vol = np.power(radl*kpc_in_cm,3.)*4./3.*np.pi
                    keden = kel*erg_in_eV/vol
                    plt.plot(radl,keden,label=dclabel,lw=2,ls=lsn,color=color)
                if wanted == 'FpFgr':
                    vol = np.power(radl*kpc_in_cm,3.)*4./3.*np.pi
                    dv = vol[1:]-vol[:-1]
                    drad = (radl[1:]-radl[:-1])*kpc_in_cm
                    dpth = (GAMMA-1.0)*(GEl[1:]-GEl[:-1])/dv/drad #cgs
                    Fg_m = Gm_in_sunl*Msun_in_g*NewtonG_in_cgs/(radl*kpc_in_cm)/(radl*kpc_in_cm)
                    dm = (Gm_in_sunl[1:]-Gm_in_sunl[:-1])*Msun_in_g
                    fg = Fg_m[:-1]*dm/dv #force per unit volume
                    plt.plot(radl[:-1],dpth/fg,label=dclabel,lw=2,ls='dashdot',color=color)
                    if havecr>0:
                        dpcr = (CRgamma-1.0)*(crel[1:]-crel[:-1])/dv/drad #cgs
                        plt.plot(radl[:-1],dpcr/fg,label=dclabel,lw=2,ls='solid',color=color)
                labcount+=1
            plt.yscale('log')
            plt.xlabel(r'$ r\; [{\rm kpc}]$',fontsize=14)
            if wanted =='gmr':
                plt.ylabel(r'enclosed $M_{\rm g} [M_{\odot}]$')
                filename=plotloc+'CRplot/gmr/gasmass_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
            if wanted =='avedenr':
                plt.ylabel(r'$\bar{n}_{\rm ISM}(<r)\; [{\rm cm^{-3}}]$',fontsize=14)
                filename=plotloc+'CRplot/avedenr/avedenr_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
                plt.legend(loc='best',fontsize=8,ncol=2)
            if wanted =='avecrdenr':
                plt.ylabel(r'$\bar{e}_{\rm cr}\; [{\rm eV/cm^{3}}]$',fontsize=14)
                filename=plotloc+'CRplot/avecrdenr/avecrdenr_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
            if wanted =='avethdenr':
                plt.ylabel(r'$\bar{e}_{\rm th}\; [{\rm eV/cm^{3}}]$',fontsize=14)
                filename=plotloc+'CRplot/avethdenr/avethdenr_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
            if wanted =='aveedenr':
                plt.ylabel(r'$\bar{e}(<r)\; [{\rm eV/cm^{3}}]$',fontsize=14)
                filename=plotloc+'CRplot/aveedenr/aveedenr_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
                plt.legend(loc='best',fontsize=14)
            if wanted=='avekedenr':
                plt.ylabel(r'$\bar{e}_{\rm KE}\; [{\rm eV/cm^{3}}]$',fontsize=14)
                filename=plotloc+'CRplot/avekedenr/avekedenr_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
            if wanted =='FpFgr':
                plt.ylabel(r'$\nabla P/\rho g$',fontsize=14)
                filename=plotloc+'CRplot/FpFgr/FpFgr_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
            plt.tight_layout()
            print 'filename', filename
            plt.title(ptitle,fontsize=16)
            plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
            plt.subplots_adjust(left  = 0.18, bottom=0.15, right=0.95, top=0.9)
            plt.savefig(filename)
            plt.clf()
        return None
