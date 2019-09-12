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
import collections


def crecumrtestinput(subdict):
        dirneed=subdict['dirneed']
        wanted=subdict['wanted']
        startno=subdict['startno']
        Nsnap=subdict['Nsnap']
        snapsep=subdict['snapsep']
        the_prefix=subdict['the_prefix']
        the_suffix=subdict['the_suffix']
        fmeat=subdict['fmeat']
        withinr=15.0
        nogrid=40
        maxlength=10.0
        print 'wanted', wanted
        print 'fmeat', fmeat
        print 'runtodo', dirneed


        if wanted=='gmz' or wanted=='crecumz' or wanted=='crecumr' or wanted=='Begyz' or wanted=='ethz' or wanted=='kez':
            plt.figure(figsize=(5,4))
            withinr=100.0
            minlength=0.5
            maxlength=100.
            nogrid=500
            normalized =1
            newlabelneed=subdict['newlabelneed']
            legendneed=subdict['legendneed']
            strlabelneed=subdict['strlabelneed']
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
                        loccen=0;
                    else:
                        loccen=1;
                    G = readsnapfromrun(runtodo,Nsnap,0,rotface=rotface,loccen=loccen)
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
                filename = plotloc+'CRplot/gmz/gasmass_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            else:
                filename = plotloc+'CRplot/gmz/gasmass_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'_normal.pdf'
        if wanted == 'crecumz':
            filename = plotloc+'CRplot/crecumz/crecumz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted == 'crecumr':
                    filename = plotloc+'CRplot/crecumr/crecumr_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted == 'Begyz':
                    filename = plotloc+'CRplot/Begyz/Begyz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted == 'ethz':
                    filename = plotloc+'CRplot/ethz/ethz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted == 'kez':
                    filename = plotloc+'CRplot/kez/kez_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        print 'filename', filename
        plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
        plt.subplots_adjust(left  = 0.18, bottom=0.15, right=0.95, top=0.9)
        plt.title(ptitle, fontsize=16)
        plt.savefig(filename)
        plt.clf()
        return None
