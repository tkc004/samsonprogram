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
import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='inferno', cmap=cmaps.inferno)
plt.register_cmap(name='plasma', cmap=cmaps.plasma)
import collections


def crdenmidplanetestinput(subdict):
        dirneed=subdict['dirneed']
        wanted=subdict['wanted']
        print 'wanted', wanted
        startno=subdict['startno']
        Nsnap=subdict['Nsnap']
        snapsep=subdict['snapsep']
        the_prefix=subdict['the_prefix']
        the_suffix=subdict['the_suffix']
        fmeat=subdict['fmeat']
        title='MW'
        ptitle='LSG'
        nolegend=0
        if wanted=='crdenmidplane' or wanted=='gasdenmidplane':
            resoneed=subdict['resoneed']
            newlabelneed=subdict['newlabelneed']
            legendneed = subdict['legendneed']
            plt.figure(figsize=(5,4))
            for runtodo in dirneed:
                print 'runtodo', runtodo
                if wanted=='crdenmidplane':
                    info=outdirname(runtodo, Nsnap)
                    havecr=info['havecr']
                    dclabel=info['dclabel']
                    haveB=info['haveB']
                    if havecr==0:
                        continue
                withinr=15.
                nogrid = 15
                maxlength=0.2 #thickness
                dr = withinr/nogrid
                radlist=np.linspace(0.001,withinr,num=nogrid)
                if wanted=='crdenmidplane':
                    credenlist=radlist*0.
                if wanted=='gasdenmidplane':
                    gasdenlist=radlist*0.
                numoftimes=0

                for i in range(startno,Nsnap,snapsep):
                    snaplist=[]
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
                    maindir=info['maindir']
                    color=info['color']
                    haveB=info['haveB']
                    M1speed=info['M1speed']
                    newlabel=info['newlabel']
                    snumadd=info['snumadd']
                    usepep=info['usepep']
                    halostr=info['halostr']
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
                    if havecr==0 and wanted=='crdenmidplane':
                        continue
                    rotface=1
                    if cosmo==1:
                        h0=1
                        datasup=0;
                    else:
                        h0=0
                        datasup=1;
                    print 'the_snapdir', the_snapdir
                    print 'Nsnapstring', Nsnapstring
                    Gextra = readsnapwcen(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix,\
                     havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
                     datasup=datasup,runtodo=runtodo,rundir=rundir,halostr=halostr)
                    Gp = Gextra['p'];
                    Gx = Gp[:,0]; Gy = Gp[:,1]; Gz = Gp[:,2];
                    Gv = Gextra['v'];
                    Gvx = Gv[:,0]; Gvy = Gv[:,1]; Gvz = Gv[:,2];
                    Grho = Gextra['rho']; Gu = Gextra['u']; Gm = Gextra['m']
                    if wanted=='crdenmidplane':
                        cregy = Gextra['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                    if wanted=='crdenmidplane':
                        cregy_in_erg = cregy*solar_mass_in_g*1e10*km_in_cm*km_in_cm
                        cregy_in_eV = cregy_in_erg*erg_in_eV
                    for irad in range(len(radlist)-1):
                        cutxy = ((Gx)*(Gx)+(Gy)*(Gy) > radlist[irad]*radlist[irad]) & ((Gx)*(Gx)+(Gy)*(Gy) < radlist[irad+1]*radlist[irad+1])
                        cutz = (Gz)*(Gz) < maxlength*maxlength/4.
                        cut = cutxy*cutz
                        shellvol_in_cm3 = np.pi*(-np.power(radlist[irad],2)+np.power(radlist[irad+1],2))*kpc_in_cm*kpc_in_cm*kpc_in_cm*maxlength
                        if wanted=='crdenmidplane':
                            cregy_in_eV_cut = np.sum(cregy_in_eV[cut])
                            creden_in_eV_per_cm3 = cregy_in_eV_cut/shellvol_in_cm3
                            credenlist[irad]+=creden_in_eV_per_cm3
                        if wanted=='gasdenmidplane':
                            Gm_in_g=Gm[cut]*1e10*Msun_in_g
                            Grho_in_g_cm_3 = np.sum(Gm_in_g)/shellvol_in_cm3
                            Gnism_in_cm_3 = 1.0/proton_mass_in_g*Grho_in_g_cm_3
                            #print 'radlist, Gnism_in_cm_3', radlist[irad], Gnism_in_cm_3
                            gasdenlist[irad]+=Gnism_in_cm_3
                    numoftimes+=1
                if haveB>0:
                    lsn='dashed'
                else:
                    lsn='solid'
                if resoneed==1:
                        if resolabel=='llr':
                                labelneed='low res'
                                lsn = 'solid'
                        if resolabel=='lr':
                                labelneed='std res'
                                lsn = 'dashed'
                        if resolabel=='mr':
                                labelneed='high res'
                                lsn = 'dashdot'
                if wanted=='crdenmidplane': 
                    plt.plot(radlist[:-1],credenlist[:-1]/numoftimes,label=labelneed,lw=2,ls=lsn,color=color)
                if wanted=='gasdenmidplane':
                    plt.plot(radlist[:-1],gasdenlist[:-1]/numoftimes,label=labelneed,lw=2,ls=lsn,color=color)
            if wanted=='crdenmidplane':
                #plt.plot(8,1.8,ls='none', marker='*', markersize=3, color=0.5)
                plt.xlabel(r'${\rm r [kpc]}$', fontsize=18)
                #plt.axhline(y=1.0,ls='--',color='k')
                plt.xlim(xmax=np.amax(radlist))
                plt.ylabel(r'$<e_{\rm CR}>[{\rm eV/cm^3}]$', fontsize=18)
                plt.yscale('log')
                #plt.legend(loc='best',frameon=False)
                filename=plotloc+'CRplot/crdenmidplane/CR_midplane_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            if wanted=='gasdenmidplane':
                plt.xlabel(r'${\rm r [kpc]}$', fontsize=18)
                plt.xlim(xmax=np.amax(radlist))
                plt.ylabel(r'$<n_{\rm ISM}>[{\rm cm^{-3}}]$',fontsize=18)
                if runtitle=='SMC':
                    plt.legend(loc='best',fontsize=10,ncol=2)
                filename=plotloc+'CRplot/gasdenmidplane/gas_midplane_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            print 'filename', filename
            plt.title(ptitle,fontsize=18)
            plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
            plt.subplots_adjust(left  = 0.18, bottom=0.15, right=0.95, top=0.9)
            plt.savefig(filename)
            plt.clf()
        return None
