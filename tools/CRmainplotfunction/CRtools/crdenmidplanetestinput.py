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
        usez1kpc = subdict['usez1kpc']
        print 'usez1kpc',usez1kpc
        title = subdict['title']
        maxlength = subdict['maxlength']
        multz = subdict['multz']
        logx = subdict['logx']
        print 'multz', multz
        #title='MW'
        ptitle='LSG'
        nolegend=0
        if wanted=='crdenmidplane' or wanted=='gasdenmidplane'\
            or wanted=='Bdenmidplane' or wanted=='B_gasdenmid'\
            or wanted=='Brmsmidplane' or wanted=='gassurden'\
            or wanted=='edenmidplane' or wanted=='thdenmidplane':
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
                if wanted=='Bdenmidplane' or wanted=='B_gasdenmid' or wanted=='Brmsmidplane':
                    info=outdirname(runtodo, Nsnap)
                    havecr=info['havecr']
                    dclabel=info['dclabel']
                    haveB=info['haveB']
                    if haveB==0:
                        continue
                withinr=15.
                nogrid = 30
                if usez1kpc==1:
                    fmeat+='_1kpc_'
                if wanted=='gassurden':
                    maxlength=6.0 #thickness
                elif usez1kpc==1:
                    maxlength=1.0
                dr = withinr/nogrid
                if multz==1:
                    zlist=[maxlength,maxlength/2.,maxlength/4.,maxlength/10.] 
                    lslist=['solid','dashed','dashdot','dotted']
                else:
                    zlist=[maxlength/2.]
                    lslist=['solid']
                for iz, zheight in enumerate(zlist): 
                    radlist=np.linspace(0.001,withinr,num=nogrid)
                    if wanted=='crdenmidplane':
                        credenlist=radlist*0.
                        credenplist=radlist*0.
                    if wanted=='thdenmidplane':
                        thdenlist=radlist*0.
                        thdenplist=radlist*0.
                    if wanted=='gasdenmidplane' or wanted=='B_gasdenmid' or wanted=='Brmsmidplane' or wanted=='gassurden':
                        gasdenlist=radlist*0.
                        gasdenplist=radlist*0.
                    if wanted=='Bdenmidplane' or wanted=='B_gasdenmid':
                        Bdenlist=radlist*0.                
                    if wanted=='Brmsmidplane':
                        Brmslist=radlist*0.
                        Brmsplist=radlist*0.
                    if wanted=='edenmidplane':
                        edenlist=radlist*0.
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
                        rotface=1
                        loccen=1
                        if cosmo==1:
                            h0=1
                        else:
                            h0=0
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        Gextra = readsnapwcen(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix,\
                         havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
                         loccen=loccen,runtodo=runtodo,rundir=rundir,halostr=halostr)
                        Gp = Gextra['p'];
                        Gx = Gp[:,0]; Gy = Gp[:,1]; Gz = Gp[:,2];
                        Gv = Gextra['v'];
                        Gvx = Gv[:,0]; Gvy = Gv[:,1]; Gvz = Gv[:,2];
                        Grho = Gextra['rho']; Gu = Gextra['u']; Gm = Gextra['m'];
                        if wanted=='edenmidplane':
                            Gne = Gextra['ne'] #(free) electron abundance per H atom
                            Gmetal = Gextra['z'] 
                            Z = Gmetal[:,0] #metal mass fraction (everything not H, He)
                            ZHe = Gmetal[:,1] # He mass fraction
                            NH = Gm*1e10*Msun_in_g*(1-(Z+ZHe))/protonmass_in_g #number of H atom
                            Ne = NH*Gne
                        if wanted=='crdenmidplane':
                            cregy = Gextra['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                            cregy_in_erg = cregy*solar_mass_in_g*1e10*km_in_cm*km_in_cm
                            cregy_in_eV = cregy_in_erg*erg_in_eV
                        if wanted=='thdenmidplane':
                            thegy_in_erg = Gu*Gm*solar_mass_in_g*1e10*km_in_cm*km_in_cm #thermal energy in erg
                            thegy_in_eV = thegy_in_erg*erg_in_eV                            
                        if wanted=='Bdenmidplane' or wanted=='B_gasdenmid':
                            GB = Gextra['B'] #B field in G
                            B2 = GB[:,0]*GB[:,0]+GB[:,1]*GB[:,1]+GB[:,2]*GB[:,2]
                            Begy_in_erg = B2/np.pi/8.*Gm/Grho*kpc_in_cm*kpc_in_cm*kpc_in_cm
                        if wanted=='Brmsmidplane':
                            GB = Gextra['B'] #B field in G
                            B2 = GB[:,0]*GB[:,0]+GB[:,1]*GB[:,1]+GB[:,2]*GB[:,2]
                            Brms_in_G = np.sqrt(B2)                      
                        for irad in range(len(radlist)-1):
                            cutxy = ((Gx)*(Gx)+(Gy)*(Gy) > radlist[irad]*radlist[irad]) & ((Gx)*(Gx)+(Gy)*(Gy) < radlist[irad+1]*radlist[irad+1])
                            cutz = (Gz)*(Gz) < zheight*zheight
                            cut = cutxy*cutz
                            shellvol_in_cm3 = np.pi*(-np.power(radlist[irad],2)+np.power(radlist[irad+1],2))*kpc_in_cm*kpc_in_cm*kpc_in_cm*zheight*2.0
                            if wanted=='edenmidplane':
                                Ne_cut = np.sum(Ne[cut])
                                ne_cm3 = Ne_cut/shellvol_in_cm3
                                edenlist[irad] += ne_cm3
                            if wanted=='crdenmidplane':
                                Gm_in_g=Gm[cut]*1e10*Msun_in_g
                                Grho_in_g_cm_3 = Grho[cut]*1e10*Msun_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm
                                Gparvol = Gm_in_g/Grho_in_g_cm_3
                                Gpsvol = np.sum(Gparvol)
                                credenplist[irad] += np.sum(cregy_in_eV[cut])/Gpsvol
                            if wanted=='thdenmidplane':
                                Gm_in_g=Gm[cut]*1e10*Msun_in_g
                                Grho_in_g_cm_3 = Grho[cut]*1e10*Msun_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm
                                Gparvol = Gm_in_g/Grho_in_g_cm_3
                                Gpsvol = np.sum(Gparvol)
                                thdenplist[irad] += np.sum(thegy_in_eV[cut])/Gpsvol
                                #credenlist[irad]+=creden_in_eV_per_cm3
                            if wanted=='gasdenmidplane' or wanted=='B_gasdenmid':
                                Gm_in_g=Gm[cut]*1e10*Msun_in_g
                                Grho_in_g_cm_3 = Grho[cut]*1e10*Msun_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm
                                Gnism_in_cm_3 = Grho_in_g_cm_3/proton_mass_in_g/1.3 #mean molecular weight in MW ~1.3
                                Gparvol = Gm_in_g/Grho_in_g_cm_3
                                Gpsvol = np.sum(Gparvol)
                                #print 'radlist, Gnism_in_cm_3', radlist[irad], Gnism_in_cm_3
                                gasdenplist[irad] += np.sum(Gnism_in_cm_3*Gparvol)/Gpsvol
                                #gasdenlist[irad] += Gnism_in_cm_3
                            if wanted=='gassurden':
                                Gm_in_Msun=Gm[cut]*1e10
                                cylarea_in_pc2 = np.pi*(-np.power(radlist[irad],2)+np.power(radlist[irad+1],2))*1e6
                                Gsurden_in_Msun_pc2 = np.sum(Gm_in_Msun)/cylarea_in_pc2
                                gasdenlist[irad] += Gsurden_in_Msun_pc2                      
                            if wanted=='Bdenmidplane' or wanted=='B_gasdenmid':
                                Begy_in_erg_cut = np.sum(Begy_in_erg[cut])
                                Begy_in_erg_per_cm3 = Begy_in_erg_cut/shellvol_in_cm3
                                Bdenlist[irad] += Begy_in_erg_per_cm3
                            if wanted=='Brmsmidplane':
                                Gm_in_g = Gm[cut]*1e10*Msun_in_g
                                Grho_in_g_cm_3 = Grho[cut]*1e10*Msun_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm
                                Gparvol = Gm_in_g/Grho_in_g_cm_3
                                Gpsvol = np.sum(Gparvol)
                                Brmsplist[irad] += np.sum(Brms_in_G[cut]*Gparvol)/Gpsvol*1e6
                                Brmslist[irad] += np.sum(Brms_in_G[cut]*Gparvol)/shellvol_in_cm3*1e6                            
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
                    lsn=lslist[iz]
                    if iz>0:
                        labelneed='_nolegend_' 
                    if wanted=='edenmidplane':
                        plt.plot(radlist[:-1],edenlist[:-1]/numoftimes,label=labelneed,lw=2,ls=lsn,color=color)
                    if wanted=='crdenmidplane': 
                        plt.plot(radlist[:-1],credenplist[:-1]/numoftimes,label=labelneed,lw=2,ls=lsn,color=color)
                    if wanted=='thdenmidplane': 
                        plt.plot(radlist[:-1],thdenplist[:-1]/numoftimes,label=labelneed,lw=2,ls=lsn,color=color)
                    if wanted=='gasdenmidplane':
                        plt.plot(radlist[:-1],gasdenplist[:-1]/numoftimes,label=labelneed,lw=2,ls=lsn,color=color)
                    if wanted=='gassurden':
                        plt.plot(radlist[:-1],gasdenlist[:-1]/numoftimes,label=labelneed,lw=2,ls=lsn,color=color)
                    if wanted=='Bdenmidplane': 
                        plt.plot(radlist[:-1],Bdenlist[:-1]/numoftimes,label=labelneed,lw=2,ls=lsn,color=color)
                    if wanted=='B_gasdenmid':
                        plt.plot(radlist[:-1],
                                 np.sqrt(Bdenlist[:-1]*8.*np.pi)/1e-6/np.power(gasdenlist[:-1],2./3.),
                                 label=labelneed,lw=2,ls=lsn,color=color)
                    if wanted=='Brmsmidplane': 
                        #plt.plot(radlist[:-1],Brmslist[:-1]/numoftimes,label=labelneed,lw=2,ls=lsn,color=color)
                        plt.plot(radlist[:-1],Brmsplist[:-1]/numoftimes,label=labelneed,lw=2,ls=lsn,color=color)
            if multz==1:
                import matplotlib.lines as mlines
                linelablist=[]
                lxlist=[]
                for icount in range(int(len(lslist))):
                    lineobj = mlines.Line2D([], [])
                    lxlist.append(lineobj)
                for ils, ls in enumerate(lslist):
                    linelabel = r'$|z| < $'+format(zlist[ils], '.2f')+' kpc'
                    lxlist[ils] = mlines.Line2D([], [], color='k', ls=ls,label=linelabel)
                    linelablist = np.append(linelablist,linelabel)
                legend1 = plt.legend(lxlist, linelablist, loc='lower left',fontsize=8,ncol=1)
                plt.gca().add_artist(legend1)
            if wanted=='edenmidplane':
                plt.xlabel(r'${\rm r [kpc]}$', fontsize=18)
                plt.xlim(xmax=np.amax(radlist))
                plt.ylabel(r'$<n_{\rm e}>[{\rm cm^{-3}}]$', fontsize=18)
                plt.yscale('log')
                filename=plotloc+'CRplot/edenmidplane/ne_midplane_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'               
            if wanted=='crdenmidplane':
                #plt.plot(8,1.8,ls='none', marker='*', markersize=3, color=0.5)
                plt.xlabel(r'${\rm r [kpc]}$', fontsize=18)
                #plt.axhline(y=1.0,ls='--',color='k')
                plt.xlim(xmax=np.amax(radlist))
                plt.ylabel(r'$<e_{\rm CR}>[{\rm eV/cm^3}]$', fontsize=18)
                plt.yscale('log')
                #plt.legend(loc='best',frameon=False)
                filename=plotloc+'CRplot/crdenmidplane/CR_midplane_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            if wanted=='thdenmidplane':
                #plt.plot(8,1.8,ls='none', marker='*', markersize=3, color=0.5)
                plt.xlabel(r'${\rm r [kpc]}$', fontsize=18)
                #plt.axhline(y=1.0,ls='--',color='k')
                plt.xlim(xmax=np.amax(radlist))
                plt.ylabel(r'$<e_{\rm th}>[{\rm eV/cm^3}]$', fontsize=18)
                plt.yscale('log')
                #plt.legend(loc='best',frameon=False)
                filename=plotloc+'CRplot/thdenmidplane/th_midplane_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            if wanted=='gasdenmidplane':
                plt.xlabel(r'${\rm r [kpc]}$', fontsize=18)
                plt.xlim(xmax=np.amax(radlist))
                plt.ylabel(r'$<n_{\rm ISM}>[{\rm cm^{-3}}]$',fontsize=18)
                if runtitle=='SMC':
                    plt.legend(loc='best',fontsize=10,ncol=2)
                plt.yscale('log')
                filename=plotloc+'CRplot/gasdenmidplane/gas_midplane_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            if wanted=='gassurden':
                plt.xlabel(r'${\rm r [kpc]}$', fontsize=18)
                plt.xlim(xmax=np.amax(radlist))
                plt.ylabel(r'$\Sigma_{\rm g}\;[{\rm M_\odot/pc^2}]$',fontsize=18)
                plt.yscale('log')
                plt.legend(loc='best',fontsize=10,ncol=2)
                filename=plotloc+'CRplot/gassurdenM_pc2/gassurdenM_pc2_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            if wanted=='Bdenmidplane':
                plt.xlabel(r'${\rm r [kpc]}$', fontsize=18)
                plt.xlim(xmax=np.amax(radlist))
                plt.ylabel(r'$<e_{\rm B}>[{\rm erg/cm^3}]$', fontsize=18)
                plt.yscale('log')
                filename=plotloc+'CRplot/Bdenmidplane/B_midplane_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            if wanted=='Brmsmidplane':
                plt.xlabel(r'${\rm r [kpc]}$', fontsize=18)
                plt.xlim(xmax=np.amax(radlist))
                plt.ylabel(r'$<B_{\rm rms}>[\mu{\rm G}]$', fontsize=18)
                plt.yscale('log')
                from matplotlib.ticker import ScalarFormatter
                ax = plt.gca()
                ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False, style='plain'))
                #ax.ticklabel_format(useOffset=False, style='plain')
                filename=plotloc+'CRplot/Brmsmidplane/Brms_midplane_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'            
            if wanted=='B_gasdenmid':
                plt.xlabel(r'${\rm r [kpc]}$', fontsize=18)
                plt.xlim(xmax=np.amax(radlist))
                plt.ylabel(r'$B_{\rm rms}/n_{\rm ISM}^{2/3}$', fontsize=18)
                plt.yscale('log')
                #from matplotlib.ticker import ScalarFormatter
                #plt.xaxis.set_major_formatter(ScalarFormatter())
                filename=plotloc+'CRplot/B_gasdenmid/B_gas_midplane_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf' 
            if cosmo==1:
                plt.legend(loc='best',fontsize=10,ncol=2)
            print 'filename', filename
            if logx==1:
                plt.xlim(xmin=0.1)
                plt.xscale('log')
            plt.title(ptitle,fontsize=18)
            plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
            plt.subplots_adjust(left  = 0.18, bottom=0.15, right=0.95, top=0.9)
            plt.savefig(filename)
            plt.clf()
        return None
