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
dirneed=['bwmwmr','bwmwmrstr','bwmwmrdc27','bwmwmrdc28', 'bwmwmrdc29']
#dirneed=['bwmwmr']
#dirneed=['bwsmclrdc28']
#dirneed=['bwmwmrdc29']
#dirneed=['m11gcr_700']
#dirneed=['bwsbclrdc27']
#dirneed=['m11dcr_b_70']
wanted='rhog'
startno=490
Nsnap=501
snapsep=10
the_prefix='snapshot'
the_suffix='.hdf5'
fmeat=''


if wanted=='crez' or wanted=='dpcrz' or wanted=='dpcrz_rhog' or wanted=='rhog'\
         or wanted=='dpcrzmrhog' or wanted=='dpcrz_dpth' or wanted=='pcr_pth'\
         or wanted=='vcr_vth' or  wanted=='dGmden':
        def func(x, a, b, c):
                return a+b*x+c*x*x
        def d10func(r, a, b, c):
                return np.power(10.,func(np.log10(r), a, b, c))*(b+2*c*np.log(r)/np.log(10.))/r
        def funcab(x, a, b):
                return a+b*x
        rcParams['figure.figsize'] = 5,4
        resoneed=0
        rotface=1
        needlogz=0
        thermalneed=1
        needfit=0
        labcount=0
        linestylabel=1
        newlabelneed=1
        twolegend=1
        ptotneed = 1
        if twolegend==1:
                import matplotlib.lines as mlines
                lxlist=[]
                lslist=[]
                dclablist=[]
                for icount in range(int(len(dirneed))):
                        lineobj = mlines.Line2D([], [])
                        lxlist.append(lineobj)
        for runtodo in dirneed:
                print 'runtodo', runtodo
                maxlength = 30.
                minlength=5.
                nogrid=8
                withinr=5.0
                diskr=10
                diskh = 1.0
                if needlogz==1:
                        zl=np.logspace(-1,np.log10(maxlength),num=nogrid)
                else:
                        zl=np.linspace(minlength,maxlength,num=nogrid)
                dz = maxlength/nogrid
                crden=zl*0.0
                thden=zl*0.0
                Gmden=zl*0.0
                mtotl=zl*0.0
                mdtot = 0.0
                numoftimes=0
                info=outdirname(runtodo, Nsnap)
                havecr=info['havecr']
                color=info['color']
                runtitle=info['runtitle']
                ptitle=title
                if runtitle=='SMC':
                        ptitle='Dwarf'
                elif runtitle=='SBC':
                        ptitle='Starburst'
                elif runtitle=='MW':
                        ptitle=r'$L\star$ Galaxy'
                if wanted=='crez':
                        if havecr==0:
                                continue
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
                        the_prefix = info['the_prefix']
                        the_suffix = info['the_suffix']
                        Nsnapstring=info['Nsnapstring']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        haveB=info['haveB']
                        newlabel=info['newlabel']
                        cosmo=info['cosmo']
                        halostr=info['halostr']
                        firever=info['firever']
                        maindir=info['maindir']
                        usepep=info['usepep']
                        snumadd=info['snumadd']
                        h0=cosmo
                        if cosmo==1:
                                loccen=0
                        else:
                                loccen=1
                        G = readsnapwcen(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix,\
                         havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
                         loccen=loccen,runtodo=runtodo,rundir=rundir,halostr=halostr)
                        Gp = G['p']; Gx = Gp[:,0]; Gy = Gp[:,1]; Gz = Gp[:,2];
                        Gv = G['v']; Gvx = Gv[:,0]; Gvy = Gv[:,1]; Gvz = Gv[:,2];
                        Grho = G['rho']; Gu = G['u']; Gm = G['m']
                        if havecr>0:
                                cregy = G['cregy']*1e10*Msun_in_g*km_in_cm*km_in_cm #cosmic ray energy in 1e10Msun km^2/sec^2

                        if haveB==0:
                                lsn='solid'
                        else:
                                lsn='dashed'
                        #if M1speed>999:
                        #       lsn='dotted'
                        #if M1speed>1999:
                        #       lsn='dashdot'
                        if wanted=='dpcrz_rhog' or wanted=='rhog' or wanted=='dpcrzmrhog' or wanted=='vcr_vth' or wanted=='dGmden':
                                mdisktot=0.0
                                print 'diskr, diskh', diskr, diskh
                                for ipd in [0,2,4]:
                                        try:
                                                Pextra = readsnapwcen(the_snapdir, Nsnapstring, ipd, snapshot_name=the_prefix, extension=the_suffix,\
                                                 havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
                                                 loccen=loccen,runtodo=runtodo,rundir=rundir,halostr=halostr)
                                                Px = Pextra['x']; Py = Pextra['y']; Pz = Pextra['z'];   
                                                Pm = Pextra['m']
                                                Prxy = np.sqrt(Px*Px+Py*Py)
                                                cutxyz = (Prxy<diskr)*(np.absolute(Pz)<diskh) #kpc
                                                mdisktot += np.sum(Pm[cutxyz]*1e10*Msun_in_g) #in g
                                                print 'mdisktot', mdisktot
                                        except KeyError:
                                                continue
                                mdtot += mdisktot
                        for iz in range(len(zl)-1):
                                cutxy = (Gx*Gx+Gy*Gy < withinr*withinr)
                                cutzu = Gz<zl[iz+1]
                                cutzd = Gz>zl[iz]
                                cutz = cutzu*cutzd
                                cut = cutxy*cutz
                                #print 'withinr', withinr
                                #print 'Gm[cut]', Gm[cut]
                                #print 'zl', zl
                                if havecr>0:
                                        crecut = cregy[cut]
                                        crden[iz] += np.sum(crecut)/(zl[iz+1]-zl[iz])/kpc_in_cm/(np.pi*withinr*withinr*kpc_in_cm*kpc_in_cm)
                                thermalcut = Gm[cut]*Gu[cut]*1e10*Msun_in_g*km_in_cm*km_in_cm
                                thden[iz] += np.sum(thermalcut)/(zl[iz+1]-zl[iz])/kpc_in_cm/(np.pi*withinr*withinr*kpc_in_cm*kpc_in_cm)
                                if wanted=='dpcrz_rhog' or wanted=='rhog' or wanted=='dpcrzmrhog' or wanted=='vcr_vth' or wanted=='dGmden':
                                        Gmcut=Gm[cut]*1e10*Msun_in_g
                                        Gmden[iz] += np.sum(Gmcut)/(zl[iz+1]-zl[iz])/kpc_in_cm/(np.pi*withinr*withinr*kpc_in_cm*kpc_in_cm)
                                        mintot = 0
                                        for ipa in [0,1,2,3,4]:
                                                try:
                                                        Pextra = readsnapwcen(the_snapdir, Nsnapstring, ipd, snapshot_name=the_prefix, extension=the_suffix,\
                                                         havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
                                                         loccen=loccen,runtodo=runtodo,rundir=rundir,halostr=halostr)
                                                        Px = Pextra['x']; Py = Pextra['y']; Pz = Pextra['z'];
                                                        Pm = Pextra['m']
                                                        Pr = np.sqrt(Px*Px+Py*Py+Pz*Pz)
                                                        cutpr = Pr<zl[iz]
                                                        if ipa==0 or ipa==2 or ipa==4:
                                                                Prxy = np.sqrt(Px*Px+Py*Py)
                                                                if zl[iz]<diskr:
                                                                        drc = zl[iz]
                                                                else:
                                                                        drc = diskr
                                                                cutdisk = (Prxy<drc)*(np.absolute(Pz)<diskh)
                                                                mdiskcut = np.sum(Pm[cutdisk]*1e10*Msun_in_g)
                                                        Pmtot = np.sum(Pm[cutpr]*1e10*Msun_in_g)-mdiskcut #in g
                                                        print 'Pmtot', Pmtot
                                                except KeyError:
                                                        continue
                                                mintot += Pmtot
                                        mtotl[iz] += mintot
                        numoftimes+=1
                if havecr>0:
                        crden = crden/numoftimes
                thden = thden/numoftimes
                Gmden = Gmden/numoftimes
                mtotl = mtotl/numoftimes
                mdtot = mdtot/numoftimes
                if havecr>0:
                        print 'crden', crden
                print 'thden', thden
                print 'mtotl', mtotl
                print 'mdtot', mdtot
                labelneed=dclabel
                if newlabelneed==1:
                        labelneed="\n".join(wrap(newlabel,17))
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
                zlm = (zl[1:]+zl[:-1])/2.
                #x0 = [0.1,-1.0,0.0]
                #Gmden1 = Gmden[:-1]
                #xdata=np.log10(zlm[~np.isnan(Gmden1)])
                #ydata=np.log10(Gmden1[~np.isnan(Gmden1)])
                #outfit=optimize.curve_fit(func, xdata, ydata, x0)
                #afit=outfit[0][0]
                #bfit=outfit[0][1]
                #cfit=outfit[0][2]
                #print 'afit,bfit,cfit', afit,bfit,cfit
                #dGmden = -d10func(zlm,afit,bfit,cfit)/kpc_in_cm
                if needfit==1:
                        x0 = [0.1,-1.0,0.0]
                        thden1 = thden[:-1]
                        xdata=np.log10(zlm[~np.isnan(thden1)])
                        ydata=np.log10(thden1[~np.isnan(thden1)])
                        outfit=optimize.curve_fit(func, xdata, ydata, x0)
                        afit=outfit[0][0]
                        bfit=outfit[0][1]
                        cfit=outfit[0][2]
                        #print 'afit,bfit,cfit', afit,bfit,cfit
                        dthcr = -d10func(zlm,afit,bfit,cfit)*(GAMMA-1.0)/kpc_in_cm
                else:
                        dthcr = -(thden[:-1]-thden[1:])/(zl[:-1]-zl[1:])*(GAMMA-1.0)/kpc_in_cm
                if havecr>0:
                        if needfit==1:
                                x0 = [0.1,-1.0,0.0]
                                crden1 = crden[:-1]
                                xdata=np.log10(zlm)
                                ydata=np.log10(crden1)
                                xdata = xdata[~np.isinf(ydata)]
                                ydata = ydata[~np.isinf(ydata)]
                                print 'xdata', xdata
                                print 'ydata', ydata
                                outfit=optimize.curve_fit(func, xdata, ydata, x0)
                                afit=outfit[0][0]
                                bfit=outfit[0][1]
                                cfit=outfit[0][2]
                                dpcr = -d10func(zlm,afit,bfit,cfit)*(CRgamma-1.0)/kpc_in_cm
                        else:
                                dpcr = -(crden[:-1]-crden[1:])/(zl[:-1]-zl[1:])*(CRgamma-1.0)/kpc_in_cm
                if wanted=='crez':
                        plt.plot(zlm, crden[:-1]*erg_in_eV, label=labelneed,lw=2,ls=lsn)
                if wanted=='dpcrz':
                        if havecr>0:
                                plt.plot(zlm, dpcr, label=labelneed,lw=2,ls=lsn,color=color)
                        plt.plot(zlm, dthcr,lw=1,color=color)
                if wanted=='dpcrz_dpth':
                        if havecr>0:
                                plt.plot(zlm, (crden[:-1]-crden[1:])/(thden[:-1]-thden[1:]),ls=lsn,color=color,marker='o')
                                plt.plot(zlm, dpcr/dthcr, label=labelneed,lw=2,ls=lsn,color=color)
                if wanted=='pcr_pth':
                        if havecr>0:
                                plt.plot(zlm, (CRgamma-1.0)*crden[:-1]/(GAMMA-1.0)/thden[:-1], label=labelneed,lw=2,ls=lsn,color=color)
                if wanted=='vcr_vth':
                        dden = -(Gmden[:-1]-Gmden[1:])/(zl[:-1]-zl[1:])/kpc_in_cm
                        gdisk = 2.0*NewtonG_cgs*mdtot/diskr/diskr*(1.0-zl/np.sqrt(diskr*diskr+zl*zl))/kpc_in_cm/kpc_in_cm
                        gsphere = mtotl*NewtonG_cgs/(zl*kpc_in_cm)/(zl*kpc_in_cm)
                        vg  = np.sqrt((gdisk+gsphere)*zl*kpc_in_cm)/km_in_cm
                        vth = np.sqrt(dthcr/dden)/km_in_cm
                        if havecr>0:
                                vcr = np.sqrt(dpcr/dden)/km_in_cm
                                vcom = np.sqrt((dthcr+dpcr)/dden)/km_in_cm
                        plt.plot(zl[:-1], vg[:-1], label=labelneed,lw=2,ls=lsn,color=color)
                        if havecr>0:
                                plt.plot(zlm, vcr,ls=lsn,marker='o',color=color,markeredgewidth=0.0,lw=1)
                        plt.plot(zlm, vth,ls=lsn,marker='^',color=color,markeredgewidth=0.0,lw=1)       
                if wanted=='dGmden':
                        dden = -(Gmden[:-1]-Gmden[1:])/(zl[:-1]-zl[1:])/kpc_in_cm
                        plt.plot(zlm,dden,ls='none',color=color,marker='o')
                        plt.plot(zlm,dGmden, label=labelneed,lw=2,ls=lsn,color=color)
                if wanted=='dpcrz_rhog':
                        print 'runtodo', runtodo
                        gdisk = 2.0*NewtonG_cgs*mdtot/diskr/diskr*(1.0-zl/np.sqrt(diskr*diskr+zl*zl))/kpc_in_cm/kpc_in_cm
                        gsphere = mtotl*NewtonG_cgs/(zl*kpc_in_cm)/(zl*kpc_in_cm)
                        rhog = Gmden*(gdisk+gsphere)
                        print 'rhog', rhog
                        if havecr>0:
                                plt.plot(zlm, dpcr/rhog[:-1], label=labelneed,lw=2,ls=lsn,color=color)
                        if thermalneed>0:
                                plt.plot(zlm, dthcr/rhog[:-1],lw=1,ls=lsn,color=color)
                if wanted=='dpcrzmrhog':
                        if havecr>0:
                                pcrden = (CRgamma-1.0)*crden
                        pthden = (GAMMA-1.0)*thden
                        print 'Gmden, mtotl', Gmden, mtotl
                        gdisk = 2.0*NewtonG_cgs*mdtot/diskr/diskr*(1.0-zl/np.sqrt(diskr*diskr+zl*zl))/kpc_in_cm/kpc_in_cm
                        gsphere = mtotl*NewtonG_cgs/(zl*kpc_in_cm)/(zl*kpc_in_cm)
                        gacc= gdisk+gsphere
                        if havecr>0:
                                print 'pcrden[:-1]', pcrden[:-1]
                                dpcr_rho = (pcrden[1:]-pcrden[:-1])/(zl[1:]-zl[:-1])/kpc_in_cm/Gmden[:-1]
                        dthcr_rho = (pthden[1:]-pthden[:-1])/(zl[1:]-zl[:-1])/kpc_in_cm/Gmden[:-1]
                        if havecr>0:
                                plt.plot(zlm[:-1], (dpcr_rho[:-1]+gacc[:-2])/km_in_cm*yr_in_sec*1e6, label=labelneed,lw=2,ls=lsn,color=color)
                        plt.plot(zlm[:-1], (dthcr_rho[:-1]+gacc[:-2])/km_in_cm*yr_in_sec*1e6,lw=1,ls=lsn,color=color)
                if wanted=='rhog':
                        gdisk = 2.0*NewtonG_cgs*mdtot/diskr/diskr*(1.0-zl/np.sqrt(diskr*diskr+zl*zl))/kpc_in_cm/kpc_in_cm
                        gsphere = mtotl*NewtonG_cgs/(zl*kpc_in_cm)/(zl*kpc_in_cm)
                        rhog = Gmden*(gdisk+gsphere)
                        if thermalneed>0:
                                if linestylabel==1 and labcount==1:
                                        plt.plot(zlm[:-1], dthcr[:-1],ms=4,ls='dashed',lw=1,color=color,marker='^',label=r'${\rm d} p_{\rm th}/{\rm d} z$')
                                else:
                                        plt.plot(zlm[:-1], dthcr[:-1],ms=4,ls='dashed',lw=1,color=color,marker='^')
                        if havecr>0:
                                if linestylabel==1 and labcount==1:
                                        plt.plot(zlm[:-1], dpcr[:-1],ms=4,ls='solid',lw=1,color=color,marker='o',label=r'${\rm d} p_{\rm CR}/{\rm d} z$')
                                        if ptotneed==1:
                                                plt.plot(zlm[:-1], dpcr[:-1]+dthcr[:-1],ms=4,ls='solid',lw=1,color=color,marker='d',label=r'${\rm d} p_{\rm tot}/{\rm d} z$')
                                else:
                                        plt.plot(zlm[:-1], dpcr[:-1],ms=4,ls='solid',lw=1,color=color,marker='o')
                                        if ptotneed==1:
                                                plt.plot(zlm[:-1], dpcr[:-1]+dthcr[:-1],ms=4,ls='solid',lw=1,color=color,marker='d')
                                #x0 = [0.1,-1.0,0.0]
                                #outfit = optimize.curve_fit(d10func,zlm[:-1],dpcr[:-1],x0)
                                #afit=outfit[0][0]
                                #bfit=outfit[0][1]
                                #cfit=outfit[0][2]
                                #print 'afit,bfit,cfit', afit,bfit,cfit
                        #plt.plot(zlm, dthcr,ms=4,ls='none',color=color,marker='s')
                        if linestylabel==1:
                                if labcount==0:
                                        plt.plot(zl[:-1], rhog[:-1],color=color,label=r'$\rho g$',ls=lsn,lw=2)
                                else:
                                        plt.plot(zl[:-1], rhog[:-1],color=color,ls=lsn,lw=2)
                        else:
                                plt.plot(zl[:-1], rhog[:-1],color=color,label=labelneed,ls=lsn,lw=2)
                        if twolegend==1:
                                lxlist[labcount] = mlines.Line2D([], [], color=color, ls='solid',label=labelneed)
                                dclablist = np.append(dclablist,labelneed)
                        labcount+=1
        if wanted=='crez':
                plt.yscale('log')
                plt.xlabel('z [kpc]')
                plt.ylabel(r'$e_{\rm cr} [{\rm eV/cm^3}]$')
                plt.legend(loc='best')
                plt.title('cylinder centered at disk with a radius ' + str(withinr) + ' kpc')
                filename=plotloc+'CRplot/crez/crez_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted=='dpcrz':
                plt.yscale('log')
                #plt.xscale('log')
                plt.xlabel(r'${\rm z\; [kpc]}$',fontsize=18)
                plt.ylabel(r'${\rm d}P/{\rm d}z [{\rm erg/cm^4}]$',fontsize=18)
                plt.legend(loc='best')
                #plt.title('cylinder centered at disk with a radius ' + str(withinr) + ' kpc')
                filename=plotloc+'CRplot/dpcrz/dpcrz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted=='dpcrz_dpth':
                plt.yscale('log')
                #plt.xscale('log')
                plt.xlabel(r'${\rm z\; [kpc]}$',fontsize=18)
                plt.ylabel(r'${\rm d}P_{\rm cr}/{\rm d}P_{\rm th}$',fontsize=18)
                #plt.legend(loc='best')
                #plt.title('cylinder centered at disk with a radius ' + str(withinr) + ' kpc')
                filename=plotloc+'CRplot/dpcrz_dpth/dpcrz_dpth_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted=='pcr_pth':
                plt.yscale('log')
                #plt.xscale('log')
                plt.xlabel(r'${\rm z\; [kpc]}$',fontsize=18)
                plt.ylabel(r'$P_{\rm cr}/P_{\rm th}$',fontsize=18)
                #plt.legend(loc='best')
                #plt.title('cylinder centered at disk with a radius ' + str(withinr) + ' kpc')
                filename=plotloc+'CRplot/pcr_pth/pcr_pth_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
                plt.title(ptitle)
        if wanted=='vcr_vth':
                #plt.yscale('log')
                #plt.xscale('log')
                plt.xlabel(r'${\rm z\; [kpc]}$',fontsize=18)
                plt.ylabel(r'$v\;{\rm [km/s]}$',fontsize=18)
                #plt.legend(loc='best')
                #plt.title('cylinder centered at disk with a radius ' + str(withinr) + ' kpc')
                filename=plotloc+'CRplot/vcr_vth/vcr_vth_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
                plt.title(ptitle,fontsize=16)
        if wanted=='dGmden':
                #plt.yscale('log')
                #plt.xscale('log')
                plt.xlabel(r'${\rm z\; [kpc]}$',fontsize=18)
                plt.ylabel(r'$\;{\rm d} \rho/{\rm d} z{\rm [g/cm^4]}$',fontsize=18)
                #plt.legend(loc='best')
                #plt.title('cylinder centered at disk with a radius ' + str(withinr) + ' kpc')
                filename=plotloc+'CRplot/dGmden/dGmden_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
                plt.title(ptitle,fontsize=16)
        if wanted=='dpcrz_rhog':
                #plt.yscale('log')
                #plt.axhspan(-2, 1, alpha=0.5,color='0.5')
                plt.ylim(ymin=-0.5)
                plt.xlabel(r'${\rm z\; [kpc]}$',fontsize=18)
                plt.ylabel(r'${\rm d}P/{\rm d}z/(\rho g)$',fontsize=18)
                #plt.legend(loc='best')
                #plt.title('cylinder centered at disk with a radius ' + str(withinr) + ' kpc')
                filename=plotloc+'CRplot/dpcrz_rhog/dpcrz_rhog_hd_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted=='dpcrzmrhog':
                plt.xlabel(r'${\rm z\; [kpc]}$',fontsize=18)
                plt.ylabel(r'$g+{\rm d}_zP/\rho\;  [{\rm km/s/Myr}]$',fontsize=18)
                #plt.legend(loc='best')
                #plt.title('cylinder centered at disk with a radius ' + str(withinr) + ' kpc')
                filename=plotloc+'CRplot/dpcrzmrhog/dpcrzmrhog_hd_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'

        if wanted=='rhog':
                if twolegend==1:
                        #print 'lxlist', lxlist
                        #print 'dclablist', dclablist
                        legend1 = plt.legend(lxlist, dclablist, loc=3,fontsize=8,ncol=2)
                        plt.gca().add_artist(legend1)
                        plt.legend(loc=3, fontsize=8)
                plt.yscale('log')
                #plt.xscale('log')
                plt.xlabel(r'${\rm z\; [kpc]}$',fontsize=18)
                plt.ylabel(r'$\rho g [{\rm erg/cm^4}]$',fontsize=18)
                if runtitle=='SMC' or linestylabel==1:
                        plt.legend(loc='best',ncol=3,fontsize=10)
                plt.title(ptitle,fontsize=16)

