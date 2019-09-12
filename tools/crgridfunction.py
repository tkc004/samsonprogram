from pathloc import *
import samson_functions as SSF



def crhisto(runtodo,i,wanted,nooftimes,Tcut,highTcut,vcut,vhcut,withinr,zup,zdown,userad,trackgas=0,rotface=1,normalized=0):
    if trackgas==1:
        vcut = -1e10 #if we track gas particles, we should consider all velocities
    if nooftimes==1:
        return None
    data = gasphase(runtodo,i,rotface=rotface,Tcut=Tcut,highTcut=highTcut,\
    vcut=vcut,vhcut=vhcut,withinr=withinr,zup=zup,zdown=zdown,userad=userad)
    partZ = data['partZ']
    partR = data['partR']
    vz = data['vz']
    vr = data['vr']
    TrueTemp = data['TrueTemp']
    Gu = data['Gu']
    rho = data['rho']
    cregy = data['cregy']
    converted_rho = data['convertedrho']
    Gmass = data['Gmass']
    vmax = data['vmax']
    vmin = data['vmin']
    Gid = data['Gid']
    if trackgas==1:
        idint0=np.in1d(Gid,Gid_t)
        if userad==1:
            partR=partR[idint0]
            vr=vr[idint0]
        else:
            partZ=partZ[idint0]
            vz=vz[idint0]
        TrueTemp=TrueTemp[idint0]
        converted_rho=converted_rho[idint0]
        Gmass = Gmass[idint0]
    if userad==1:
        x = partR
    else:
        x = partZ
    if wanted=='Tvz':
        if userad==1:
            x = np.log10(vr)
        else:
            x = np.log10(vz)
    if wanted=='Tvznolog':
        if userad==1:
            x = vr
        else:
            x = vz    
    if wanted=='pcrpth':
        pth = (GAMMA-1.0)*Gu*u_codetocgs*(rho*rho_codetocgs) #cgs
        x = np.log10(pth)
        print 'x', x
    if wanted=='Tz' or wanted=='Tvz' or wanted=='Tvznolog':
        y = np.log10(TrueTemp)
    if wanted=='rhoz':
        y = np.log10(converted_rho)
    if wanted=='vzz':
        if userad==1:
            y = np.log10(vr)
        else:
            y = np.log10(vz)
            
        if uselog==0:
            if userad==1:
                y = vr
            else:
                y = vz
    if wanted=='pcrpth':
        pcr = (CRgamma-1.0)*cregy*cregy_codetocgs*(rho*rho_codetocgs)/(Gmass*m_codetocgs) #cgs
        y = np.log10(pcr)
        print 'y', y
    gridxx = np.linspace(extent[0],extent[1],nobin)
    gridyy = np.linspace(extent[2],extent[3],nobin)
    if normalized==1:
        weights = Gmass*1e10
    else:
        Gmtot = np.sum(Gmass)
        weights = Gmass/Gmtot
    H, xedges, yedges = np.histogram2d(x, y, bins=[gridxx, gridyy],weights=weights)
    return H


def crgridfunction(gridag,runtodo,wanted,idir,plotupt,plotit,startno,Nsnap,snapsep,cbcolor='hot',
                   rotface=1,
                   Tcut=-1.0, #K
                   highTcut=1e10,
                   vcut=-1e10, #outflow velocity cut
                   vhcut=1e10, #upper outflow velocity cut
                   withinr=100.0, 
                   withoutr=0.01,
                   zup=100.0,
                   zdown=-100.0,
                   userad=0,
                   trackgas=0,
                   fmeat='',
                   usesmallB=0,
                   normalized=0,
                   usedatamaxmin=1,
                   addTcr=1,
                   extendedrange=0):
    if wanted=='Tz' or wanted=='rhoz' or wanted=='vzz' or wanted=='TcrT' or wanted=='rhoTcr' or wanted=='rhoPcr'\
       or wanted=='Tvz' or wanted=='pcrpth' or wanted=='rhoT' or wanted=='rhoB' or wanted=='rhoPth'\
       or wanted=='Tvznolog':
        newlabelneed=1
        if fmeat=='': fmeat=runtodo
        from crtestfunction import *
        from textwrap import wrap
        tlabel=''
        usehalfz=0
        if trackgas==1:
            snaptrack = Nsnap
            Tcut_t=1.0 #K
            highTcut_t = 1e10
            vcut_t=0.0 #outflow velocity cut
            vhcut_t = 1e10 #upper outflow velocity cut
            withinr_t=20.0
            withoutr_t=0.01
            zup_t=30.0
            zdown_t=20.0
            tlabel='track'
        if userad==1:
            tlabel+='rad'
        if wanted=='Tz':
            extent = [zdown,zup,1,9]
        if wanted=='rhoz':
            extent = [zdown,zup,-7,4]
        if wanted=='vzz':
            uselog=0
            extent = [zdown,zup,0,3.2]
            if uselog==0:
                extent = [zdown,zup,-1000.0,1000.0]
        if wanted=='Tvz':
            extent = [0,3.2,1,9]
            #extent = [0,9,0,9]
        if wanted=='Tvznolog':
            extent = [-5e2,5e2,3,7]
            usehalfz=1
        if wanted=='pcrpth':
            extent = [-17,-9,-17,-9]
        if wanted=='rhoT':
            if extendedrange==1:
                extent = [-6,4,0,7.5]
            else:
                extent = [-3.5,4,0,9]
        if wanted=='rhoTcr':
            extent = [-3.5,4,0,9]
        if wanted=='rhoPcr':
            extent = [-3.5,4,-17,-9]
        if wanted=='rhoPth':
            extent = [-3.5,4,-17,-9] 
        if wanted=='TcrT':
            extent = [1,7,1,7]
        if wanted=='rhoB':
            if usesmallB==1:
                extent = [-4,4,-10,-3]
            else:
                extent = [-4,4,-10,-3]
        nobin=51
        needcontour=1
        info=SSF.outdirname(runtodo, Nsnap)
        cosmo=info['cosmo']
        runtitle=info['runtitle']
        runlabel=info['runlabel']
        havecr=info['havecr']
        dclabel=info['dclabel']
        newlabel=info['newlabel']
        if wanted=='pcrpth' and havecr==0:
            return None
        ptitle=runlabel
        if runtitle=='SMC':
            ptitle='Dwarf'
        elif runtitle=='SBC':
            ptitle='Starburst'
        elif runtitle=='MW':
            ptitle=r'$L\star$ Galaxy'
        labelneed=dclabel
        if newlabelneed==1:
            labelneed="\n".join(wrap(newlabel,17))
        plotupt=np.append(plotupt,ptitle)
        plotit =np.append(plotit,labelneed)
        Hadd = np.zeros((nobin-1,nobin-1))
        if trackgas==1:
            data = gasphase(runtodo,snaptrack,rotface=rotface,Tcut=Tcut_t,highTcut=highTcut_t,\
            vcut=vcut_t,vhcut=vhcut_t,withinr=withinr_t,withoutr=withoutr_t,zup=zup_t,zdown=zdown_t,userad=userad)
            Gid_t = data['Gid']
        nooftimes=0
        for i in range(startno,Nsnap,snapsep):
            if trackgas==1:
                vcut = -1e10 #if we track gas particles, we should consider all velocities
            if nooftimes==1:
                return None
            data = gasphase(runtodo,i,rotface=rotface,Tcut=Tcut,highTcut=highTcut,\
            vcut=vcut,vhcut=vhcut,withinr=withinr,withoutr=withoutr,zup=zup,zdown=zdown,userad=userad,usehalfz=usehalfz)
            partZ = data['partZ']
            partR = data['partR']
            vz = data['vz']
            vr = data['vr']
            TrueTemp = data['TrueTemp']
            Gu = data['Gu']
            rho = data['rho']
            cregy = data['cregy']
            Brms = data['Brms']
            converted_rho = data['convertedrho']
            Gmass = data['Gmass']
            if normalized==1:
                vmax = 1.0
                vmin = 0
            else:
                vmax = data['vmax']
                vmin = data['vmin']
            Gid_n = data['Gid']
            if trackgas==1:
                from gastracklib import multapply
                Gid_tu, unindex_t = np.unique(Gid_t, return_index=True)
                Gid_nu, repindex_n = np.unique(Gid_n, return_counts=True)
                unindex_n = np.in1d(Gid_n,Gid_nu[repindex_n<2]) #consider only unique id, drop all repeated
                Gidnr=Gid_n[unindex_n]
                idint0=np.in1d(Gidnr,Gid_tu)

                if userad==1:
                    partR=multapply(partR,unindex_n,idint0)
                    vr=multapply(vr,unindex_n,idint0)
                else:
                    partZ=multapply(partZ,unindex_n,idint0)
                    vz=multapply(vz,unindex_n,idint0)
                TrueTemp=multapply(TrueTemp,unindex_n,idint0)
                converted_rho=multapply(converted_rho,unindex_n,idint0)
                Gmass=multapply(Gmass,unindex_n,idint0)
            if userad==1:
                x = partR
            else:
                x = partZ
            if wanted=='Tvz':
                if userad==1:
                    x = np.log10(vr)
                else:
                    x = np.log10(vz)
            if wanted=='Tvznolog':
                if userad==1:
                    x = vr
                else:
                    x = vz                
            if wanted=='pcrpth':
                pth = (GAMMA-1.0)*Gu*u_codetocgs*(rho*rho_codetocgs) #cgs
                x = np.log10(pth)
                print 'x', x
            if wanted=='rhoT':
                x = np.log10(converted_rho)
            if wanted=='TcrT':
                Tcr = (CRgamma-1.0)*cregy*cregy_codetocgs*(rho*rho_codetocgs)/(Gmass*m_codetocgs)/BOLTZMANN/converted_rho #K
                x = np.log10(Tcr)
            if wanted=='rhoTcr':
                x = np.log10(converted_rho)                
            if wanted=='rhoPcr':
                x = np.log10(converted_rho)
            if wanted=='rhoPth':
                x = np.log10(converted_rho)                
            if wanted=='rhoB':
                x = np.log10(converted_rho)
            if wanted=='rhoTcr':
                Tcr = (CRgamma-1.0)*cregy*cregy_codetocgs*(rho*rho_codetocgs)/(Gmass*m_codetocgs)/BOLTZMANN/converted_rho #K
                y = np.log10(Tcr)  
            if wanted=='rhoPcr':
                Pcr = (CRgamma-1.0)*cregy*cregy_codetocgs*(rho*rho_codetocgs)/(Gmass*m_codetocgs) #cgs
                y = np.log10(Pcr)
            if wanted=='rhoPth':
                pth = (GAMMA-1.0)*Gu*u_codetocgs*(rho*rho_codetocgs) #cgs
                y = np.log10(pth)
            if wanted=='rhoT' or wanted=='TcrT':
                y = np.log10(TrueTemp)
            if wanted=='rhoB':
                y = np.log10(Brms)
            if wanted=='Tz' or wanted=='Tvz' or wanted=='Tvznolog':
                y = np.log10(TrueTemp)
            if wanted=='rhoz':
                y = np.log10(converted_rho)
            if wanted=='vzz':
                if userad==1:
                    y = np.log10(vr)
                else:
                    y = np.log10(vz)
                if uselog==0:
                    if userad==1:
                        y = vr
                    else:
                        y = vz                    
            if wanted=='pcrpth':
                pcr = (CRgamma-1.0)*cregy*cregy_codetocgs*(rho*rho_codetocgs)/(Gmass*m_codetocgs) #cgs
                y = np.log10(pcr)
                print 'y', y
            gridxx = np.linspace(extent[0],extent[1],nobin)
            gridyy = np.linspace(extent[2],extent[3],nobin)
            if normalized==1:
                Gmtot = np.sum(Gmass)            
                weights = Gmass/Gmtot
                #needdensity = True
            else:
                weights = Gmass*1e10
                #needdensity = False
            #, density=needdensity
            H, xedges, yedges = np.histogram2d(x, y, bins=[gridxx, gridyy],weights=weights)
            print 'H', H
            H=H.T
            Hadd += H
            nooftimes += 1
        Hadd=Hadd/nooftimes
        if normalized==1:
            vmax = -1.5
            vmin = -7.0
        elif usedatamaxmin==1:
            vmax = np.amax(np.log10(Hadd))
            vmin = vmax-6.0
            print 'usedatamaxmin'
            print 'vmin', vmin
            print 'vmax', vmax
        if needcontour==1:
            numaspect = (extent[0]-extent[1])/(extent[2]-extent[3])
            gridag[idir].set_aspect(numaspect)
            print 'np.sum(Hadd)', np.sum(Hadd)
            print 'xedges', xedges
            print 'extent', extent
            levels = np.linspace(vmin,vmax,num=9)
            im = gridag[idir].contourf(xedges[:-1], yedges[:-1], np.log10(Hadd),\
            aspect='auto',cmap=cbcolor,levels=levels, extent=extent,origin='lower')
        else:
            im = gridag[idir].imshow(np.log10(Hadd), extent=extent,\
            aspect='auto',interpolation='nearest', cmap=cbcolor, origin='lower',vmax=vmax,vmin=vmin)
        if normalized==1:
            cblabel = r'$\mathrm{Log}_{10} ({\rm probability})$'
        else:
            cblabel = r'$\mathrm{Log}_{10} (M {\rm [M_\odot\;per\;grid]})$'
        if userad==1:
            gridag[idir].set_xlabel(r"$r {\rm [kpc]}$")
        else:
            gridag[idir].set_xlabel(r"$z {\rm [kpc]}$")
        if wanted=='Tvz':
            if userad==1:
                gridag[idir].set_xlabel(r"$\mathrm{Log}_{10} (v_r {\rm [km/s]})$")
            else:
                gridag[idir].set_xlabel(r"$\mathrm{Log}_{10} (v_z {\rm [km/s]})$")
        if wanted=='Tvznolog':
            if userad==1:
                gridag[idir].set_xlabel(r"$v_r {\rm [km/s]}$")
            else:
                gridag[idir].set_xlabel(r"$v_z {\rm [km/s]}$")
        if wanted=='TcrT':
            if addTcr==1:
                Ttest = SSF.calstraightline(gridxx, point1=(4.,4.), slope= 1.0)
                gridag[idir].plot(gridxx,Ttest,ls='dashed',color='k',lw=1)
                #Ttest = SSF.calstraightline(gridxx, point1=(0.,logTcr-1.0), slope= -1.0)
                #gridag[idir].plot(gridxx,Ttest,color='k',lw=1,ls='dotted',label=r'$T_{\rm cr}({\rm 10\;eV/cm^3})$')          
            gridag[idir].set_xlabel(r"$\mathrm{Log}_{10} (T_{\rm cr} {\rm [K]})$")
            gridag[idir].set_ylabel(r"$\mathrm{Log}_{10} (T {\rm [K]})$")
        if wanted=='rhoT':
            if addTcr==1:
                logTcr = np.log10(1.0*eV_in_erg*(CRgamma-1.0)/BOLTZMANN) #effective temperature from CR pressure at 1eV/cm^3 and nism=1cm^-3
                Ttest = SSF.calstraightline(gridxx, point1=(0.,logTcr), slope= -1.0)
                gridag[idir].plot(gridxx,Ttest,ls='dashed',color='k',lw=1,label=r'$T_{\rm cr}({\rm 1\;eV/cm^3})$')
                Ttest = SSF.calstraightline(gridxx, point1=(0.,logTcr-1), slope= -1.0)
                gridag[idir].plot(gridxx,Ttest,ls='dotted',color='k',lw=1,label=r'$T_{\rm cr}({\rm 0.1\;eV/cm^3})$')
                if havecr>0:
                    gridag[idir].legend(loc='best',fontsize=12,frameon=False)
                #Ttest = SSF.calstraightline(gridxx, point1=(0.,logTcr-1.0), slope= -1.0)
                #gridag[idir].plot(gridxx,Ttest,color='k',lw=1,ls='dotted',label=r'$T_{\rm cr}({\rm 10\;eV/cm^3})$')
            gridag[idir].set_xlim((extent[0],extent[1]))
            gridag[idir].set_ylim((extent[2],extent[3]))            
            gridag[idir].set_xlabel(r"$\mathrm{Log}_{10} (n {\rm [cm^{-3}]})$")
            gridag[idir].set_ylabel(r"$\mathrm{Log}_{10} (T {\rm [K]})$")
        if wanted=='rhoPcr':        
            gridag[idir].set_xlabel(r"$\mathrm{Log}_{10} (n {\rm [cm^{-3}]})$")
            gridag[idir].set_ylabel(r"$\mathrm{Log}_{10} (P_{\rm cr} {\rm [erg/cm^3]})$")
        if wanted=='rhoPth':        
            gridag[idir].set_xlabel(r"$\mathrm{Log}_{10} (n {\rm [cm^{-3}]})$")
            gridag[idir].set_ylabel(r"$\mathrm{Log}_{10} (P_{\rm th} {\rm [erg/cm^3]})$")
        if wanted=='rhoTcr':        
            gridag[idir].set_xlabel(r"$\mathrm{Log}_{10} (n {\rm [cm^{-3}]})$")
            gridag[idir].set_ylabel(r"$\mathrm{Log}_{10} (T_{\rm cr} {\rm [K]})$")
        if wanted=='rhoB':
            Btest = SSF.calstraightline(gridxx, point1=(0.,-5.22), slope= 2./3.)
            gridag[idir].plot(gridxx,Btest,ls='dashed',color='k')
            gridxobsup = np.linspace(np.log10(300.0),extent[1],nobin)
            gridxobsdown = np.linspace(np.log10(10.0),np.log10(300.0),nobin)
            gridyobsup = np.log10(49e-6*np.power(np.power(10.0,gridxobsup-4.0),0.65))
            gridyobsdown = np.ones(nobin)*np.log10(5e-6)
            gridag[idir].plot(gridxobsup,gridyobsup,color='k',lw=1,label='Obs')
            gridag[idir].plot(gridxobsdown,gridyobsdown,color='k',lw=1)
            gridag[idir].set_xlabel(r"$\mathrm{Log}_{10} (n {\rm [cm^{-3}]})$")
            gridag[idir].set_ylabel(r"$\mathrm{Log}_{10} (B_{\rm rms} {\rm [G]})$")            
        if wanted=='pcrpth':
            gridag[idir].set_xlabel(r"$\mathrm{Log}_{10} (P_{\rm th} {\rm [erg/cm^3]})$")
            gridag[idir].set_ylabel(r"$\mathrm{Log}_{10} (P_{\rm cr} {\rm [erg/cm^3]})$")
        if wanted=='rhoz':
            gridag[idir].set_ylabel(r"$\mathrm{Log}_{10} (n {\rm [cm^{-3}]})$")
        if wanted=='Tz' or wanted=='Tvz' or wanted=='Tvznolog':
            gridag[idir].set_ylabel(r"$\mathrm{Log}_{10} (T {\rm [K]})$")
        if wanted=='Tvznolog':
            gridag[idir].axvline(x=0.0,color='k')
        if wanted=='vzz':
            gridag[idir].axhline(y=0.0,color='k')
            if userad==1:
                gridag[idir].set_ylabel(r"$\mathrm{Log}_{10} (v_r {\rm [km/s]})$")
            else:
                gridag[idir].set_ylabel(r"$\mathrm{Log}_{10} (v_z {\rm [km/s]})$")
            if uselog==0:
                if userad==1:
                    gridag[idir].set_ylabel(r"$v_r {\rm [km/s]}$")
                else:
                    gridag[idir].set_ylabel(r"$v_z {\rm [km/s]}$")                
            #gridag[idir].set_yscale('symlog')
        if wanted == 'rhoT':
            if needcontour==1:
                totalname = plotloc+'/CRplot/rhoT/rhoT_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'_contour'+tlabel+'ag.pdf'
            else:
                totalname = plotloc+'/CRplot/rhoT/rhoT_'+fmeat+'sn'+str(startno)+'_'+str(Nsnap)+tlabel+'ag.pdf'
        if wanted == 'rhoTcr':
            if needcontour==1:
                totalname = plotloc+'/CRplot/rhoTcr/rhoTcr_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'_contour'+tlabel+'ag.pdf'
            else:
                totalname = plotloc+'/CRplot/rhoTcr/rhoTcr_'+fmeat+'sn'+str(startno)+'_'+str(Nsnap)+tlabel+'ag.pdf'
        if wanted == 'rhoPcr':
            if needcontour==1:
                totalname = plotloc+'/CRplot/rhoPcr/rhoPcr_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'_contour'+tlabel+'ag.pdf'
            else:
                totalname = plotloc+'/CRplot/rhoPcr/rhoPcr_'+fmeat+'sn'+str(startno)+'_'+str(Nsnap)+tlabel+'ag.pdf'
        if wanted == 'rhoPth':
            if needcontour==1:
                totalname = plotloc+'/CRplot/rhoPth/rhoPth_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'_contour'+tlabel+'ag.pdf'
            else:
                totalname = plotloc+'/CRplot/rhoPth/rhoPth_'+fmeat+'sn'+str(startno)+'_'+str(Nsnap)+tlabel+'ag.pdf'
        if wanted == 'TcrT':
            if needcontour==1:
                totalname = plotloc+'/CRplot/TcrT/TcrT_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'_contour'+tlabel+'ag.pdf'
            else:
                totalname = plotloc+'/CRplot/TcrT/TcrT_'+fmeat+'sn'+str(startno)+'_'+str(Nsnap)+tlabel+'ag.pdf'                
        if wanted == 'rhoB':
            if needcontour==1:
                totalname = plotloc+'/CRplot/rhoB/rhoB_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'_contour'+tlabel+'ag.pdf'
            else:
                totalname = plotloc+'/CRplot/rhoB/rhoB_'+fmeat+'sn'+str(startno)+'_'+str(Nsnap)+tlabel+'ag.pdf'
        if wanted == 'Tz':
            if needcontour==1:
                totalname = plotloc+'/CRplot/Tz/Tz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'_contour'+tlabel+'ag.pdf'
            else:
                totalname = plotloc+'/CRplot/Tz/Tz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+tlabel+'ag.pdf'
        if wanted == 'rhoz':
            if needcontour==1:
                totalname = plotloc+'/CRplot/rhoz/rhoz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'_contour'+tlabel+'ag.pdf'
            else:
                totalname = plotloc+'/CRplot/rhoz/rhoz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+tlabel+'ag.pdf'
        if wanted == 'vzz':
            totalname = plotloc+'/CRplot/vzz/vzz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+tlabel+'ag.pdf'
        if wanted == 'Tvz':
            totalname = plotloc+'/CRplot/Tvz/Tvz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+tlabel+'ag.pdf'
        if wanted == 'Tvznolog':
            totalname = plotloc+'/CRplot/Tvz/Tvznolog_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+tlabel+'ag.pdf'
        if wanted == 'pcrpth':
            totalname = plotloc+'/CRplot/pcrpth/pcrpth_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+tlabel+'ag.pdf'
        del Hadd
    return im, plotupt,plotit,cblabel,totalname
