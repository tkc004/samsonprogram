from pathloc import *



def crhisto(runtodo,i,wanted,nooftimes,Tcut,highTcut,vcut,vhcut,withinr,zup,zdown,userad,trackgas=0):
    if trackgas==1:
        vcut = -1e10 #if we track gas particles, we should consider all velocities
    if nooftimes==1:
        return None
    data = gasphase(runtodo,i,Tcut=Tcut,highTcut=highTcut,\
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
    if wanted=='pcrpth':
        pth = (GAMMA-1.0)*Gu*u_codetocgs*(rho*rho_codetocgs) #cgs
        x = np.log10(pth)
        print 'x', x
    if wanted=='Tz' or wanted=='Tvz':
        y = np.log10(TrueTemp)
    if wanted=='rhoz':
        y = np.log10(converted_rho)
    if wanted=='vzz':
        if userad==1:
            y = np.log10(vr)
        else:
            y = np.log10(vz)
    if wanted=='pcrpth':
        pcr = (CRgamma-1.0)*cregy*cregy_codetocgs*(rho*rho_codetocgs)/(Gmass*m_codetocgs) #cgs
        y = np.log10(pcr)
        print 'y', y
    gridxx = np.linspace(extent[0],extent[1],nobin)
    gridyy = np.linspace(extent[2],extent[3],nobin)
    H, xedges, yedges = np.histogram2d(x, y, bins=[gridxx, gridyy],weights=Gmass*1e10)
    return H


def crgridfunction(gridag,runtodo,wanted,idir,plotupt,plotit,startno,Nsnap,snapsep,cbcolor='hot'):
    if wanted=='Tz' or wanted=='rhoz' or wanted=='vzz' or wanted=='Tvz' or wanted=='pcrpth' or wanted=='rhoT':
        newlabelneed=1
        from crtestfunction import *
        from textwrap import wrap
        Tcut=-1.0 #K
        highTcut = 1e10
        vcut=-1000.0 #outflow velocity cut
        vhcut = 1e10 #upper outflow velocity cut
        withinr=8.0
        withoutr=0.0
        zup=0.5
        zdown=-0.5
        trackgas=0
        tlabel=''
        userad=0
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
            extent = [zdown,zup,0,3.2]
        if wanted=='Tvz':
            extent = [0,3.2,1,9]
            #extent = [0,9,0,9]
        if wanted=='pcrpth':
            extent = [-20,-11,-20,-11]
        if wanted=='rhoT':
            extent = [-7,4,1,9]
        nobin=51
        needcontour=1
        info=outdirname(runtodo, Nsnap)
        cosmo=info['cosmo']
        runtitle=info['runtitle']
        havecr=info['havecr']
        dclabel=info['dclabel']
        newlabel=info['newlabel']
        if wanted=='pcrpth' and havecr==0:
            return None
        ptitle=title
        if runtitle=='SMC':
            ptitle='Dwarf'
        elif runtitle=='SBC':
            ptitle='Starburst'
        elif runtitle=='MW':
            ptitle=r'$L\star$ Galaxy'
        ptitle=''
        labelneed=dclabel
        if newlabelneed==1:
            labelneed="\n".join(wrap(newlabel,17))
        plotupt=np.append(plotupt,ptitle)
        plotit =np.append(plotit,labelneed)
        Hadd = np.zeros((nobin-1,nobin-1))
        if trackgas==1:
            data = gasphase(runtodo,snaptrack,Tcut=Tcut_t,highTcut=highTcut_t,\
            vcut=vcut_t,vhcut=vhcut_t,withinr=withinr_t,withoutr=withoutr_t,zup=zup_t,zdown=zdown_t,userad=userad)
            Gid_t = data['Gid']
        nooftimes=0
        for i in range(startno,Nsnap,snapsep):
            if trackgas==1:
                vcut = -1e10 #if we track gas particles, we should consider all velocities
            if nooftimes==1:
                return None
            data = gasphase(runtodo,i,Tcut=Tcut,highTcut=highTcut,\
            vcut=vcut,vhcut=vhcut,withinr=withinr,withoutr=withoutr,zup=zup,zdown=zdown,userad=userad)
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
            if wanted=='pcrpth':
                pth = (GAMMA-1.0)*Gu*u_codetocgs*(rho*rho_codetocgs) #cgs
                x = np.log10(pth)
                print 'x', x
            if wanted=='rhoT':
                x = np.log10(converted_rho)
            if wanted=='rhoT':
                y = np.log10(TrueTemp)            
            if wanted=='Tz' or wanted=='Tvz':
                y = np.log10(TrueTemp)
            if wanted=='rhoz':
                y = np.log10(converted_rho)
            if wanted=='vzz':
                if userad==1:
                    y = np.log10(vr)
                else:
                    y = np.log10(vz)
            if wanted=='pcrpth':
                pcr = (CRgamma-1.0)*cregy*cregy_codetocgs*(rho*rho_codetocgs)/(Gmass*m_codetocgs) #cgs
                y = np.log10(pcr)
                print 'y', y
            gridxx = np.linspace(extent[0],extent[1],nobin)
            gridyy = np.linspace(extent[2],extent[3],nobin)
            H, xedges, yedges = np.histogram2d(x, y, bins=[gridxx, gridyy],weights=Gmass*1e10)
            H=H.T
            Hadd += H
            nooftimes += 1
        Hadd=Hadd/nooftimes
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
        if wanted=='rhoT':
            gridag[idir].set_xlabel(r"$\mathrm{Log}_{10} (n {\rm [cm^{-3}]})$")
        if wanted=='rhoT':
            gridag[idir].set_ylabel(r"$\mathrm{Log}_{10} (T {\rm [K]})$")
        if wanted=='pcrpth':
            gridag[idir].set_xlabel(r"$\mathrm{Log}_{10} (P_{\rm th} {\rm [erg/cm^3]})$")
            gridag[idir].set_ylabel(r"$\mathrm{Log}_{10} (P_{\rm cr} {\rm [erg/cm^3]})$")
        if wanted=='rhoz':
            gridag[idir].set_ylabel(r"$\mathrm{Log}_{10} (n {\rm [cm^{-3}]})$")
        if wanted=='Tz' or wanted=='Tvz':
            gridag[idir].set_ylabel(r"$\mathrm{Log}_{10} (T {\rm [K]})$")
        if wanted=='vzz':
            if userad==1:
                gridag[idir].set_ylabel(r"$\mathrm{Log}_{10} (v_r {\rm [km/s]})$")
            else:
                gridag[idir].set_ylabel(r"$\mathrm{Log}_{10} (v_z {\rm [km/s]})$")
        if wanted == 'rhoT':
            if needcontour==1:
                totalname = plotloc+'/CRplot/rhoT/rhoT_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+'_contour'+tlabel+'ag.pdf'
            else:
                totalname = plotloc+'/CRplot/rhoT/rhoT_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+tlabel+'ag.pdf'
        if wanted == 'Tz':
            if needcontour==1:
                totalname = plotloc+'/CRplot/Tz/Tz_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+'_contour'+tlabel+'ag.pdf'
            else:
                totalname = plotloc+'/CRplot/Tz/Tz_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+tlabel+'ag.pdf'
        if wanted == 'rhoz':
            if needcontour==1:
                totalname = plotloc+'/CRplot/rhoz/rhoz_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+'_contour'+tlabel+'ag.pdf'
            else:
                totalname = plotloc+'/CRplot/rhoz/rhoz_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+tlabel+'ag.pdf'
        if wanted == 'vzz':
            totalname = plotloc+'/CRplot/vzz/vzz_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+tlabel+'ag.pdf'
        if wanted == 'Tvz':
            totalname = plotloc+'/CRplot/Tvz/Tvz_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+tlabel+'ag.pdf'
        if wanted == 'pcrpth':
            totalname = plotloc+'/CRplot/pcrpth/pcrpth_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+tlabel+'ag.pdf'
        del Hadd
    return im, plotupt,plotit,cblabel,totalname
