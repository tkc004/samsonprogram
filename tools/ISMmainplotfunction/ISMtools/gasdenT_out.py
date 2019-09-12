from stdmodandoption import *
import cameron_functions as CF
import collections

def gasdenT_out(ssdict):    
    nested_dict = lambda: collections.defaultdict(nested_dict)
    plotdict = nested_dict()
    runtodo=ssdict['runtodo']
    wanted=ssdict['wanted']
    print 'wanted', wanted
    startno=ssdict['startno']
    Nsnap=ssdict['Nsnap']
    snapsep=ssdict['snapsep']
    the_prefix=ssdict['the_prefix']
    the_suffix=ssdict['the_suffix']
    fmeat=ssdict['fmeat']
    maxlength=ssdict['maxlength'] #thickness
    withinr=ssdict['withinr']
    withoutr = ssdict['withoutr']
    title = ssdict['title']
    phicut = ssdict['phicut']
    usehz = ssdict['usehz']
    titleneed=title
    nogrid=20
    xaxis_snapno=0
    if not (wanted=='gasdenTz' or wanted=='gasdenTr'):
        print 'wrong wanted'
        return None
    resoneed=0
    rotface=1
    loccen=1
    newlabelneed=1
    print 'runtodo', runtodo
    if wanted=='gasdenTz':
        dr = withinr/nogrid
        #zlist = np.linspace(-maxlength/2.,maxlength/2.,num=nogrid)
        #zlist = SSF.symlogspace(maxlength/2.,nogrid/2,base=50.0)
        zlist = SSF.symlogspace(maxlength/2.,nogrid,base=30.0,halfz=usehz)
        mgzlist=zlist*0.
        rhogzlist=zlist*0.
        rhommgzlist=zlist*0.
        rhocoldgzlist=zlist*0.
        rhowarmgzlist=zlist*0.
        rhowimgzlist=zlist*0.
        rhownmgzlist=zlist*0.
        rhohotgzlist=zlist*0.
    if wanted=='gasdenTr':
        dr = withinr/nogrid
        rlist = np.linspace(withoutr,withinr,num=nogrid)
        mgrlist=rlist*0.
        rhogrlist=rlist*0.
        rhocoldgrlist=rlist*0.
        rhowarmgrlist=rlist*0.
        rhowimgrlist=rlist*0.
        rhownmgrlist=rlist*0.
        rhohotgrlist=rlist*0.        
    numoftimes=0
    for i in range(startno,Nsnap,snapsep):
        snaplist=[]
        info=SSF.outdirname(runtodo, i)
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
        labelneed=dclabel
        if newlabelneed==1:
            labelneed="\n".join(wrap(newlabel,17))
        if cosmo==1:
            h0=1
        else:
            h0=0
        Gextra = SSF.readsnapwcen(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix,\
         havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
         loccen=loccen,runtodo=runtodo,rundir=rundir,halostr=halostr)
        Gx = Gextra['p'][:,0]; Gy = Gextra['p'][:,1]; Gz = Gextra['p'][:,2];
        Grxy = np.sqrt(Gx*Gx+Gy*Gy)
        Gvx = Gextra['v'][:,0]; Gvy = Gextra['v'][:,1]; Gvz = Gextra['v'][:,2];
        Grho = Gextra['rho']; Gu = Gextra['u']; Gm = Gextra['m']; Neb = Gextra['ne'];
        Nnh = Gextra['nh']; Gh = Gextra['h']; Gmetal = Gextra['z'];
        Gh = Gextra['h']; Gmetal = Gextra['z'];
        TrueTemp, converted_rho  = SF.convertTemp(Gu, Neb, Grho)
        datanH = CF.calnH(Gm,Nnh,Gh,Grho,Gmetal)
        mHII = datanH['mHII_in_g']; mH2=datanH['mH2_in_g']; mHI=datanH['mHI_in_g']
        Zmetal = Gmetal[:,0] #metal mass fraction (everything not H, He)
        ZHe = Gmetal[:,1] # He mass fraction
        Hmfrac = 1.0-Zmetal-ZHe #H mass fraction
        if wanted=='gasdenTz':
            for ir in range(len(zlist)-1):
                cutxy = (Grxy < withinr) & (Grxy > withoutr)
                cutz = (Gz < zlist[ir+1]) & (Gz > zlist[ir])
                cutxyz = cutxy*cutz             
                if phicut==1:
                    Gxcutxy=Gx[cutxy]; Gycutxy=Gy[cutxy];
                    phi = np.arctan2(Gycutxy,Gxcutxy)
                    print 'phi', phi
                    cutphi = np.absolute(phi) < np.pi/10.0
                    cutz=cutz*cutphi
                Gmcutz = Gm[cutxyz]*m_codetocgs
                Grhocutz = Grho[cutxyz]*rho_codetocgs
                Tcutz = TrueTemp[cutxyz]
                mH2cutz = mH2[cutxyz]
                Nnhcutz = Nnh[cutxyz]
                Hfcutz = Hmfrac[cutxyz]
                GHmcutz = Gmcutz*Hfcutz
                mH2cutz = mH2[cutxyz]
                hotcut = Tcutz>1e5 #K
                warmcut = (Tcutz>1e3)&(Tcutz<1e5)
                coldcut = Tcutz<1e3
                mgzlist[ir] += np.sum(GHmcutz)
                rhogzlist[ir] += np.sum(GHmcutz)/np.sum(Gmcutz/Grhocutz) #in g/cm^3
                rhommgzlist[ir] += np.sum(mH2cutz)/np.sum(Gmcutz/Grhocutz) #in g/cm^3
                rhocoldgzlist[ir] += (np.sum(GHmcutz[coldcut])-np.sum(mH2cutz))/np.sum(Gmcutz/Grhocutz) #in g/cm^3
                rhowarmgzlist[ir] += np.sum(GHmcutz[warmcut])/np.sum(Gmcutz/Grhocutz) #in g/cm^3
                rhowimgzlist[ir] += np.sum(GHmcutz[warmcut]*(1.0-Nnhcutz[warmcut]))/np.sum(Gmcutz/Grhocutz) #in g/cm^3
                rhownmgzlist[ir] += np.sum(GHmcutz[warmcut]*Nnhcutz[warmcut])/np.sum(Gmcutz/Grhocutz) #in g/cm^3
                rhohotgzlist[ir] += np.sum(GHmcutz[hotcut])/np.sum(Gmcutz/Grhocutz) #in g/cm^3
        if wanted=='gasdenTr':
            for ir in range(len(rlist)-1):
                cutxy = (Grxy < rlist[ir+1]) & (Grxy > rlist[ir])
                Gmcutxy = Gm[cutxy]*m_codetocgs #in g
                Grhocutxy = Grho[cutxy]*rho_codetocgs #in g/cm^3
                Gzcutxy = Gz[cutxy]
                Tcutxy = TrueTemp[cutxy]
                Nnhcutxy = Nnh[cutxy]
                cutz = np.absolute(Gzcutxy) < maxlength/2.
                Gmcutz = Gmcutxy[cutz]
                Grhocutz = Grhocutxy[cutz]
                Tcutz = Tcutxy[cutz]
                Nnhcutz = Nnhcutxy[cutz]
                if phicut==1:
                    Gxcutxy=Gx[cutxy]; Gycutxy=Gy[cutxy];
                    phi = np.arctan2(Gycutxy,Gxcutxy)
                    cutphi = np.absolute(phi) < np.pi/10.0
                    cutz=cutz*cutphi
                hotcut = Tcutz>1e5 #K
                warmcut = (Tcutz>1e3)&(Tcutz<1e5)
                coldcut = Tcutz<1e3
                mgrlist[ir] += np.sum(Gmcutz)
                rhogrlist[ir] += np.sum(Gmcutz)/np.sum(Gmcutz/Grhocutz) #in g/cm^3
                rhocoldgrlist[ir] += np.sum(Gmcutz[coldcut])/np.sum(Gmcutz/Grhocutz) #in g/cm^3
                rhowarmgrlist[ir] += np.sum(Gmcutz[warmcut])/np.sum(Gmcutz/Grhocutz) #in g/cm^3
                rhowimgrlist[ir] += np.sum(Gmcutz[warmcut]*(1.0-Nnhcutz[warmcut]))/np.sum(Gmcutz/Grhocutz) #in g/cm^3
                rhownmgrlist[ir] += np.sum(Gmcutz[warmcut]*Nnhcutz[warmcut])/np.sum(Gmcutz/Grhocutz) #in g/cm^3
                rhohotgrlist[ir] += np.sum(Gmcutz[hotcut])/np.sum(Gmcutz/Grhocutz) #in g/cm^3
        numoftimes+=1
        linelist = []
        if wanted=='gasdenTz':
            plotdict[wanted]['xlab'] = r'${\rm z\;[kpc]}$'
            zmlist = (zlist[:-1]+zlist[1:])/2.
            plotdict[wanted]['ylab'] = r'$\rho_{\rm gas}\;[{\rm g/cm^3}]$'
            plotdict[wanted]['ynl']['tot'] = rhogzlist[:-1]/numoftimes
            plotdict[wanted]['xnl']['tot'] = zmlist
            plotdict[wanted]['lsn']['tot'] = 'solid'
            plotdict[wanted]['linelab']['tot'] = 'Total'
            linelist.append('tot') 
            plotdict[wanted]['ynl']['mm'] = rhommgzlist[:-1]/numoftimes
            plotdict[wanted]['xnl']['mm'] = zmlist
            plotdict[wanted]['lsn']['mm'] = 'solid'
            plotdict[wanted]['linelab']['mm'] = 'MM'
            linelist.append('mm') 
            plotdict[wanted]['ynl']['cold'] = rhocoldgzlist[:-1]/numoftimes
            plotdict[wanted]['xnl']['cold'] = zmlist
            plotdict[wanted]['lsn']['cold'] = 'dashed'
            plotdict[wanted]['linelab']['cold'] = 'CNM'
            linelist.append('cold') 
            #plotdict[wanted]['linelab']['cold'] = r'${\rm T < 1e3K}$'
            #plotdict[wanted]['ynl']['warm'] = rhowarmgzlist[:-1]/numoftimes
            #plotdict[wanted]['xnl']['warm'] = zmlist
            #plotdict[wanted]['lsn']['warm'] = 'dashdot'
            #plotdict[wanted]['linelab']['warm'] = r'${\rm 1e3K < T < 1e5K}$'
            plotdict[wanted]['ynl']['wnm'] = rhownmgzlist[:-1]/numoftimes
            plotdict[wanted]['xnl']['wnm'] = zmlist
            plotdict[wanted]['lsn']['wnm'] = 'dashdot'
            #plotdict[wanted]['linelab']['wnm'] = r'${\rm 1e3K < T < 1e5K}$; neutral'
            plotdict[wanted]['linelab']['wnm'] = 'WNM'
            linelist.append('wnm') 
            plotdict[wanted]['ynl']['wim'] = rhowimgzlist[:-1]/numoftimes
            plotdict[wanted]['xnl']['wim'] = zmlist
            plotdict[wanted]['lsn']['wim'] = (0, (3, 5, 1, 5, 1, 5))
            #plotdict[wanted]['linelab']['wim'] = r'${\rm 1e3K < T < 1e5K}$; ionized'
            plotdict[wanted]['linelab']['wim'] = 'WIM'
            linelist.append('wim')             
            plotdict[wanted]['ynl']['hot'] = rhohotgzlist[:-1]/numoftimes
            plotdict[wanted]['xnl']['hot'] = zmlist
            plotdict[wanted]['lsn']['hot'] = 'dotted'
            #plotdict[wanted]['linelab']['hot'] = r'${\rm T > 1e5K}$'
            plotdict[wanted]['linelab']['hot'] = 'HIM'
            linelist.append('hot') 
        if wanted=='gasdenTr':
            plotdict[wanted]['xlab'] = r'${\rm r\;[kpc]}$'
            rmlist = (rlist[:-1]+rlist[1:])/2.
            plotdict[wanted]['ylab'] = r'$\rho_{\rm gas}\;[{\rm g/cm^3}]$'
            plotdict[wanted]['ynl']['tot'] = rhogrlist[:-1]/numoftimes
            plotdict[wanted]['xnl']['tot'] = rmlist
            plotdict[wanted]['lsn']['tot'] = 'solid'
            plotdict[wanted]['linelab']['tot'] = 'Total'
            linelist.append('tot') 
            plotdict[wanted]['ynl']['cold'] = rhocoldgrlist[:-1]/numoftimes
            plotdict[wanted]['xnl']['cold'] = rmlist
            plotdict[wanted]['lsn']['cold'] = 'dashed'
            #plotdict[wanted]['linelab']['cold'] = r'${\rm T < 1e3K}$'
            plotdict[wanted]['linelab']['cold'] = 'CNM'
            linelist.append('cold') 
            #plotdict[wanted]['ynl']['warm'] = rhowarmgrlist[:-1]/numoftimes
            #plotdict[wanted]['xnl']['warm'] = rmlist
            #plotdict[wanted]['lsn']['warm'] = 'dashdot'
            #plotdict[wanted]['linelab']['warm'] = r'${\rm 1e3K < T < 1e5K}$' 
            plotdict[wanted]['ynl']['wnm'] = rhownmgrlist[:-1]/numoftimes
            plotdict[wanted]['xnl']['wnm'] = rmlist
            plotdict[wanted]['lsn']['wnm'] = 'dashdot'
            #plotdict[wanted]['linelab']['wnm'] = r'${\rm 1e3K < T < 1e5K}$; neutral'
            plotdict[wanted]['linelab']['wnm'] = 'WNM'
            linelist.append('wnm') 
            plotdict[wanted]['ynl']['wim'] = rhowimgrlist[:-1]/numoftimes
            plotdict[wanted]['xnl']['wim'] = rmlist
            plotdict[wanted]['lsn']['wim'] = (0, (3, 5, 1, 5, 1, 5))
            #plotdict[wanted]['linelab']['wim'] = r'${\rm 1e3K < T < 1e5K}$; ionized'
            plotdict[wanted]['linelab']['wim'] = 'WIM'
            linelist.append('wim')
            plotdict[wanted]['ynl']['hot'] = rhohotgrlist[:-1]/numoftimes
            plotdict[wanted]['xnl']['hot'] = rmlist
            plotdict[wanted]['lsn']['hot'] = 'dotted'
            #plotdict[wanted]['linelab']['hot'] = r'${\rm T > 1e5K}$'
            plotdict[wanted]['linelab']['hot'] = 'HIM'
            linelist.append('hot') 
        plotdict[wanted]['runtodo'] = runtodo
        plotdict[wanted]['labelneed'] = labelneed
        plotdict[wanted]['lw'] = 2
        plotdict[wanted]['marker'] = 'None'
        plotdict[wanted]['color'] = color
        plotdict[wanted]['runtitle'] = runtitle
        plotdict[wanted]['ptitle'] = ptitle
        plotdict[wanted]['linelist']=linelist
        if wanted=='gasdenTz':
            filename=plotloc+'CRplot/gasdenTz/gasdenTz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'_wim.pdf'
            plotdict[wanted]['filename'] = filename      
        if wanted=='gasdenTr':
            filename=plotloc+'CRplot/gasdenTr/gasdenTr_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'_wim.pdf'
            plotdict[wanted]['filename'] = filename    
        return plotdict
