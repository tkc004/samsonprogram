from stdmodandoption import *
import plot_setup as PS
import collections


def CR_rz_energy_distribution_volave(ssdict):
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
    HIonly=ssdict['HIonly']
    titleneed=title
    ptitle=title
    needlog=0
    dclabelneed=1
    correctIa=0
    useM1=0
    resoneed=0
    rotface=1
    newlabelneed=1
    findradiusatnism = 1 #find the radius that has that ISM density (negative to turn off)
    print 'runtodo', runtodo
    info=SSF.outdirname(runtodo, Nsnap)
    havecr=info['havecr']
    dclabel=info['dclabel']
    haveB=info['haveB']
    withinr = ssdict['withinr']
    withoutr = ssdict['withoutr']
    nogrid = ssdict['nogrid']
    maxlength=ssdict['maxlength'] #thickness 
    usehalfz=ssdict['usehalfz']
    addsurden=ssdict['addsurden']
    if wanted=='crdenr' or wanted=='crpr' or wanted=='crvr':
        rlist=np.linspace(withoutr,withinr,num=nogrid)
        if havecr>0: credenlist=rlist*0.
        if haveB>0:  Bdenlist=rlist*0.
        thdenlist = rlist*0.
        denlist = rlist*0.
        turdenlist = rlist*0.
        gsur=rlist*0.
    if wanted=='crdenz' or wanted=='crpz' or wanted=='crdpz' or wanted=='crdpzcox' or wanted=='crpzcox': 
        if usehalfz==1:
            zlist=np.linspace(0.0,maxlength/2.,num=nogrid)
        else:
            zlist=np.linspace(-maxlength/2.,maxlength/2.,num=nogrid)
        if havecr>0: credenlist=zlist*0.
        if haveB>0:  Bdenlist=zlist*0.
        thdenlist = zlist*0.
        turdenlist = zlist*0.        
    numoftimes=0
    snaplist=[]
    if wanted=='crdenr' or wanted=='crpr' or wanted=='crdenz' or wanted=='crpz'\
    or wanted=='crvr'\
    or wanted=='crdpz' or wanted=='crdpzcox' or wanted=='crpzcox':
        for isnap in range(startno,Nsnap+1,snapsep):
            info=SSF.outdirname(runtodo, isnap)
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
            labelneed=dclabel
            if newlabelneed==1:
                labelneed="\n".join(wrap(newlabel,17))
            ptitle=labelneed
            if cosmo==1:
                h0=1
            else:
                h0=0
            egydata=CRTF.energyoutput(runtodo,isnap,
                                 shiftz=1,usesolarcircle=0,
                                 rotface=1,usecentral=0,cylr=withinr,
                                 cylrin=withoutr,cylz=maxlength,nbin = 7, HIonly=HIonly)   
            turl = egydata['turl']; xcell=egydata['xcell']; ycell=egydata['ycell']; zcell=egydata['zcell'];
            gvcell = egydata['gvcell']
            Gcrxy = np.sqrt((xcell)*(xcell)+(ycell)*(ycell))
            if wanted=='crdenr' or wanted=='crpr' or wanted=='crvr': 
                cutcz = np.absolute(zcell) < maxlength/2.
                for i in range(len(rlist)-1):
                    cutcrxy = (Gcrxy < rlist[i+1]) & (Gcrxy > rlist[i])
                    cutc = cutcrxy*cutcz
                    turl_in_erg_cut = np.sum(turl[cutc])
                    gvcellcut = np.sum(gvcell[cutc])
                    turl_in_erg_per_cm3 = turl_in_erg_cut/gvcellcut
                    turdenlist[i] += turl_in_erg_per_cm3
            if wanted=='crdenz' or wanted=='crpz' or wanted=='crdpz' or wanted=='crpzcox' or wanted=='crdpzcox':
                cutcrxy = np.logical_and(Gcrxy<withinr,Gcrxy>withoutr)
                for i in range(len(zlist)-1):
                    cutcz = (zcell < zlist[i+1]) & (zcell > zlist[i])
                    cutc = cutcrxy*cutcz
                    turl_in_erg_cut = np.sum(turl[cutc])
                    gvcellcut = np.sum(gvcell[cutc])
                    turl_in_erg_per_cm3 = turl_in_erg_cut/gvcellcut
                    turdenlist[i] += turl_in_erg_per_cm3
            if wanted=='crpzcox' or wanted=='crdpzcox':
                Evzl=egydata['Evzl']; Begy_2Bzl=egydata['Begy_2Bzl'];
            therml=egydata['therml']; Begyl=egydata['Begyl']; cregyl=egydata['cregyl'];
            Gxpl=egydata['Gxpl']; Gypl=egydata['Gypl']; Gzpl=egydata['Gzpl']; 
            gminl=egydata['gminl']; voll=egydata['voll']
            Grxy = np.sqrt((Gxpl)*(Gxpl)+(Gypl)*(Gypl))
            if wanted=='crdenr' or wanted=='crpr' or wanted=='crvr':         
                cutz = np.absolute(Gzpl) < maxlength/2.
                for i in range(len(rlist)-1):
                    cutrxy = (Grxy < rlist[i+1]) & (Grxy > rlist[i])
                    cut = cutrxy*cutz
                    vollcut = np.sum(voll[cut])
                    if havecr>0: 
                        cregy_in_erg_cut = np.sum(cregyl[cut])
                        creden_in_erg_per_cm3 = cregy_in_erg_cut/vollcut
                        credenlist[i]+=creden_in_erg_per_cm3
                    if haveB>0: 
                        Begy_in_erg_cut = np.sum(Begyl[cut])
                        Bden_in_erg_per_cm3 = Begy_in_erg_cut/vollcut
                        Bdenlist[i]+=Bden_in_erg_per_cm3
                    Gm_in_g_cut = np.sum(gminl[cut])
                    den_in_g_per_cm3 = Gm_in_g_cut/vollcut
                    denlist[i]+=den_in_g_per_cm3
                    thden_in_erg_cut = np.sum(therml[cut])
                    thden_in_erg_per_cm3 = thden_in_erg_cut/vollcut
                    thdenlist[i]+=thden_in_erg_per_cm3
            if wanted=='crdenz' or wanted=='crpz' or wanted=='crdpz' or wanted=='crdpzcox' or wanted=='crpzcox':        
                cutrxy = np.logical_and(Grxy<withinr,Grxy>withoutr)
                for i in range(len(zlist)-1):
                    cutz = (Gzpl < zlist[i+1]) & (Gzpl > zlist[i])
                    cut = cutrxy*cutz
                    vollcut = np.sum(voll[cut])
                    if havecr>0: 
                        cregy_in_erg_cut = np.sum(cregyl[cut])
                        creden_in_erg_per_cm3 = cregy_in_erg_cut/vollcut
                        credenlist[i]+=creden_in_erg_per_cm3
                    if haveB>0:
                        if wanted=='crpzcox' or wanted=='crdpzcox':
                            Begy_in_erg_cut = np.sum(Begy_2Bzl[cut])
                        else:
                            Begy_in_erg_cut = np.sum(Begyl[cut])
                        Bden_in_erg_per_cm3 = Begy_in_erg_cut/vollcut
                        Bdenlist[i]+=Bden_in_erg_per_cm3
                    thden_in_erg_cut = np.sum(therml[cut])
                    thden_in_erg_per_cm3 = thden_in_erg_cut/vollcut
                    thdenlist[i]+=thden_in_erg_per_cm3
            if wanted=='crpr' and addsurden==1:
                G=egydata['G']
                Gsurdenlist = SSF.calsurden(G,rlist,maxlength)
                print 'Gsurdenlist', Gsurdenlist
                Gsurdenlist *= Msun_in_g/pc_in_cm/pc_in_cm #in cgs
                S = SSF.readsnapwcen(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix,\
                 havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
                 loccen=1,runtodo=runtodo,rundir=rundir,halostr=halostr,importLcen=1,\
                 angLin=G['angL'],cenin=G['cen'],vcenin=G['vcen'])
                Ssurdenlist = SSF.calsurden(S,rlist,maxlength)
                print 'Ssurdenlist', Ssurdenlist
                Ssurdenlist *= Msun_in_g/pc_in_cm/pc_in_cm #in cgs
                gsur += np.pi*NewtonG_in_cgs*Gsurdenlist*(Ssurdenlist+Gsurdenlist) #in cgs
                print 'Gsurdenlist,Ssurdenlist', Gsurdenlist,Ssurdenlist
            numoftimes+=1
        if haveB>0:
                lsn='dashed'
        else:
                lsn='solid'
        if wanted=='crdenr' or wanted=='crpr' or wanted=='crvr': 
            rlistm = (rlist[:-1]+rlist[1:])/2.
        if wanted=='crdenz' or wanted=='crpz' or wanted=='crdpz' or wanted=='crpzcox' or wanted=='crdpzcox':         
            zlistm = (zlist[:-1]+zlist[1:])/2.
            zlistmm = (zlistm[:-1]+zlistm[1:])/2.
        if wanted=='crdenr' or wanted=='crpr' or wanted=='crvr':  
            plotdict[wanted]['xlab'] = r'${\rm r\;[kpc]}$'        
        if wanted=='crdenz' or wanted=='crpz' or wanted=='crdpz' or wanted=='crpzcox' or wanted=='crdpzcox':  
            plotdict[wanted]['xlab'] = r'${\rm z\;[kpc]}$'
        if wanted=='crdenr' or wanted=='crdenz':
            plotdict[wanted]['ylab'] = r'$e\;[{\rm erg/cm^3}]$'
        elif wanted=='crpr' or wanted=='crpz' or wanted=='crpzcox':
            plotdict[wanted]['ylab'] = r'$P\;[{\rm dyne/cm^2}]$'
        elif wanted=='crvr':
            plotdict[wanted]['ylab'] = r'$v\;[{\rm km/s}]$'
        elif wanted=='crdpz' or wanted=='crdpzcox':
            plotdict[wanted]['ylab'] = r'${\rm d} p/{\rm d} z\;[{\rm dyne/cm^3}]$'
        if wanted=='crdenr' or wanted=='crpr' or wanted=='crvr':   
            plotdict[wanted]['xnl']['etur'] = rlistm
        if wanted=='crdenz' or wanted=='crpz' or wanted=='crpzcox':   
            plotdict[wanted]['xnl']['etur'] = zlistm
        if wanted=='crdpz' or wanted=='crdpzcox':
            plotdict[wanted]['xnl']['etur'] = zlistmm
        if wanted=='crdenr' or wanted=='crdenz':
            plotdict[wanted]['ynl']['etur'] = turdenlist[:-1]/numoftimes
        elif wanted=='crpr' or wanted=='crpz' or wanted=='crpzcox':
            #0.5*rho*sigmaz^2 is the turbulence energy, rho*sigmaz^2 is the vertical turbulence pressure
            plotdict[wanted]['ynl']['etur'] = turdenlist[:-1]/numoftimes*2.0
        elif wanted=='crvr':
            plotdict[wanted]['ynl']['etur'] = np.sqrt(turdenlist[:-1]/denlist[:-1]*2.0)/km_in_cm
        elif wanted=='crdpz':
            plotdict[wanted]['ynl']['etur'] = np.absolute((turdenlist[1:-1]-turdenlist[:-2])*2.0/(zlist[1:-1]-zlist[:-2])/kpc_in_cm)/numoftimes
        elif wanted=='crdpzcox':
            plotdict[wanted]['ynl']['etur'] = np.absolute((turdenlist[1:-1]-turdenlist[:-2])*2.0/(zlist[1:-1]-zlist[:-2])/kpc_in_cm)/3.0/numoftimes
        plotdict[wanted]['lw']['etur'] = 2    
        plotdict[wanted]['marker']['etur'] = 'd'    
        plotdict[wanted]['lsn']['etur'] = 'solid'
        if wanted=='crpzcox' or wanted=='crdpzcox':
            plotdict[wanted]['linelab']['etur'] = r'$\rho \sigma_z^2$'
        else:
            plotdict[wanted]['linelab']['etur'] = 'turbulence'            
        if wanted=='crdenr' or wanted=='crpr' or wanted=='crvr':        
            plotdict[wanted]['xnl']['eth'] = rlistm
        if wanted=='crdenz' or wanted=='crpz' or wanted=='crpzcox':        
            plotdict[wanted]['xnl']['eth'] = zlistm
        if wanted=='crdpz' or wanted=='crdpzcox':
            plotdict[wanted]['xnl']['eth'] = zlistmm
        if wanted=='crdenr' or wanted=='crdenz':
            plotdict[wanted]['ynl']['eth'] = thdenlist[:-1]/numoftimes
        elif wanted=='crpr' or wanted=='crpz' or wanted=='crpzcox':
            plotdict[wanted]['ynl']['eth'] = thdenlist[:-1]*(GAMMA-1.)/numoftimes
        elif wanted=='crvr':
            plotdict[wanted]['ynl']['eth'] = np.sqrt(thdenlist[:-1]*(GAMMA-1.)*GAMMA/denlist[:-1])/km_in_cm
        elif wanted=='crdpz' or wanted=='crdpzcox':
            plotdict[wanted]['ynl']['eth'] = np.absolute((thdenlist[1:-1]-thdenlist[:-2])*(GAMMA-1.)/(zlist[1:-1]-zlist[:-2])/kpc_in_cm)/numoftimes
        plotdict[wanted]['lw']['eth'] = 2    
        plotdict[wanted]['marker']['eth'] = '^'
        plotdict[wanted]['lsn']['eth'] = 'dashed'
        plotdict[wanted]['linelab']['eth'] = 'thermal'
        if haveB>0:
            if wanted=='crdenr' or wanted=='crpr':
                plotdict[wanted]['xnl']['eB'] = rlistm
            if wanted=='crdenz' or wanted=='crpz' or wanted=='crpzcox':
                plotdict[wanted]['xnl']['eB'] = zlistm
            if wanted=='crdpz' or wanted=='crdpzcox':
                plotdict[wanted]['xnl']['eB'] = zlistmm
            if wanted=='crdenr' or wanted=='crdenz':
                plotdict[wanted]['ynl']['eB'] = Bdenlist[:-1]/numoftimes
            elif wanted=='crpr' or wanted=='crpz' or wanted=='crpzcox':
                plotdict[wanted]['ynl']['eB'] = Bdenlist[:-1]/numoftimes
            elif wanted=='crdpz' or wanted=='crdpzcox':
                plotdict[wanted]['ynl']['eB'] = np.absolute((Bdenlist[1:-1]-Bdenlist[:-2])/(zlist[1:-1]-zlist[:-2])/kpc_in_cm)/numoftimes           
            plotdict[wanted]['lw']['eB'] = 1
            plotdict[wanted]['lsn']['eB'] = 'solid'        
            plotdict[wanted]['marker']['eB'] = 's'
            if wanted=='crpzcox' or wanted=='crdpzcox':
                plotdict[wanted]['linelab']['eB'] = r'$(B^2-2B_z^2)/(8\pi)$'
            else:
                plotdict[wanted]['linelab']['eB'] = 'Bfield'
        if havecr>0: 
            if wanted=='crdenr' or wanted=='crpr' or wanted=='crvr':
                plotdict[wanted]['xnl']['ecr'] = rlistm
            if wanted=='crdenz' or wanted=='crpz' or wanted=='crpzcox':
                plotdict[wanted]['xnl']['ecr'] = zlistm
            if wanted=='crdpz' or wanted=='crdpzcox':
                plotdict[wanted]['xnl']['ecr'] = zlistmm
            if wanted=='crdenr' or wanted=='crdenz':
                plotdict[wanted]['ynl']['ecr'] = credenlist[:-1]/numoftimes
            elif wanted=='crpr' or wanted=='crpz' or wanted=='crpzcox':
                plotdict[wanted]['ynl']['ecr'] = credenlist[:-1]*(CRgamma-1.)/numoftimes
            elif wanted=='crvr':
                plotdict[wanted]['ynl']['ecr'] = np.sqrt(credenlist[:-1]*(CRgamma-1.)*CRgamma/denlist[:-1])/km_in_cm/numoftimes
            elif wanted=='crdpz' or wanted=='crdpzcox':
                plotdict[wanted]['ynl']['ecr'] = np.absolute((credenlist[1:-1]-credenlist[:-2])*(CRgamma-1.)/(zlist[1:-1]-zlist[:-2])/kpc_in_cm)/numoftimes 
            plotdict[wanted]['lw']['ecr'] = 1  
            plotdict[wanted]['lsn']['ecr'] = 'dashdot' 
            plotdict[wanted]['linelab']['ecr'] = 'CR' 
            plotdict[wanted]['marker']['ecr'] = 'o'
        if wanted=='crpr' and addsurden==1:
            plotdict[wanted]['xnl']['gsur'] = rlistm
            plotdict[wanted]['ynl']['gsur'] = gsur[:-1]/numoftimes
            plotdict[wanted]['linelab']['gsur'] = 'gravity' 
            plotdict[wanted]['lw']['gsur'] = 3    
            plotdict[wanted]['marker']['gsur'] = 'None'    
            plotdict[wanted]['lsn']['gsur'] = 'dotted'            
        plotdict[wanted]['runtodo'] = runtodo
        plotdict[wanted]['labelneed'] = labelneed
        #plotdict[wanted]['lsn'] = lsn
        plotdict[wanted]['color'] = color
        plotdict[wanted]['runtitle'] = runtitle
        plotdict[wanted]['ptitle'] = ptitle
        if wanted=='crdenz':
            filename=plotloc+'CRplot/crdenz/CRz_volave_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        elif wanted=='crpz':
            filename=plotloc+'CRplot/crdenz/CRpz_volave_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'        
        elif wanted=='crdenr':
            filename=plotloc+'CRplot/crdenr/CRr_volave_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        elif wanted=='crpr':
            filename=plotloc+'CRplot/crdenr/CRpr_volave_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        elif wanted=='crvr':
            filename=plotloc+'CRplot/crdenr/CRvr_volave_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        elif wanted=='crdpz':
            filename=plotloc+'CRplot/crdenz/CRdpz_volave_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        elif wanted=='crpzcox':
            filename=plotloc+'CRplot/crpzcox/CRpzcox_volave_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        elif wanted=='crdpzcox':
            filename=plotloc+'CRplot/crdenz/CRdpzcox_volave_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        plotdict[wanted]['filename'] = filename
        return plotdict