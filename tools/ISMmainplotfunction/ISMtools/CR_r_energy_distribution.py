from stdmodandoption import *
import plot_setup as PS
import collections


def CR_r_energy_distribution(ssdict):
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
    rlist=np.linspace(withoutr,withinr,num=nogrid)
    if havecr>0: credenlist=rlist*0.
    if haveB>0:  Bdenlist=rlist*0.
    thdenlist = rlist*0.
    turdenlist = rlist*0.
    numoftimes=0
    snaplist=[]
    if wanted=='crdenr' or wanted=='crpr':
        info=SSF.outdirname(runtodo, Nsnap)
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
        if cosmo==1:
            datasup=0;
        else:
            datasup=1;
        egydata=CRTF.energyoutput(runtodo,Nsnap,
                             shiftz=1,usesolarcircle=0,
                             rotface=1,usecentral=0,cylr=withinr,
                             cylrin=withoutr,cylz=maxlength,nbin = 5)   
        turl = egydata['turl']; xcell=egydata['xcell']; ycell=egydata['ycell']; zcell=egydata['zcell'];    
        cutcz = np.absolute(zcell) < maxlength/2.
        Gcrxy = np.sqrt((xcell)*(xcell)+(ycell)*(ycell))
        for i in range(len(rlist)-1):
            shellcvol = np.pi*(np.square(rlist[i+1])-np.square(rlist[i]))*maxlength
            shellcvol_in_cm3 = shellcvol*kpc_in_cm*kpc_in_cm*kpc_in_cm
            cutcrxy = (Gcrxy < rlist[i+1]) & (Gcrxy > rlist[i])
            cutc = cutcrxy*cutcz
            turl_in_erg_cut = np.sum(turl[cutc])
            turl_in_erg_per_cm3 = turl_in_erg_cut/shellcvol_in_cm3
            turdenlist[i] = turl_in_erg_per_cm3

        therml=egydata['therml']; Begyl=egydata['Begyl']; cregyl=egydata['cregyl'];
        Gxpl=egydata['Gxpl']; Gypl=egydata['Gypl']; Gzpl=egydata['Gzpl'];        
        cutz = np.absolute(Gzpl) < maxlength/2.
        Grxy = np.sqrt((Gxpl)*(Gxpl)+(Gypl)*(Gypl))
        for i in range(len(rlist)-1):
            shellvol = np.pi*(np.square(rlist[i+1])-np.square(rlist[i]))*maxlength
            shellvol_in_cm3 = shellvol*kpc_in_cm*kpc_in_cm*kpc_in_cm
            cutrxy = (Grxy < rlist[i+1]) & (Grxy > rlist[i])
            cut = cutrxy*cutz
            if havecr>0: cregy_in_erg_cut = np.sum(cregyl[cut])
            if havecr>0: creden_in_erg_per_cm3 = cregy_in_erg_cut/shellvol_in_cm3
            if havecr>0: credenlist[i]=creden_in_erg_per_cm3
            if haveB>0: Begy_in_erg_cut = np.sum(Begyl[cut])
            if haveB>0: Bden_in_erg_per_cm3 = Begy_in_erg_cut/shellvol_in_cm3
            if haveB>0: Bdenlist[i]=Bden_in_erg_per_cm3
            thden_in_erg_cut = np.sum(therml[cut])
            thden_in_erg_per_cm3 = thden_in_erg_cut/shellvol_in_cm3
            thdenlist[i]=thden_in_erg_per_cm3        
        if haveB>0:
                lsn='dashed'
        else:
                lsn='solid'
        plotdict[wanted]['xlab'] = r'${\rm r\;[kpc]}$'
        if wanted=='crdenr':
            plotdict[wanted]['ylab'] = r'$e\;[{\rm erg/cm^3}]$'
        elif wanted=='crpr':
            plotdict[wanted]['ylab'] = r'$p\;[{\rm dyne/cm^2}]$'
        plotdict[wanted]['xnl']['etur'] = rlist[:-1]
        if wanted=='crdenz':
            plotdict[wanted]['ynl']['etur'] = turdenlist[:-1]
        elif wanted=='crpz':
            # https://doi.org/10.1111/j.1365-2966.2011.18550.x not used
            plotdict[wanted]['ynl']['etur'] = turdenlist[:-1]
        plotdict[wanted]['lw']['etur'] = 2    
        plotdict[wanted]['marker']['etur'] = 'd'    
        plotdict[wanted]['lsn']['etur'] = 'solid'
        plotdict[wanted]['linelab']['etur'] = 'turbulence'
        plotdict[wanted]['xnl']['eth'] = rlist[:-1]
        if wanted=='crdenr':
            plotdict[wanted]['ynl']['eth'] = thdenlist[:-1]
        elif wanted=='crpr':
            plotdict[wanted]['ynl']['eth'] = thdenlist[:-1]*(GAMMA-1.)
        plotdict[wanted]['lw']['eth'] = 2    
        plotdict[wanted]['marker']['eth'] = '^'
        plotdict[wanted]['lsn']['eth'] = 'dashed'
        plotdict[wanted]['linelab']['eth'] = 'thermal'
        if haveB>0:
            if wanted=='crdenr':
                plotdict[wanted]['ynl']['eB'] = Bdenlist[:-1]
            elif wanted=='crpr':
                plotdict[wanted]['ynl']['eB'] = Bdenlist[:-1]
            plotdict[wanted]['lw']['eB'] = 1
            plotdict[wanted]['lsn']['eB'] = 'solid'        
            plotdict[wanted]['marker']['eB'] = 's'
            plotdict[wanted]['linelab']['eB'] = 'Bfield'    
        if havecr>0: 
            plotdict[wanted]['xnl']['ecr'] = rlist[:-1]
            if wanted=='crdenr':
                plotdict[wanted]['ynl']['ecr'] = credenlist[:-1]
            elif wanted=='crpr':
                plotdict[wanted]['ynl']['ecr'] = credenlist[:-1]*(CRgamma-1.)
            plotdict[wanted]['lw']['ecr'] = 1  
            plotdict[wanted]['lsn']['ecr'] = 'dashdot' 
            plotdict[wanted]['linelab']['ecr'] = 'CR' 
            plotdict[wanted]['marker']['ecr'] = 'o'
        plotdict[wanted]['runtodo'] = runtodo
        plotdict[wanted]['labelneed'] = labelneed
        #plotdict[wanted]['lsn'] = lsn
        plotdict[wanted]['color'] = color
        plotdict[wanted]['runtitle'] = runtitle
        plotdict[wanted]['ptitle'] = ptitle
        if wanted=='crdenr':
            filename=plotloc+'CRplot/crdenr/CRr_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        elif wanted=='crpr':
            filename=plotloc+'CRplot/crdenr/CRpr_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        plotdict[wanted]['filename'] = filename
        return plotdict