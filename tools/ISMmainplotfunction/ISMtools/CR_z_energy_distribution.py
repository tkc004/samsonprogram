from stdmodandoption import *
import plot_setup as PS
import collections


def CR_z_energy_distribution(ssdict):
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
    title='MW'
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
    withinr=5. #withinr will change according to the above
    dr=2.
    
    
    print 'runtodo', runtodo
    info=SSF.outdirname(runtodo, Nsnap)
    havecr=info['havecr']
    dclabel=info['dclabel']
    haveB=info['haveB']

    nogrid = 25
    maxlength=10. #thickness
    zlist=np.linspace(0.001,1,num=nogrid)
    if havecr>0: credenlist=zlist*0.
    if haveB>0:  Bdenlist=zlist*0.
    thdenlist = zlist*0.
    numoftimes=0
    snaplist=[]
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
    if cosmo==1:
        h0=1
    else:
        h0=0
    if cosmo==1:
        datasup=0;
    else:
        datasup=1;
    Gextra = SSF.readsnapwcen(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix,\
     havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
     datasup=datasup,runtodo=runtodo,rundir=rundir,halostr=halostr)
    Gx = Gextra['p'][:,0]; Gy = Gextra['p'][:,1]; Gz = Gextra['p'][:,2];
    Gvx = Gextra['v'][:,0]; Gvy = Gextra['v'][:,1]; Gvz = Gextra['v'][:,2];
    Grho = Gextra['rho']; Gu = Gextra['u']; Gm = Gextra['m']; 
    if haveB>0: GB = Gextra['B']; GBmag2  = GB[:,0]*GB[:,0]+GB[:,1]*GB[:,1]+GB[:,2]*GB[:,2]
    if haveB>0: Begy = GBmag2/8./np.pi/Grho*Gm*kpc_in_cm*kpc_in_cm*kpc_in_cm
    if havecr>0: cregy = Gextra['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
    if havecr>0: cregy_in_erg = cregy*solar_mass_in_g*1e10*km_in_cm*km_in_cm
    thden_in_erg = Gm*1e10*Msun_in_g*Gu*km_in_cm*km_in_cm

    if findradiusatnism>0.: withinr = SSF.findradwnism(Gextra, findradiusatnism)
    
    cutxy = (np.sqrt((Gx)*(Gx)+(Gy)*(Gy)) > withinr) & (np.sqrt((Gx)*(Gx)+(Gy)*(Gy)) < withinr+dr)
    zlist = np.linspace(-maxlength,maxlength,num=nogrid)
    for i in range(len(zlist)-1):
        shellvol = np.pi*2.*((withinr+dr)*(withinr+dr)-withinr*withinr)*(zlist[i+1]-zlist[i])
        shellvol_in_cm3 = shellvol*kpc_in_cm*kpc_in_cm*kpc_in_cm
        cutz = (Gz < zlist[i+1]) & (Gz > zlist[i])
        cut = cutxy*cutz
        if havecr>0: cregy_in_erg_cut = np.sum(cregy_in_erg[cut])
        if havecr>0: creden_in_erg_per_cm3 = cregy_in_erg_cut/shellvol_in_cm3
        if havecr>0: credenlist[i]=creden_in_erg_per_cm3
        if haveB>0: Begy_in_erg_cut = np.sum(Begy[cut])
        if haveB>0: Bden_in_erg_per_cm3 = Begy_in_erg_cut/shellvol_in_cm3
        if haveB>0: Bdenlist[i]=Bden_in_erg_per_cm3
        thden_in_erg_cut = np.sum(thden_in_erg[cut])
        thden_in_erg_per_cm3 = thden_in_erg_cut/shellvol_in_cm3
        thdenlist[i]=thden_in_erg_per_cm3        
    if haveB>0:
            lsn='dashed'
    else:
            lsn='solid'
    plotdict[wanted]['xlab'] = r'${\rm z [kpc]}$'
    plotdict[wanted]['ylab'] = r'$e\;[{\rm erg/cm^3}]$'
    plotdict[wanted]['xnl']['eth'] = zlist[:-1]
    plotdict[wanted]['ynl']['eth'] = thdenlist[:-1]
    plotdict[wanted]['lw']['eth'] = 1    
    plotdict[wanted]['marker']['eth'] = '^'    
    plotdict[wanted]['linelab']['eth'] = 'thermal'
    if haveB>0: plotdict[wanted]['xnl']['eB'] = zlist[:-1]
    if haveB>0: plotdict[wanted]['ynl']['eB'] = Bdenlist[:-1]
    if haveB>0: plotdict[wanted]['lw']['eB'] = 1
    if haveB>0: plotdict[wanted]['marker']['eB'] = 's'
    if haveB>0: plotdict[wanted]['linelab']['eB'] = 'Bfield'
    if havecr>0: plotdict[wanted]['xnl']['ecr'] = zlist[:-1]
    if havecr>0: plotdict[wanted]['ynl']['ecr'] = credenlist[:-1]
    if havecr>0: plotdict[wanted]['lw']['ecr'] = 1   
    if havecr>0: plotdict[wanted]['linelab']['ecr'] = 'CR' 
    if havecr>0: plotdict[wanted]['marker']['ecr'] = 'o'
    plotdict[wanted]['runtodo'] = runtodo
    plotdict[wanted]['labelneed'] = labelneed
    plotdict[wanted]['lsn'] = lsn
    plotdict[wanted]['color'] = color
    plotdict[wanted]['runtitle'] = runtitle
    plotdict[wanted]['ptitle'] = ptitle
    filename=homedir+'CRplot/crdenz/CRz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    plotdict[wanted]['filename'] = filename
    return plotdict