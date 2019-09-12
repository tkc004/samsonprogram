from stdmodandoption import *
import cameron_functions as CF
import collections

def star_midplane_out(ssdict):
    printsfrclustering =1
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
    maxlength=6.0 #thickness
    withinr=10.0
    withoutr = 0.01
    title='MW'
    titleneed=title
    dclabelneed=1
    useM1=1
    nogrid=10
    rotface=1
    loccen=1
    newlabelneed=1
    print 'runtodo', runtodo
    if not wanted=='SFRsurden':
        print 'not correct key'
        return None
    numoftimes=0
    dr = withinr/nogrid
    radlist=np.linspace(withoutr,withinr,num=nogrid)
    SFRdenlist=0.0*radlist
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
        agecount=10 #Myr
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
        haveB=info['haveB']
        Fcal=info['Fcal']
        iavesfr=info['iavesfr']
        timestep=info['timestep']
        color=info['color']
        ptitle=title
        if runtitle=='SMC':
                ptitle='Dwarf'
        elif runtitle=='SBC':
                ptitle='Starburst'
        elif runtitle=='MW':
                ptitle=r'$L\star$ Galaxy'
        S = SSF.readsnapwcen(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix,\
         havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
         loccen=loccen,runtodo=runtodo,rundir=rundir,halostr=halostr)
        Sage = S['age']
        Spos = S['p']
        Sx = Spos[:,0]
        Sy = Spos[:,1]
        Sz = Spos[:,2]
        Sm = S['m']
        maxage = np.amax(Sage)        
        if cosmo==1:
            readtimelist=SSF.readtime(firever=2)
            snap2list=readtimelist['snaplist']
            time2list=readtimelist['timelist']
            a2list=readtimelist['alist']
            maxage_in_Gyr = np.interp(maxage,a2list,time2list)
            oldtime = np.interp(maxage_in_Gyr-agecount/Gyr_in_Myr,time2list,a2list)
        else:
            oldtime = maxage-agecount/Gyr_in_Myr
        agecutlow = oldtime
        agecuthigh = maxage
        Smdenlist = []
        Smlist = []
        cuta = (Sage> agecutlow) & (Sage <agecuthigh)
        Sx=Sx[cuta]; Sy=Sy[cuta]; Sz=Sz[cuta]; Sm=Sm[cuta];
        SFRp = Sm*m_codetoMsun/agecount/Myr_in_yr #in Msun/yr
        print 'SFRp', SFRp
        for irad in range(len(radlist)-1):
            cutxy = ((Sx)*(Sx)+(Sy)*(Sy) > radlist[irad]*radlist[irad]) & ((Sx)*(Sx)+(Sy)*(Sy) < radlist[irad+1]*radlist[irad+1])
            cutz = (Sz)*(Sz) < maxlength*maxlength/4.
            cut = cutxy*cutz
            cylarea_in_pc2 = np.pi*(-np.power(radlist[irad],2)+np.power(radlist[irad+1],2))
            SFRsurden_in_Msun_yr_kpc2 = np.sum(SFRp[cut])/cylarea_in_pc2
            print 'SFRsurden_in_Msun_yr_kpc2', SFRsurden_in_Msun_yr_kpc2
            SFRdenlist[irad] += SFRsurden_in_Msun_yr_kpc2
        if printsfrclustering==1:
            ncutxy = (Sx)*(Sx)+(Sy)*(Sy) < radlist[-1]*radlist[-1]
            ncutz = (Sz)*(Sz) < maxlength*maxlength/4.
            ncut = ncutxy*ncutz
            SFRpncut = SFRp[ncut]
            Sxncut = Sx[ncut]; Syncut = Sy[ncut]; 
            xlist = np.linspace(-radlist[-1], radlist[-1], num=80)
            ylist = xlist
            clustering = SSF.cal2dclustering(Sxncut,Syncut,SFRpncut,xlist,ylist)
            print 'clustering', clustering
        numoftimes+=1
    if haveB>0:
            lsn='dashed'
    else:
            lsn='solid'
    rmlist = (radlist[:-1]+radlist[1:])/2.
    plotdict[wanted]['xlab'] = r'${\rm r\;[kpc]}$'
    plotdict[wanted]['xnl'] = rmlist          
    plotdict[wanted]['ylab'] = r'${\rm \dot{\Sigma}_\star}\;[{\rm M_{\odot}/yr/kpc^2}]$'
    plotdict[wanted]['ynl'] = SFRdenlist[:-1]/numoftimes
    plotdict[wanted]['runtodo'] = runtodo
    plotdict[wanted]['labelneed'] = labelneed
    plotdict[wanted]['lsn'] = '-'
    plotdict[wanted]['lw'] = 2
    plotdict[wanted]['marker'] = 'None'
    plotdict[wanted]['color'] = color
    plotdict[wanted]['runtitle'] = runtitle
    plotdict[wanted]['ptitle'] = ptitle
    filename=plotloc+'CRplot/SFRsurden/SFRsurden_midplane_'+fmeat+'.pdf'
    plotdict[wanted]['filename'] = filename        
    return plotdict
