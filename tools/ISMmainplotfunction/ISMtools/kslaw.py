from stdmodandoption import *
import cameron_functions as CF
import collections

def kslaw(ssdict):
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
    maxlength=ssdict['maxlength']
    withinr=ssdict['withinr']
    withoutr=ssdict['withoutr']
    nogrid=ssdict['nogrid']
    title='MW'
    titleneed=title
    dclabelneed=1
    useM1=1
    nogrid=10
    rotface=1
    loccen=1
    newlabelneed=1
    print 'runtodo', runtodo
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
        highres = info['highres']
        ptitle=title
        if runtitle=='SMC':
                ptitle='Dwarf'
        elif runtitle=='SBC':
                ptitle='Starburst'
        elif runtitle=='MW':
                ptitle=r'$L\star$ Galaxy'
        G=SSF.readsnapfromrun(runtodo,Nsnap,0,rotface=1,loccen=1)
        cenin = G['cen']; vcenin = G['vcen']; angLin = G['angL'];
        S = SSF.readsnapfromrun(runtodo,Nsnap,4,rotface=1,loccen=1,importLcen=1,angLin=angLin,cenin=cenin,vcenin=vcenin)
        xlist=ylist=np.linspace(withoutr,withinr,num=nogrid)
        gasdenlist = SSF.calsurdenxy(G,xlist,ylist,maxlength)      
        Sdata = SSF.calsfr(S,tintval=0.01,cosmo=1,withinr=20)
        sfrl = Sdata['sfrl']; Sxl = Sdata['Sxl']; Syl = Sdata['Syl']; Szl = Sdata['Szl']
        SFRdenlist = SSF.calSFRsurdenxy(sfrl,Sxl,Syl,Szl,xlist,ylist,maxlength)
    plotdict[wanted]['xlab'] = r'${\rm \Sigma_{gas}\;[M_{\odot}/pc^2]}$'
    plotdict[wanted]['xnl'] =  np.ravel(gasdenlist)        
    plotdict[wanted]['ylab'] = r'${\rm \dot{\Sigma}_\star}\;[{\rm M_{\odot}/yr/kpc^2}]$'
    plotdict[wanted]['ynl'] = np.ravel(SFRdenlist)
    plotdict[wanted]['runtodo'] = runtodo
    plotdict[wanted]['labelneed'] = labelneed
    plotdict[wanted]['lsn'] = 'None'
    plotdict[wanted]['lw'] = 2
    plotdict[wanted]['marker'] = 's'
    plotdict[wanted]['color'] = color
    plotdict[wanted]['runtitle'] = runtitle
    plotdict[wanted]['ptitle'] = ptitle
    filename=plotloc+'CRplot/SFRsurden/kslaw_'+fmeat+'.pdf'
    plotdict[wanted]['filename'] = filename        
    return plotdict
