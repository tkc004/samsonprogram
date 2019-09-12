from stdmodandoption import *
import collections

def hgashalf_out(ssdict):    
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
    titleneed=title
    nogrid=10
    xaxis_snapno=0
    if not (wanted=='hgas' or wanted=='mgasz' or wanted=='hcr' or wanted=='hth'):
        print 'wrong wanted'
        return None
    resoneed=0
    rotface=1
    loccen=1
    newlabelneed=1
    print 'runtodo', runtodo
    if wanted=='hgas':
        nogrid = 15
        dr = withinr/nogrid
        rlist = np.linspace(withoutr,withinr,num=nogrid)
        zlist = np.linspace(0.01,maxlength,num=50)
        hglist=rlist*0. 
    if wanted=='mgasz':
        nogrid = 15
        dr = withinr/nogrid
        rlist = np.linspace(withoutr,withinr,num=nogrid)
        mgzlist=rlist*0.
    if wanted=='hcr':
        nogrid = 15
        dr = withinr/nogrid
        rlist = np.linspace(withoutr,withinr,num=nogrid)
        zlist = np.linspace(0.01,maxlength,num=50)
        hcrlist=rlist*0.
    if wanted=='hth':
        nogrid = 15
        dr = withinr/nogrid
        rlist = np.linspace(withoutr,withinr,num=nogrid)
        zlist = np.linspace(0.01,maxlength,num=50)
        hthlist=rlist*0.
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
        #if runtitle=='SMC':
        #    ptitle='Dwarf'
        #elif runtitle=='SBC':
        #    ptitle='Starburst'
        #elif runtitle=='MW':
        #    ptitle=r'$L\star$ Galaxy'
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
        Grho = Gextra['rho']; Gu = Gextra['u']; Gm = Gextra['m'];
        if wanted=='hgas':
            for ir in range(len(rlist)-1):
                cutxy = (Grxy > rlist[ir]) & (Grxy < rlist[ir+1])
                Gmcutxy = Gm[cutxy]*1e10
                Gzcutxy = Gz[cutxy]
                gmzlist = []
                for iz in range(len(zlist)-1):
                    cutz = (np.absolute(Gzcutxy) < zlist[iz+1])
                    Gmcutz = Gmcutxy[cutz]
                    gmzlist = np.append(gmzlist,np.sum(Gmcutz))
                print 'zlist', zlist
                print 'gmzlist', gmzlist
                hgas = np.interp(gmzlist[-1]/2.0,gmzlist,zlist[1:])
                hglist[ir] += hgas
        if wanted=='mgasz':
            for ir in range(len(rlist)-1):
                cutxy = (Grxy > rlist[ir]) & (Grxy < rlist[ir+1])
                Gmcutxy = Gm[cutxy]*1e10 #in Msun
                Gzcutxy = Gz[cutxy]
                cutz = (np.absolute(Gzcutxy) < maxlength/2.0) & (np.absolute(Gzcutxy) > 1.0)
                Gmcutz = Gmcutxy[cutz]
                mgzlist[ir] += np.sum(Gmcutz)
        if wanted=='hcr':
            cregy = Gextra['cregy'] 
            for ir in range(len(rlist)-1):
                cutxy = (Grxy > rlist[ir]) & (Grxy < rlist[ir+1])
                cregycutxy = cregy[cutxy]
                Gzcutxy = Gz[cutxy]
                cregyzlist = []
                for iz in range(len(zlist)-1):
                    cutz = (np.absolute(Gzcutxy) < zlist[iz+1])
                    cregycutz = cregycutxy[cutz]
                    cregyzlist = np.append(cregyzlist,np.sum(cregycutz))
                hcr = np.interp(cregyzlist[-1]/2.0,cregyzlist,zlist[1:])
                hcrlist[ir] += hcr
        if wanted=='hth':
            for ir in range(len(rlist)-1):
                cutxy = (Grxy > rlist[ir]) & (Grxy < rlist[ir+1])
                Ethcutxy = Gu[cutxy]*Gm[cutxy]
                Gzcutxy = Gz[cutxy]
                thzlist = []
                for iz in range(len(zlist)-1):
                    cutz = (np.absolute(Gzcutxy) < zlist[iz+1])
                    Ethcutz = Ethcutxy[cutz]
                    thzlist = np.append(thzlist,np.sum(Ethcutz))
                hth = np.interp(thzlist[-1]/2.0,thzlist,zlist[1:])
                hthlist[ir] += hth                
        numoftimes+=1
        rmlist = (rlist[:-1]+rlist[1:])/2.
        plotdict[wanted]['xlab'] = r'${\rm r\;[kpc]}$'
        if wanted=='hgas':
            plotdict[wanted]['ylab'] = r'$h_{\rm gas,1/2}\;[{\rm kpc}]$'
            plotdict[wanted]['ynl'] = hglist[:-1]/numoftimes
            plotdict[wanted]['xnl'] = rmlist
        if wanted=='mgasz':
            plotdict[wanted]['ylab'] = r'$M_{\rm gas}\;[{\rm M_\odot}]$'
            plotdict[wanted]['ynl'] = mgzlist[:-1]/numoftimes
            plotdict[wanted]['xnl'] = rmlist
        if wanted=='hcr':
            plotdict[wanted]['ylab'] = r'$h_{\rm cr,1/2}\;[{\rm kpc}]$'
            plotdict[wanted]['ynl'] = hcrlist[:-1]/numoftimes
            plotdict[wanted]['xnl'] = rmlist
        if wanted=='hth':
            plotdict[wanted]['ylab'] = r'$h_{\rm th,1/2}\;[{\rm kpc}]$'
            plotdict[wanted]['ynl'] = hthlist[:-1]/numoftimes
            plotdict[wanted]['xnl'] = rmlist
        plotdict[wanted]['runtodo'] = runtodo
        plotdict[wanted]['labelneed'] = labelneed
        plotdict[wanted]['lw'] = 2
        plotdict[wanted]['lsn'] = 'solid'
        plotdict[wanted]['marker'] = 'None'
        plotdict[wanted]['color'] = color
        plotdict[wanted]['runtitle'] = runtitle
        plotdict[wanted]['ptitle'] = ptitle
        if wanted=='hgas':
            filename=plotloc+'CRplot/hgas/hgas_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            plotdict[wanted]['filename'] = filename  
        if wanted=='mgasz':
            filename=plotloc+'CRplot/mgasz/mgasz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            plotdict[wanted]['filename'] = filename 
        if wanted=='hcr':
            filename=plotloc+'CRplot/hcr/hcr_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            plotdict[wanted]['filename'] = filename 
        if wanted=='hth':
            filename=plotloc+'CRplot/hth/hth_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            plotdict[wanted]['filename'] = filename             
        return plotdict
