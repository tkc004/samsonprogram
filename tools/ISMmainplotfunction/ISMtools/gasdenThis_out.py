from stdmodandoption import *
import collections

def gasdenThis_out(ssdict):    
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
    outz=ssdict['outz'] #use z>zdown and z<zup 
    zdown=ssdict['zdown']
    zup=ssdict['zup']
    titleneed=title
    nogrid=10
    xaxis_snapno=0
    if not (wanted=='rhoTvol' or wanted=='rhoTmass'):
        print 'wrong wanted'
        return None
    resoneed=0
    rotface=1
    loccen=1
    newlabelneed=1
    print 'runtodo', runtodo
    nogrid = 15
    dr = withinr/nogrid
    lognlist = np.linspace(-4.0,4.0,num=nogrid)      
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
        Grho = Gextra['rho']; Gu = Gextra['u']; Gm = Gextra['m']; Neb = Gextra['ne']
        TrueTemp, converted_rho = SF.convertTemp(Gu, Neb, Grho)
        cutxy = (Grxy < withinr) & (Grxy > withoutr)
        Gmcutxy = Gm[cutxy]*m_codetocgs #in g
        Grhocutxy = Grho[cutxy]*rho_codetocgs #in g/cm^3
        Gzcutxy = Gz[cutxy]
        Tcutxy = TrueTemp[cutxy]
        Gncutxy = converted_rho[cutxy]
        if outz==0:
            cutz = np.absolute(Gzcutxy) < maxlength/2.0
        elif outz==1:
            cutz = (Gzcutxy > zdown) & (Gzcutxy < zup)
        Gmcutz = Gmcutxy[cutz]
        Grhocutz = Grhocutxy[cutz]
        Tcutz = Tcutxy[cutz]
        Gncutz = Gncutxy[cutz]
        hotcut = Tcutz>1e5 #K
        warmcut = (Tcutz>1e3)&(Tcutz<1e5)
        coldcut = Tcutz<1e3
        if wanted=='rhoTvol':
            Gvolcutz = Gmcutz/Grhocutz
            nweight = Gvolcutz/np.sum(Gvolcutz)
        if wanted=='rhoTmass':
            nweight = Gmcutz/np.sum(Gmcutz)
        histhot, bin_edges= np.histogram(np.log10(Gncutz[hotcut]), bins=lognlist,weights=nweight[hotcut])
        histwarm, bin_edges= np.histogram(np.log10(Gncutz[warmcut]), bins=lognlist,weights=nweight[warmcut])
        histcold, bin_edges= np.histogram(np.log10(Gncutz[coldcut]), bins=lognlist,weights=nweight[coldcut])
    lognmlist = (lognlist[:-1]+lognlist[1:])/2.0
    if wanted=='rhoTvol':
        plotdict[wanted]['ylab'] = r'${\rm d}V_/{\rm d}\log_{10}n_{\rm ISM}$'            
    if wanted=='rhoTmass':
        plotdict[wanted]['ylab'] = r'${\rm d}m/{\rm d}\log_{10}n_{\rm ISM}$'
    plotdict[wanted]['xlab'] = r'${\log_{10}n_{\rm ISM}\;[{\rm cm^{-3}}]}$'
    plotdict[wanted]['ynl']['cold'] = histcold
    plotdict[wanted]['xnl']['cold'] = lognmlist
    plotdict[wanted]['lsn']['cold'] = 'dashed'
    plotdict[wanted]['linelab']['cold'] = r'${\rm T < 1e3K}$'
    plotdict[wanted]['ynl']['warm'] = histwarm
    plotdict[wanted]['xnl']['warm'] = lognmlist
    plotdict[wanted]['lsn']['warm'] = 'dashdot'
    plotdict[wanted]['linelab']['warm'] = r'${\rm 1e3K < T < 1e5K}$'
    plotdict[wanted]['ynl']['hot'] = histhot
    plotdict[wanted]['xnl']['hot'] = lognmlist
    plotdict[wanted]['lsn']['hot'] = 'dotted'
    plotdict[wanted]['linelab']['hot'] = r'${\rm T > 1e5K}$'
    plotdict[wanted]['runtodo'] = runtodo
    plotdict[wanted]['labelneed'] = labelneed
    plotdict[wanted]['lw'] = 2
    plotdict[wanted]['marker'] = 'None'
    plotdict[wanted]['color'] = color
    plotdict[wanted]['runtitle'] = runtitle
    plotdict[wanted]['ptitle'] = ptitle 
    if wanted=='rhoTvol':
        filename=plotloc+'CRplot/rhoTvol/rhoTvol_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if outz==1:
            filename=plotloc+'CRplot/rhoTvol/rhoTvoloutz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    if wanted=='rhoTmass':
        filename=plotloc+'CRplot/rhoTmass/rhoTmass_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if outz==1:
            filename=plotloc+'CRplot/rhoTmass/rhoTmassoutz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'            
    plotdict[wanted]['filename'] = filename      
        
    return plotdict
