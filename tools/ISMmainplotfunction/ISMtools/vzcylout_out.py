from stdmodandoption import *
import collections

def vzcylout_out(ssdict):    
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
    nogrid=40
    xaxis_snapno=0
    if not (wanted=='absvzout' or wanted=='gmz'):
        print 'wrong wanted'
        return None
    resoneed=0
    rotface=1
    loccen=1
    newlabelneed=1
    print 'runtodo', runtodo
    if wanted=='absvzout':
        nogrid = 15
        dr = withinr/nogrid
        zlist=np.linspace(1.0,maxlength,num=nogrid)
        absvzlist=zlist*0.
    if wanted=='gmz':
        nogrid = 15
        dr = withinr/nogrid
        zlist=np.linspace(1.0,maxlength,num=nogrid)
        gmzlist=zlist*0. 
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
        if wanted=='absvzout':
            for iz in range(len(zlist)-1):
                cutxy = (Grxy > withoutr) & (Grxy < withinr)
                cutz = (Gz > zlist[iz]) & (Gz < zlist[iz+1])
                cut = cutxy*cutz
                Gvzc_km_s = Gvz[cut] #km/s
                Gmc_in_g = Gm[cut]*1e10*Msun_in_g
                #print 'Gmc_in_g', Gmc_in_g
                #print 'Gvzc_km_s', Gvzc_km_s
                absvz = np.sqrt(np.sum(Gmc_in_g*Gvzc_km_s*Gvzc_km_s)/np.sum(Gmc_in_g)) 
                absvzlist[iz] += absvz       
        if wanted=='gmz':
            for iz in range(len(zlist)-1):
                cutxy = (Grxy > withoutr) & (Grxy < withinr)
                cutz = (Gz < zlist[iz+1])
                cut = cutxy*cutz
                Gmc_in_Msun = Gm[cut]*1e10
                gmzlist[iz] += np.sum(Gmc_in_Msun)
        numoftimes+=1
        zmlist = (zlist[:-1]+zlist[1:])/2.
        plotdict[wanted]['xlab'] = r'${\rm z\;[kpc]}$'
        if wanted=='absvzout':
            plotdict[wanted]['ylab'] = r'$<|v_{\rm z}|\;[{\rm km/s}]>$'
            plotdict[wanted]['ynl'] = absvzlist[:-1]/numoftimes
            plotdict[wanted]['xnl'] = zmlist
        if wanted=='gmz':
            plotdict[wanted]['ylab'] = r'$M_{\rm gas}(<z)\;[{\rm M_\odot}]$'
            plotdict[wanted]['ynl'] = gmzlist[:-1]/numoftimes
            plotdict[wanted]['xnl'] = zlist[1:]
        plotdict[wanted]['runtodo'] = runtodo
        plotdict[wanted]['labelneed'] = labelneed
        plotdict[wanted]['lw'] = 2
        plotdict[wanted]['lsn'] = 'solid'
        plotdict[wanted]['marker'] = 'None'
        plotdict[wanted]['color'] = color
        plotdict[wanted]['runtitle'] = runtitle
        plotdict[wanted]['ptitle'] = ptitle
        if wanted=='absvzout':
            filename=plotloc+'CRplot/absvzout/absvzout_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            plotdict[wanted]['filename'] = filename
        if wanted=='gmz':
            filename=plotloc+'CRplot/gmz/gmz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            plotdict[wanted]['filename'] = filename        
        return plotdict
