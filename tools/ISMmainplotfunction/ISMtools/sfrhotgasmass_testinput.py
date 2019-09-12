from stdmodandoption import *
import plot_setup as PS
import collections



def sfrhotgasmass_testinput(subdict):
    nested_dict = lambda: collections.defaultdict(nested_dict)
    plotdict = nested_dict()
    startno=subdict['startno']
    Nsnap=subdict['Nsnap']
    snapsep=subdict['snapsep']
    wanted=subdict['wanted']
    dirneed=subdict['dirneed']
    fmeat=subdict['fmeat']
    withinr=subdict['withinr']
    withoutr=subdict['withoutr']
    maxlength=subdict['maxlength']
    minsfr=subdict['minsfr']
    maxsfr=subdict['maxsfr']
    title=''
    titleneed=title
    ptitle=title
    needlog=0
    dclabelneed=1
    correctIa=0
    useM1=0
    resoneed=0
    rotface=1
    loccen=1
    newlabelneed=1
    
    findradiusatnism = 1 #find the radius that has that ISM density (negative to turn off)
    #withinr=5. #withinr will change according to the above
    dr=2.
    fig, ax = PS.setupfig(nrows=1, ncols=1,sharex=True,sharey=False)    
    if wanted=='sfrhotgasmass':
        ncount=0
        ncrcount=0
        nbcount=0
        energylabel=1
        usesolarcircle=0
        usecentral=0
        if usesolarcircle==1:
            fmeat+='_usesolarcircle'
        elif usecentral==1:
            fmeat+='_usecentral'
        sfrl=[]
        gmhiml=[]
        coloril=[]
        runtodol=[]
        for runtodo in dirneed:
            for i in range(startno,Nsnap,snapsep):
                info=SSF.outdirname(runtodo, i)
                rundir=info['rundir']
                runtitle=info['runtitle']
                slabel=info['slabel']
                snlabel=info['snlabel']
                dclabel=info['dclabel']
                resolabel=info['resolabel']
                the_snapdir=info['the_snapdir']
                the_prefix = info['the_prefix']
                the_suffix = info['the_suffix']
                Nsnapstring=info['Nsnapstring']
                havecr=info['havecr']
                Fcal=info['Fcal']
                iavesfr=info['iavesfr']
                timestep=info['timestep']
                maindir=info['maindir']
                haveB=info['haveB']
                cosmo=info['cosmo']
                newlabel=info['newlabel']
                color = info['color']
                usepep = info['usepep']
                snumadd = info['snumadd']
                halostr = info['halostr']
                h0 = info['h0']
                highres = info['highres']
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
                try:
                    usesnipshot=1
                    if usesnipshot==1:
                        print 'usesnipshot'
                        commonpath='/home/tkc004/scratch/snipshot/philruns/'
                        fname=commonpath+rundir+'/output/withinr20GS/snipshot_'+Nsnapstring+'.hdf5'
                        data = RSS.readsnipshot(fname,ptypelist = [0,4])
                        G = data[0];
                        S = data[1];
                        oheader = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo,header_only=0,oheader_only=1)
                        S['header']=oheader
                        #print 'oheader', oheader
                    else:
                        G = SSF.readsnapwcen(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix,\
                         havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
                         loccen=loccen,runtodo=runtodo,rundir=rundir,halostr=halostr)
                        S = SSF.readsnapwcen(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix,\
                         havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
                         loccen=loccen,runtodo=runtodo,rundir=rundir,halostr=halostr)
                    Gmdata = SSF.calGmatT(G,withinr=withinr,nbin=2)
                    gmhim = np.sum(Gmdata['mgrhot'])
                    gmhiml = np.append(gmhiml,gmhim)
                    sfrdata = SSF.calsfr(S,cosmo=cosmo)
                    del G, S
                    sfr = sfrdata['sfr']
                    sfrl=np.append(sfrl,sfr)
                    colori = highres
                    coloril = np.append(coloril, colori)
                    runtodol = np.append(runtodol, runtodo)
                except (ZeroDivisionError,IndexError):
                    continue
                
        
        if wanted=='sfrhotgasmass':
            neederrbar=1
            if neederrbar==1:
                for runtodo in dirneed:
                    Sds = sfrl[runtodol==runtodo]
                    sHs = gmhiml[runtodol==runtodo]
                    cls = coloril[runtodol==runtodo]
                    Sdmed = np.median(Sds)
                    sHmed = np.median(sHs)
                    clmed = np.median(cls)
                    Sdu = np.percentile(Sds,84)-Sdmed
                    Sdd = Sdmed-np.percentile(Sds,15.8)
                    sHu = np.percentile(sHs,84)-sHmed
                    sHd = sHmed-np.percentile(sHs,15.8)
                    if clmed==0:
                        color='k'; marker='s';
                    if clmed==1:
                        color='r'; marker='^';
                    if clmed==3:
                        color='g'; marker='o';
                    ax.errorbar(Sdmed,sHmed, xerr=[[Sdd],[Sdu]], yerr=[[sHd],[sHu]],color=color,fmt=marker,markersize=7)
            else:
                ax.plot(sfrl[coloril==1],gmhiml[coloril==1],color='r',ls='None',marker='^')
                ax.plot(sfrl[coloril==3],gmhiml[coloril==3],color='g',ls='None',marker='s')
        xlab = r'${\rm SFR\;[M_\odot/yr]}$'
        if wanted=='sfrhotgasmass':
            ylab = r'$M(T>10^5\;{\rm K})\;[{\rm M_\odot}]$'
        ax.set_title(ptitle, fontsize=16)
        PS.miscsetup(ax,logx=0,logy=0,xlab=xlab,ylab=ylab,legendneed=0,labfs=22,legfs=12)
        import matplotlib.lines as mlines
        potr = mlines.Line2D([], [], color='r', ls='None',marker='^')
        potg = mlines.Line2D([], [], color='g', ls='None',marker='s')
        lxlist=[potr,potg]
        lablist=['MHD+', 'CR+']
        legend1 = plt.legend(lxlist, lablist, loc='upper left',fontsize=14,ncol=1)
        plt.gca().add_artist(legend1)
        if wanted=='sfrhotgasmass':
            filename=plotloc+'CRplot/sfrhotgasmass/sfrhotgasmass_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        PS.finishsave(plt,filename)
                
                
                
