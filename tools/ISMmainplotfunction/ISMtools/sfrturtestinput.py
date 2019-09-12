from stdmodandoption import *
import plot_setup as PS
import collections



def sfrturtestinput(subdict):
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
    withinr=5. #withinr will change according to the above
    dr=2.
    fig, ax = PS.setupfig(nrows=1, ncols=1,sharex=True,sharey=False)    
    if wanted=='sfrtur' or wanted=='sfrtur_m' or wanted=='sigmatur_mweighted'\
    or wanted=='vth_mweighted' or wanted=='vcr_mweighted'\
    or wanted=='vA_mweighted' or wanted=='v_mweighted':
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
        turt=[]
        thermt=[]
        crt=[]
        Begyt=[]
        sfrl=[]
        gml=[]
        coloril=[]
        for runtodo in dirneed:
            info=SSF.outdirname(runtodo, Nsnap)
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
            for i in range(startno,Nsnap,snapsep):
                try:
                    info=SSF.outdirname(runtodo, i)
                    Nsnapstring=info['Nsnapstring']
                    egydata=CRTF.energyoutput(runtodo,i,usesolarcircle=usesolarcircle,usecentral=usecentral,\
                                              cylr=withinr,cylrin=withoutr,cylz=maxlength)
                    turl=egydata['turl']
                    therml=egydata['therml']
                    cregyl=egydata['cregyl']
                    Begyl=egydata['Begyl']
                    gm = np.sum(egydata['gminl']) #in g
                    turt=np.append(turt,np.sum(turl))
                    thermt=np.append(thermt,np.sum(therml))
                    crt=np.append(crt,np.sum(cregyl))
                    Begyt=np.append(Begyt,np.sum(Begyl))
                    gml = np.append(gml,gm)
                    S = SSF.readsnapwcen(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix,\
                     havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
                     loccen=loccen,runtodo=runtodo,rundir=rundir,halostr=halostr)
                    sfrdata = SSF.calsfr(S,cosmo=cosmo)
                    sfr = sfrdata['sfr']
                    sfrl=np.append(sfrl,sfr)
                    colori = highres
                    coloril = np.append(coloril, colori)
                except (ZeroDivisionError,IndexError):
                    continue

        turt=np.array(turt)
        thermt=np.array(thermt)
        crt=np.array(crt)
        Begyt=np.array(Begyt)
        sigmaz = np.sqrt(turt/gml)/km_in_cm
        vth = np.sqrt(thermt*(GAMMA-1)*(GAMMA)/gml)/km_in_cm
        vcr = np.sqrt(crt*(CRgamma-1)*(CRgamma)/gml)/km_in_cm
        vA = np.sqrt(Begyt/gml*2.0)/km_in_cm
        print 'sfrl', sfrl
        print 'turt', turt
        print 'sigmaz', sigmaz
        print 'vth', vth
        print 'vA', vA
        print 'coloril', coloril
        
        if wanted=='sfrtur':
            ax.plot(sfrl[coloril==0],turt[coloril==0],color='k',ls='None',marker='s')
            ax.plot(sfrl[coloril==1],turt[coloril==1],color='r',ls='None',marker='s')
            ax.plot(sfrl[coloril==3],turt[coloril==3],color='g',ls='None',marker='s')
        if wanted=='sfrtur_m':
            ax.plot(sfrl[coloril==0],turt[coloril==0]/gml[coloril==0],color='k',ls='None',marker='s')
            ax.plot(sfrl[coloril==1],turt[coloril==1]/gml[coloril==1],color='r',ls='None',marker='s')
            ax.plot(sfrl[coloril==3],turt[coloril==3]/gml[coloril==3],color='g',ls='None',marker='s')
        if wanted=='sigmatur_mweighted' or wanted=='v_mweighted':
            ax.plot(sfrl[coloril==0],sigmaz[coloril==0],color='k',ls='None',marker='d')
            ax.plot(sfrl[coloril==1],sigmaz[coloril==1],color='r',ls='None',marker='d')
            ax.plot(sfrl[coloril==3],sigmaz[coloril==3],color='g',ls='None',marker='d')
            ax.set_ylim(ymin=1,ymax=30)
            ax.set_xlim(xmin=minsfr,xmax=maxsfr)
        if wanted=='vth_mweighted' or wanted=='v_mweighted':
            ax.plot(sfrl[coloril==0],vth[coloril==0],color='k',ls='None',marker='^')
            ax.plot(sfrl[coloril==1],vth[coloril==1],color='r',ls='None',marker='^')
            ax.plot(sfrl[coloril==3],vth[coloril==3],color='g',ls='None',marker='^')
            ax.set_ylim(ymin=1,ymax=30)
            ax.set_xlim(xmin=minsfr,xmax=maxsfr)
        if wanted=='vcr_mweighted' or wanted=='v_mweighted':
            ax.plot(sfrl[coloril==0],vcr[coloril==0],color='k',ls='None',marker='o')
            ax.plot(sfrl[coloril==1],vcr[coloril==1],color='r',ls='None',marker='o')
            ax.plot(sfrl[coloril==3],vcr[coloril==3],color='g',ls='None',marker='o')
            ax.set_ylim(ymin=1,ymax=30)
            ax.set_xlim(xmin=minsfr,xmax=maxsfr)
        if wanted=='vA_mweighted' or wanted=='v_mweighted':
            ax.plot(sfrl[coloril==0],vA[coloril==0],color='k',ls='None',marker='s')
            ax.plot(sfrl[coloril==1],vA[coloril==1],color='r',ls='None',marker='s')
            ax.plot(sfrl[coloril==3],vA[coloril==3],color='g',ls='None',marker='s')
            ax.set_ylim(ymin=1,ymax=30)
            ax.set_xlim(xmin=minsfr,xmax=maxsfr)
        xlab = r'${\rm SFR\;[M_\odot/yr]}$'
        if wanted=='sfrtur':
            ylab = r'${\rm E_{tur}\;[erg]}$'
        if wanted=='sfrtur_m':
            ylab = r'${\rm E_{tur}/m\;[erg/g]}$'
        if wanted=='sigmatur_mweighted':
            ylab = r'${\rm  \langle\sigma_{tur}\rangle\;[km/s]}$'            
        if wanted=='vth_mweighted':
            ylab = r'${\rm  \langle c_{th}\rangle\;[km/s]}$'
        if wanted=='vcr_mweighted':
            ylab = r'${\rm  \langle c_{cr}\rangle\;[km/s]}$'
        if wanted=='vA_mweighted':
            ylab = r'${\rm  \langle v_{A}\rangle\;[km/s]}$'            
        if wanted=='v_mweighted':
            ylab = r'${\rm  \langle c\rangle\;[km/s]}$'  
        ax.set_title(ptitle, fontsize=16)
        PS.miscsetup(ax,logx=1,logy=1,xlab=xlab,ylab=ylab,legendneed=0,labfs=22,legfs=12)
        import matplotlib.lines as mlines
        potk = mlines.Line2D([], [], color='k', ls='None',marker='s')
        potr = mlines.Line2D([], [], color='r', ls='None',marker='s')
        potg = mlines.Line2D([], [], color='g', ls='None',marker='s')
        #potk = plt.plot([],[], ls='None',marker='s', color='k')
        #potr = plt.plot([],[], ls='None',marker='s', color='r')
        #potg = plt.plot([],[], ls='None',marker='s', color='g')
        lxlist=[potk,potr,potg]
        lablist=['no CR', 'MHD no CR', r'MHD $\kappa$=3e29']
        legend1 = plt.legend(lxlist, lablist, loc='lower left',fontsize=14,ncol=1)
        plt.gca().add_artist(legend1)
        if wanted=='v_mweighted':
            pottur = mlines.Line2D([], [], color='r', ls='None',marker='d')
            potth = mlines.Line2D([], [], color='r', ls='None',marker='^')
            potcr = mlines.Line2D([], [], color='r', ls='None',marker='o')
            potvA = mlines.Line2D([], [], color='r', ls='None',marker='s')
            lxvlist=[pottur,potth,potcr,potvA]
            labvlist=[r'${\rm \langle\sigma_{tur}\rangle}$', r'${\rm \langle c_{th}\rangle}$', r'${\rm \langle c_{cr}\rangle}$',r'${\rm \langle v_{A}\rangle}$']
            legend1v = plt.legend(lxvlist, labvlist, loc='lower right',fontsize=14,ncol=1)
            plt.gca().add_artist(legend1v)        
        if wanted=='sfrtur':
            filename=plotloc+'CRplot/sfrtur/sfrtur_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted=='sfrtur_m':
            filename=plotloc+'CRplot/sfrtur/sfrtur_m_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted=='sigmatur_mweighted':
            filename=plotloc+'CRplot/sigmatur_mweighted/sigmatur_mweighted_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted=='vth_mweighted':
            filename=plotloc+'CRplot/vth_mweighted/vth_mweighted_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted=='vcr_mweighted':
            filename=plotloc+'CRplot/vcr_mweighted/vcr_mweighted_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted=='vA_mweighted':
            filename=plotloc+'CRplot/vA_mweighted/vA_mweighted_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted=='v_mweighted':
            filename=plotloc+'CRplot/v_mweighted/v_mweighted_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        PS.finishsave(plt,filename)
                
                
                
