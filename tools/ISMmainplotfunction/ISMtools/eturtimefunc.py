from stdmodandoption import *
import plot_setup as PS
import collections



def eturtimefunc(ssdict):
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
    
    if wanted=='etur_mtime' or wanted=='eturtime':
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
        turt=[]
        thermt=[]
        crt=[]
        bt=[]
        timel=[]
        Gmt=[]
        for i in range(startno,Nsnap,snapsep):
            #try:
            egydata=CRTF.energyoutput(runtodo,i,usesolarcircle=usesolarcircle,usecentral=usecentral)
            Begycutz=egydata['Begyl']
            cregycutz=egydata['cregyl']
            GEintcutz=egydata['therml']
            turl=egydata['turl']
            Gxpl=egydata['Gxpl']
            Gypl=egydata['Gypl']
            Gzpl=egydata['Gzpl']
            Gxcl=egydata['xcell']
            Gycl=egydata['ycell']
            Gzcl=egydata['zcell']
            gminl=egydata['gminl']
            timen=egydata['timen']
            cylrin=egydata['cylrin']
            cylr=egydata['cylr']
            cylz=egydata['cylz']
            nbin=egydata['nbin']
            print 'gml', np.sum(egydata['gml'])
            print 'gminl', np.sum(gminl)
            Gmt = np.append(Gmt, np.sum(gminl))
            turt=np.append(turt, np.sum(turl))
            thermt=np.append(thermt,np.sum(GEintcutz))
            if havecr>0:
                crt = np.append(crt,np.sum(cregycutz))
            if haveB>0:
                bt = np.append(bt,np.sum(Begycutz))
            timel = np.append(timel,timen/1e6) #Myr 
            #except (ZeroDivisionError,IndexError):
            #    continue
            turt=np.array(turt)
            thermt=np.array(thermt)
            crt=np.array(crt)
            bt=np.array(bt)
            timel=np.array(timel)
            Gmt=np.array(Gmt)       
        print 'Gmt', np.sum(Gmt)
        ls='solid'
        tlb = 'turb'
        thlb = 'therm'
        clb = 'CR'
        Blb = 'B'

        if wanted =='etur_mtime': 
            plotdict[wanted]['xlab'] = r'$t\;{\rm Myr}$'
            plotdict[wanted]['ylab'] = r'$E/m {\rm [erg/g]}$'
            plotdict[wanted]['labelneed'] = labelneed
            plotdict[wanted]['xnl']['turt'] = timel
            print 'turt', turt
            print 'Gmt', Gmt
            plotdict[wanted]['ynl']['turt'] = turt/Gmt
            plotdict[wanted]['lw']['turt'] = 2
            plotdict[wanted]['lsn']['turt'] = ls
            plotdict[wanted]['marker']['turt'] = 'None'    
            plotdict[wanted]['linelab']['turt'] = tlb
            plotdict[wanted]['color']['turt'] = 'g'
            plotdict[wanted]['xnl']['thermt'] = timel
            plotdict[wanted]['ynl']['thermt'] = thermt/Gmt
            plotdict[wanted]['lw']['thermt'] = 2
            plotdict[wanted]['lsn']['thermt'] = ls
            plotdict[wanted]['marker']['thermt'] = 'None'    
            plotdict[wanted]['linelab']['thermt'] = thlb
            plotdict[wanted]['color']['thermt'] = 'r'            
            if havecr>0:
                plotdict[wanted]['xnl']['crt'] = timel
                plotdict[wanted]['ynl']['crt'] = crt/Gmt
                plotdict[wanted]['lw']['crt'] = 2
                plotdict[wanted]['lsn']['crt'] = ls
                plotdict[wanted]['marker']['crt'] = 'None'    
                plotdict[wanted]['linelab']['crt'] = clb
                plotdict[wanted]['color']['crt'] = 'y'
            if haveB>0:
                plotdict[wanted]['xnl']['bt'] = timel
                plotdict[wanted]['ynl']['bt'] = bt/Gmt
                plotdict[wanted]['lw']['bt'] = 2
                plotdict[wanted]['lsn']['bt'] = ls
                plotdict[wanted]['marker']['bt'] = 'None'    
                plotdict[wanted]['linelab']['bt'] = Blb
                plotdict[wanted]['color']['bt'] = 'b'
            plotdict[wanted]['runtodo'] = runtodo
            plotdict[wanted]['labelneed'] = labelneed
            plotdict[wanted]['runtitle'] = runtitle
            plotdict[wanted]['ptitle'] = labelneed
            figname=plotloc+'CRplot/etur_mtime/etur_mtime_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            plotdict[wanted]['filename'] = figname
        if wanted =='eturtime': 
            plotdict[wanted]['xlab'] = r'$t\;{\rm Myr}$'
            plotdict[wanted]['ylab'] = r'$E {\rm [erg]}$'
            plotdict[wanted]['labelneed'] = labelneed
            plotdict[wanted]['xnl']['turt'] = timel
            plotdict[wanted]['ynl']['turt'] = turt
            plotdict[wanted]['lw']['turt'] = 2
            plotdict[wanted]['lsn']['turt'] = ls
            plotdict[wanted]['marker']['turt'] = 'None'    
            plotdict[wanted]['linelab']['turt'] = tlb
            plotdict[wanted]['color']['turt'] = 'g'
            plotdict[wanted]['xnl']['thermt'] = timel
            plotdict[wanted]['ynl']['thermt'] = thermt
            plotdict[wanted]['lw']['thermt'] = 2
            plotdict[wanted]['lsn']['thermt'] = ls
            plotdict[wanted]['marker']['thermt'] = 'None'    
            plotdict[wanted]['linelab']['thermt'] = thlb
            plotdict[wanted]['color']['thermt'] = 'r'            
            if havecr>0:
                plotdict[wanted]['xnl']['crt'] = timel
                plotdict[wanted]['ynl']['crt'] = crt
                plotdict[wanted]['lw']['crt'] = 2
                plotdict[wanted]['lsn']['crt'] = ls
                plotdict[wanted]['marker']['crt'] = 'None'    
                plotdict[wanted]['linelab']['crt'] = clb
                plotdict[wanted]['color']['crt'] = 'y'
            if haveB>0:
                plotdict[wanted]['xnl']['bt'] = timel
                plotdict[wanted]['ynl']['bt'] = bt
                plotdict[wanted]['lw']['bt'] = 2
                plotdict[wanted]['lsn']['bt'] = ls
                plotdict[wanted]['marker']['bt'] = 'None'    
                plotdict[wanted]['linelab']['bt'] = Blb
                plotdict[wanted]['color']['bt'] = 'b'
            plotdict[wanted]['runtodo'] = runtodo
            plotdict[wanted]['labelneed'] = labelneed
            plotdict[wanted]['runtitle'] = runtitle
            plotdict[wanted]['ptitle'] = labelneed
            figname=plotloc+'CRplot/etur_mtime/eturtime_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            plotdict[wanted]['filename'] = figname
    return plotdict