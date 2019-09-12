from stdmodandoption import *
import plot_setup as PS
import collections



def sfrtestinput(subdict):
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
    if wanted=='sfr':
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
        runtodol=[]
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
                    S = SSF.readsnapwcen(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix,\
                     havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
                     loccen=loccen,runtodo=runtodo,rundir=rundir,halostr=halostr)
                    sfrdata = SSF.calsfr(S,cosmo=cosmo)
                    sfr = sfrdata['sfr']
                    runtodol=np.append(runtodol,runtodo)
                    sfrl=np.append(sfrl,sfr)
                    colori = highres
                    coloril = np.append(coloril, colori)
                except (ZeroDivisionError,IndexError):
                    continue
        print 'runtodol', runtodol
        print 'sfrl', sfrl
                
                
