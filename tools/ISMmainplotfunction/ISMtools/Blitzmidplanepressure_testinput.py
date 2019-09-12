from stdmodandoption import *
import plot_setup as PS
import collections


def Blitzmidplanepressure_testinput(subdict):
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
    maxlength=subdict['maxlength']
    usesnipshot=subdict['usesnipshot']
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
      

    the_prefix = 'snapshot'
    the_suffix = '.hdf5'
    sigmaHIlosgridl=np.array([])
    SFRdengridl=np.array([])
    coloril=np.array([])
    runtodol=np.array([])
    for runtodo in dirneed:
        fig, ax = PS.setupfig(nrows=1, ncols=1,sharex=True,sharey=False)  
        i=Nsnap
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
        h0=info['h0']
        cosmo=info['cosmo']
        usepep=info['usepep']
        maindir=info['maindir']
        snumadd=info['snumadd']
        halostr=info['halostr']
        highres=info['highres']
        rotface = 1
        loccen = 1
        if usesnipshot==1:
            commonpath='/home/tkc004/scratch/snipshot/philruns/'
            fname=commonpath+rundir+'/output/withinr20GS/snipshot_'+Nsnapstring+'.hdf5'
            data = RSS.readsnipshot(fname,ptypelist = [0,4])
            G = data[0];
            S = data[1];
            oheader = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix,\
                                 extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo,\
                                 header_only=0,oheader_only=1)
            S['header']=oheader
            #print 'oheader', oheader
        else:
            S = SSF.readsnapwcen(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix,\
             havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
             loccen=loccen,runtodo=runtodo,rundir=rundir,halostr=halostr)
            G = SSF.readsnapwcen(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix,\
             havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
             loccen=loccen,runtodo=runtodo,rundir=rundir,halostr=halostr)           
        rlist = np.linspace(0,10,num=10)
        gasdenlist = SSF.calsurden(G,rlist,maxlength)
        stardenlist = SSF.calsurden(S,rlist,maxlength)
        propconst = Msun_in_g/pc_in_cm/pc_in_cm
        gasdenlist = gasdenlist*propconst
        stardenlist = stardenlist*propconst
        Sscaledata = SSF.findscalelength(S, halfmassscale=0, usecylr=0, usecylz=1, cutrout =10., cutzout = 5.)
        hstar = Sscaledata['rs']
        hstar =hstar*kpc_in_cm
        Gscaledata = SSF.findscalelength(G, halfmassscale=0, usecylr=0, usecylz=1, cutrout =10., cutzout = 5.)
        hgas = Gscaledata['rs']
        hgas = hgas*kpc_in_cm
        vertical=0; horizontal=1;
        dr = 0.5 #kpc
        griddir = 'grid0_5kpc'
        usehalfz=0
        data = SSF.calrhogfrompreexist(runtodo,i,withinr,maxlength,\
                                       usehalfz=usehalfz,griddir=griddir,\
                                       cutcold=0,vertical=vertical,\
                                       horizontal=horizontal,withoutr=withoutr,\
                                      dr=dr,outHI=1)
        pthl = data['pthl']; rhol = data['rhol'];
        pturHIl = data['pturHIl']; rhoHIl = data['rhoHIl'];
        volll = data['volll']
        massl = volll/rhol
        vturHI2l = pturHIl/rhoHIl*massl
        vth2l = pthl/rhol*GAMMA*massl
        massl[~np.isfinite(massl)]=0.0
        vturHI2l[~np.isfinite(vturHI2l)]=0.0
        vth2l[~np.isfinite(vth2l)]=0.0
        vHIdis = np.sqrt(np.sum(vturHI2l+vth2l)/np.sum(massl))
        #finally Blitz & Rosolowsky 2004 formula:
        pext = np.sqrt(2.0*NewtonG_in_cgs*stardenlist/hstar)*gasdenlist*vHIdis
        #with gas gravity:
        pextg = np.sqrt(2.0*NewtonG_in_cgs)*gasdenlist*vHIdis*\
        (np.sqrt(stardenlist/hstar)+np.sqrt(np.pi/4.0*gasdenlist/hgas))
        ax.plot(rlist, pext)
        ax.plot(rlist, pextg,ls='dashed')
        xlab = r'$ r\;[{\rm kpc]}$'
        ylab = r'$  P_{\rm ext}[{\rm dyne/cm^2]}$'             
        ax.set_title(ptitle, fontsize=16)
        PS.miscsetup(ax,logx=0,logy=1,xlab=xlab,ylab=ylab,legendneed=0,labfs=22,legfs=12)
        filename=plotloc+'CRplot/Pext/Pext_Blitz_'+runtodo+'_sn'+str(Nsnap)+'.pdf'
        PS.finishsave(plt,filename)
