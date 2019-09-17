from stdmodandoption import *
import plot_setup as PS
import collections
from mockvelocitydispersion import calsigmaHIlos


def sfrmockturtestinput(subdict):
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
    
    fig, ax = PS.setupfig(nrows=1, ncols=1,sharex=True,sharey=False)    

    the_prefix = 'snapshot'
    the_suffix = '.hdf5'
    sigmaHIlosgridl=np.array([])
    SFRdengridl=np.array([])
    coloril=np.array([])
    runtodol=np.array([])
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
                oheader = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo,header_only=0,oheader_only=1)
                S['header']=oheader
                #print 'oheader', oheader
            else:
                S = SSF.readsnapwcen(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix,\
                 havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
                 loccen=loccen,runtodo=runtodo,rundir=rundir,halostr=halostr)
                G = SSF.readsnapwcen(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix,\
                 havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
                 loccen=loccen,runtodo=runtodo,rundir=rundir,halostr=halostr)
            Sdata = SSF.calsfr(S,tintval=0.01,cosmo=1,withinr=20)
            sfrl = Sdata['sfrl']; Sxl = Sdata['Sxl']; Syl = Sdata['Syl']; Szl = Sdata['Szl']
            
            if wanted=='SFRturrad':
                xlist = np.linspace(-10,10,num=40)
                ylist = np.linspace(-10,10,num=40)
                rlist = np.linspace(0,15,num=10)
                data = calsigmaHIlos(G,xlist,ylist,maxlength)
                newdata = SSF.calSFRsurdenxynew(sfrl,Sxl,Syl,Szl,xlist,ylist,maxlength)
                sigmaHIlosgrid = data['sigmaHIlosgrid']
                HIm = data['HImgrid']
                Hpos = data['pos']
                sHIr = SSF.fromxytorad(sigmaHIlosgrid,Hpos,HIm,rlist)
                SFRdengrid = newdata['SFRdengrid']
                sarea = newdata['areagrid']
                spos = newdata['pos']
                SFRr = SSF.fromxytorad(SFRdengrid,spos,sarea,rlist)
                sigmaHIlosgridl=np.append(sigmaHIlosgridl,sHIr)
                SFRdengridl=np.append(SFRdengridl,SFRr)
            if wanted=='SFRturtot':
                Sscaledata = SSF.findscalelength(S, halfmassscale=0, usecylr=1, usecylz=0, cutrout =10., cutzout = 5.)
                rs = Sscaledata['rs']
                Sscalearea = 4.0*np.pi*rs*rs #in kpc2
                SFRtot = np.sum(sfrl) #in Msun/yr
                SFRperarea = SFRtot/Sscalearea
                SFRdengridl=np.append(SFRdengridl,SFRperarea)
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
                print 'rhol', rhol
                massl = volll/rhol
                vturHI2l = pturHIl/rhoHIl*massl
                vth2l = pthl/rhol*GAMMA*massl
                massl[~np.isfinite(massl)]=0.0
                vturHI2l[~np.isfinite(vturHI2l)]=0.0
                vth2l[~np.isfinite(vth2l)]=0.0
                vHIdis = np.sqrt(np.sum(vturHI2l+vth2l)/np.sum(massl))/km_in_cm
                print 'runtodo', runtodo
                print 'np.sqrt(np.sum(vturHI2l))',  np.sqrt(np.sum(vturHI2l))
                print 'np.sqrt(np.sum(vth2l))', np.sqrt(np.sum(vth2l))
                sigmaHIlosgridl = np.append(sigmaHIlosgridl,vHIdis)
                runtodol = np.append(runtodol, runtodo)
                #print 'vHIdis', vHIdis
                #print 'runtodo', runtodo
                #print 'i', i
            colori = highres
            coloril = np.append(coloril, colori)
    #print 'coloril', coloril
    #print 'sigmaHIlosgridl', sigmaHIlosgridl
    #print 'SFRdengridl', SFRdengridl
    if wanted=='SFRturrad':
        ax.plot(SFRdengridl,sigmaHIlosgridl,color='k',ls='None',marker='s')
    if wanted=='SFRturtot':
        neederrbar=1
        if neederrbar==1:
            for runtodo in dirneed:
                Sds = SFRdengridl[runtodol==runtodo]
                sHs = sigmaHIlosgridl[runtodol==runtodo]
                cls = coloril[runtodol==runtodo]
                Sdmed = np.median(Sds)
                sHmed = np.median(sHs)
                clmed = np.median(cls)
                Sdu = np.percentile(Sds,84)-Sdmed
                Sdd = Sdmed-np.percentile(Sds,15.8)
                sHu = np.percentile(sHs,84)-sHmed
                sHd = sHmed-np.percentile(sHs,15.8)
                print 'Sdu,Sdmed,Sdd', Sdu,Sdmed,Sdd
                print 'sHu,sHmed,sHd', sHu,sHmed,sHd
                if clmed==0:
                    color='k'; marker='s';
                if clmed==1:
                    color='r'; marker='^';
                if clmed==3:
                    color='g'; marker='o';
                    
                ax.errorbar(Sdmed,sHmed, xerr=[[Sdd],[Sdu]], yerr=[[sHd],[sHu]],color=color,fmt=marker,markersize=7)
        else:
            ax.plot(SFRdengridl[coloril==0],sigmaHIlosgridl[coloril==0],color='k',ls='None',marker='s')
            ax.plot(SFRdengridl[coloril==1],sigmaHIlosgridl[coloril==1],color='r',ls='None',marker='^')
            ax.plot(SFRdengridl[coloril==3],sigmaHIlosgridl[coloril==3],color='g',ls='None',marker='o')
        plotobserved=1
        if plotobserved==1:
            ftxt = open(programdir+'/data/DibsigmaSFRA.txt', 'r')
            ftxt.readline()
            dars = ftxt.readlines()
            ftxt.close()
            obslogSFRA=[]
            obssigma=[]
            for line in dars:
                    xsd = line.split()
                    obslogSFRA=np.append(obslogSFRA, float(xsd[0]))
                    obssigma=np.append(obssigma, float(xsd[1]))
            ax.plot(np.power(10,obslogSFRA),obssigma,markeredgecolor='k',markerfacecolor="None",ls='None',marker='D')
            #ax.plot(np.power(10,-1.59832), 34.8944,markeredgecolor='k',markerfacecolor="k",ls='None',marker='D')
    ax.set_ylim(ymin=1,ymax=50)
    ax.set_xlim(xmin=1e-5,xmax=0.5)
    xlab = r'${\rm \Sigma_{SFR}\;[M_\odot/yr/kpc^2]}$'
    ylab = r'${\rm  \sigma_{HI}[km/s]}$'             
    ax.set_title(ptitle, fontsize=16)
    PS.miscsetup(ax,logx=1,logy=0,xlab=xlab,ylab=ylab,legendneed=0,labfs=22,legfs=12)
    import matplotlib.lines as mlines
    potk = mlines.Line2D([], [], color='k', ls='None',marker='s')
    potr = mlines.Line2D([], [], color='r', ls='None',marker='^')
    potg = mlines.Line2D([], [], color='g', ls='None',marker='o')
    poto = mlines.Line2D([], [],markeredgecolor='k',markerfacecolor="None",ls='None',marker='D')
    lxlist=[
        #potk,
        potr,
        potg,
        poto
        ]
    lablist=[
        #'no CR', 
        #'MHD no CR',
        'MHD+',
        'CR+',
        #r'MHD $\kappa$=3e29',
        'Obs'
        ]
    legend1 = plt.legend(lxlist, lablist, loc='upper left',fontsize=14,ncol=1)
    plt.gca().add_artist(legend1)
    if wanted=='SFRturrad':
        filename=plotloc+'CRplot/sfrtur/Sigmasfr_sigmaHI_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    if wanted=='SFRturtot':
        filename=plotloc+'CRplot/sfrtur/Sigmasfr_sigmaHItot_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    PS.finishsave(plt,filename)
                
                
                
