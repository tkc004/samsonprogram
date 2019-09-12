from stdmodandoption import *
import cameron_functions as CF
import collections

def Faradayrotation_out(ssdict):    
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
    title='MW'
    titleneed=title
    needlog=0
    dclabelneed=1
    correctIa=0
    useM1=1
    nogrid=40
    med = -0.1 
    wholeboxgas=1
    diskgas=1
    Rfrac = 0.5
    nosum=0
    xaxis_snapno=0
    if not (wanted=='Faradayrotation'):
        print 'wrong wanted'
        return None
    resoneed=0
    rotface=1
    loccen=1
    newlabelneed=1
    print 'runtodo', runtodo
    if wanted=='Faradayrotation':
        nogrid = 15
        dr = withinr/nogrid
        radlist=np.linspace(withoutr,withinr,num=nogrid)
        xlist = np.linspace(-withinr,withinr,num=withinr*4.0)
        ylist = np.linspace(-withinr,withinr,num=withinr*4.0)
        dx = xlist[1]-xlist[0]; dy = ylist[1]-ylist[0];
        frmlist=radlist*0.
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
        Gextra = SSF.readsnapwcen(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix,\
         havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
         loccen=loccen,runtodo=runtodo,rundir=rundir,halostr=halostr)
        Gx = Gextra['p'][:,0]; Gy = Gextra['p'][:,1]; Gz = Gextra['p'][:,2];
        Gr = np.sqrt(Gx*Gx+Gy*Gy+Gz*Gz)
        Grho = Gextra['rho']; Gu = Gextra['u']; Gm = Gextra['m']; Gnh = Gextra['nh'];
        ne = Gextra['ne'];
        Gh = Gextra['h']; Gmetal = Gextra['z']; Gcregy = Gextra['cregy']
        density = Grho*1e10*Msun_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm
        Zmetal = Gmetal[:,0]; ZHe = Gmetal[:,1]
        edensity = density/protonmass_in_g*(1.0-Zmetal-ZHe)*ne
        GBx = Gextra['B'][:,0]; GBy = Gextra['B'][:,1]; GBz = Gextra['B'][:,2];
        RMprefactor = 2.64e-17 #unit in 1/Gauss
        RMdensity = RMprefactor*edensity*GBz
        GV = Gm/Grho #kpc^3
        
        
        if wanted=='Faradayrotation':
            cutz = np.absolute(Gz)<maxlength/2.
            x = Gx[cutz]; y = Gy[cutz]; weights=RMdensity[cutz]*GV[cutz]/dx/dy*kpc_in_cm # in rad/cm^2
            meshx,meshy,meshweights = SSF.depositandcor(x,y,weights,xlist,ylist)
            meshr = np.sqrt(meshx*meshx+meshy*meshy)
            for irad in range(len(radlist)-1):
                cutr = (meshr < radlist[irad+1])*(meshr > radlist[irad]) 
                cut = cutr
                fRM = np.sum(np.absolute(meshweights[cut]))/len(meshweights[cut])
                frmlist[irad] += fRM*m_in_cm*m_in_cm
        print 'np.amax(meshweights)', np.amax(meshweights)
        numoftimes+=1
        if haveB>0:
                lsn='dashed'
        else:
                lsn='solid'
        if resoneed==1:
                if resolabel=='llr':
                    labelneed='low res'
                    lsn = 'solid'
                if resolabel=='lr':
                    labelneed='std res'
                    lsn = 'dashed'
                if resolabel=='mr':
                    labelneed='high res'
                    lsn = 'dashdot'
        if wanted=='Faradayrotation':
            rmlist = (radlist[:-1]+radlist[1:])/2.
            plotdict[wanted]['xlab'] = r'${\rm r\;[kpc]}$'
            plotdict[wanted]['xnl'] = rmlist
            if wanted=='Faradayrotation':
                plotdict[wanted]['ylab'] = r'${\rm RM}\;{\rm [rad/m^2]}$'
                plotdict[wanted]['ynl'] = frmlist[:-1]/numoftimes    
        plotdict[wanted]['runtodo'] = runtodo
        plotdict[wanted]['labelneed'] = labelneed
        plotdict[wanted]['lsn'] = '--'
        plotdict[wanted]['lw'] = 2
        plotdict[wanted]['marker'] = 'None'
        plotdict[wanted]['color'] = color
        plotdict[wanted]['runtitle'] = runtitle
        plotdict[wanted]['ptitle'] = ptitle
        if wanted=='Faradayrotation':
            filename=plotloc+'CRplot/Faradayrotation/Faradayrotation_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        plotdict[wanted]['filename'] = filename        
        return plotdict
