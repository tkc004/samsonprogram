from stdmodandoption import *
import cameron_functions as CF
import collections

def synchrotron_out(ssdict):    
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
    if not (wanted=='syncr' or wanted=='Bmin'):
        print 'wrong wanted'
        return None
    resoneed=0
    rotface=1
    loccen=1
    newlabelneed=1
    print 'runtodo', runtodo
    if wanted=='syncr' or wanted=='Bmin':
        nogrid = 15
        dr = withinr/nogrid
        radlist=np.linspace(withoutr,withinr,num=nogrid)
        synlist=radlist*0.
        Bminlist=radlist*0.
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
        Gh = Gextra['h']; Gmetal = Gextra['z']; Gcregy = Gextra['cregy']
        GBx = Gextra['B'][:,0]; GBy = Gextra['B'][:,1]; GBz = Gextra['B'][:,2];
        B_in_T = np.sqrt(GBx*GBx+GBy*GBy+GBz*GBz)/Tesla_in_Gauss
        GV = Gm/Grho #kpc^3
        GV_in_m3 = GV*kpc_in_m*kpc_in_m*kpc_in_m
        beta = 100. #fraction between CR proton and electron energy
        eta = beta+1
        nuneed = 1e9 #synchrotron frequency
        Gecre = Gcregy/beta*cregy_codetocgs*erg_in_GeV/GV_in_m3 #CR electron density in GeV/m^3
        
        if wanted=='syncr' or wanted=='Bmin':
            for irad in range(len(radlist)-1):
                cutr = (Gr < radlist[irad+1])*(Gr > radlist[irad])
                cutz = (np.absolute(Gz)<maxlength/2.) 
                cut = cutr*cutz
                vol_in_cm3 = np.pi*(np.power(radlist[irad+1],2)-np.power(radlist[irad],2))*kpc_in_cm*kpc_in_cm*kpc_in_cm*maxlength
                vol_in_m3 = vol_in_cm3/m_in_cm/m_in_cm/m_in_cm
                Gecrecut = Gecre[cut]; B_in_Tcut = B_in_T[cut]; GV_in_m3cut = GV_in_m3[cut];
                radiolum_in_W_Hz = CRTF.calsyndefault(B_in_Tcut,Gecrecut,GV_in_m3cut,nu=nuneed)
                synlist[irad] += np.sum(radiolum_in_W_Hz)
                Bmin= CRTF.calBmin(eta,np.sum(radiolum_in_W_Hz),vol_in_m3,nuneed)
                Bminlist[irad] += Bmin
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
        if wanted=='syncr' or wanted=='Bmin':
            rmlist = (radlist[:-1]+radlist[1:])/2.
            plotdict[wanted]['xlab'] = r'${\rm r\;[kpc]}$'
            plotdict[wanted]['xnl'] = rmlist
            if wanted=='syncr':
                plotdict[wanted]['ylab'] = r'$L_{\rm syn}\;{\rm [W\,Hz^{-1}]}$'
                plotdict[wanted]['ynl'] = synlist[:-1]/numoftimes
            if wanted=='Bmin':
                plotdict[wanted]['ylab'] = r'$B_{\rm min}\;{\rm [G]}$'
                plotdict[wanted]['ynl'] = Bminlist[:-1]/numoftimes        
        plotdict[wanted]['runtodo'] = runtodo
        plotdict[wanted]['labelneed'] = labelneed
        plotdict[wanted]['lsn'] = '--'
        plotdict[wanted]['lw'] = 2
        plotdict[wanted]['marker'] = 'None'
        plotdict[wanted]['color'] = color
        plotdict[wanted]['runtitle'] = runtitle
        plotdict[wanted]['ptitle'] = ptitle
        if wanted=='syncr':
            filename=plotloc+'CRplot/synr/synr_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted=='Bmin':
            filename=plotloc+'CRplot/Bmin/Bmin_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        plotdict[wanted]['filename'] = filename        
        return plotdict
