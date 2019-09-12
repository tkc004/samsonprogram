from stdmodandoption import *
import cameron_functions as CF
import collections

def nH_midplane_out(ssdict):    
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
    phicut=ssdict['phicut']
    useobs=ssdict['useobs']
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
    if not (wanted=='nHmidplane' or wanted=='nHz'):
        print 'wrong wanted'
        return None
    resoneed=0
    rotface=1
    loccen=1
    newlabelneed=1
    print 'runtodo', runtodo
    if wanted=='nHmidplane':
        nogrid = 30
        dr = withinr/nogrid
        radlist=np.linspace(withoutr,withinr,num=nogrid)
        nHIlist=radlist*0.
        nHIIlist=radlist*0.
        nH2list=radlist*0.
    if wanted=='nHz':
        nogrid = 40
        dr = withinr/nogrid
        zlist=SSF.symlogspace(maxlength/2.,nogrid,base=30.0,halfz=0)
        #zlist=np.linspace(-maxlength/2.,maxlength/2.,num=nogrid)
        nHIlist=zlist*0.
        nHIIlist=zlist*0.
        nH2list=zlist*0.            
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
        Grxy = np.sqrt(Gx*Gx+Gy*Gy)
        Gvx = Gextra['v'][:,0]; Gvy = Gextra['v'][:,1]; Gvz = Gextra['v'][:,2];
        Grho = Gextra['rho']; Gu = Gextra['u']; Gm = Gextra['m']; Gnh = Gextra['nh'];
        Gh = Gextra['h']; Gmetal = Gextra['z'];
        #print 'Gm', Gm
        #print 'np.amax(Gnh)', np.amax(Gnh)
        #print 'np.amin(Gnh)', np.amin(Gnh)
        datanH = CF.calnH(Gm,Gnh,Gh,Grho,Gmetal)
        #print 'Gm', Gm
        NHI = datanH['NHI']; NHII = datanH['NHII']; NH2 = datanH['NH2'];
        if wanted=='nHmidplane':
            for irad in range(len(radlist)-1):
                cutxy = (Grxy > radlist[irad]) & (Grxy < radlist[irad+1])
                cutz = (Gz)*(Gz) < maxlength*maxlength/4.
                cut = cutxy*cutz
                #shellvol_in_cm3 = np.pi*(-np.power(radlist[irad],2)+np.power(radlist[irad+1],2))*kpc_in_cm*kpc_in_cm*kpc_in_cm*maxlength
                if phicut==1:
                    phi = np.arctan2(Gy,Gx)
                    cutphi = np.absolute(phi) < np.pi/10.0
                    cut=cut*cutphi
                Gmcut = Gm[cut]; Grhocut=Grho[cut];
                
                shellvol_in_cm3 = np.sum(Gmcut/Grhocut)*kpc_in_cm*kpc_in_cm*kpc_in_cm
                nHI = np.sum(NHI[cut])/shellvol_in_cm3
                nHII = np.sum(NHII[cut])/shellvol_in_cm3
                nH2 = np.sum(NH2[cut])/shellvol_in_cm3
                nHIlist[irad] += nHI
                nHIIlist[irad] += nHII
                nH2list[irad] += nH2
        if wanted=='nHz':
            for iz in range(len(zlist)-1):
                cutxy = (Grxy > withoutr) & (Grxy < withinr)
                cutz = (Gz > zlist[iz]) & (Gz < zlist[iz+1])
                cut = cutxy*cutz
                #shellvol_in_cm3 = np.pi*(-withoutr*withoutr+withinr*withinr)*kpc_in_cm*kpc_in_cm*kpc_in_cm*(zlist[iz+1]-zlist[iz])
                if phicut==1:
                    phi = np.arctan2(Gy,Gx)
                    cutphi = np.absolute(phi) < np.pi/10.0
                    cut=cut*cutphi
                Gmcut = Gm[cut]; Grhocut=Grho[cut];
                #print 'Gmcut', Gmcut
                shellvol_in_cm3 = np.sum(Gmcut/Grhocut)*kpc_in_cm*kpc_in_cm*kpc_in_cm
                nHI = np.sum(NHI[cut])/shellvol_in_cm3
                nHII = np.sum(NHII[cut])/shellvol_in_cm3
                nH2 = np.sum(NH2[cut])/shellvol_in_cm3
                nHIlist[iz] += nHI
                nHIIlist[iz] += nHII
                nH2list[iz] += nH2                
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
    linelist = []            
    if wanted=='nHmidplane':
        rmlist = (radlist[:-1]+radlist[1:])/2.
        plotdict[wanted]['xlab'] = r'${\rm r\;[kpc]}$'
        plotdict[wanted]['xnl']['nH2'] = rmlist
        plotdict[wanted]['xnl']['nHI'] = rmlist
        plotdict[wanted]['xnl']['nHII'] = rmlist
    if wanted=='nHz':
        zmlist = (zlist[:-1]+zlist[1:])/2.
        plotdict[wanted]['xlab'] = r'${\rm z\;[kpc]}$'
        plotdict[wanted]['xnl']['nH2'] = zmlist
        plotdict[wanted]['xnl']['nHI'] = zmlist
        plotdict[wanted]['xnl']['nHII'] = zmlist
        if useobs==1:
            #HI distribution from Dickey Lockman 1990; parameterized in Wood 2010
            nHIobs = 0.4*np.exp(-0.5*zmlist*zmlist/0.09/0.09)+0.11*np.exp(-0.5*zmlist*zmlist/0.225/0.225)+0.06*np.exp(-np.absolute(zmlist)/0.4)
            #HII distribution from Reynolds 1984
            nHIIobs = 0.015*np.exp(-np.absolute(zmlist)/0.07)+0.025*np.exp(-np.absolute(zmlist)/0.9)
            #TKtest
            #nHIobs =5.0*nHIobs; nHIIobs = 5.0*nHIIobs;
            #testend
            plotdict[wanted]['xnl']['nHIobs'] = zmlist
            plotdict[wanted]['xnl']['nHIIobs'] = zmlist
            plotdict[wanted]['linelab']['nHIobs'] = 'HI obs'
            plotdict[wanted]['linelab']['nHIIobs'] = 'HII obs'
            plotdict[wanted]['ynl']['nHIobs'] = nHIobs
            plotdict[wanted]['ynl']['nHIIobs'] = nHIIobs
            plotdict[wanted]['lsn']['nHIobs'] = '--'
            plotdict[wanted]['lsn']['nHIIobs'] = '-.'
    plotdict[wanted]['ylab'] = r'$<{\rm nH}>\;[{\rm atom/cm^3}]$'
    plotdict[wanted]['ynl']['nHI'] = nHIlist[:-1]/numoftimes
    plotdict[wanted]['ynl']['nH2'] = nH2list[:-1]/numoftimes
    plotdict[wanted]['ynl']['nHII'] = nHIIlist[:-1]/numoftimes
    plotdict[wanted]['runtodo'] = runtodo
    plotdict[wanted]['labelneed'] = labelneed
    plotdict[wanted]['lsn']['nHI'] = '--'
    plotdict[wanted]['lsn']['nH2'] = '-'
    plotdict[wanted]['lsn']['nHII'] = '-.'
    plotdict[wanted]['linelab']['nHI'] = 'HI'
    plotdict[wanted]['linelab']['nH2'] = 'H2'
    plotdict[wanted]['linelab']['nHII'] = 'HII'
    linelist.append('nHI')
    linelist.append('nH2')
    linelist.append('nHII')
    for inkey in linelist:
        plotdict[wanted]['lw'][inkey] = 2
        plotdict[wanted]['marker'][inkey] = 'None'
        plotdict[wanted]['color'][inkey] = color
        plotdict[wanted]['runtitle'][inkey] = runtitle

    if useobs==1 and wanted=='nHz':
        for inkey in ['nHIobs','nHIIobs']:
            plotdict[wanted]['lw'][inkey] = 1
            plotdict[wanted]['marker'][inkey] = 'None'
            plotdict[wanted]['color'][inkey] = 'k'
            plotdict[wanted]['runtitle'][inkey] = runtitle               
        linelist.append('nHIobs')
        linelist.append('nHIIobs')
    plotdict[wanted]['linelist']=linelist
    plotdict[wanted]['ptitle'] = labelneed
    if wanted=='nHmidplane':
        filename=plotloc+'CRplot/nHmidplane/nH_midplane_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        plotdict[wanted]['filename'] = filename
    if wanted=='nHz':
        filename=plotloc+'CRplot/nHz/nHz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        plotdict[wanted]['filename'] = filename        
    return plotdict
