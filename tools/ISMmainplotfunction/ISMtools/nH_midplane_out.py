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
    title='MW'
    titleneed=title
    needlog=0
    dclabelneed=1
    correctIa=0
    useM1=1
    withinr=20.0
    nogrid=40
    med = -0.1 
    wholeboxgas=1
    diskgas=1
    Rfrac = 0.5
    nosum=0
    xaxis_snapno=0
    if wanted=='nHmidplane':
        resoneed=0
        rotface=1
        newlabelneed=1
        print 'runtodo', runtodo
        withinr=15.
        nogrid = 15
        maxlength=0.2 #thickness
        dr = withinr/nogrid
        radlist=np.linspace(0.001,withinr,num=nogrid)
        if wanted=='nHmidplane':
            nHIlist=radlist*0.
            nHIIlist=radlist*0.
            nH2list=radlist*0.
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
            if cosmo==1:
                datasup=0;
            else:
                datasup=1;
            Gextra = SSF.readsnapwcen(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix,\
             havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
             datasup=datasup,runtodo=runtodo,rundir=rundir,halostr=halostr)
            Gx = Gextra['p'][:,0]; Gy = Gextra['p'][:,1]; Gz = Gextra['p'][:,2];
            Gvx = Gextra['v'][:,0]; Gvy = Gextra['v'][:,1]; Gvz = Gextra['v'][:,2];
            Grho = Gextra['rho']; Gu = Gextra['u']; Gm = Gextra['m']; Gnh = Gextra['nh'];
            Gh = Gextra['h']; Gmetal = Gextra['z'];
            datanH = CF.calnH(Gm,Gnh,Gh,Grho,Gmetal)
            NHI = datanH['NHI']; NHII = datanH['NHII']; NH2 = datanH['NH2'];
            for irad in range(len(radlist)-1):
                cutxy = ((Gx)*(Gx)+(Gy)*(Gy) > radlist[irad]*radlist[irad]) & ((Gx)*(Gx)+(Gy)*(Gy) < radlist[irad+1]*radlist[irad+1])
                cutz = (Gz)*(Gz) < maxlength*maxlength/4.
                cut = cutxy*cutz
                shellvol_in_cm3 = np.pi*(-np.power(radlist[irad],2)+np.power(radlist[irad+1],2))*kpc_in_cm*kpc_in_cm*kpc_in_cm*maxlength
                if wanted=='nHmidplane':
                    nHI = np.sum(NHI[cut])/shellvol_in_cm3
                    nHII = np.sum(NHII[cut])/shellvol_in_cm3
                    nH2 = np.sum(NH2[cut])/shellvol_in_cm3
                    nHIlist[irad] += nHI
                    nHIIlist[irad] += nHII
                    nH2list[irad] += nH2
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
            if wanted=='nHmidplane':
                plotdict[wanted]['xlab'] = r'${\rm r\;[kpc]}$'
                plotdict[wanted]['xnl']['nH2'] = radlist[:-1]
                plotdict[wanted]['xnl']['nHI'] = radlist[:-1]
                plotdict[wanted]['xnl']['nHII'] = radlist[:-1]
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
                plotdict[wanted]['lw'] = 2
                plotdict[wanted]['marker'] = 'None'
                plotdict[wanted]['color'] = color
                plotdict[wanted]['runtitle'] = runtitle
                plotdict[wanted]['ptitle'] = ptitle
        if wanted=='nHmidplane':
            filename=plotloc+'CRplot/nHmidplane/nH_midplane_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            plotdict[wanted]['filename'] = filename
        return plotdict
