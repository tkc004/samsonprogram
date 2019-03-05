from stdmodandoption import *
import collections

def CR_midplane_out(runtodo,wanted,startno,Nsnap,snapsep,fmeat):
    title='MW'
    titleneed=title
    needlog=0
    dclabelneed=1
    correctIa=0
    useM1=1
    the_prefix='snapshot'
    the_suffix='.hdf5'
    withinr=15.0
    nogrid=40
    maxlength=10.0
    med = -0.1 
    wholeboxgas=1
    diskgas=1
    Rfrac = 0.5
    nosum=0
    xaxis_snapno=0

    plotdict = collections.defaultdict(dict)
    if wanted=='crdenmidplane' or wanted=='gasdenmidplane' or wanted=='stardenmidplane': #cosmic ray energy density at solar circle (observation ~ 1 eV/cm^3)
        resoneed=0
        rotface=1
        newlabelneed=1
        print 'runtodo', runtodo
        if wanted=='crdenmidplane':
            info=SSF.outdirname(runtodo, Nsnap)
            havecr=info['havecr']
            dclabel=info['dclabel']
            haveB=info['haveB']
            if havecr==0:
                return None
        withinr=15.
        nogrid = 15
        maxlength=0.1 #thickness
        dr = withinr/nogrid
        radlist=np.linspace(0.001,withinr,num=nogrid)
        if wanted=='crdenmidplane':
            credenlist=radlist*0.
        if wanted=='gasdenmidplane':
            gasdenlist=radlist*0.
        if wanted=='stardenmidplane':
            stardenlist=radlist*0.            
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
            if havecr==0 and wanted=='crdenmidplane':
                continue
            if cosmo==1:
                h0=1
            else:
                h0=0
            if cosmo==1:
                datasup=0;
            else:
                datasup=1;
            if wanted=='crdenmidplane' or wanted == 'gasdenmidplane': ptype=0;
            if wanted == 'stardenmidplane': ptype=4;
            Gextra = SSF.readsnapwcen(the_snapdir, Nsnapstring, ptype, snapshot_name=the_prefix, extension=the_suffix,\
             havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
             datasup=datasup,runtodo=runtodo,rundir=rundir,halostr=halostr)
            Gx = Gextra['p'][:,0]; Gy = Gextra['p'][:,1]; Gz = Gextra['p'][:,2];
            Gvx = Gextra['v'][:,0]; Gvy = Gextra['v'][:,1]; Gvz = Gextra['v'][:,2];
            Gm = Gextra['m']
            if wanted=='crdenmidplane':
                cregy = Gextra['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                cregy_in_erg = cregy*solar_mass_in_g*1e10*km_in_cm*km_in_cm
                cregy_in_eV = cregy_in_erg*erg_in_eV
            for irad in range(len(radlist)-1):
                cutxy = ((Gx)*(Gx)+(Gy)*(Gy) > radlist[irad]*radlist[irad]) & ((Gx)*(Gx)+(Gy)*(Gy) < radlist[irad+1]*radlist[irad+1])
                cutz = (Gz)*(Gz) < maxlength*maxlength/4.
                cut = cutxy*cutz
                shellvol_in_cm3 = np.pi*(-np.power(radlist[irad],2)+np.power(radlist[irad+1],2))*kpc_in_cm*kpc_in_cm*kpc_in_cm*maxlength
                if wanted=='crdenmidplane':
                    cregy_in_eV_cut = np.sum(cregy_in_eV[cut])
                    creden_in_eV_per_cm3 = cregy_in_eV_cut/shellvol_in_cm3
                    credenlist[irad]+=creden_in_eV_per_cm3
                if wanted=='gasdenmidplane':
                    Gm_in_g=Gm[cut]*1e10*Msun_in_g
                    Grho_in_g_cm_3 = np.sum(Gm_in_g)/shellvol_in_cm3
                    Gnism_in_cm_3 = 1.0/proton_mass_in_g*Grho_in_g_cm_3
                    gasdenlist[irad]+=Gnism_in_cm_3
                if wanted=='stardenmidplane':
                    Gm_in_g=Gm[cut]*1e10*Msun_in_g
                    Grho_in_g_cm_3 = np.sum(Gm_in_g)/shellvol_in_cm3
                    stardenlist[irad]+=Grho_in_g_cm_3
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
            #if M1speed>999:
            #   lsn='dotted'
            #if M1speed>1999:
        #       lsn='dashdot'
            if wanted=='crdenmidplane': 
                plotdict[wanted]['xlab'] = r'${\rm r [kpc]}$'
                plotdict[wanted]['xnl'] = radlist[:-1]
                plotdict[wanted]['ylab'] = r'$<e_{\rm CR}>[{\rm eV/cm^3}]$'
                plotdict[wanted]['ynl'] = credenlist[:-1]/numoftimes
                plotdict[wanted]['runtodo'] = runtodo
                plotdict[wanted]['labelneed'] = labelneed
                plotdict[wanted]['lsn'] = lsn
                plotdict[wanted]['color'] = color
                plotdict[wanted]['runtitle'] = runtitle
                plotdict[wanted]['ptitle'] = ptitle
            if wanted=='gasdenmidplane':
                plotdict[wanted]['xlab'] = r'${\rm r [kpc]}$'
                plotdict[wanted]['xnl'] = radlist[:-1]
                plotdict[wanted]['ylab'] = r'$<n_{\rm ISM}>[{\rm cm^{-3}}]$'
                plotdict[wanted]['ynl'] = gasdenlist[:-1]/numoftimes
                plotdict[wanted]['runtodo'] = runtodo
                plotdict[wanted]['labelneed'] = labelneed
                plotdict[wanted]['lsn'] = lsn
                plotdict[wanted]['color'] = color
                plotdict[wanted]['runtitle'] = runtitle
                plotdict[wanted]['ptitle'] = ptitle
            if wanted=='stardenmidplane':
                plotdict[wanted]['xlab'] = r'${\rm r [kpc]}$'
                plotdict[wanted]['xnl'] = radlist[:-1]
                plotdict[wanted]['ylab'] = r'$<\rho_*>[{\rm g/cm^{3}}]$'
                plotdict[wanted]['ynl'] = stardenlist[:-1]/numoftimes
                plotdict[wanted]['runtodo'] = runtodo
                plotdict[wanted]['labelneed'] = labelneed
                plotdict[wanted]['lsn'] = lsn
                plotdict[wanted]['color'] = color
                plotdict[wanted]['runtitle'] = runtitle
                plotdict[wanted]['ptitle'] = ptitle
        if wanted=='crdenmidplane':
            filename=homedir+'/CRplot/crdenmidplane/CR_midplane_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            plotdict[wanted]['filename'] = filename
        if wanted=='gasdenmidplane':
            filename=homedir+'/CRplot/gasdenmidplane/gas_midplane_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            plotdict[wanted]['filename'] = filename
        if wanted=='stardenmidplane':
            filename=homedir+'/CRplot/stardenmidplane/star_midplane_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            plotdict[wanted]['filename'] = filename
        return plotdict
