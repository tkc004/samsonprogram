from stdmodandoption import *
#matplotlib.use('agg')
def calmassloading(dirneed,Nsnaplist,avesnap=1):
    tsep=10.0 #time step in Myr to calculate SFR
    sfrl=np.array([])
    mll =np.array([]) #massloading factor
    mvl =np.array([])
    mwl =np.array([]) #outflow rate
    msl =np.array([])
    runtodol = np.array([])
    dcl = np.array([])
    nsml = np.array([])
    for runtodo in dirneed:
        for Nsnap in Nsnaplist:
            #haloinfo=cosmichalo(runtodo)
            haloinfo=SSF.outdirname(runtodo)
            rundir=haloinfo['rundir']
            maindir=haloinfo['maindir']
            print 'maindir', maindir
            subdir=haloinfo['subdir']
            halostr=haloinfo['halostr']
            snumadd=haloinfo['snumadd']
            usepep=haloinfo['usepep']
            hubble=haloinfo['h0']
            firever=haloinfo['firever']
            highres=haloinfo['highres']
            Rvir=haloinfo['Rvir']
            cosmo=haloinfo['cosmo']
            Ms = haloinfo['Msini']
            Mv = haloinfo['Mhini']
            if cosmo==1:
                try:
                    if usepep==0:
                        halosA = SF.read_halo_history(rundir, halonostr=halostr,comoving=0,hubble=hubble,maindir=maindir,snumadd=snumadd)
                    else:
                        halosA = SF.read_halo_history_pep(rundir, Nsnap, singlesnap=1, firever=firever, halonostr=halostr, comoving=0, maindir=maindir)
                    redlist = halosA['redshift']
                    Rvirlist = halosA['R']
                    Mslist = halosA['Ms']
                    Mvlist = halosA['M']
                    #print 'Mslist', Mslist
                    #print 'Mvlist', Mvlist
                except (IOError,KeyError):
                    continue
                try:
                    Rvir = Rvirlist[-1]
                    Ms = Mslist[-1]
                    Mv = Mvlist[-1]
                except TypeError:
                    Rvir=Rvirlist; Ms = Mslist; Mv = Mvlist;
                
            Rup = Rvir*0.3
            Rdown = Rvir*0.2
            try:
                emdata=CRTF.gaswindphase(runtodo,Nsnap,rotface=0,withinr=100.0,zup=Rup,zdown=Rdown,userad=1)
                sfrdata=CRTF.outsfr(runtodo, Nsnap,tsep=tsep)
                vr =  emdata['vr'] #in km/s
                Gmass = emdata['Gmass'] #in 1e10Msun
                print 'np.median(vr)', np.median(vr)
                print 'np.sum(Gmass)', np.sum(Gmass)
                SFR = sfrdata['SFR']# in Msun/yr
                Nsm = sfrdata['Nsm']# in solar mass
                print 'runtodo', runtodo
                print 'dc', highres
                print 'Gmass', np.sum(Gmass)
                print 'SFR', SFR
                print 'Rvir', Rvir
                mwind = np.sum(Gmass*1e10*vr*km_in_cm)/(Rup-Rdown)/kpc_in_cm*yr_in_sec #Msun/yr
                ml = mwind/SFR
                msl=np.append(msl,Ms)
                mvl=np.append(mvl,Mv)
                nsml = np.append(nsml,Nsm)
                mwl = np.append(mwl,mwind)
                mll=np.append(mll,ml)
                sfrl=np.append(sfrl,SFR)
                runtodol= np.append(runtodol,runtodo)
                dcl=np.append(dcl,highres)
            except (IOError,KeyError):
                continue
            
    print 'dcl', dcl
    print 'dcl==0', dcl==0




    if avesnap==1:
        print 'runtodol', runtodol
        print 'mvl', mvl
        print 'sfrl', sfrl
        cutinf = np.isfinite(mvl)
        
        mvlold=mvl[cutinf];mslold=msl[cutinf];mllold=mll[cutinf];sfrlold=sfrl[cutinf];
        runtodolold=runtodol[cutinf];mwlold=mwl[cutinf];dclold=dcl[cutinf];
        nsmoldl=nsml[cutinf]
        mlld=mllu=mvl=msl=mll=runtodol=mwl=dcl=sfrl=avensml=mwmedl=mwdownl=mwupl=nsmmedl=np.array([])
        for runtodo in dirneed:
            print 'runtodo', runtodo 
            print 'mvlold', mvlold[runtodolold==runtodo], 
            print 'ave', np.average(mvlold[runtodolold==runtodo])
            numoftimes = len(mvlold[runtodolold==runtodo])
            mvl=np.append(mvl,np.average(mvlold[runtodolold==runtodo]))
            msl=np.append(msl,np.average(mslold[runtodolold==runtodo]))
            mlmedian = np.median(mllold[runtodolold==runtodo])
            mlup = np.percentile(mllold[runtodolold==runtodo],84)-mlmedian
            mldown = mlmedian-np.percentile(mllold[runtodolold==runtodo],15.8)
            mll=np.append(mll,mlmedian)
            mllu=np.append(mllu,mlup)
            mlld=np.append(mlld,mldown)
            mwmed = np.median(mwlold[runtodolold==runtodo])
            mwup = np.percentile(mwlold[runtodolold==runtodo],84)-mwmed
            mwdown = mwmed-np.percentile(mwlold[runtodolold==runtodo],15.8)
            mwdownl = np.append(mwdownl,mwdown)
            mwupl = np.append(mwupl,mwup)
            dcl=np.append(dcl,np.median(dclold[runtodolold==runtodo]))
            sfrl=np.append(sfrl,np.average(sfrlold[runtodolold==runtodo]))
            runtodol=np.append(runtodol,np.unique(runtodolold[runtodolold==runtodo]))
            mwmed = np.median(mwlold[runtodolold==runtodo])
            mwup = np.percentile(mwlold[runtodolold==runtodo],84)-mwmed
            mwdown = mwmed-np.percentile(mwlold[runtodolold==runtodo],15.8)
            mwl=np.append(mwl,np.average(mwlold[runtodolold==runtodo]))
            mwmedl = np.append(mwmedl,mwmed)
            mwdownl = np.append(mwdownl,mwdown)
            mwupl = np.append(mwupl,mwup)
            if avesnap==1: avensml = np.append(avensml,np.average(nsmoldl[runtodolold==runtodo]))
                
            nsmmedl = np.append(nsmmedl,np.median(nsmoldl[runtodolold==runtodo]))
    if avesnap==1: mlla = mwl/(avensml/numoftimes/tsep/1e6)
    if avesnap==0:
        return {'sfrl':sfrl,'mll':mll,'mvl':mvl,'mwl':mwl,'msl':msl,'runtodol':runtodol,'dcl':dcl,'nsml':nsml}
    if avesnap==1:
        return {'sfrl':sfrl,'mll':mll,'mvl':mvl,'mwl':mwl,'msl':msl,'runtodol':runtodol,'dcl':dcl,'nsml':nsml,
                'mlld':mlld,'mllu':mllu,'dcl':dcl,'avensml':avensml,'mwmedl':mwmedl,'mwdownl':mwdownl,'mwupl':mwupl,
                'nsmmedl':nsmmedl}
