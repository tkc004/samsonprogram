def energyoutputold(runtodo,Nsnap,shiftz=1,usesolarcircle=0,rotface=1,usecentral=0,cylr=9.0,cylrin=0.01,cylz=0.5,nbin = 5):
        #see Su 2016 for the definition of turbulent energy
        #cylr=10. #kpc
        #cylz=2. #kpc #thickness of the cylinder we consider
        
        if usesolarcircle==1:
                cylrin=7.0; cylr=9.0
        if usecentral==1:
                cylrin=0.01; cylr=3.0
        fraccut = 68 #in percentage
        turl = []
        gml = []
        gminl = []
        cregyl = []
        therml = []
        Begyl = []
        xcell = []
        ycell = []
        zcell = []
        info=outdirname(runtodo, Nsnap)
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
        maindir=info['maindir']
        haveB=info['haveB']
        the_prefix=info['the_prefix']
        the_suffix=info['the_suffix']
        usepep=info['usepep']
        print 'the_snapdir', the_snapdir
        print 'Nsnapstring', Nsnapstring
        print 'havecr', havecr
        cosmo=info['cosmo']
        if cosmo==1:
                h0=1
        else:
                h0=0
        header = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo, header_only=1)
        ascale = header['time']
        #print 'this time', ascale
        thisred = 1./ascale-1.
        hubble = header['hubble']
        #print 'hubble', hubble
        G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo)
        if cosmo==1:
                if usepep==0:
                        halosingle = SF.read_halo_history(rundir, halonostr='00',hubble=hubble, comoving=0, maindir=maindir, singlesnap=1, atime=ascale)
                        afactor=1.0
                else:
                        halosingle = SF.read_halo_history_pep(rundir, Nsnap, singlesnap=1,comoving=1, halonostr='00', maindir=maindir,firever=2)
                        afactor=atime
                xcen = halosingle['x']*afactor
                ycen = halosingle['y']*afactor
                zcen = halosingle['z']*afactor
                xvcen = halosingle['xv']
                yvcen = halosingle['yv']
                zvcen = halosingle['zv']
                Rvirnow = halosingle['R']*afactor
                MgAHF = halosingle['Mg']
                Lxstar = halosingle['Lxstar']
                Lystar = halosingle['Lystar']
                Lzstar = halosingle['Lzstar']
                print 'xvcen', xvcen
                print 'xcen', xcen
        #       print 'xcenl[0]', xcenl[0]
        #       print 'thisred', thisred
                #print 'cen', xcen, ycen, zcen
                #print 'MgAHF', MgAHF
                #print 'Rvir', Rvirnow
        else:
                xcen=0
                ycen=0
                if shiftz==1:
                        zcen=findcenz(runtodo,Nsnap)
                else:
                        zcen=0.
                xvcen=0
                yvcen=0
                zvcen=0
        Gpos = G['p']
        Gvel = G['v']*km_in_cm #now in cm/s

        Gx = Gpos[:,0]-xcen
        Gy = Gpos[:,1]-ycen
        Gz = Gpos[:,2]-zcen
        Gvx = Gvel[:,0]-xvcen
        Gvy = Gvel[:,1]-yvcen
        Gvz = Gvel[:,2]-zvcen
        Grho = G['rho']
        if havecr>0:
                cregy  = G['cregy']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
        Gm = G['m']*1e10*solar_mass_in_g #in g
        #print 'np.average(Gvx)', np.average(Gvx)
        #print 'np.average(Gx)', np.average(Gx)
        Gtotm = np.sum(Gm)
        GEint = G['u']*km_in_cm*km_in_cm*Gm #in erg
        if haveB>0:
                GB = G['B']
                Bx = GB[:,0]
                By = GB[:,1]
                Bz = GB[:,2]
                B2 = Bx*Bx+By*By+Bz*Bz
                Begy = B2/np.pi/8.*Gm/(Grho*1e10*solar_mass_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm)

        if rotface==1:
                Gr =  np.sqrt(Gx*Gx+Gy*Gy+Gz*Gz)
                cutr = Gr < 5. #kpc
                Gxcutr = Gx[cutr]; Gycutr = Gy[cutr]; Gzcutr = Gz[cutr];
                Gvxcutr = Gvx[cutr]; Gvycutr = Gvy[cutr]; Gvzcutr = Gvz[cutr];
                Lang = [0.,0.,0.]
                for i in range(len(Gxcutr)):
                        Lang += np.cross([Gxcutr[i],Gycutr[i],Gzcutr[i]],[Gvxcutr[i],Gvycutr[i],Gvzcutr[i]])
                Gx, Gy, Gz = SF.rotateL_to_z(Gx,Gy,Gz,Lang[0],Lang[1],Lang[2])
                Gvx, Gvy, Gvz = SF.rotateL_to_z(Gvx,Gvy,Gvz,Lang[0],Lang[1],Lang[2])
                
        Grxy = np.sqrt(Gx*Gx+Gy*Gy)
        dr = (cylr-cylrin)/float(nbin)
        cutz = np.absolute(Gz)<cylz/2.
        radl=[]
        for nrad in range(nbin+1):
                radl = np.append(radl,cylrin+nrad*dr)
        for nrad in range(nbin):
                cutr = (Grxy<radl[nrad+1]) & (Grxy>radl[nrad])
                #print 'radl[nrad+1]', radl[nrad+1]
                #print 'radl[nrad]', radl[nrad]
                cut=cutr*cutz
                Gthc=GEint[cut]
                Gmc=Gm[cut]
                Grc=Grxy[cut]
                Gxc=Gx[cut]
                Gyc=Gy[cut]
                Gzc=Gz[cut]
                Gvxc = Gvx[cut]
                Gvyc = Gvy[cut]
                Gvzc = Gvz[cut]
                #print 'len(Gmc)', len(Gmc)
                Gvzcave = np.average(Gvz[cut],weights=Gmc)
                #calculate average rotational velocity:
                #theta hat dot v = vtheta
                Grxyc = np.sqrt(Gxc*Gxc+Gyc*Gyc)
                vth = -Gyc/Grxyc*Gvxc+Gxc/Grxyc*Gvyc
                vthave = np.average(vth,weights=Gmc)
                #print 'Gvxc', Gvxc
                Gvxc = Gvxc+Gyc/Grxyc*vthave
                Gvyc = Gvyc-Gxc/Grxyc*vthave
                vzlog = np.log10(np.absolute(Gvzc-Gvzcave))
                logvzsigma = weighted_quantile(vzlog,fraccut,sample_weight=Gmc)
                cutvz = np.log10(vzlog)<logvzsigma
                Gthcz=Gthc[cutvz]
                Grcz = Grc[cutvz]
                Gmcz = Gmc[cutvz]
                Gzcz = Gzc[cutvz]
                Gxcz = Gxc[cutvz]
                Gycz = Gyc[cutvz]
                Gvzcz = Gvzc[cutvz]
                Gvxcz = Gvxc[cutvz]
                Gvycz = Gvyc[cutvz]
                totnum = len(Gmcz)
                print 'totnum', totnum
                vol = np.pi*(radl[nrad+1]*radl[nrad+1]-radl[nrad]*radl[nrad])*cylz
                numden = totnum/vol
                vol15 = 15.0/numden #TK test
                dl = np.power(vol15,1./3.)
                print 'dl', dl
                #print 'int(dr/dl)+1', int(cylz/dl)+1
                numcount=0
                rscl = []
                for nrs in range(int(dr/dl)+2):
                        prsc = 100.0/(int(dr/dl)+1)*nrs
                        if prsc>100.0:
                                prsc =100.0
                        rsc = np.percentile(Grcz,prsc)
                        rscl = np.append(rscl,rsc)
                for nrs in range(int(dr/dl)+1):
                        rcyls = (nrad)*dr+nrs*dl
                        cutrs = (Grcz<rscl[nrs+1]) & (Grcz>rscl[nrs])
                        Gthcrs = Gthcz[cutrs]
                        Gmcrs = Gmcz[cutrs]
                        Gzcrs = Gzcz[cutrs]
                        Gxcrs = Gxcz[cutrs]
                        Gycrs = Gycz[cutrs]
                        Gvzcrs = Gvzcz[cutrs]
                        Gvxcrs = Gvxcz[cutrs]
                        Gvycrs = Gvycz[cutrs]
                        rzscl = []
                        for nzs in range(int(cylz/dl)+2):
                                przsc = 100.0/(int(cylz/dl)+1)*nzs
                                if przsc>100.0:
                                        przsc=100.0
                                rzsc = np.percentile(Gzcrs,przsc)
                                rzscl = np.append(rzscl,rzsc)
                        for nzs in range(int(cylz/dl)+1):
                                cutzs = (Gzcrs<rzscl[nzs+1]) & (Gzcrs>rzscl[nzs])
                                Gthcrzs = Gthcrs[cutzs]
                                Gmcrzs = Gmcrs[cutzs]
                                Gxcrzs = Gxcrs[cutzs]
                                Gycrzs = Gycrs[cutzs]
                                Gzcrzs = Gzcrs[cutzs]
                                Gvxcrzs = Gvxcrs[cutzs]
                                Gvycrzs = Gvycrs[cutzs]
                                Gvzcrzs = Gvzcrs[cutzs]
                                totnphi=2.0*np.pi*rcyls/dl
                                phi = np.arctan2(Gycrzs,Gxcrzs)
                                phi[phi<0] += 2.0*np.pi
                                phil=[]
                                for nphi in range(int(totnphi)+2):
                                        pphi = 100.0*nphi/(int(totnphi)+1)
                                        if pphi>100.0:
                                                pphi=100.0
                                        phineed = np.percentile(phi,pphi)
                                        phil = np.append(phil, phineed)
                                #print 'phil', phil
                                for nphi in range(int(totnphi)+1):
                                        phis = nphi*2.0*np.pi/totnphi
                                        cutphi = (phi>phil[nphi]) & (phi<phil[nphi+1])
                                        Gths = Gthcrzs[cutphi]
                                        Gms = Gmcrzs[cutphi]
                                        Gvxs = Gvxcrzs[cutphi]
                                        Gvys = Gvycrzs[cutphi]
                                        Gvzs = Gvzcrzs[cutphi]
                                        Gvxsrel = Gvxs-np.average(Gvxs,weights=Gms)
                                        Gvysrel = Gvys-np.average(Gvys,weights=Gms)
                                        Gvzsrel = Gvzs-np.average(Gvzs,weights=Gms)
                                        logGv2s = np.log10(Gvxsrel*Gvxsrel+Gvysrel*Gvysrel+Gvzsrel*Gvzsrel)
                                        if len(Gms)<2:
                                                Gmcell = Gms
                                                Gvxcell = Gvxsrel
                                                Gvycell = Gvysrel
                                                Gvzcell = Gvzsrel
                                        else:
                                                logv2sigma = weighted_quantile(logGv2s,fraccut,sample_weight=Gms)
                                                cutv2 = logGv2s<logv2sigma
                                                Gmcell = Gms[cutv2]
                                                Gvxcell = Gvxsrel[cutv2]
                                                Gvycell = Gvysrel[cutv2]
                                                Gvzcell = Gvzsrel[cutv2]
                                        turenergy = np.sum(0.5*Gmcell*(Gvxcell*Gvxcell+Gvycell*Gvycell+Gvzcell*Gvzcell)) #in erg
                                        Gthcelll = np.append(Gthcelll,np.sum(Gths))
                                        turl = np.append(turl,turenergy)
                                        gml = np.append(gml,np.sum(Gmcell))
                                        xcell = np.append(xcell,0.5*(rscl[nrs+1]+rscl[nrs])*np.cos(0.5*(phil[nphi]+phil[nphi+1])))
                                        ycell = np.append(ycell,0.5*(rscl[nrs+1]+rscl[nrs])*np.sin(0.5*(phil[nphi]+phil[nphi+1])))
                                        zcell = np.append(zcell,0.5*(rzscl[nzs+1]+rzscl[nzs]))
        #print 'xcell', xcell
        #print 'ycell', ycell
        #print 'zcell', zcell
        Gxcutz = Gx[cutz]
        Gycutz = Gy[cutz]
        Gzcutz = Gz[cutz]
        if havecr>0:
                cregyl = cregy[cutz]
        therml = GEint[cutz]
        if haveB>0:
                Begyl = Begy[cutz]
        Gmcutz = Gm[cutz]
        Grxycz = np.sqrt(Gxcutz*Gxcutz+Gycutz*Gycutz)
        rxycell = np.sqrt(xcell*xcell+ycell*ycell)
        print 'np.amax(radl)', np.amax(radl)
        print 'np.amax(rscl)', np.amax(rscl)
        print 'np.amax(xcell)', np.amax(xcell)
        print 'np.amax(rxycell)', np.amax(rxycell)
        cutr = (Grxycz<cylr) & (Grxycz>cylrin)
        Gmcutzr = Gmcutz[cutr]
        print 'Gmcutzr', np.sum(Gmcutzr)
        if haveB>0:
                Begyl=Begyl[cutr]
        therml=therml[cutr]
        if havecr>0:
                cregyl=cregyl[cutr]
        Gxpl=Gxcutz[cutr]
        Gypl=Gycutz[cutr]
        Gzpl=Gzcutz[cutr]
        if cosmo==1:
                readtimelist=readtime(firever=2)
                snap2list=readtimelist['snaplist']
                time2list=readtimelist['timelist']
                a2list=readtimelist['alist']
                tnow = np.interp(ascale,a2list,time2list)*1e9
        if cosmo==1:
                timen = tnow
        else:
                timen = float(Nsnap)*0.98*1e6
        return {'turl':turl,'therml':therml,'Begyl':Begyl,'cregyl':cregyl,\
        'cylrin':cylrin,'cylr':cylr,'cylz':cylz,'nbin':nbin,\
        'Gxpl':Gxpl, 'Gypl':Gypl, 'Gzpl':Gzpl, 'xcell':xcell,\
         'ycell':ycell, 'zcell':zcell, 'gminl':Gmcutzr, 'Gthcelll':Gthcelll, 'timen':timen}
