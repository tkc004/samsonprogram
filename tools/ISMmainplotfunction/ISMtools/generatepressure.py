from stdmodandoption import *



def pressureXYZlocal(G, xyz, dx, dy, dz,havecr,haveB,cutcold,outHI,i):
    #print 'in local', xyz
    #print 'i, data', i, data
    data = CRTF.pressureXYZ(G, xyz, dx, dy, dz,havecr=havecr,haveB=haveB,cutcold=cutcold,outHI=outHI)
    return [[i,data]]



def pparallel(pos,G, dx, dy, dz,havecr,haveB,cutcold,outHI):
    import multiprocessing as mp
    import paralleltool as PT
    def collect_results(result):
        results.extend(result)
    nocpu =  16
    lenpos = len(pos)
    chunksize = lenpos/nocpu
    listpos = PT.chunks(pos, chunksize)
    pool = mp.Pool(nocpu)
    results=[]
    #print 'pos', pos
    #print 'listpos[0]', listpos[0]
    #print 'len(pos)', len(pos)
    lno = 0
    #for l in listpos:
    #   print 'len(l)', len(l)
    #    lno += len(l)
    #print 'lno', lno    
    # Step 2: `pool.apply` the `howmany_within_range()
    #for i, xyz in enumerate(listpos):
        #print 'i', i
        #pressureXYZlocal(G, xyz, dx, dy, dz,havecr,haveB,cutcold,i)
        #resultp = pool.apply_async(pressureXYZlocal, args=(G, xyz, dx, dy, dz,havecr,haveB,cutcold,i),\
        #                 callback=collect_results)
        #print resultp.get()
    pxyz = [pool.apply_async(pressureXYZlocal, args=(G, xyz, dx, dy, dz,havecr,haveB,cutcold,outHI,i),\
                         callback=collect_results)  for i, xyz in enumerate(listpos)]
    # Step 3: Don't forget to close
    pool.close()
    pool.join()
    #for item in results:
    #    print 'item[0]', item[0]
    comdict = PT.joindata(results)
    for key in comdict:
        print 'key', key
    return comdict



def generatepressure(subdict):
    startno=subdict['startno']
    Nsnap=subdict['Nsnap']
    snapsep=subdict['snapsep']
    wanted=subdict['wanted']
    dirneed=subdict['dirneed']
    fmeat=subdict['fmeat']
    nogrid=subdict['nogrid']
    maxlength=subdict['maxlength']
    griddir=subdict['griddir']
    cutcold=subdict['cutcold']
    parallel=subdict['parallel']
    outHI=subdict['outHI']
    print 'snapsep', snapsep
    print 'parallel', parallel

    for runtodo in dirneed:
        for i in range(startno,Nsnap+1,snapsep):
            info=SSF.outdirname(runtodo, i)
            rundir=info['rundir']
            Nsnapstring=info['Nsnapstring']
            havecr=info['havecr']
            haveB=info['haveB']
            commonpath='/home/tkc004/scratch/snipshot/philruns/'
            #SSF.mkdir_p(commonpath+rundir+'/output/withinr200spno100')
            fname=commonpath+rundir+'/output/withinr20G/snipshot_'+Nsnapstring+'.hdf5'
            data = RSS.readsnipshot(fname,ptypelist = [0])
            G = data[0];
            zmax=maxlength/2.0
            xlist = ylist = zlist = np.linspace(-zmax,zmax,num=nogrid)
            xl,yl,zl = np.meshgrid(xlist,ylist,zlist)
            xl = np.ravel(xl); yl = np.ravel(yl); zl = np.ravel(zl);
            #dx=dy=dz=1
            dx=np.absolute(zlist[1]-zlist[0])
            dy=np.absolute(zlist[1]-zlist[0])
            dz=np.absolute(zlist[1]-zlist[0])
            #pos = [[x,y,z] for x,y,z in zip(xl,yl,zl)]
            pos=[]
            for x,y,z in zip(xl,yl,zl):
                pos.append([x,y,z])
            
            if parallel==0:
                dendata = CRTF.pressureXYZ(G, pos, dx, dy, dz,havecr=havecr,haveB=haveB,cutcold=cutcold,outHI=outHI)
            elif parallel==1:
                dendata = pparallel(pos, G, dx, dy, dz,havecr,haveB,cutcold,outHI)          
            rhol = dendata['rhol']; pthl = dendata['pthl'];
            #print 'len(rhol)', len(rhol)
            if outHI==1:
                rhoHIl = dendata['rhoHIl'];
                pturHIl = dendata['pturHIl'];
                pthHIl = dendata['pthHIl'];
                rhohotl = dendata['rhohotl'];
                rhocoldl = dendata['rhocoldl'];
            pturl = dendata['pturl']; pcrl = dendata['pcrl'];
            pturhotl = dendata['pturhotl'];
            pturcoldl = dendata['pturcoldl'];
            pBl = dendata['pBl']; voll = dendata['voll'];
            vzavel = dendata['vzavel']; kezl = dendata['kezl'];
            pBtl = dendata['pBtl']; 
            print 'rhol in gp', rhol
            #print 'pthHIl', pthHIl
            #print 'pthl', pthl
            rhol[~np.isfinite(rhol)] = 0; pthl[~np.isfinite(pthl)] = 0;
            if outHI==1:
                rhoHIl[~np.isfinite(rhoHIl)] = 0;
                pturHIl[~np.isfinite(pturHIl)] = 0;
                pthHIl[~np.isfinite(pthHIl)] = 0;
                rhohotl[~np.isfinite(rhohotl)] = 0;
                rhocoldl[~np.isfinite(rhocoldl)] = 0;
            pturl[~np.isfinite(pturl)] = 0; pcrl[~np.isfinite(pcrl)] = 0;
            pturhotl[~np.isfinite(pturhotl)] = 0;
            pturcoldl[~np.isfinite(pturcoldl)] = 0;            
            pBl[~np.isfinite(pBl)] = 0; rhol[~np.isfinite(rhol)] = 0;
            vzavel[~np.isfinite(vzavel)] = 0; pBtl[~np.isfinite(pBtl)] = 0;
            Gpre = {}
            Gpre['xlist']=xlist; Gpre['ylist']=ylist; Gpre['zlist']=zlist;
            Gpre['xl']=xl; Gpre['yl']=yl; Gpre['zl']=zl; 
            Gpre['rhol']=rhol; Gpre['pthl']=pthl;
            #print 'len(rhol)', len(rhol)
            Gpre['rhohotl']=rhohotl;
            Gpre['rhocoldl']=rhocoldl;
            Gpre['rhoHIl']=rhoHIl;
            Gpre['pturl']=pturl; Gpre['voll']=voll;
            Gpre['pturhotl']=pturhotl;
            Gpre['pturcoldl']=pturcoldl;
            Gpre['pturHIl']=pturHIl;
            Gpre['pcrl']=pcrl; Gpre['pBl']=pBl;
            Gpre['pBtl']=pBtl;
            Gpre['vzavel']=vzavel; Gpre['kezl']=kezl;
            Gpre['pthHIl']=pthHIl;
            #print 'pturl', pturl
            dirpath=commonpath+rundir+'/deriveddata/'+griddir+'/pressure/'
            SSF.mkdir_p(dirpath)
            fname=dirpath+'/snipshot_'+Nsnapstring+'.hdf5'
            SSF.ssrm(fname)
            fhdf5 = h5py.File(fname, 'w')
            grp0 = fhdf5.create_group("pressure/")
            RSS.write_layer(grp0, Gpre)
            fhdf5.close()
            del fhdf5