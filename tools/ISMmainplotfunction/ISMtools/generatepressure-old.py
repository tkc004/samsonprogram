from stdmodandoption import *



def pressureXYZlocal(G, xyz, dx, dy, dz,havecr,haveB,cutcold,i):
    #print 'in local', xyz
    #print 'i, data', i, data
    data = CRTF.pressureXYZ(G, xyz, dx, dy, dz,havecr=havecr,haveB=haveB,cutcold=cutcold)
    return [[i,data]]



def pparallel(pos,G, dx, dy, dz,havecr,haveB,cutcold):
    import multiprocessing as mp
    import paralleltool as PT
    results=[]
    def collect_results(result):
        results.extend(result)
    nocpu =  2
    lenpos = len(pos)
    chunksize = lenpos/nocpu
    listpos = PT.chunks(pos, chunksize)
    pool = mp.Pool(nocpu)
    #print 'pos', pos
    #print 'listpos[0]', listpos[0]
    #print 'len(pos)', len(pos)
    # Step 2: `pool.apply` the `howmany_within_range()
    #for i, xyz in enumerate(listpos):
        #print 'i', i
        #pressureXYZlocal(G, xyz, dx, dy, dz,havecr,haveB,cutcold,i)
        #resultp = pool.apply_async(pressureXYZlocal, args=(G, xyz, dx, dy, dz,havecr,haveB,cutcold,i),\
        #                 callback=collect_results)
        #print resultp.get()
    pxyz = [pool.apply_async(pressureXYZlocal, args=(G, xyz, dx, dy, dz,havecr,haveB,cutcold,i),\
                         callback=collect_results)  for i, xyz in enumerate(listpos)]
    # Step 3: Don't forget to close
    pool.close()
    pool.join()
    #print 'pxyz', pxyz
    #for p in pxyz:
    #    print 'p.get()', p.get()
    print 'results', results
    comdict = PT.joindata(results)
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
            fname=commonpath+rundir+'/output/withinr200GS/snipshot_'+Nsnapstring+'.hdf5'
            data = RSS.readsnipshot(fname,ptypelist = [0,4])
            G = data[0]; S = data[1];
            zmax=maxlength/2.0
            xlist = ylist = zlist = np.linspace(-zmax,zmax,num=nogrid)
            xl,yl,zl = np.meshgrid(xlist,ylist,zlist)
            xl = np.ravel(xl); yl = np.ravel(yl); zl = np.ravel(zl);
            dx=dy=dz=1
            #dx=np.absolute(zlist[1]-zlist[0])
            #dy=np.absolute(zlist[1]-zlist[0])
            #dz=np.absolute(zlist[1]-zlist[0])
            #pos = [[x,y,z] for x,y,z in zip(xl,yl,zl)]
            pos=[]
            for x,y,z in zip(xl,yl,zl):
                pos.append([x,y,z])
            
            if parallel==0:
                dendata = CRTF.pressureXYZ(G, pos, dx, dy, dz,havecr=havecr,haveB=haveB,cutcold=cutcold)
            elif parallel==1:
                dendata = pparallel(pos, G, dx, dy, dz,havecr,haveB,cutcold)
            
            rhol = dendata['rhol']; pthl = dendata['pthl'];
            print 'rhol', rhol
            pturl = dendata['pturl']; pcrl = dendata['pcrl'];
            pBl = dendata['pBl']; voll = dendata['voll'];
            vzavel = dendata['vzavel']; kezl = dendata['kezl'];
            rhol[~np.isfinite(rhol)] = 0; pthl[~np.isfinite(pthl)] = 0;
            pturl[~np.isfinite(pturl)] = 0; pcrl[~np.isfinite(pcrl)] = 0;
            pBl[~np.isfinite(pBl)] = 0; rhol[~np.isfinite(rhol)] = 0;
            vzavel[~np.isfinite(vzavel)] = 0;
            Gpre = {}
            Gpre['xlist']=xlist; Gpre['ylist']=ylist; Gpre['zlist']=zlist;
            Gpre['xl']=xl; Gpre['yl']=yl; Gpre['zl']=zl; 
            Gpre['rhol']=rhol; Gpre['pthl']=pthl;
            Gpre['pturl']=pturl; Gpre['voll']=voll;
            Gpre['pcrl']=pcrl; Gpre['pBl']=pBl;
            Gpre['vzavel']=vzavel; Gpre['kezl']=kezl;
            dirpath=commonpath+rundir+'/deriveddata/'+griddir+'/pressure/'
            SSF.mkdir_p(dirpath)
            fname=dirpath+'/snipshot_'+Nsnapstring+'.hdf5'
            SSF.ssrm(fname)
            fhdf5 = h5py.File(fname, 'w')
            grp0 = fhdf5.create_group("pressure/")
            RSS.write_layer(grp0, Gpre)
            fhdf5.close()
            del fhdf5