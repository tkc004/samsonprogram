from stdmodandoption import *


def calgfromparlocal(Glist,pos,withinr,spno,i):
    gxl = np.zeros(len(pos)); gyl = np.zeros(len(pos)); gzl = np.zeros(len(pos))
    for G in Glist:
        data = SSF.calgfromparlist(G,pos=pos,withinr=withinr,spno=spno)
        gxl += data['gx']; gyl += data['gy']; gzl += data['gz'];
    data = {}
    data['gx']=gxl; data['gy']=gyl; data['gz']=gzl;
    return [[i,data]]


def gserial(pos,G,S,DM,withinr,spno):
    gxl=np.array([]); gyl=np.array([]); gzl=np.array([]);
    for inpos in pos:
        Gdata = SSF.calgfrompar(G,pos=inpos,withinr=withinr,spno=spno) 
        Ggx = Gdata['gx']; Ggy = Gdata['gy']; Ggz = Gdata['gz']
        Sdata = SSF.calgfrompar(S,pos=inpos,withinr=withinr,spno=spno) 
        Sgx = Sdata['gx']; Sgy = Sdata['gy']; Sgz = Sdata['gz']
        DMdata = SSF.calgfrompar(DM,pos=inpos,withinr=withinr,spno=spno) 
        DMgx = DMdata['gx']; DMgy = DMdata['gy'];DMgz = DMdata['gz']
        gx = Ggx+Sgx+DMgx
        gy = Ggy+Sgy+DMgy
        gz = Ggz+Sgz+DMgz
        gxl=np.append(gxl,gx); gyl=np.append(gyl,gy); gzl=np.append(gzl,gz);
    return {'gx':gxl, 'gy':gyl, 'gz':gzl}


def gparallel(pos,G,S,DM,withinr,spno):
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
    # Step 2: `pool.apply` the `howmany_within_range()
    gxyz = [pool.apply_async(calgfromparlocal, args=([G,S,DM],xyz,withinr,spno,i),\
                             callback=collect_results)  for i, xyz in enumerate(listpos)]
    # Step 3: Don't forget to close
    pool.close()
    pool.join()
    comdict = PT.joindata(results)
    return comdict


def generategfield(subdict):
    startno=subdict['startno']
    Nsnap=subdict['Nsnap']
    snapsep=subdict['snapsep']
    wanted=subdict['wanted']
    dirneed=subdict['dirneed']
    fmeat=subdict['fmeat']
    nogrid=subdict['nogrid']
    maxlength=subdict['maxlength']
    griddir=subdict['griddir']
    parallel=subdict['parallel']

    for runtodo in dirneed:
        for i in range(startno,Nsnap+1,snapsep):
            info=SSF.outdirname(runtodo, i)
            rundir=info['rundir']
            Nsnapstring=info['Nsnapstring']
            havecr=info['havecr']
            haveB=info['haveB']
            commonpath='/home/tkc004/scratch/snipshot/philruns/'
            #SSF.mkdir_p(commonpath+rundir+'/output/withinr200spno100')
            fname=commonpath+rundir+'/output/withinr200spno100/snipshot_'+Nsnapstring+'.hdf5'
            data = RSS.readsnipshot(fname,ptypelist = [0,1,4])

            G = data[0]; DM = data[1]; S = data[2];

            zmax=maxlength/2.0
            xlist = ylist = zlist = np.linspace(-zmax,zmax,num=nogrid)
            xl,yl,zl = np.meshgrid(xlist,ylist,zlist)
            xl = np.ravel(xl); yl = np.ravel(yl); zl = np.ravel(zl);
            pos = zip(xl,yl,zl)
            print 'parallel', parallel
            if parallel==0:
                gdata = gserial(pos,G,S,DM,200,1)
            elif parallel==1:
                gdata = gparallel(pos,G,S,DM,200,1)
           
            gxl = gdata['gx']; gyl = gdata['gy']; gzl = gdata['gz']
            
            print 'gxl', gxl
            print 'len(gxl)', len(gxl)
            Gf = {}
            Gf['xlist']=xlist; Gf['ylist']=ylist; Gf['zlist']=zlist;
            Gf['xl']=xl; Gf['yl']=yl; Gf['zl']=zl;
            Gf['gxl']=gxl; Gf['gyl']=gyl; Gf['gzl']=gzl;
            dirpath=commonpath+rundir+'/deriveddata/'+griddir+'/gfield/'
            SSF.mkdir_p(dirpath)
            fname=dirpath+'/snipshot_'+Nsnapstring+'.hdf5'
            SSF.ssrm(fname)
            fhdf5 = h5py.File(fname, 'w')
            grp0 = fhdf5.create_group("gfield/")
            RSS.write_layer(grp0, Gf)
            fhdf5.close()
            del fhdf5