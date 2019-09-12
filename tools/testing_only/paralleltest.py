from stdmodandoption import *
import time

runtodo='m12imhdcvhr'
i=580
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

def calgfromparlocal(Glist,pos,withinr,spno,i):
    gxl = np.zeros(len(pos)); gyl = np.zeros(len(pos)); gzl = np.zeros(len(pos))
    for G in Glist:
        data = SSF.calgfromparlist(G,pos=pos,withinr=withinr,spno=spno)
        gxl += data['gx']; gyl += data['gy']; gzl += data['gz'];
    return [[i,gxl,gyl,gzl]]



def chunks(l, n):
    n = max(1, n)
    return [l[i:i+n] for i in xrange(0, len(l), n)]

def collect_results(result):
    results.extend(result)
    
maxlength=20
nogrid=10
withinr=200; spno=1
compareserial=1
compareparallel=1



zmax=maxlength/2.0
xlist = ylist = zlist = np.linspace(-zmax,zmax,num=nogrid)
xl,yl,zl = np.meshgrid(xlist,ylist,zlist)
xl = np.ravel(xl); yl = np.ravel(yl); zl = np.ravel(zl);
gxl=np.array([]); gyl=np.array([]); gzl=np.array([]);

print 'len(xl)', len(xl)


if compareserial==1:
    start = time.time()
    for x,y,z in zip(xl,yl,zl):
        Gdata = SSF.calgfrompar(G,pos=[x,y,z],withinr=withinr,spno=spno) 
        Ggx = Gdata['gx']; Ggy = Gdata['gy']; Ggz = Gdata['gz']
        Sdata = SSF.calgfrompar(S,pos=[x,y,z],withinr=withinr,spno=spno) 
        Sgx = Sdata['gx']; Sgy = Sdata['gy']; Sgz = Sdata['gz']
        DMdata = SSF.calgfrompar(DM,pos=[x,y,z],withinr=withinr,spno=spno) 
        DMgx = DMdata['gx']; DMgy = DMdata['gy'];DMgz = DMdata['gz']
        gx = Ggx+Sgx+DMgx
        gy = Ggy+Sgy+DMgy
        gz = Ggz+Sgz+DMgz
        gxl=np.append(gxl,gx); gyl=np.append(gyl,gy); gzl=np.append(gzl,gz);
    end = time.time()
    print 'series time'
    print(end - start)
    print 'gxl', gxl


if compareparallel==1:        
    import multiprocessing as mp
    pos=[]
    for x,y,z in zip(xl,yl,zl):
        pos.append([x,y,z])

    nocpu = mp.cpu_count()

    lenpos = len(pos)

    chunksize = lenpos/nocpu

    listpos = chunks(pos, chunksize)

    start = time.time()
    # Step 1: Init multiprocessing.Pool()

    pool = mp.Pool(nocpu)
    results=[]

    # Step 2: `pool.apply` the `howmany_within_range()`
    gxyz = [pool.apply_async(calgfromparlocal, args=([G,S,DM],xyz,withinr,spno,i), callback=collect_results)  for i, xyz in enumerate(listpos)]
    #for i, xyz in enumerate(listpos):
    #    pool.apply_async(calgfromparlocal, args=([G,S,DM],xyz,withinr,spno,i), callback=collect_results)
    #gxyz = [pool.apply_async(calgfromparlocal, args=([G,S,DM],xyz,withinr,spno,i),\
    #                         callback=collect_results)  for i, xyz in enumerate(listpos)]
    # Step 3: Don't forget to close
    pool.close()   
    pool.join()
    
    argsortlist=np.argsort([results[i][0] for i in range(len(results))])
    gxlist =np.array([results[i][1] for i in range(len(results))])
    gxlnew = gxlist[argsortlist].flatten()
    gylist =np.array([results[i][2] for i in range(len(results))])
    gylnew = gylist[argsortlist].flatten()
    gzlist =np.array([results[i][3] for i in range(len(results))])
    gzlnew = gzlist[argsortlist].flatten()

    end = time.time()
    print 'parallel time'
    print(end - start)
    print 'gxlnew.ravel()', gxlnew.ravel()
    #print 'gxlist', gxlist
