from stdmodandoption import *
import time
import multiprocessing as mp

runtodo='m12fmhdcvhr'
i=580
info=SSF.outdirname(runtodo, i)
rundir=info['rundir']
Nsnapstring=info['Nsnapstring']
havecr=info['havecr']
haveB=info['haveB']
cutcold=0
dx=dy=dz=1
commonpath='/home/tkc004/scratch/snipshot/philruns/'
#SSF.mkdir_p(commonpath+rundir+'/output/withinr200spno100')
fname=commonpath+rundir+'/output/withinr20G/snipshot_'+Nsnapstring+'.hdf5'
data = RSS.readsnipshot(fname,ptypelist = [0])

G = data[0];

def pressureXYZlocal(G, pos, dx, dy, dz,havecr,haveB,cutcold,i):
    data = CRTF.pressureXYZ(G, pos, dx, dy, dz,havecr=havecr,haveB=haveB,cutcold=cutcold)
    return [[i,data]]

def chunks(l, n):
    n = max(1, n)
    lenl = len(l)
    chucklist = [l[i:i+n] for i in xrange(0, lenl, n)]
    return chucklist
        
def collect_results(result):
    results.extend(result)
    
def joindata(results):
    argsortlist=np.argsort([results[i][0] for i in range(len(results))])
    totdict = np.array([results[i][1] for i in range(len(results))])
    totdict=totdict[argsortlist]
    comdict={}
    for key in totdict[0]:
        totarray=np.array([])
        for indict in totdict:
            totarray = np.append(totarray,indict[key])
        comdict[key] = totarray.flatten()
    return comdict




start = time.time()


withinr=200; spno=1
# Step 1: Init multiprocessing.Pool()

maxlength=20
nogrid=40

zmax=maxlength/2.0
xlist = ylist = zlist = np.linspace(-zmax,zmax,num=nogrid)
xl,yl,zl = np.meshgrid(xlist,ylist,zlist)
xl = np.ravel(xl); yl = np.ravel(yl); zl = np.ravel(zl);

pos=[]
for x,y,z in zip(xl,yl,zl):
    pos.append([x,y,z])

nocpu = 10

lenpos = len(pos)

chunksize = lenpos/nocpu

listpos = chunks(pos, chunksize)



pool = mp.Pool(nocpu)
results=[]
# Step 2: `pool.apply` the `howmany_within_range()
pxyz = [pool.apply_async(pressureXYZlocal, args=(G, xyz, dx, dy, dz,havecr,haveB,cutcold,i),\
                         callback=collect_results)  for i, xyz in enumerate(listpos)]
# Step 3: Don't forget to close
pool.close()
pool.join()

#for p in pxyz:
#    print 'p.get()', p.get()

end = time.time()
print(end - start)

comdict=joindata(results)

#print 'comdict', comdict
