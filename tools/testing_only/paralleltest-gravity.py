from stdmodandoption import *
import time
import multiprocessing as mp
from paralleltool import *

        
def collect_results(result):
    results.extend(result)

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
    data = {}
    data['gx']=gxl; data['gy']=gyl; data['gz']=gzl;
    return [[i,data]]

maxlength=20
nogrid=10

zmax=maxlength/2.0
xlist = ylist = zlist = np.linspace(-zmax,zmax,num=nogrid)
xl,yl,zl = np.meshgrid(xlist,ylist,zlist)
xl = np.ravel(xl); yl = np.ravel(yl); zl = np.ravel(zl);
gxl=np.array([]); gyl=np.array([]); gzl=np.array([]);
pos=[]
for x,y,z in zip(xl,yl,zl):
    pos.append([x,y,z])

nocpu = 20

lenpos = len(pos)

chunksize = lenpos/nocpu

listpos = chunks(pos, chunksize)


start = time.time()


withinr=200; spno=1
# Step 1: Init multiprocessing.Pool()


pool = mp.Pool(nocpu)
results=[]

# Step 2: `pool.apply` the `howmany_within_range()
gxyz = [pool.apply_async(calgfromparlocal, args=([G,S,DM],xyz,withinr,spno,i), callback=collect_results)  for i, xyz in enumerate(listpos)]

# Step 3: Don't forget to close
pool.close()
pool.join()

end = time.time()
print(end - start)

comdict=joindata(results)

print comdict['gx']
