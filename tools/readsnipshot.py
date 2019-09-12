import numpy as np
import h5py

def write_layer(gp, adict,hdf5keyexceptions=[]):
    for k,v in adict.items():
        if k not in hdf5keyexceptions:
            if isinstance(v, dict):
                gp1 = gp.create_group(k)
                write_layer(gp1, v)
            else:
                gp.create_dataset(k, data=np.atleast_1d(v))
    return None                
                
def savesnipshot(parlist,fname,ptypelist = [0]):
    hdf5keyexceptions=['header']
    fhdf5 = h5py.File(fname, 'w')
    for ptype,pdata in zip(ptypelist,parlist):
        bname="PartType"+str(ptype)+"/"
        grp = fhdf5.create_group(bname)
        write_layer(grp, pdata,hdf5keyexceptions=hdf5keyexceptions)
    fhdf5.close()
    del fhdf5
    return None
    
def readsnipshot(fname,ptypelist = [0]):
    fhdf5 = h5py.File(fname, 'r')
    parlist=[]
    bnamelist=[]
    for ptype in ptypelist:
        parlist.append(dict())
        bnamelist.append("PartType"+str(ptype)+"/")
    for i, bname in enumerate(bnamelist):
        fgroup = fhdf5[bname]
        for name in fgroup:
            parlist[i][name] = fgroup[name][:]
    fhdf5.close()
    del fhdf5
    return parlist

def readgeneralhdf5(fname):
    fhdf5 = h5py.File(fname, 'r')
    data = {}
    for gp in fhdf5:
        group = fhdf5[gp]
        for name in group:
            #print 'name', name
            data[name] = group[name][:]
    fhdf5.close()
    del fhdf5
    return data

def Gspcut(G,spno=1):
    Gnew={}
    Gm = G['m']
    Glen = len(Gm)
    if spno>1:
        for key in G:
            Gshape = np.array(G[key]).shape
            if len(Gshape)==0: 
                Gnew[key]=G[key]
            elif Gshape[0]<Glen-1:
                Gnew[key]=G[key]
            else:
                if len(Gshape)==1:
                    Gnew[key]=G[key][::spno]
                else:
                    Glist = []
                    for iy in range(Gshape[1]):
                        Gx = np.array(G[key][:,iy]);
                        Glist.append(Gx[::spno])
                    Gtri = np.vstack((Glist)).T    
                    Gnew[key]=Gtri
        Gnew['m']=spno*Gnew['m']
        return Gnew
    else:
        return G
        


def Gwithinrcut(Gn,withinr=-1):
    Gnew={}
    if withinr>0:
        Gp = Gn['p']; Gx = Gp[:,0]; Gy = Gp[:,1]; Gz = Gp[:,2];
        Glen = len(Gx)
        #consider only particles within withinr
        Gr = np.sqrt(Gx*Gx+Gy*Gy+Gz*Gz);
        cutr = Gr<withinr;
        for key in Gn:
            Gshape = np.array(Gn[key]).shape
            if len(Gshape)==0:
                Gnew[key]=Gn[key]
            elif Gshape[0]<Glen-1:
                Gnew[key]=Gn[key]
            else:
                if len(Gshape)==1:
                    Gnew[key]=Gn[key][cutr]
                else:
                    Glist = []
                    for iy in range(Gshape[1]):
                        Gx = np.array(Gn[key][:,iy]);
                        Glist.append(Gx[cutr])
                    Gtri = np.vstack((Glist)).T    
                    Gnew[key]=Gtri
        return Gnew
    else:
        return Gn

    
def compresswithinr(parlist,ptypelist,fname,withinr=200):
    parGcutlist=[]
    for parG, ptype in zip(parlist,ptypelist):
        parGcut = Gwithinrcut(parG,withinr=withinr)
        parGcutlist.append(parGcut)
    savesnipshot(parGcutlist,fname,ptypelist = ptypelist)
    return None

def compresswithinrandspno(parlist,ptypelist,fname,withinr=200,spno=100):
    parGcutlist=[]
    for parG, ptype in zip(parlist,ptypelist):
        parGmid = Gwithinrcut(parG,withinr=withinr)
        parGcut = Gspcut(parGmid,spno=spno)
        parGcutlist.append(parGcut)
    savesnipshot(parGcutlist,fname,ptypelist = ptypelist)
    return None