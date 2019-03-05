from stdmodandoption import *
import collections

def calreff(S,rlist=np.linspace(0.1,50,num=50)):        
        mlist=np.array([])
        Spos = S['p'][:,:]
        Smass = S['m'][:]
        Sx = Spos[:,0]
        Sy = Spos[:,1]
        Sz = Spos[:,2]
        Sr = np.sqrt(np.square(Sx)+np.square(Sy)+np.square(Sz))
        for i in range(len(rlist)):
            withinrlist = Sr<rlist[i]
            withinrm = np.sum(Smass[withinrlist]*1e10)
            mlist = np.append(mlist,withinrm)
        r12 = np.interp(mlist[-1]*0.5,mlist,rlist)
        reff = r12*3.0/4.0
        return {'r12':r12, 'reff':reff, 'mlist':mlist, 'rlist':rlist} # mlist in solar sun and rlist in kpc


def calveldis(S,r12=-1.,withinkpc=20.):
        rlist=np.linspace(0.1,20,num=200)
        mlist=[]
        Spos = S['p'][:,:]
        Smass = S['m'][:]
        Svel = S['v'][:,:]
        Sx = Spos[:,0]
        Sy = Spos[:,1]
        Sz = Spos[:,2]
        Svx=Svel[:,0]
        Svy=Svel[:,1]
        Svz=Svel[:,2]
        Sr = np.sqrt(np.square(Sx)+np.square(Sy)+np.square(Sz))
        Srx = np.sqrt(np.square(Sy)+np.square(Sz))
        Sry = np.sqrt(np.square(Sx)+np.square(Sz))
        Srz = np.sqrt(np.square(Sx)+np.square(Sy))
        if r12<0:
            for i in range(len(rlist)):
                withinrlist = Sr<rlist[i]
                withinrm = np.sum(Smass[withinrlist]*1e10)
                mlist = np.append(mlist,withinrm)
            r12 = np.interp(mlist[-1]*0.5,mlist,rlist)
        reff = r12*3.0/4.0
        withinr = Sr<r12
        withinx = Srx<withinkpc
        withiny = Sry<withinkpc
        withinz = Srz<withinkpc
        Swithinm = np.sum(Smass[withinr]*1e10)
        xdis=np.std(Svx[withinx])
        ydis=np.std(Svy[withiny])
        zdis=np.std(Svz[withinz])
        return {'withinm':Swithinm, 'xdis':xdis, 'ydis':ydis, 'zdis': zdis, 'r12':r12, 'reff':reff}
    
def stellarveldis(runtodo,Nsnap):
        gasveldisneed=0
        if (int(Nsnap) < 10):
                Nsnapstring = '00'+str(Nsnap)
        elif (int(Nsnap) < 100):
                Nsnapstring = '0'+str(Nsnap)
        else:
                Nsnapstring = str(Nsnap)
        haloinfo=SSF.outdirname(runtodo)
        rundir=haloinfo['rundir']
        subdir=haloinfo['subdir']
        maindir=haloinfo['maindir']
        multifile=haloinfo['multifile']
        halostr=haloinfo['halostr']
        cosmo=haloinfo['cosmo']
        the_snapdir = '/home/tkc004/'+maindir+'/'+rundir+'/'+subdir
        the_prefix ='snapshot'
        the_suffix = '.hdf5'
        rotface=0;
        if cosmo==1:
            datasup=0;
        else:
            datasup=1;
        S = SSF.readsnapfromrun(runtodo,Nsnap,4,rotface=rotface,datasup=datasup)
        dataS = calveldis(S,r12=-1.)
        Swithinm = dataS['withinm']; xdis = dataS['xdis']; ydis = dataS['ydis']; zdis = dataS['zdis'];
        r12 = dataS['r12']; reff = dataS['reff']
        print 'x dispersion', xdis
        print 'y dispersion', ydis
        print 'z dispersion', zdis
        del S

        G = SSF.readsnapfromrun(runtodo,Nsnap,0,rotface=rotface,datasup=datasup)
        dataG = calveldis(G,r12=r12)
        xgdis = dataG['xdis']; ygdis = dataG['ydis']; zgdis = dataG['zdis']; Gwithinm = dataG['withinm']
        print 'xgdis', xgdis
        print 'ygdis', ygdis
        print 'zgdis', zgdis
        del G
        
        DM = SSF.readsnapfromrun(runtodo,Nsnap,1,rotface=rotface,datasup=datasup)
        dataDM = calveldis(DM,r12=r12)
        DMwithinm = dataDM['withinm']
        masswithin = Gwithinm+Swithinm+DMwithinm
        del DM
        print 'masswithin', masswithin
        print 'log(masswithin)', np.log10(masswithin)
        print 'expected dispersion', np.sqrt(masswithin/9.3e5/reff)
        return reff, xdis, ydis, zdis, masswithin, Gwithinm
    
