from stdmodandoption import *



def outputsnipshot(subdict):
    startno=subdict['startno']
    Nsnap=subdict['Nsnap']
    snapsep=subdict['snapsep']
    wanted=subdict['wanted']
    dirneed=subdict['dirneed']
    fmeat=subdict['fmeat']


    for runtodo in dirneed:
        for i in range(startno,Nsnap+1,snapsep):
            info=SSF.outdirname(runtodo, i)
            rundir=info['rundir']
            Nsnapstring=info['Nsnapstring']
            commonpath='/home/tkc004/scratch/snipshot/philruns/'
            SSF.mkdir_p(commonpath+rundir+'/output/withinr200spno100')
            fname=commonpath+rundir+'/output/withinr200spno100/snipshot_'+Nsnapstring+'.hdf5'
            S=SSF.readsnapfromrun(runtodo,i,4,rotface=1,loccen=1)
            cenin = S['cen']; vcenin = S['vcen']; angLin = S['angL'];
            G = SSF.readsnapfromrun(runtodo,i,0,rotface=1,loccen=1,\
                                   importLcen=1,angLin=angLin,cenin=cenin,vcenin=vcenin)
            DM = SSF.readsnapfromrun(runtodo,i,1,rotface=1,\
                                   importLcen=1,angLin=angLin,cenin=cenin,vcenin=vcenin)
            parlist=[G,DM,S]
            ptypelist=[0,1,4]
            SSF.ssrm(fname)
            RSS.compresswithinrandspno(parlist,ptypelist,fname,withinr=200,spno=100)
            SSF.ssmkdir(commonpath+rundir+'/output/withinr20G')
            fname=commonpath+rundir+'/output/withinr20G/snipshot_'+Nsnapstring+'.hdf5'
            parlist=[G]
            ptypelist=[0]
            SSF.ssrm(fname)
            RSS.compresswithinrandspno(parlist,ptypelist,fname,withinr=20,spno=1)
            SSF.ssmkdir(commonpath+rundir+'/output/withinr20GS')
            fname=commonpath+rundir+'/output/withinr20GS/snipshot_'+Nsnapstring+'.hdf5'
            parlist=[G,S]
            ptypelist=[0,4]
            SSF.ssrm(fname)
            RSS.compresswithinrandspno(parlist,ptypelist,fname,withinr=20,spno=1)            
    return None
            