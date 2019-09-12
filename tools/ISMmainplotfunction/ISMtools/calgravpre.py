from stdmodandoption import *
import plot_setup as PS
import collections
from CR_rz_energy_distribution import CR_rz_energy_distribution

def calrhogfromfit(runtodo,Nsnap,withinr,withoutr,maxlength,nogrid):
    rmax=withinr; zmax=maxlength/2.;
    rneed = (withinr+withoutr)/2.
    zlist = np.linspace(0.01,zmax,num=10)
    G=SSF.readsnapfromrun(runtodo,Nsnap,0,rotface=1,loccen=1)
    cenin = G['cen']; vcenin = G['vcen']; angLin = G['angL'];
    fitGdata = SSF.fitMNpot(G, rmax=rmax,zmax=zmax,nogrid=nogrid)
    aG = fitGdata['afit']; bG = fitGdata['bfit']; MG = fitGdata['Mfit'];
    rhoreal = fitGdata['rhoreal']; rlistr = fitGdata['rlist']; zlistr = fitGdata['zlist'];
    rhoreal = rhoreal.reshape((len(rlistr),len(zlistr)))
    rhoatrl = np.array([])
    for irho in range(len(zlistr)):
        rhoatr = np.interp(rneed, rlistr,rhoreal[:,irho])
        rhoatrl = np.append(rhoatrl,rhoatr)
    rhorzl = np.interp(zlist, zlistr,rhoatrl)
    Gdata = SSF.calMNpot(rneed,zlist,MG,aG,bG)
    GdPhidz = Gdata['dPhidz']
    S = SSF.readsnapfromrun(runtodo,Nsnap,4,rotface=1,loccen=1,\
                           importLcen=1,angLin=angLin,cenin=cenin,vcenin=vcenin)
    fitSdata = SSF.fitMNpot(S, rmax=rmax,zmax=zmax,nogrid=nogrid)
    aS = fitSdata['afit']; bS = fitSdata['bfit']; MS = fitSdata['Mfit'];
    Sdata = SSF.calMNpot(rneed,zlist,MS,aS,bS)
    SdPhidz = Sdata['dPhidz'] 
    DM = SSF.readsnapfromrun(runtodo,Nsnap,1,rotface=1,\
                           importLcen=1,angLin=angLin,cenin=cenin,vcenin=vcenin)
    DMdata = SSF.DMtosphpot(DM,cutrout=20., nogrid=20.0)
    DMdPhidr = DMdata['dPhidr']; rsph = DMdata['rsph']
    rsphneed = np.sqrt(rneed*rneed+zlist*zlist)
    DMdPhidrneed = np.interp(rsphneed,rsph,DMdPhidr)
    DMdPhidz = DMdPhidrneed/rsphneed*zlist
    dPhidz=SdPhidz+GdPhidz+DMdPhidz
    rhog = np.absolute(rhorzl*dPhidz)
    return {'zlist':zlist, 'rhog':rhog, 'dPhidz':dPhidz, 'rhorzl':rhorzl }


def calgravpre(ssdict):
    nested_dict = lambda: collections.defaultdict(nested_dict)
    plotdict = nested_dict()
    runtodo=ssdict['runtodo']
    wanted=ssdict['wanted']
    print 'wanted', wanted
    startno=ssdict['startno']
    Nsnap=ssdict['Nsnap']
    snapsep=ssdict['snapsep']
    the_prefix=ssdict['the_prefix']
    the_suffix=ssdict['the_suffix']
    fmeat=ssdict['fmeat']
    titleneed=title
    ptitle=title
    needlog=0
    dclabelneed=1
    correctIa=0
    useM1=0
    resoneed=0
    rotface=1
    newlabelneed=1
    findradiusatnism = 1 #find the radius that has that ISM density (negative to turn off)
    print 'runtodo', runtodo
    info=SSF.outdirname(runtodo, Nsnap)
    havecr=info['havecr']
    dclabel=info['dclabel']
    haveB=info['haveB']
    withinr = ssdict['withinr']
    withoutr = ssdict['withoutr']
    nogrid = ssdict['nogrid']
    maxlength=ssdict['maxlength'] #thickness 
    data = calrhogfromfit(runtodo,Nsnap,withinr,withoutr,maxlength,nogrid)
    zlist=data['zlist']; rhog=data['rhog'];
    plotdict = CR_rz_energy_distribution(ssdict)
    print 'zlist', zlist
    print 'rhog', rhog
    plotdict[wanted]['xnl']['rhog'] = zlist
    plotdict[wanted]['ynl']['rhog'] = rhog
    plotdict[wanted]['lw']['rhog'] = 2    
    plotdict[wanted]['marker']['rhog'] = 'None'    
    plotdict[wanted]['lsn']['rhog'] = 'dashed'
    plotdict[wanted]['linelab']['rhog'] = r'$\rho g$'
    filename=plotloc+'CRplot/crdenz/rhogdpz_'+wanted+'_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    plotdict[wanted]['filename'] = filename
    return plotdict
    
    
    