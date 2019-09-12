from stdmodandoption import *
import plot_setup as PS
import collections




def calgravprefrompar(ssdict):
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
    withinr = ssdict['withinr']
    withoutrnew = ssdict['withoutrnew']
    nogrid = ssdict['nogrid']
    maxlength=ssdict['maxlength'] #thickness
    usehalfz=ssdict['usehalfz']
    usepredata=ssdict['usepredata']
    griddir=ssdict['griddir']
    cutcold=ssdict['cutcold']
    usekez=ssdict['usekez']
    if not (wanted=='dpdz' or wanted=='pz'):
        print 'wrong wanted'
        return None
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
    rhogl=np.zeros(nogrid); intrhogl=np.zeros(nogrid)
    pthl=np.zeros(nogrid);pturl=np.zeros(nogrid);
    pcrl=np.zeros(nogrid);pBl=np.zeros(nogrid);
    pkezl=np.zeros(nogrid); pvzl=np.zeros(nogrid);
    dpthl=np.zeros(nogrid-1);dpturl=np.zeros(nogrid-1);
    dptotl=np.zeros(nogrid-1);dpcrl=np.zeros(nogrid-1);
    dpBl=np.zeros(nogrid-1); dpkezl=np.zeros(nogrid-1);
    nooftimes=0
    for i in range(startno,Nsnap+1,snapsep):
        info=SSF.outdirname(runtodo, i)
        havecr=info['havecr']
        dclabel=info['dclabel']
        haveB=info['haveB']
        newlabel=info['newlabel']
        color=info['color']
        runlabel=info['runlabel']
        if newlabelneed==1:
            labelneed="\n".join(wrap(newlabel,17))
        ptitle=labelneed
        if usepredata==1:
            data = SSF.calrhogfrompreexist(runtodo,Nsnap,withinr,maxlength,\
                                           usehalfz=usehalfz,griddir=griddir,cutcold=cutcold,withoutr=withoutrnew)            
        else:
            data = SSF.calrhogfrompar(runtodo,i,withinr,maxlength,nogrid,\
                                  havecr=havecr,haveB=haveB,usehalfz=usehalfz)
        zlist=data['zlist']; rhog=data['rhog'];
        pthl += data['pthl']; pturl += data['pturl'];
        gzlist = data['gzlist']; pkezl += data['pkezl'];
        pvzl += data['pvzl'];
        #print 'pkezl', pkezl
        if havecr>0: pcrl += data['pcrl']; 
        if haveB>0: pBl += data['pBl'];
        dzlist = (zlist[1:]-zlist[:-1])*kpc_in_cm
        dpthl += np.absolute(pthl[1:]-pthl[:-1])/dzlist;
        dpturl += np.absolute(pturl[1:]-pturl[:-1])/dzlist;
        dpkezl += np.absolute(pkezl[1:]-pkezl[:-1])/dzlist;
        if havecr>0:
            dpcrl += np.absolute(pcrl[1:]-pcrl[:-1])/dzlist;
        if haveB>0:
            dpBl += np.absolute(pBl[1:]-pBl[:-1])/dzlist;
        rhogl += rhog
        dz = np.absolute(zlist[1]-zlist[0])*kpc_in_cm
        poscut =  gzlist>0
        negrhog = rhog[~poscut]*dz
        posrhog = rhog[poscut]*dz
        negintrhog = np.cumsum(negrhog)
        revarr = np.cumsum(posrhog[::-1])
        posintrhog = revarr[::-1]
        intrhog = np.concatenate([negintrhog,posintrhog])
        intrhogl += np.absolute(intrhog)
        nooftimes+=1
    zlistm = (zlist[1:]+zlist[:-1])/2.0
    rhogl = rhogl/nooftimes
    intrhogl = intrhogl/nooftimes
    dpthl = dpthl/nooftimes; dpturl = dpturl/nooftimes;
    pthl = pthl/nooftimes; pturl = pturl/nooftimes;
    pkezl = pkezl/nooftimes; pvzl=pvzl/nooftimes;
    dpcrl = dpcrl/nooftimes; dpBl = dpBl/nooftimes;
    pcrl = pcrl/nooftimes; pBl=pBl/nooftimes;
    dptotl = dpthl
    print 'pthl', pthl
    smoothing=1
    if smoothing==1:
        smoothdict={}
        smoothlist = [rhogl,intrhogl,dpthl,dpturl,pthl,pturl,pkezl,pvzl,pkezl,pvzl,dpcrl,dpBl,pcrl,pBl,dptotl]
        namelist=['rhogl','intrhogl','dpthl','dpturl','pthl','pturl','pkezl','pvzl','pkezl','pvzl','dpcrl','dpBl','pcrl','pBl','dptotl']
        box_pts = 100
        for i, namel in enumerate(namelist):
            smoothdict[namel] = smoothlist[i]
            smoothdict[namel] = SSF.smooth_convolve(smoothdict[namel], box_pts)
    print 'pthl', smoothdict['pthl']
    if usekez==0:
        dptotl += dpturl
    else:
        dptotl += dpkezl        
    ptotl = pthl+pturl
    if usekez==0:
        ptotl += pturl
    else:
        ptotl += pkezl
    if havecr>0:
        dptotl += dpcrl
        ptotl += pcrl
    if haveB>0:
        dptotl += dpBl
        ptotl += pBl

    #print 'zlist', zlist
    #print 'rhog', rhog
    #print 'zip(zlist,rhog)', zip(zlist,rhog)
    
    if wanted=='dpdz':
        plotdict[wanted]['xlab'] = r'$z\;{\rm[kpc]}$'
        plotdict[wanted]['ylab'] = r'${\rm d} P/{\rm d} z\; {\rm [dyne/cm^3]}$'
        plotdict[wanted]['ptitle'] = ptitle
        plotdict[wanted]['runtitle'] = runlabel
        plotdict[wanted]['xnl']['rhog'] = zlist
        plotdict[wanted]['ynl']['rhog'] = np.absolute(rhogl)
        plotdict[wanted]['lw']['rhog'] = 1    
        plotdict[wanted]['marker']['rhog'] = 'None'    
        plotdict[wanted]['lsn']['rhog'] = 'dashed'
        #plotdict[wanted]['linelab']['rhog'] = r'$\rho g$'
        plotdict[wanted]['linelab']['rhog'] = 'gravity'
        plotdict[wanted]['color']['rhog'] = 'k'
        plotdict[wanted]['xnl']['eth'] = zlistm
        plotdict[wanted]['ynl']['eth'] = dpthl    
        #plotdict[wanted]['linelab']['eth'] = r'$p_{\rm th}$'
        plotdict[wanted]['linelab']['eth'] = 'thermal'
        plotdict[wanted]['lw']['eth'] = 1    
        plotdict[wanted]['marker']['eth'] = 'none'
        plotdict[wanted]['lsn']['eth'] = 'dotted'
        plotdict[wanted]['color']['eth'] = 'r'
        plotdict[wanted]['xnl']['etur'] = zlistm
        if usekez==0:
            plotdict[wanted]['ynl']['etur'] = dpturl 
            plotdict[wanted]['linelab']['etur'] = 'turbulence'
        else:
            plotdict[wanted]['ynl']['etur'] = dpkezl 
            plotdict[wanted]['linelab']['etur'] = 'kinetic'
        plotdict[wanted]['lw']['etur'] = 2    
        plotdict[wanted]['marker']['etur'] = 'none'    
        plotdict[wanted]['lsn']['etur'] = 'solid'
        plotdict[wanted]['color']['etur'] = 'g'        
        if haveB>0:
            plotdict[wanted]['xnl']['eB'] = zlistm
            plotdict[wanted]['ynl']['eB'] = dpBl
            plotdict[wanted]['linelab']['eB'] = 'Bfield'
            #plotdict[wanted]['linelab']['eB'] = r'$(B^2-2B_z^2)/(8\pi)$'
            plotdict[wanted]['lw']['eB'] = 1
            plotdict[wanted]['lsn']['eB'] = 'solid'        
            plotdict[wanted]['marker']['eB'] = 'none'
            plotdict[wanted]['color']['eB'] = 'b'            
        if havecr>0:
            plotdict[wanted]['xnl']['ecr'] = zlistm
            plotdict[wanted]['ynl']['ecr'] = dpcrl 
            plotdict[wanted]['lw']['ecr'] = 1  
            plotdict[wanted]['lsn']['ecr'] = 'dashdot' 
            plotdict[wanted]['linelab']['ecr'] = 'CR' 
            plotdict[wanted]['marker']['ecr'] = 'none'
            plotdict[wanted]['color']['ecr'] = 'y'
        plotdict[wanted]['xnl']['etot'] = zlistm
        plotdict[wanted]['ynl']['etot'] = dptotl 
        plotdict[wanted]['linelab']['etot'] = 'total'
        plotdict[wanted]['lw']['etot'] = 1    
        plotdict[wanted]['marker']['etot'] = 'None'    
        plotdict[wanted]['lsn']['etot'] = 'solid'
        plotdict[wanted]['color']['etot'] = '0.5'   
        filename=plotloc+'CRplot/crdenz/rhogdpzfrompar_'+wanted+'_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        plotdict[wanted]['filename'] = filename
    if wanted=='pz':
        delpthl = np.absolute(pthl-pthl[-1])
        delpturl = np.absolute(pturl-pturl[-1])
        delpvzl = np.absolute(pvzl-pvzl[-1])
        if np.absolute(pkezl[-1])<1e-20 or np.absolute(pkezl[-1])>1e-12:
            delpkezl = np.absolute(pkezl-pkezl[-2])
        else:
            delpkezl = np.absolute(pkezl-pkezl[-1])            
        delptotl = delpthl
        if usekez==0:
            delptotl = delptotl+delpturl
        else:
            delptotl = delptotl+delpkezl
            #print 'usekez'
        if haveB>0: 
            delpBl = np.absolute(pBl-pBl[-1])
            delptotl = delptotl+delpBl
            #print 'haveB'
        if havecr>0:
            if np.absolute(pcrl[-1])<1e-20:
                delpcrl = np.absolute(pcrl-pcrl[-2])
            else:
                delpcrl = np.absolute(pcrl-pcrl[-1])
            delptotl = delptotl+delpcrl
        linelist = []
        plotdict[wanted]['xlab'] = r'$z\;{\rm[kpc]}$'
        plotdict[wanted]['ylab'] = r'$\Delta\Pi_z {\rm [dyne/cm^2]}$'
        plotdict[wanted]['ptitle'] = ptitle
        plotdict[wanted]['runtitle'] = runlabel
        plotdict[wanted]['xnl']['rhog'] = zlist
        plotdict[wanted]['ynl']['rhog'] = np.absolute(intrhogl)
        plotdict[wanted]['lw']['rhog'] = 2    
        plotdict[wanted]['marker']['rhog'] = ''    
        plotdict[wanted]['lsn']['rhog'] = 'dashed'
        #plotdict[wanted]['linelab']['rhog'] = r'$\rho g$'
        plotdict[wanted]['linelab']['rhog'] = 'gravity'
        plotdict[wanted]['color']['rhog'] = '0.5'
        linelist.append('rhog')
        plotdict[wanted]['xnl']['eth'] = zlist
        plotdict[wanted]['ynl']['eth'] = delpthl    
        #plotdict[wanted]['linelab']['eth'] = r'$p_{\rm th}$'
        plotdict[wanted]['linelab']['eth'] = 'thermal'
        plotdict[wanted]['lw']['eth'] = 1    
        plotdict[wanted]['marker']['eth'] = ''
        plotdict[wanted]['lsn']['eth'] = 'dashed'
        plotdict[wanted]['color']['eth'] = 'r'
        linelist.append('eth')
        plotdict[wanted]['xnl']['etur'] = zlist
        plotdict[wanted]['ynl']['etur'] = delpturl
        plotdict[wanted]['linelab']['etur'] = 'turbulence'
        plotdict[wanted]['lw']['etur'] = 1    
        plotdict[wanted]['marker']['etur'] = ''    
        plotdict[wanted]['lsn']['etur'] = 'dashed'
        plotdict[wanted]['color']['etur'] = 'g'
        linelist.append('etur')
        if usekez==1:
            plotdict[wanted]['xnl']['ekez'] = zlist
            plotdict[wanted]['ynl']['ekez'] = delpkezl 
            plotdict[wanted]['linelab']['ekez'] = 'kinetic'
            plotdict[wanted]['lw']['ekez'] = 2    
            plotdict[wanted]['marker']['ekez'] = ''    
            plotdict[wanted]['lsn']['ekez'] = 'dotted'
            plotdict[wanted]['color']['ekez'] = 'g'
            linelist.append('ekez')
            plotdict[wanted]['xnl']['evz'] = zlist
            plotdict[wanted]['ynl']['evz'] = delpvzl 
            plotdict[wanted]['linelab']['evz'] = 'vz'
            plotdict[wanted]['lw']['evz'] = 1    
            plotdict[wanted]['marker']['evz'] = ''    
            plotdict[wanted]['lsn']['evz'] = 'dashed'
            plotdict[wanted]['color']['evz'] = 'c'
            #linelist.append('evz')
        if haveB>0:
            plotdict[wanted]['xnl']['eB'] = zlist
            plotdict[wanted]['ynl']['eB'] = delpBl
            plotdict[wanted]['linelab']['eB'] = 'Bfield'
            #plotdict[wanted]['linelab']['eB'] = r'$(B^2-2B_z^2)/(8\pi)$'
            plotdict[wanted]['lw']['eB'] = 1
            plotdict[wanted]['lsn']['eB'] = 'solid'        
            plotdict[wanted]['marker']['eB'] = ''
            plotdict[wanted]['color']['eB'] = 'b'
            linelist.append('eB')
        if havecr>0:
            plotdict[wanted]['xnl']['ecr'] = zlist
            plotdict[wanted]['ynl']['ecr'] = delpcrl
            plotdict[wanted]['lw']['ecr'] = 1  
            plotdict[wanted]['lsn']['ecr'] = 'dashdot' 
            plotdict[wanted]['linelab']['ecr'] = 'CR' 
            plotdict[wanted]['marker']['ecr'] = ''
            plotdict[wanted]['color']['ecr'] = 'y'
            linelist.append('ecr')
        plotdict[wanted]['xnl']['etot'] = zlist
        plotdict[wanted]['ynl']['etot'] = delptotl
        plotdict[wanted]['linelab']['etot'] = 'total'
        plotdict[wanted]['lw']['etot'] = 1    
        plotdict[wanted]['marker']['etot'] = ''    
        plotdict[wanted]['lsn']['etot'] = 'solid'
        plotdict[wanted]['color']['etot'] = 'k'
        linelist.append('etot')
        plotdict[wanted]['linelist']=linelist
        filename=plotloc+'CRplot/crdenz/rhogpzfrompar_'+wanted+'_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        plotdict[wanted]['filename'] = filename    
    return plotdict
    
    
    