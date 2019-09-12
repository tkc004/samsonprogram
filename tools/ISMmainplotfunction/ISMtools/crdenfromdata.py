from stdmodandoption import *
import plot_setup as PS
import collections
from intrhogfunc import intrhogfunc, intrhogtorlist



def crdenfromdata(ssdict):
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
    withoutr = ssdict['withoutr']
    nogrid = ssdict['nogrid']
    maxlength=ssdict['maxlength'] #thickness
    usehalfz=ssdict['usehalfz']
    usepredata=ssdict['usepredata']
    griddir=ssdict['griddir']
    cutcold=ssdict['cutcold']
    usekez=ssdict['usekez']
    withoutr=ssdict['withoutr']
    outHI=ssdict['outHI']
    if not (wanted=='realpr' or wanted=='er' or wanted=='pz' or wanted=='pr' or wanted=='vz'\
            or wanted=='vr' or wanted=='vturr' or wanted=='vturz'):
        print 'wrong wanted'
        return None
    titleneed=title
    ptitle=title
    dclabelneed=1
    newlabelneed=1
    print 'runtodo', runtodo
    dr = 1.0 #kpc
    if wanted=='er' or wanted=='realpr' or wanted=='pr' or wanted=='vr' or wanted=='vturr':
        nogrid=int(withinr/dr)
    rhogl=np.zeros(nogrid); intrhogl=np.zeros(nogrid)
    rhol=np.zeros(nogrid);
    pthl=np.zeros(nogrid);pturl=np.zeros(nogrid);
    pthHIl=np.zeros(nogrid);
    pturhotl=np.zeros(nogrid);
    pturcoldl=np.zeros(nogrid);
    pturHIl=np.zeros(nogrid);
    rhohotl=np.zeros(nogrid);
    rhocoldl=np.zeros(nogrid);
    rhoHIl=np.zeros(nogrid);
    pcrl=np.zeros(nogrid);pBl=np.zeros(nogrid);
    pBtl=np.zeros(nogrid);
    pkezl=np.zeros(nogrid); pvzl=np.zeros(nogrid);
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
        if wanted=='pz' or wanted=='vz' or wanted=='vturz':
            vertical=1; horizontal=0;
        if wanted=='er' or wanted=='realpr' or wanted=='pr' or wanted=='vr' or wanted=='vturr':
            vertical=0; horizontal=1;            
        if usepredata==1:
            predata = SSF.readpreexist(runtodo,Nsnap,griddir,cutcold=cutcold,outHI=outHI)
            data = SSF.calrhofrompredata(predata,withinr,maxlength,\
                                usehalfz=usehalfz,cutcold=cutcold,vertical=vertical,\
                              horizontal=horizontal,withoutr=withoutr,dr=dr,outHI=outHI)           
        else:
            data = SSF.calrhogfrompar(runtodo,i,withinr,maxlength,nogrid,\
                                  havecr=havecr,haveB=haveB,usehalfz=usehalfz)
        if wanted=='pz' or wanted=='vz' or wanted=='vturz':
            zlist=data['zlist']; 
        elif wanted=='er' or wanted=='realpr' or wanted=='pr' or wanted=='vr' or wanted=='vturr':
            zlist=data['rlist'];

        rhol += data['rhol'];
        pthl += data['pthl']; pturl += data['pturl'];
        if outHI==1:
            pturhotl += data['pturhotl'];
            pturcoldl += data['pturcoldl'];
            pturHIl += data['pturHIl'];
            pthHIl += data['pthHIl'];
            rhohotl += data['rhohotl'];
            rhocoldl += data['rhocoldl'];
            rhoHIl += data['rhoHIl'];
        gzlist = data['gzlist']; pkezl += data['pkezl'];
        pvzl += data['pvzl'];
        #print 'pkezl', pkezl
        if havecr>0: pcrl += data['pcrl']; 
        if haveB>0: pBl += data['pBl'];
        if haveB>0: pBtl += data['pBtl'];    
        if wanted=='pz' or wanted=='vz' or wanted=='vturz':
            rhog=data['rhog'];            
            rhogl += rhog
            intrhogl += intrhogfunc(rhog,gzlist,zlist)
        if wanted=='pr':
            intrhogl += intrhogtorlist(predata,zlist,10.0)
        nooftimes+=1
    zlistm = (zlist[1:]+zlist[:-1])/2.0
    intrhogl = intrhogl/nooftimes
    rhol = rhol/nooftimes
    pthl = pthl/nooftimes; pturl = pturl/nooftimes;
    pkezl = pkezl/nooftimes; pvzl=pvzl/nooftimes;
    pcrl = pcrl/nooftimes; pBl=pBl/nooftimes;
    pBtl=pBtl/nooftimes;
    vthl = np.sqrt(pthl*GAMMA/rhol)/km_in_cm
    vturl = np.sqrt(pturl/rhol)/km_in_cm
    if outHI==1:
        pturcoldl = pturcoldl/nooftimes;
        pturhotl = pturhotl/nooftimes;
        pturHIl = pturHIl/nooftimes;
        pthHIl = pthHIl/nooftimes;
        rhocoldl = rhocoldl/nooftimes;
        rhohotl = rhohotl/nooftimes;
        rhoHIl = rhoHIl/nooftimes;
        vturhotl = np.sqrt(pturhotl/rhohotl)/km_in_cm
        vturcoldl = np.sqrt(pturcoldl/rhocoldl)/km_in_cm
        vturHIl = np.sqrt(pturHIl/rhoHIl)/km_in_cm
        vthHIl = np.sqrt(pthHIl*GAMMA/rhoHIl)/km_in_cm
    vkezl = np.sqrt(pkezl/rhol)/km_in_cm
    vvzl = np.sqrt(pvzl/rhol)/km_in_cm
    vcrl = np.sqrt(pcrl*CRgamma/rhol)/km_in_cm
    vBl =  np.sqrt(pBtl*2.0/rhol)/km_in_cm
    ptotl = pthl+pturl
    if usekez==0:
        ptotl += pturl
    else:
        ptotl += pkezl
    if havecr>0:
        ptotl += pcrl
    if haveB>0:
        ptotl += pBl
    vtotl = np.sqrt(ptotl/rhol)/km_in_cm
    if outHI==1:
        print 'vturHIl', vturHIl
        print 'vturhotl', vturhotl
        print 'vturcoldl', vturcoldl
        print 'vthHIl', vthHIl
    linelist = []
    if wanted=='pz' or wanted=='vz' or wanted=='vturz':
        plotdict[wanted]['xlab'] = r'$z\;{\rm[kpc]}$'
    elif wanted=='er' or wanted=='realpr' or wanted=='pr' or wanted=='vr' or wanted=='vturr':
        plotdict[wanted]['xlab'] = r'$r\;{\rm[kpc]}$'
    if wanted=='pz':
        filename=plotloc+'CRplot/crdenz/crpzfrompar_'+wanted+'_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    elif wanted=='pr':
        filename=plotloc+'CRplot/crdenr/crprfrompar_'+wanted+'_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    elif wanted=='realpr':
        filename=plotloc+'CRplot/crdenr/crrealprfrompar_'+wanted+'_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'        
    elif wanted=='er':
        filename=plotloc+'CRplot/crdenr/crerfrompar_'+wanted+'_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'        
    elif wanted=='vz':
        filename=plotloc+'CRplot/crdenz/crvzfrompar_'+wanted+'_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    elif wanted=='vr':
        filename=plotloc+'CRplot/crdenr/crvrfrompar_'+wanted+'_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    elif wanted=='vturr':
        filename=plotloc+'CRplot/crdenr/crvturrfrompar_'+wanted+'_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    elif wanted=='vturz':
        filename=plotloc+'CRplot/crdenz/crvturzfrompar_'+wanted+'_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    if wanted=='pz' or wanted=='pr':
        plotdict[wanted]['ylab'] = r'$\Pi_z {\rm [dyne/cm^2]}$'
    elif wanted=='vz' or wanted=='vr' or wanted=='vturr':
        plotdict[wanted]['ylab'] = r'$\sqrt{\left \langle v_z^2 \right \rangle} {\rm [km/s]}$'
    elif wanted=='er':
        plotdict[wanted]['ylab'] = r'$e {\rm [eV/cm^3]}$'
    elif wanted=='realpr':
        plotdict[wanted]['ylab'] = r'$P {\rm [dyne/cm^2]}$'
    plotdict[wanted]['ptitle'] = ptitle
    plotdict[wanted]['runtitle'] = runlabel
    plotdict[wanted]['xnl']['rhog'] = zlist[:-1]
    plotdict[wanted]['ynl']['rhog'] = np.absolute(intrhogl[:-1])
    plotdict[wanted]['lw']['rhog'] = 2    
    plotdict[wanted]['marker']['rhog'] = ''    
    plotdict[wanted]['lsn']['rhog'] = 'dashed'
    plotdict[wanted]['linelab']['rhog'] = 'gravity'
    plotdict[wanted]['color']['rhog'] = '0.5'
    if wanted=='pz' or wanted=='pr':
        linelist.append('rhog')
    plotdict[wanted]['xnl']['eth'] = zlist[:-1]
    if wanted=='pz' or wanted=='pr' or wanted=='realpr':
        plotdict[wanted]['ynl']['eth'] = pthl[:-1] 
    elif wanted=='vz' or wanted=='vr':
        plotdict[wanted]['ynl']['eth'] = vthl[:-1] 
    elif wanted=='er':
        plotdict[wanted]['ynl']['eth'] = pthl[:-1]*erg_in_eV/(GAMMA-1)
    plotdict[wanted]['linelab']['eth'] = 'thermal'
    plotdict[wanted]['lw']['eth'] = 1    
    plotdict[wanted]['marker']['eth'] = ''
    plotdict[wanted]['lsn']['eth'] = 'dashed'
    plotdict[wanted]['color']['eth'] = 'r'
    if wanted=='pz' or wanted=='pr' or wanted=='realpr' or wanted=='vz' or wanted=='vr' or wanted=='er':
        linelist.append('eth')
    plotdict[wanted]['xnl']['etur'] = zlist[:-1]
    if wanted=='pz' or wanted=='pr' or wanted=='realpr':
        plotdict[wanted]['ynl']['etur'] = pturl[:-1]
    elif wanted=='vz' or wanted=='vr' or wanted=='vturr':
        plotdict[wanted]['ynl']['etur'] = vturl[:-1]
    elif wanted=='er':
        plotdict[wanted]['ynl']['etur'] = pturl[:-1]*erg_in_eV/2.0 #kinetic energy density
    plotdict[wanted]['linelab']['etur'] = 'turbulence'
    plotdict[wanted]['lw']['etur'] = 2    
    plotdict[wanted]['marker']['etur'] = ''    
    plotdict[wanted]['lsn']['etur'] = 'dashed'
    plotdict[wanted]['color']['etur'] = 'g'
    if wanted=='pz' or wanted=='realpr' or wanted=='pr' or wanted=='vturr' or wanted=='er':
        linelist.append('etur')
    if outHI==1:
        if wanted=='vturr' or wanted=='vturr':
            plotdict[wanted]['xnl']['eturh'] = zlist[:-1]
            plotdict[wanted]['ynl']['eturh'] = vturhotl[:-1]
            plotdict[wanted]['linelab']['eturh'] = 'hot turbulence'
            plotdict[wanted]['lw']['eturh'] = 2    
            plotdict[wanted]['marker']['eturh'] = ''    
            plotdict[wanted]['lsn']['eturh'] = 'dotted'
            plotdict[wanted]['color']['eturh'] = 'r'
            linelist.append('eturh')
            plotdict[wanted]['xnl']['eturc'] = zlist[:-1]
            plotdict[wanted]['ynl']['eturc'] = vturcoldl[:-1]
            plotdict[wanted]['linelab']['eturc'] = 'cold turbulence'
            plotdict[wanted]['lw']['eturc'] = 2    
            plotdict[wanted]['marker']['eturc'] = ''    
            plotdict[wanted]['lsn']['eturc'] = 'solid'
            plotdict[wanted]['color']['eturc'] = 'c'
            linelist.append('eturc')
            plotdict[wanted]['xnl']['eturHI'] = zlist[:-1]
            plotdict[wanted]['ynl']['eturHI'] = vturHIl[:-1]
            plotdict[wanted]['linelab']['eturHI'] = 'HI turbulence'
            plotdict[wanted]['lw']['eturHI'] = 2    
            plotdict[wanted]['marker']['eturHI'] = ''    
            plotdict[wanted]['lsn']['eturHI'] = 'dashdot'
            plotdict[wanted]['color']['eturHI'] = 'b'
            linelist.append('eturHI')
            plotdict[wanted]['xnl']['ethHI'] = zlist[:-1]
            plotdict[wanted]['ynl']['ethHI'] = vthHIl[:-1]
            plotdict[wanted]['linelab']['ethHI'] = 'HI thermal'
            plotdict[wanted]['lw']['ethHI'] = 2    
            plotdict[wanted]['marker']['ethHI'] = ''    
            plotdict[wanted]['lsn']['ethHI'] = 'dashdot'
            plotdict[wanted]['color']['ethHI'] = 'm'
            linelist.append('ethHI')
            plotdict[wanted]['xnl']['etotHI'] = zlist[:-1]
            plotdict[wanted]['ynl']['etotHI'] = np.sqrt(np.square(vthHIl[:-1])+np.square(vturHIl[:-1]))
            plotdict[wanted]['linelab']['etotHI'] = 'HI total'
            plotdict[wanted]['lw']['etotHI'] = 1    
            plotdict[wanted]['marker']['etotHI'] = ''    
            plotdict[wanted]['lsn']['etotHI'] = 'solid'
            plotdict[wanted]['color']['etotHI'] = 'k'
            linelist.append('etotHI')
    if usekez==1:
        plotdict[wanted]['xnl']['ekez'] = zlist[:-1]
        if wanted=='pz' or wanted=='pr' or wanted=='realpr':
            plotdict[wanted]['ynl']['ekez'] = pkezl[:-1]
        elif wanted=='vz' or wanted=='vr':
            plotdict[wanted]['ynl']['ekez'] = vkezl[:-1]
        elif wanted=='er':
            plotdict[wanted]['ynl']['ekez'] = pkezl[:-1]*erg_in_eV/2.0
        plotdict[wanted]['linelab']['ekez'] = 'kinetic'
        plotdict[wanted]['lw']['ekez'] = 1    
        plotdict[wanted]['marker']['ekez'] = ''    
        plotdict[wanted]['lsn']['ekez'] = 'solid'
        plotdict[wanted]['color']['ekez'] = 'g'
        if wanted=='pz' or wanted=='pr' or wanted=='vz' or wanted=='vr' or wanted=='er':
            linelist.append('ekez')
        plotdict[wanted]['xnl']['evz'] = zlist[:-1]
        if wanted=='pz' or wanted=='pr':
            plotdict[wanted]['ynl']['evz'] = pvzl[:-1]
        elif wanted=='vz' or wanted=='vr':
            plotdict[wanted]['ynl']['evz'] = vvzl[:-1]
        plotdict[wanted]['linelab']['evz'] = 'vz'
        plotdict[wanted]['lw']['evz'] = 1    
        plotdict[wanted]['marker']['evz'] = ''    
        plotdict[wanted]['lsn']['evz'] = 'dashed'
        plotdict[wanted]['color']['evz'] = 'c'
        if wanted=='pz':
            linelist.append('evz')
    if haveB>0:
        plotdict[wanted]['xnl']['eB'] = zlist[:-1]
        if wanted=='pz' or wanted=='pr':
            plotdict[wanted]['ynl']['eB'] = pBl[:-1]
        elif wanted=='vz' or wanted=='vr':
            plotdict[wanted]['ynl']['eB'] = vBl[:-1]
        elif wanted=='er':
            plotdict[wanted]['ynl']['eB'] = pBtl[:-1]*erg_in_eV
        elif wanted=='realpr':
            plotdict[wanted]['ynl']['eB'] = pBtl[:-1]
        plotdict[wanted]['linelab']['eB'] = 'Bfield'
        plotdict[wanted]['lw']['eB'] = 1
        plotdict[wanted]['lsn']['eB'] = 'solid'        
        plotdict[wanted]['marker']['eB'] = ''
        plotdict[wanted]['color']['eB'] = 'b'
        if wanted=='pz' or wanted=='pr' or wanted=='realpr' or wanted=='er' or wanted=='vz' or wanted=='vr':
            linelist.append('eB')
    if havecr>0:
        plotdict[wanted]['xnl']['ecr'] = zlist[:-1]
        if wanted=='pz' or wanted=='pr' or wanted=='realpr':
            plotdict[wanted]['ynl']['ecr'] = pcrl[:-1]
        elif wanted=='vz' or wanted=='vr':
            plotdict[wanted]['ynl']['ecr'] = vcrl[:-1]
        elif wanted=='er':
            plotdict[wanted]['ynl']['ecr'] = pcrl[:-1]*erg_in_eV/(CRgamma-1.0)
        plotdict[wanted]['lw']['ecr'] = 1  
        plotdict[wanted]['lsn']['ecr'] = 'dashdot' 
        plotdict[wanted]['linelab']['ecr'] = 'CR' 
        plotdict[wanted]['marker']['ecr'] = ''
        plotdict[wanted]['color']['ecr'] = 'y'
        if wanted=='pz' or wanted=='pr' or wanted=='realpr' or wanted=='vz' or wanted=='vr' or wanted=='er':
            linelist.append('ecr')
    plotdict[wanted]['xnl']['etot'] = zlist[:-1]
    if wanted=='pz' or wanted=='pr':
        plotdict[wanted]['ynl']['etot'] = ptotl[:-1]
    elif wanted=='vz' or wanted=='vr':
        plotdict[wanted]['ynl']['etot'] = vtotl[:-1]
    plotdict[wanted]['linelab']['etot'] = 'total'
    plotdict[wanted]['lw']['etot'] = 1    
    plotdict[wanted]['marker']['etot'] = ''    
    plotdict[wanted]['lsn']['etot'] = 'solid'
    plotdict[wanted]['color']['etot'] = 'k'
    if wanted=='pz' or wanted=='pr':
        linelist.append('etot')
    plotdict[wanted]['linelist']=linelist
    plotdict[wanted]['filename'] = filename        
    return plotdict
    
    
    