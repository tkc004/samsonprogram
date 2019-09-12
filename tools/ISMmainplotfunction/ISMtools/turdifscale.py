from stdmodandoption import *
import plot_setup as PS
import collections


def outvtur(runtodo,wanted,startno,Nsnap,snapsep,nogrid,griddir,withinr,maxlength,\
                                usehalfz=0,cutcold=0,vertical=1,\
                              horizontal=0,withoutr=0,dr=1,outHI=0,usepredata=1):
    if vertical==1:
        if griddir=='grid0_125kpc':
            nogrid=160
        if griddir=='grid0_25kpc':
            nogrid=80            
        if griddir=='grid0_5kpc':
            nogrid=40
        if griddir=='grid1kpc':
            nogrid=21            
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
    pkezl=np.zeros(nogrid); pvzl=np.zeros(nogrid);
    nooftimes=0
    for i in range(startno,Nsnap+1,snapsep):
        info=SSF.outdirname(runtodo, i)
        havecr=info['havecr']
        dclabel=info['dclabel']
        haveB=info['haveB']
        color=info['color']
        runlabel=info['runlabel']
        if wanted=='pz' or wanted=='vz' or wanted=='vturz':
            vertical=1; horizontal=0;
        if wanted=='pr' or wanted=='vr' or wanted=='vturr':
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
        elif wanted=='pr' or wanted=='vr' or wanted=='vturr':
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
        nooftimes+=1
    zlistm = (zlist[1:]+zlist[:-1])/2.0
    intrhogl = intrhogl/nooftimes
    rhol = rhol/nooftimes
    pthl = pthl/nooftimes; pturl = pturl/nooftimes;
    pkezl = pkezl/nooftimes; pvzl=pvzl/nooftimes;
    pcrl = pcrl/nooftimes; pBl=pBl/nooftimes; 
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
    vBl =  np.sqrt(pBl*2.0/rhol)/km_in_cm
    ptotl = pthl+pturl
    ptotl += pkezl
    if havecr>0:
        ptotl += pcrl
    if haveB>0:
        ptotl += pBl
    vtotl = np.sqrt(ptotl/rhol)/km_in_cm
    outdict = {'zlist':zlist,'zlistm':zlistm, 'intrhogl':intrhogl, 'rhol':rhol, 'pthl':pthl, 'pturl':pturl,\
               'pkezl':pkezl, 'pvzl':pvzl, 'pcrl':pcrl, 'pBl':pBl, 'vthl':vthl, 'vturl':vturl,\
               'vkezl':vkezl, 'vvzl':vvzl, 'vcrl':vcrl, 'vBl':vBl, 'ptotl':ptotl, 'vtotl':vtotl}
    if outHI==1:
        outdict['pturcoldl'] = pturcoldl;
        outdict['pturhotl'] = pturhotl;
        outdict['pturHIl'] = pturHIl;
        outdict['pthHIl'] = pthHIl;
        outdict['rhocoldl'] = rhocoldl;
        outdict['rhohotl'] = rhohotl;
        outdict['rhoHIl'] = rhoHIl;
        outdict['vturhotl'] = vturhotl
        outdict['vturcoldl'] = vturcoldl
        outdict['vturHIl'] = vturHIl
        outdict['vthHIl'] = vthHIl
    return outdict



def turdifscale(ssdict):
    nested_dict = lambda: collections.defaultdict(nested_dict)
    plotdict = nested_dict()
    outgriddict = nested_dict()
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
    cutcold=ssdict['cutcold']
    usekez=ssdict['usekez']
    withoutr=ssdict['withoutr']
    outHI=ssdict['outHI']
    info=SSF.outdirname(runtodo, Nsnap)
    newlabel=info['newlabel']
    runlabel=info['runlabel']
    newlabelneed=1
    if newlabelneed==1:
        labelneed="\n".join(wrap(newlabel,17))
    ptitle=labelneed
    if not (wanted=='pz' or wanted=='pr' or wanted=='vz' or wanted=='vr' or wanted=='vturr' or wanted=='vturz'):
        print 'wrong wanted'
        return None
    dclabelneed=1
    newlabelneed=1
    print 'runtodo', runtodo
    dr = 1.0 #kpc
    if wanted=='pr' or wanted=='vr' or wanted=='vturr':
        nogrid=int(withinr/dr)
        vertical=0
        horizontal=1
    if wanted=='pz' or wanted=='vz' or wanted=='vturz':
        
        vertical=1
        horizontal=0

    griddirlist = ['grid0_125kpc','grid0_25kpc','grid0_5kpc','grid1kpc']
    lsnlist = ['solid','dashed','dashdot','dotted']
    colorlist = ['y','c','m','k']
    labellist=['l = 0.125 kpc', 'l = 0.25 kpc','l = 0.5 kpc','l = 1 kpc']
    for griddir in griddirlist:
        outdict = outvtur(runtodo,wanted,startno,Nsnap,snapsep,nogrid,griddir,withinr,maxlength,\
                         usehalfz=usehalfz,cutcold=cutcold,vertical=vertical,\
                         horizontal=horizontal,withoutr=withoutr,dr=dr,outHI=outHI,\
                         usepredata=usepredata)
        outgriddict[griddir]=outdict
    #for key in outgriddict:
    #    print 'key', key
    #    for inkey in outgriddict[key]:
    #        print 'inkey', inkey
    linelist = []
    if wanted=='pz' or wanted=='vz' or wanted=='vturz':
        plotdict[wanted]['xlab'] = r'$z\;{\rm[kpc]}$'
        plotdir = 'crdenz'
    elif wanted=='pr' or wanted=='vr' or wanted=='vturr':
        plotdict[wanted]['xlab'] = r'$r\;{\rm[kpc]}$'
        plotdir = 'crdenr'
    filename=plotloc+'CRplot/'+plotdir+'/turdifscale_'+wanted+'_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    if wanted=='pz' or wanted=='pr':
        plotdict[wanted]['ylab'] = r'$\Pi_z {\rm [dyne/cm^2]}$'
    elif wanted=='vz' or wanted=='vr' or wanted=='vturr':
        plotdict[wanted]['ylab'] = r'$\sqrt{\left \langle v_z^2 \right \rangle} {\rm [km/s]}$'    
    plotdict[wanted]['ptitle'] = ptitle
    plotdict[wanted]['runtitle'] = runlabel
    for i, griddir in enumerate(griddirlist):
        plotdict[wanted]['xnl'][griddir] = outgriddict[griddir]['zlist'][:-1]
        if wanted=='pz' or wanted=='pr':    
            plotdict[wanted]['ynl'][griddir] = outgriddict[griddir]['pturl'][:-1]
        if wanted=='vz' or wanted=='vr' or wanted=='vturr':
            plotdict[wanted]['ynl'][griddir] = outgriddict[griddir]['vturl'][:-1]
        plotdict[wanted]['linelab'][griddir] = labellist[i]
        plotdict[wanted]['lw'][griddir] = 2    
        plotdict[wanted]['marker'][griddir] = ''    
        plotdict[wanted]['lsn'][griddir] = lsnlist[i]
        plotdict[wanted]['color'][griddir] = colorlist[i]
        linelist.append(griddir)
    plotdict[wanted]['linelist']=linelist
    plotdict[wanted]['filename'] = filename        
    return plotdict
    
    
    