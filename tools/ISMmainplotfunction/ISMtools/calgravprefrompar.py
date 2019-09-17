from stdmodandoption import *
import plot_setup as PS
import collections
from intrhogfunc import intrhogfunc, intdpfunc




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
    if not (wanted=='dpdz' or wanted=='pz' or wanted=='vg' or wanted=='ag'):
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
    vardict={}
    namelist=['rhogl','intrhogl','pthl','pturl',\
              'pkezl','pvzl','pcrl','pBl','intgl','gl']
    dnamelist=['dpthl','dpturl','dpcrl','dpBl','dptotl','dpkezl','intdpthl',\
              'intdpturl','intdpkezl','intdpcrl','intdpBl',\
              'intdpthrhol','intdpturrhol','intdpkezrhol',\
              'intdpcrrhol','intdpBrhol',\
              'dpthrhol','dpturrhol','dpkezrhol',\
              'dpcrrhol','dpBrhol',\
              ]
    for namel in namelist:
        vardict[namel]=np.zeros(nogrid);
    for dnamel in dnamelist:
        vardict[dnamel]=np.zeros(nogrid-1);
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
            data = SSF.calrhogfrompreexist(runtodo,i,withinr,maxlength,\
                                           usehalfz=usehalfz,griddir=griddir,cutcold=cutcold,withoutr=withoutrnew,outHI=1)            
        else:
            data = SSF.calrhogfrompar(runtodo,i,withinr,maxlength,nogrid,\
                                  havecr=havecr,haveB=haveB,usehalfz=usehalfz)
        #datalist=['rhog','rhol','pthl','pturl','gzlist','pkezl','pvzl','pcrl','pBl']
        zlist=data['zlist'];
        #cutzmax = np.absolute(zlist)
        rhog=data['rhog']; rhol=data['rhol'];
        vardict['pthl'] += data['pthl']; vardict['pturl'] += data['pturl'];
        gzlist = data['gzlist']; vardict['pkezl'] += data['pkezl'];
        vardict['pvzl'] += data['pvzl'];
        if havecr>0: vardict['pcrl'] += data['pcrl']; 
        if haveB>0: vardict['pBl'] += data['pBl'];
        dzlist = (zlist[1:]-zlist[:-1])*kpc_in_cm
        zlistm = (zlist[1:]+zlist[:-1])/2.0
        rholm = np.interp(zlistm,zlist,rhol)
        dpthl = (data['pthl'][1:]-data['pthl'][:-1])/dzlist;
        dpturl = (data['pturl'][1:]-data['pturl'][:-1])/dzlist;
        dpkezl = (data['pkezl'][1:]-data['pkezl'][:-1])/dzlist;
        vardict['dpthl'] += dpthl;
        vardict['dpturl'] += dpturl;
        vardict['dpkezl'] += dpkezl;
        if havecr>0:
            dpcrl = (data['pcrl'][1:]-data['pcrl'][:-1])/dzlist;
            vardict['dpcrl'] += dpcrl;
            #print 'pcrl', vardict['pcrl']
        if haveB>0:
            dpBl = (data['pBl'][1:]-data['pBl'][:-1])/dzlist;
            vardict['dpBl'] += dpBl;
        vardict['rhogl'] += rhog;
        zeropoint = np.interp(0.0,rhog,zlist)
        #print 'zeropoint', zeropoint
        vardict['intrhogl'] += intdpfunc(rhog,zeropoint,zlist)
        vardict['intdpthl'] += intdpfunc(dpthl,zeropoint,zlistm)
        vardict['intdpturl'] += intdpfunc(dpturl,zeropoint,zlistm)
        vardict['intdpkezl'] += intdpfunc(dpkezl,zeropoint,zlistm)
        if havecr>0:
            intdpcrl=intdpfunc(dpcrl,zeropoint,zlistm)
            vardict['intdpcrl'] += intdpcrl
            #print 'intdpcrl', intdpcrl
        if haveB>0:
            vardict['intdpBl'] += intdpfunc(dpBl,zeropoint,zlistm)
        if wanted=='vg' or wanted=='ag':
            dpthrhol = dpthl/rholm;
            dpturrhol = dpturl/rholm;
            dpkezrhol = dpkezl/rholm;
            if havecr>0:
                dpcrrhol = dpcrl/rholm;
            if haveB>0:
                dpBrhol  = dpBl/rholm;
            gl = rhog/rhol;
            vardict['gl'] += gl;
            vardict['dpthrhol'] += dpthrhol;
            vardict['dpturrhol'] += dpturrhol;
            vardict['dpkezrhol'] += dpkezrhol;
            if havecr>0:                
                vardict['dpcrrhol'] += dpcrrhol;
            if haveB>0:
                vardict['dpBrhol'] += dpBrhol;
            vardict['intgl'] += intdpfunc(gl,zeropoint,zlist,sumfromcen=1)
            vardict['intdpthrhol'] += intdpfunc(dpthrhol,zeropoint,zlistm,sumfromcen=1)
            vardict['intdpturrhol'] += intdpfunc(dpturrhol,zeropoint,zlistm,sumfromcen=1)
            vardict['intdpkezrhol'] += intdpfunc(dpkezrhol,zeropoint,zlistm,sumfromcen=1)
            if havecr>0:                
                vardict['intdpcrrhol'] += intdpfunc(dpcrrhol,zeropoint,zlistm,sumfromcen=1)
            if haveB>0:
                vardict['intdpBrhol'] += intdpfunc(dpBrhol,zeropoint,zlistm,sumfromcen=1)            
        nooftimes+=1
    for namel in namelist:
        vardict[namel]=vardict[namel]/nooftimes
    for dnamel in dnamelist:
        vardict[dnamel]=vardict[dnamel]/nooftimes
    vardict['dptotl'] = vardict['dpthl']
    smoothing=0
    if smoothing==1:
        box_pts = 4
        for i, namel in enumerate(namelist):
            vardict[namel] = np.log10(SSF.moving_average(np.power(10,vardict[namel]), box_pts))
    if usekez==0:
         vardict['dptotl'] +=  vardict['dpturl']
    else:
         vardict['dptotl'] +=  vardict['dpkezl']        
    vardict['ptotl'] =  vardict['pthl']+ vardict['pturl']
    if usekez==0:
        vardict['ptotl'] +=  vardict['pturl']
    else:
        vardict['ptotl'] +=  vardict['pkezl']
    if havecr>0:
        vardict['dptotl'] +=  vardict['dpcrl']
        vardict['ptotl'] +=  vardict['pcrl']
    if haveB>0:
        vardict['dptotl'] +=  vardict['dpBl']
        vardict['ptotl'] +=  vardict['pBl']

    #print 'zlist', zlist
    #print 'rhog', rhog
    #print 'zip(zlist,rhog)', zip(zlist,rhog)

    if wanted=='ag':
        linelist = []
        plotdict[wanted]['xlab'] = r'$z\;{\rm[kpc]}$'
        plotdict[wanted]['ylab'] = r'$a\;{\rm [cm/s^2]}$'
        plotdict[wanted]['ptitle'] = ptitle
        plotdict[wanted]['runtitle'] = runlabel
        plotdict[wanted]['xnl']['rhog'] = zlist
        ag = np.absolute(vardict['gl'])
        print 'ag', ag
        plotdict[wanted]['ynl']['rhog'] = ag
        plotdict[wanted]['lw']['rhog'] = 1    
        plotdict[wanted]['marker']['rhog'] = 'None'    
        plotdict[wanted]['lsn']['rhog'] = 'dashed'
        plotdict[wanted]['linelab']['rhog'] = 'gravity'
        plotdict[wanted]['color']['rhog'] = 'k'
        atot = np.interp(zlistm,zlist,ag)
        linelist.append('rhog')
        plotdict[wanted]['xnl']['eth'] = zlistm
        ath = np.absolute(vardict['dpthrhol'])
        plotdict[wanted]['ynl']['eth'] = ath   
        plotdict[wanted]['linelab']['eth'] = 'thermal'
        plotdict[wanted]['lw']['eth'] = 1    
        plotdict[wanted]['marker']['eth'] = 'None'
        plotdict[wanted]['lsn']['eth'] = 'dotted'
        plotdict[wanted]['color']['eth'] = 'r'
        linelist.append('eth')
        atot = atot-ath
        if havecr>0:
            plotdict[wanted]['xnl']['ecr'] = zlistm
            acr = np.absolute(vardict['dpcrrhol'])
            print 'acr', acr
            plotdict[wanted]['ynl']['ecr'] = acr
            plotdict[wanted]['lw']['ecr'] = 1  
            plotdict[wanted]['lsn']['ecr'] = 'dashdot' 
            plotdict[wanted]['linelab']['ecr'] = 'CR' 
            plotdict[wanted]['marker']['ecr'] = 'None'
            plotdict[wanted]['color']['ecr'] = 'y'
            atot = atot-acr
            linelist.append('ecr')
        plotdict[wanted]['xnl']['etot'] = zlistm
        plotdict[wanted]['ynl']['etot'] = atot
        plotdict[wanted]['linelab']['etot'] = 'total'
        plotdict[wanted]['lw']['etot'] = 1    
        plotdict[wanted]['marker']['etot'] = 'None'    
        plotdict[wanted]['lsn']['etot'] = 'solid'
        plotdict[wanted]['color']['etot'] = '0.5'
        linelist.append('etot')
        plotdict[wanted]['linelist']=linelist
        filename=plotloc+'CRplot/crdenz/agfrompar_'+wanted+'_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        plotdict[wanted]['filename'] = filename    
    
    
    
    if wanted=='vg':
        linelist = []
        plotdict[wanted]['xlab'] = r'$z\;{\rm[kpc]}$'
        plotdict[wanted]['ylab'] = r'$v_{\rm turnover}\;{\rm [km/s]}$'
        plotdict[wanted]['ptitle'] = ptitle
        plotdict[wanted]['runtitle'] = runlabel
        plotdict[wanted]['xnl']['rhog'] = zlist
        vg = np.sqrt(np.absolute(vardict['intgl']))/km_in_cm
        plotdict[wanted]['ynl']['rhog'] = vg
        plotdict[wanted]['lw']['rhog'] = 1    
        plotdict[wanted]['marker']['rhog'] = 'None'    
        plotdict[wanted]['lsn']['rhog'] = 'dashed'
        plotdict[wanted]['linelab']['rhog'] = 'gravity'
        plotdict[wanted]['color']['rhog'] = 'k'
        vtot = np.interp(zlistm,zlist,vg)
        linelist.append('rhog')
        #plotdict[wanted]['xnl']['eth'] = zlistm
        #vth = np.sqrt(np.absolute(vardict['intdpthrhol']))/km_in_cm
        #print 'np.amax(intdpthrhol)',np.amax(vardict['intdpthrhol'])
        #plotdict[wanted]['ynl']['eth'] = vth   
        #plotdict[wanted]['linelab']['eth'] = 'thermal'
        #plotdict[wanted]['lw']['eth'] = 1    
        #plotdict[wanted]['marker']['eth'] = 'None'
        #plotdict[wanted]['lsn']['eth'] = 'dotted'
        #plotdict[wanted]['color']['eth'] = 'r'
        #linelist.append('eth')
        #vtot = vtot-vth
        if havecr>0:
            plotdict[wanted]['xnl']['ecr'] = zlistm
            vcr = np.sqrt(np.absolute(vardict['intdpcrrhol']))/km_in_cm
            plotdict[wanted]['ynl']['ecr'] = vcr
            plotdict[wanted]['lw']['ecr'] = 1  
            plotdict[wanted]['lsn']['ecr'] = 'dashdot' 
            plotdict[wanted]['linelab']['ecr'] = 'CR' 
            plotdict[wanted]['marker']['ecr'] = 'None'
            plotdict[wanted]['color']['ecr'] = 'y'
            vtot = vtot-vcr
            linelist.append('ecr')
            print 'np.amax(intdpcrrhol)',np.amax(vardict['intdpcrrhol'])
        plotdict[wanted]['xnl']['etot'] = zlistm
        plotdict[wanted]['ynl']['etot'] = vtot
        plotdict[wanted]['linelab']['etot'] = 'total (- CR)'
        plotdict[wanted]['lw']['etot'] = 1    
        plotdict[wanted]['marker']['etot'] = 'None'    
        plotdict[wanted]['lsn']['etot'] = 'solid'
        plotdict[wanted]['color']['etot'] = '0.5'
        linelist.append('etot')
        plotdict[wanted]['linelist']=linelist
        filename=plotloc+'CRplot/crdenz/vgfrompar_'+wanted+'_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        plotdict[wanted]['filename'] = filename
        
    if wanted=='dpdz':
        plotdict[wanted]['xlab'] = r'$z\;{\rm[kpc]}$'
        plotdict[wanted]['ylab'] = r'${\rm d} P/{\rm d} z\; {\rm [dyne/cm^3]}$'
        plotdict[wanted]['ptitle'] = ptitle
        plotdict[wanted]['runtitle'] = runlabel
        plotdict[wanted]['xnl']['rhog'] = zlist
        plotdict[wanted]['ynl']['rhog'] = np.absolute(vardict['rhogl'])
        plotdict[wanted]['lw']['rhog'] = 1    
        plotdict[wanted]['marker']['rhog'] = 'None'    
        plotdict[wanted]['lsn']['rhog'] = 'dashed'
        #plotdict[wanted]['linelab']['rhog'] = r'$\rho g$'
        plotdict[wanted]['linelab']['rhog'] = 'gravity'
        plotdict[wanted]['color']['rhog'] = 'k'
        plotdict[wanted]['xnl']['eth'] = zlistm
        plotdict[wanted]['ynl']['eth'] = vardict['dpthl']    
        #plotdict[wanted]['linelab']['eth'] = r'$p_{\rm th}$'
        plotdict[wanted]['linelab']['eth'] = 'thermal'
        plotdict[wanted]['lw']['eth'] = 1    
        plotdict[wanted]['marker']['eth'] = 'none'
        plotdict[wanted]['lsn']['eth'] = 'dotted'
        plotdict[wanted]['color']['eth'] = 'r'
        plotdict[wanted]['xnl']['etur'] = zlistm
        if usekez==0:
            plotdict[wanted]['ynl']['etur'] = vardict['dpturl'] 
            plotdict[wanted]['linelab']['etur'] = 'turbulence'
        else:
            plotdict[wanted]['ynl']['etur'] = vardict['dpkezl'] 
            plotdict[wanted]['linelab']['etur'] = 'kinetic'
        plotdict[wanted]['lw']['etur'] = 2    
        plotdict[wanted]['marker']['etur'] = 'none'    
        plotdict[wanted]['lsn']['etur'] = 'solid'
        plotdict[wanted]['color']['etur'] = 'g'        
        if haveB>0:
            plotdict[wanted]['xnl']['eB'] = zlistm
            plotdict[wanted]['ynl']['eB'] = vardict['dpBl']
            plotdict[wanted]['linelab']['eB'] = 'Bfield'
            #plotdict[wanted]['linelab']['eB'] = r'$(B^2-2B_z^2)/(8\pi)$'
            plotdict[wanted]['lw']['eB'] = 1
            plotdict[wanted]['lsn']['eB'] = 'solid'        
            plotdict[wanted]['marker']['eB'] = 'none'
            plotdict[wanted]['color']['eB'] = 'b'            
        if havecr>0:
            plotdict[wanted]['xnl']['ecr'] = zlistm
            plotdict[wanted]['ynl']['ecr'] = vardict['dpcrl'] 
            plotdict[wanted]['lw']['ecr'] = 1  
            plotdict[wanted]['lsn']['ecr'] = 'dashdot' 
            plotdict[wanted]['linelab']['ecr'] = 'CR' 
            plotdict[wanted]['marker']['ecr'] = 'none'
            plotdict[wanted]['color']['ecr'] = 'y'
        plotdict[wanted]['xnl']['etot'] = zlistm
        plotdict[wanted]['ynl']['etot'] = vardict['dptotl'] 
        plotdict[wanted]['linelab']['etot'] = 'total'
        plotdict[wanted]['lw']['etot'] = 1    
        plotdict[wanted]['marker']['etot'] = 'None'    
        plotdict[wanted]['lsn']['etot'] = 'solid'
        plotdict[wanted]['color']['etot'] = '0.5'
        filename=plotloc+'CRplot/crdenz/rhogdpzfrompar_'+wanted+'_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        plotdict[wanted]['filename'] = filename
    if wanted=='pz':
        delpthl = vardict['intdpthl']
        delpturl = vardict['intdpturl']
        delpkezl = vardict['intdpkezl']
        delptotl = delpthl
        if usekez==0:
            delptotl = delptotl+delpturl
        else:
            delptotl = delptotl+delpkezl
        if haveB>0: 
            delpBl = vardict['intdpBl']
            delptotl = delptotl+delpBl
        if havecr>0:
            delpcrl = vardict['intdpcrl']
            delptotl = delptotl+delpcrl
            print 'delpcrl', delpcrl 
        linelist = []
        plotdict[wanted]['xlab'] = r'$z\;{\rm[kpc]}$'
        plotdict[wanted]['ylab'] = r'$\Delta\Pi_z {\rm [dyne/cm^2]}$'
        plotdict[wanted]['ptitle'] = ptitle
        plotdict[wanted]['runtitle'] = runlabel
        plotdict[wanted]['xnl']['rhog'] = zlist
        plotdict[wanted]['ynl']['rhog'] = np.absolute(vardict['intrhogl'])
        plotdict[wanted]['lw']['rhog'] = 2    
        plotdict[wanted]['marker']['rhog'] = ''    
        plotdict[wanted]['lsn']['rhog'] = 'dashed'
        #plotdict[wanted]['linelab']['rhog'] = r'$\rho g$'
        plotdict[wanted]['linelab']['rhog'] = 'gravity'
        plotdict[wanted]['color']['rhog'] = '0.5'
        linelist.append('rhog')
        plotdict[wanted]['xnl']['eth'] = zlistm
        plotdict[wanted]['ynl']['eth'] = delpthl    
        #plotdict[wanted]['linelab']['eth'] = r'$p_{\rm th}$'
        plotdict[wanted]['linelab']['eth'] = 'thermal'
        plotdict[wanted]['lw']['eth'] = 1    
        plotdict[wanted]['marker']['eth'] = ''
        plotdict[wanted]['lsn']['eth'] = 'dashed'
        plotdict[wanted]['color']['eth'] = 'r'
        linelist.append('eth')
        plotdict[wanted]['xnl']['ethend'] = zlist[-3:-1]
        plotdict[wanted]['ynl']['ethend'] = np.amax(vardict['pthl'])*np.ones(len(zlist[-3:-1]))
        plotdict[wanted]['lw']['ethend'] = 1  
        plotdict[wanted]['lsn']['ethend'] = 'dashdot' 
        plotdict[wanted]['linelab']['ethend'] = '' 
        plotdict[wanted]['marker']['ethend'] = '>'
        plotdict[wanted]['color']['ethend'] = 'r'            
        linelist.append('ethend') 
        plotdict[wanted]['xnl']['etur'] = zlistm
        plotdict[wanted]['ynl']['etur'] = delpturl
        plotdict[wanted]['linelab']['etur'] = 'turbulence'
        plotdict[wanted]['lw']['etur'] = 1    
        plotdict[wanted]['marker']['etur'] = ''    
        plotdict[wanted]['lsn']['etur'] = 'dashed'
        plotdict[wanted]['color']['etur'] = 'g'
        linelist.append('etur')
        if usekez==1:
            plotdict[wanted]['xnl']['ekez'] = zlistm
            plotdict[wanted]['ynl']['ekez'] = delpkezl 
            plotdict[wanted]['linelab']['ekez'] = 'kinetic'
            plotdict[wanted]['lw']['ekez'] = 2    
            plotdict[wanted]['marker']['ekez'] = ''    
            plotdict[wanted]['lsn']['ekez'] = 'dotted'
            plotdict[wanted]['color']['ekez'] = 'g'
            linelist.append('ekez')
            plotdict[wanted]['xnl']['ekezend'] = zlist[-3:-1]
            plotdict[wanted]['ynl']['ekezend'] = np.amax(vardict['pkezl'])*np.ones(len(zlist[-3:-1]))
            plotdict[wanted]['lw']['ekezend'] = 2  
            plotdict[wanted]['lsn']['ekezend'] = 'dotted' 
            plotdict[wanted]['linelab']['ekezend'] = '' 
            plotdict[wanted]['marker']['ekezend'] = '>'
            plotdict[wanted]['color']['ekezend'] = 'g'            
            linelist.append('ekezend')
        if haveB>0:
            plotdict[wanted]['xnl']['eB'] = zlistm
            plotdict[wanted]['ynl']['eB'] = delpBl
            plotdict[wanted]['linelab']['eB'] = 'B field'
            #plotdict[wanted]['linelab']['eB'] = r'$(B^2-2B_z^2)/(8\pi)$'
            plotdict[wanted]['lw']['eB'] = 1
            plotdict[wanted]['lsn']['eB'] = 'solid'        
            plotdict[wanted]['marker']['eB'] = ''
            plotdict[wanted]['color']['eB'] = 'b'
            linelist.append('eB')
        if havecr>0:
            plotdict[wanted]['xnl']['ecr'] = zlistm
            plotdict[wanted]['ynl']['ecr'] = delpcrl
            plotdict[wanted]['lw']['ecr'] = 1  
            plotdict[wanted]['lsn']['ecr'] = 'dashdot' 
            plotdict[wanted]['linelab']['ecr'] = 'CR' 
            plotdict[wanted]['marker']['ecr'] = ''
            plotdict[wanted]['color']['ecr'] = 'y'            
            linelist.append('ecr')
            plotdict[wanted]['xnl']['ecrend'] = zlist[-3:-1]
            plotdict[wanted]['ynl']['ecrend'] = np.amax(vardict['pcrl'])*np.ones(len(zlist[-3:-1]))
            plotdict[wanted]['lw']['ecrend'] = 1  
            plotdict[wanted]['lsn']['ecrend'] = 'dashdot' 
            plotdict[wanted]['linelab']['ecrend'] = '' 
            plotdict[wanted]['marker']['ecrend'] = '>'
            plotdict[wanted]['color']['ecrend'] = 'y'            
            linelist.append('ecrend')         
        plotdict[wanted]['xnl']['etot'] = zlistm
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
    
    
    