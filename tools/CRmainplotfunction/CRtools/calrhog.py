from stdmodandoption import *
import plot_setup as PS
import collections


def calrhog(runtodo,wanted,startno,Nsnap,snapsep,fmeat):
    nested_dict = lambda: collections.defaultdict(nested_dict)
    plotdict = nested_dict()
    if wanted=='crez' or wanted=='rhog':
        def func(x, a, b, c):
                return a+b*x+c*x*x
        def d10func(r, a, b, c):
                return np.power(10.,func(np.log10(r), a, b, c))*(b+2*c*np.log(r)/np.log(10.))/r
        def funcab(x, a, b):
                return a+b*x
        rcParams['figure.figsize'] = 5,4
        resoneed=0
        rotface=1
        needlogz=0
        thermalneed=1
        needfit=0
        labcount=0
        linestylabel=0
        newlabelneed=1
        twolegend=1
        ptotneed = 1
        print 'runtodo', runtodo
        maxlength = 30.
        minlength=5.
        nogrid=8
        withinr=5.0
        diskr=10
        diskh = 1.0
        if needlogz==1:
                zl=np.logspace(-1,np.log10(maxlength),num=nogrid)
        else:
                zl=np.linspace(minlength,maxlength,num=nogrid)
        dz = maxlength/nogrid
        crden=zl*0.0
        thden=zl*0.0
        Gmden=zl*0.0
        mtotl=zl*0.0
        mdtot = 0.0
        numoftimes=0
        info=SSF.outdirname(runtodo, Nsnap)
        havecr=info['havecr']
        color=info['color']
        runtitle=info['runtitle']
        ptitle=title
        if runtitle=='SMC':
                ptitle='Dwarf'
        elif runtitle=='SBC':
                ptitle='Starburst'
        elif runtitle=='MW':
                ptitle=r'$L\star$ Galaxy'
        if wanted=='crez':
                if havecr==0:
                        return None
        for i in range(startno,Nsnap,snapsep):
                info=SSF.outdirname(runtodo, i)
                M1speed=info['M1speed']
                rundir=info['rundir']
                runtitle=info['runtitle']
                slabel=info['slabel']
                snlabel=info['snlabel']
                dclabel=info['dclabel']
                resolabel=info['resolabel']
                the_snapdir=info['the_snapdir']
                snapshot_name = info['the_prefix']
                extension = info['the_suffix']
                Nsnapstring=info['Nsnapstring']
                Fcal=info['Fcal']
                iavesfr=info['iavesfr']
                timestep=info['timestep']
                haveB=info['haveB']
                newlabel=info['newlabel']
                cosmo=info['cosmo']
                halostr=info['halostr']
                firever=info['firever']
                maindir=info['maindir']
                usepep=info['usepep']
                snumadd=info['snumadd']
                h0=cosmo
                datasup=0
                Mdict = collections.defaultdict(dict)
                ptlist = np.array([])
                ptypelist=[1,2,3,4] #gas is special
                G = SSF.readsnapwcen(the_snapdir, Nsnapstring, 0, snapshot_name=snapshot_name, extension=extension,\
                 havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
                 datasup=datasup,runtodo=runtodo,rundir=rundir,halostr=halostr)
                Mdict[0]=G
                ptlist=np.append(ptlist,0)
                Gp = G['p']; Gx = Gp[:,0]; Gy = Gp[:,1]; Gz = Gp[:,2];
                Gv = G['v']; Gvx = Gv[:,0]; Gvy = Gv[:,1]; Gvz = Gv[:,2];
                Grho = G['rho']; Gu = G['u']; Gm = G['m']
                for ptype in ptypelist:
                    try:
                        Pt = SSF.readsnapwcen(the_snapdir, Nsnapstring, ptype, snapshot_name=snapshot_name, extension=extension,\
                         havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
                         datasup=datasup,runtodo=runtodo,rundir=rundir,halostr=halostr)
                        Mdict[ptype]=Pt
                        ptlist=np.append(ptlist,ptype)
                    except KeyError:
                        continue
                if havecr>0:
                        cregy = G['cregy']*1e10*Msun_in_g*km_in_cm*km_in_cm #cosmic ray energy in 1e10Msun km^2/sec^2
                if haveB==0:
                        lsn='solid'
                else:
                        lsn='dashed'
                #if M1speed>999:
                #       lsn='dotted'
                #if M1speed>1999:
                #       lsn='dashdot'
                if wanted=='rhog':
                        mdisktot=0.0
                        print 'diskr, diskh', diskr, diskh
                        for ipd in [0,2,4]:
                                try:
                                        Pe = Mdict[ipd]
                                        Px = Pe['p'][:,0]; Py = Pe['p'][:,1]; Pz = Pe['p'][:,2];  
                                        Pm = Pe['m']
                                        Prxy = np.sqrt(Px*Px+Py*Py)
                                        cutxyz = (Prxy<diskr)*(np.absolute(Pz)<diskh) #kpc
                                        mdisktot += np.sum(Pm[cutxyz]*1e10*Msun_in_g) #in g
                                        print 'mdisktot', mdisktot
                                except KeyError:
                                        continue
                        mdtot += mdisktot
                for iz in range(len(zl)-1):
                        cutxy = (Gx*Gx+Gy*Gy < withinr*withinr)
                        cutzu = Gz<zl[iz+1]
                        cutzd = Gz>zl[iz]
                        cutz = cutzu*cutzd
                        cut = cutxy*cutz
                        #print 'withinr', withinr
                        #print 'Gm[cut]', Gm[cut]
                        #print 'zl', zl
                        if havecr>0:
                                crecut = cregy[cut]
                                crden[iz] += np.sum(crecut)/(zl[iz+1]-zl[iz])/kpc_in_cm/(np.pi*withinr*withinr*kpc_in_cm*kpc_in_cm)
                        thermalcut = Gm[cut]*Gu[cut]*1e10*Msun_in_g*km_in_cm*km_in_cm
                        thden[iz] += np.sum(thermalcut)/(zl[iz+1]-zl[iz])/kpc_in_cm/(np.pi*withinr*withinr*kpc_in_cm*kpc_in_cm)
                        if wanted=='rhog':
                                Gmcut=Gm[cut]*1e10*Msun_in_g
                                Gmden[iz] += np.sum(Gmcut)/(zl[iz+1]-zl[iz])/kpc_in_cm/(np.pi*withinr*withinr*kpc_in_cm*kpc_in_cm)
                                mintot = 0
                                for ipa in ptlist:
                                        try:
                                                Pe = Mdict[ipa]
                                                Px = Pe['p'][:,0]; Py = Pe['p'][:,1]; Pz = Pe['p'][:,2];
                                                Pm = Pe['m']
                                                Pr = np.sqrt(Px*Px+Py*Py+Pz*Pz)
                                                cutpr = Pr<zl[iz]
                                                if ipa==0 or ipa==2 or ipa==4:
                                                        Prxy = np.sqrt(Px*Px+Py*Py)
                                                        if zl[iz]<diskr:
                                                                drc = zl[iz]
                                                        else:
                                                                drc = diskr
                                                        cutdisk = (Prxy<drc)*(np.absolute(Pz)<diskh)
                                                        mdiskcut = np.sum(Pm[cutdisk]*1e10*Msun_in_g)
                                                Pmtot = np.sum(Pm[cutpr]*1e10*Msun_in_g)-mdiskcut #in g
                                                print 'Pmtot', Pmtot
                                        except KeyError:
                                                continue
                                        mintot += Pmtot
                                mtotl[iz] += mintot
                numoftimes+=1
        if havecr>0:
                crden = crden/numoftimes
        thden = thden/numoftimes
        Gmden = Gmden/numoftimes
        mtotl = mtotl/numoftimes
        mdtot = mdtot/numoftimes
        if havecr>0:
                print 'crden', crden
        print 'thden', thden
        print 'mtotl', mtotl
        print 'mdtot', mdtot
        labelneed=dclabel
        if newlabelneed==1:
                labelneed="\n".join(wrap(newlabel,17))
        if resoneed==1:
                if resolabel=='llr':
                        labelneed='low res'
                        lsn = 'solid'
                if resolabel=='lr':
                        labelneed='std res'
                        lsn = 'dashed'
                if resolabel=='mr':
                        labelneed='high res'
                        lsn = 'dashdot'
        zlm = (zl[1:]+zl[:-1])/2.
        if needfit==1:
                x0 = [0.1,-1.0,0.0]
                thden1 = thden[:-1]
                xdata=np.log10(zlm[~np.isnan(thden1)])
                ydata=np.log10(thden1[~np.isnan(thden1)])
                outfit=optimize.curve_fit(func, xdata, ydata, x0)
                afit=outfit[0][0]
                bfit=outfit[0][1]
                cfit=outfit[0][2]
                dthcr = -d10func(zlm,afit,bfit,cfit)*(GAMMA-1.0)/kpc_in_cm
        else:
                dthcr = -(thden[:-1]-thden[1:])/(zl[:-1]-zl[1:])*(GAMMA-1.0)/kpc_in_cm
        if havecr>0:
                if needfit==1:
                        x0 = [0.1,-1.0,0.0]
                        crden1 = crden[:-1]
                        xdata=np.log10(zlm)
                        ydata=np.log10(crden1)
                        xdata = xdata[~np.isinf(ydata)]
                        ydata = ydata[~np.isinf(ydata)]
                        print 'xdata', xdata
                        print 'ydata', ydata
                        outfit=optimize.curve_fit(func, xdata, ydata, x0)
                        afit=outfit[0][0]
                        bfit=outfit[0][1]
                        cfit=outfit[0][2]
                        dpcr = -d10func(zlm,afit,bfit,cfit)*(CRgamma-1.0)/kpc_in_cm
                else:
                        dpcr = -(crden[:-1]-crden[1:])/(zl[:-1]-zl[1:])*(CRgamma-1.0)/kpc_in_cm
        if wanted=='crez':
                plt.plot(zlm, crden[:-1]*erg_in_eV, label=labelneed,lw=2,ls=lsn)    
        if wanted=='rhog':
                gdisk = 2.0*NewtonG_cgs*mdtot/diskr/diskr*(1.0-zl/np.sqrt(diskr*diskr+zl*zl))/kpc_in_cm/kpc_in_cm
                gsphere = mtotl*NewtonG_cgs/(zl*kpc_in_cm)/(zl*kpc_in_cm)
                rhog = Gmden*(gdisk+gsphere)
    if wanted=='crez':
        plotdict[wanted]['xlab'] = 'z [kpc]'
        plotdict[wanted]['xnl'] = zl[:-1]
        plotdict[wanted]['ylab'] = r'$e_{\rm cr} [{\rm eV/cm^3}]$'
        plotdict[wanted]['ynl'] = crden[:-1]*erg_in_eV
        plotdict[wanted]['runtodo'] = runtodo
        plotdict[wanted]['labelneed'] = labelneed
        plotdict[wanted]['lsn'] = lsn
        plotdict[wanted]['color'] = color
        plotdict[wanted]['runtitle'] = runtitle
        plotdict[wanted]['ptitle'] = ptitle
        filename=plotloc+'CRplot/crez/crez_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        plotdict[wanted]['filename'] = filename
    if wanted=='rhog':
        plotdict[wanted]['xlab'] = 'z [kpc]'
        plotdict[wanted]['xnl']['rhog'] = zl[:-1]
        plotdict[wanted]['xnl']['dpth'] = zlm[:-1]
        if havecr>0: plotdict[wanted]['xnl']['dpcr'] = zlm[:-1]
        plotdict[wanted]['ylab'] = r'$\frac{\mathrm{d} P}{\mathrm{d} z}\;[{\rm erg/cm^4}]$'
        plotdict[wanted]['ynl']['rhog'] = rhog[:-1]
        plotdict[wanted]['ynl']['dpth'] = dthcr[:-1]
        if havecr>0: plotdict[wanted]['ynl']['dpcr'] = dpcr[:-1]
        plotdict[wanted]['runtodo'] = runtodo
        plotdict[wanted]['labelneed'] = labelneed
        plotdict[wanted]['linelab']['rhog'] = r'$\rho g$'
        plotdict[wanted]['linelab']['dpth'] = r'$\mathrm{d} P_{\rm th}/\mathrm{d} z$'
        if havecr>0: plotdict[wanted]['linelab']['dpcr'] = r'$\mathrm{d} P_{\rm cr}/\mathrm{d} z$'
        plotdict[wanted]['lsn'] = lsn
        plotdict[wanted]['lw']['rhog'] = 2
        plotdict[wanted]['lw']['dpth'] = 1
        if havecr>0: plotdict[wanted]['lw']['dpcr'] = 1
        plotdict[wanted]['marker']['rhog'] = 'None'
        plotdict[wanted]['marker']['dpth'] = '^'
        if havecr>0: plotdict[wanted]['marker']['dpcr'] = 'o'
        plotdict[wanted]['color'] = color
        plotdict[wanted]['runtitle'] = runtitle
        plotdict[wanted]['ptitle'] = ptitle
        filename=plotloc+'CRplot/rhog/rhog_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        plotdict[wanted]['filename'] = filename
        plotdict[wanted]['legendneed'] = 0                
        if runtitle=='SMC' or linestylabel==1: plotdict[wanted]['legendneed'] = 1
    return plotdict