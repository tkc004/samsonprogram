def calrhogfrompreexist(runtodo,Nsnap,withinr,maxlength,griddir='grid1kpc',\
                        usehalfz=1,cutcold=0,vertical=1,horizontal=0,withoutr=-1.0,dr=1.0,outHI=0):
    zmax=maxlength/2.;
    predata = readpreexist(runtodo,Nsnap,griddir,cutcold=cutcold,outHI=outHI)
    #for key in predata:
    #    print 'key', key
    xlist = predata['xlist'];ylist = predata['ylist'];
    zlist = predata['zlist'];gzgrid = predata['gzgrid'];
    rhogrid = predata['rhogrid']; pthgrid = predata['pthgrid'];
    if outHI==1:
        rhohotgrid = predata['rhohotgrid'];
        rhocoldgrid = predata['rhocoldgrid'];
        rhoHIgrid = predata['rhoHIgrid'];
        pturhotgrid = predata['pturhotgrid'];
        pturcoldgrid = predata['pturcoldgrid'];
        pturHIgrid = predata['pturHIgrid'];
        pthHIgrid = predata['pthHIgrid'];
    pBgrid = predata['pBgrid']; pturgrid = predata['pturgrid'];
    pkezgrid = predata['pkezgrid']; vzgrid = predata['vzgrid'];
    pcrgrid = predata['pcrgrid']; volgrid = predata['volgrid'];
    havecr=predata['havecr']; haveB=predata['haveB'];
    #print 'SSF pturgrid', pturgrid
    #print 'SSF volgrid', volgrid
    #print 'SSF rhogrid', rhogrid
    if vertical==1:
        gzl=zlist*0.0; rholl=zlist*0.0;
        if outHI==1:
            rhohotll=zlist*0.0;
            rhocoldll=zlist*0.0;
            rhoHIll=zlist*0.0;
            pturhotll=zlist*0.0;
            pturcoldll=zlist*0.0;
            pturHIll=zlist*0.0;
            pthHIll=zlist*0.0;
        pthll=zlist*0.0; pturll=zlist*0.0;
        pcrll=zlist*0.0; pBll=zlist*0.0;
        volll=zlist*0.0; pkezll=zlist*0.0;
        pvzll=zlist*0.0;
        dr = np.absolute(zlist[1]-zlist[0])
        for i in range(len(xlist)):
            for j in range(len(ylist)):
                rxy = np.sqrt(xlist[i]*xlist[i]+ylist[j]*ylist[j])
                if withoutr>0:
                    smallr = withoutr; bigr = withinr
                else:
                    smallr = withinr-dr; bigr = withinr+dr;
                if ((rxy<bigr) and (rxy>smallr)):
                    gzl+=gzgrid[i,j,:]*volgrid[i,j,:]
                    rholl+=rhogrid[i,j,:]*volgrid[i,j,:]
                    pthll+=pthgrid[i,j,:]*volgrid[i,j,:]
                    pturll+=pturgrid[i,j,:]*volgrid[i,j,:]
                    if outHI==1:
                        rhohotll+=rhohotgrid[i,j,:]*volgrid[i,j,:]
                        rhocoldll+=rhocoldgrid[i,j,:]*volgrid[i,j,:]
                        rhoHIll+=rhoHIgrid[i,j,:]*volgrid[i,j,:]
                        pturhotll+=pturhotgrid[i,j,:]*volgrid[i,j,:]
                        pturcoldll+=pturcoldgrid[i,j,:]*volgrid[i,j,:]
                        pturHIll+=pturHIgrid[i,j,:]*volgrid[i,j,:]
                        pthHIll+=pthHIgrid[i,j,:]*volgrid[i,j,:]
                    pkezll+=pkezgrid[i,j,:]*volgrid[i,j,:]
                    pvzll+=vzgrid[i,j,:]*vzgrid[i,j,:]*rhogrid[i,j,:]*volgrid[i,j,:]
                    if havecr>0:
                        pcrll+=pcrgrid[i,j,:]*volgrid[i,j,:]
                    if haveB>0:
                        pBll+=pBgrid[i,j,:]*volgrid[i,j,:]
                    volll+=volgrid[i,j,:]
        gz=gzlist=gzl/volll
        rhol=rholl/volll
        pthl=pthll/volll
        pturl=pturll/volll
        if outHI==1:
            rhohotl=rhohotll/volll
            rhocoldl=rhocoldll/volll
            rhoHIl=rhoHIll/volll
            pturhotl=pturhotll/volll
            pturcoldl=pturcoldll/volll
            pturHIl=pturHIll/volll
            pthHIl=pthHIll/volll
        pkezl=pkezll/volll
        pcrl=pcrll/volll
        pBl=pBll/volll
        pvzl=pvzll/volll
        rhog = gz*rhol
        rhol[~np.isfinite(rhol)] = 0;
        if outHI==1:
            rhohotl[~np.isfinite(rhohotl)] = 0; 
            rhocoldl[~np.isfinite(rhocoldl)] = 0;
            rhoHIl[~np.isfinite(rhoHIl)] = 0;
            pturhotl[~np.isfinite(pturhotl)] = 0;
            pturcoldl[~np.isfinite(pturcoldl)] = 0;
            pturHIl[~np.isfinite(pturHIl)] = 0;
            pthHIl[~np.isfinite(pthHIl)] = 0;
        pthl[~np.isfinite(pthl)] = 0;
        pturl[~np.isfinite(pturl)] = 0;
        pkezl[~np.isfinite(pkezl)] = 0;
        pvzl[~np.isfinite(pvzl)] = 0;
        if havecr>0:
            pcrl[~np.isfinite(pcrl)] = 0;
        if haveB>0:
            pBl[~np.isfinite(pBl)] = 0; 
        rhol[~np.isfinite(rhol)] = 0;
        rhog[~np.isfinite(rhog)] = 0;
        if usehalfz==1:
            rhog=rhog[zlist>0];
            rhol=rhol[zlist>0]; 
            gz=gz[zlist>0];
            pthl=pthl[zlist>0];
            pturl=pturl[zlist>0];
            if outHI==1:
                rhohotl=rhohotl[zlist>0]; 
                rhocoldl=rhocoldl[zlist>0]; 
                rhoHIl=rhol[zlist>0]; 
                pturhotl=pturhotl[zlist>0];
                pturcoldl=pturcoldl[zlist>0];
                pturHIl=pturHIl[zlist>0];
                pthHIl=pthHIl[zlist>0];
            pkezl=pkezl[zlist>0];
            pvzl=pvzl[zlist>0];
            if havecr>0:
                pcrl=pcrl[zlist>0];
            if haveB>0:
                pBl=pBl[zlist>0];
            gzlist=gzlist[zlist>0];
            zlist=zlist[zlist>0];
        outdict = {'zlist':zlist, 'xlist':xlist,'ylist':ylist,'zlist':zlist,'rhog':rhog,'dPhidz':gz,'rhol':rhol,\
               'pthl':pthl, 'pturl':pturl, 'pkezl':pkezl, 'pcrl':pcrl, 'pBl':pBl,'gzlist':gzlist, 'pvzl':pvzl, 'volll':volll} 
        if outHI==1:
            outdict['rhohotl']=rhohotl
            outdict['rhocoldl']=rhocoldl
            outdict['rhoHIl']=rhoHIl
            outdict['pturhotl']=pturhotl
            outdict['pturcoldl']=pturcoldl
            outdict['pturHIl']=pturHIl
            outdict['pthHIl']=pthHIl
        return outdict
    elif horizontal==1:
        if withoutr>0:
            smallr=withoutr
        else:
            smallr=0.1
        rlist = np.linspace(smallr,withinr+dr,num=withinr/dr)
        gzl=rlist*0.0; rholl=rlist*0.0;
        pthll=rlist*0.0; pturll=rlist*0.0;
        if outHI==1:
            rhohotll=rlist*0.0;
            rhocoldll=rlist*0.0;
            rhoHIll=rlist*0.0;
            pturhotll=rlist*0.0;
            pturcoldll=rlist*0.0;
            pturHIll=rlist*0.0;
            pthHIll=rlist*0.0;   
        pcrll=rlist*0.0; pBll=rlist*0.0;
        volll=rlist*0.0; pkezll=rlist*0.0;
        pvzll=rlist*0.0;
        for i in range(len(xlist)):
            for j in range(len(ylist)):
                rxy = np.sqrt(xlist[i]*xlist[i]+ylist[j]*ylist[j])
                for k in range(len(zlist)):
                    for l in range(len(rlist)-1):
                        if ((rxy>rlist[l]) and (rxy<rlist[l+1]) and np.absolute(zlist[k])<zmax):
                            gzl[l]+=gzgrid[i,j,k]*volgrid[i,j,k]
                            rholl[l]+=rhogrid[i,j,k]*volgrid[i,j,k]
                            pthll[l]+=pthgrid[i,j,k]*volgrid[i,j,k]
                            pturll[l]+=pturgrid[i,j,k]*volgrid[i,j,k]
                            if outHI==1:
                                rhohotll[l]+=rhohotgrid[i,j,k]*volgrid[i,j,k]
                                rhocoldll[l]+=rhocoldgrid[i,j,k]*volgrid[i,j,k]
                                rhoHIll[l]+=rhoHIgrid[i,j,k]*volgrid[i,j,k] 
                                pturhotll[l]+=pturhotgrid[i,j,k]*volgrid[i,j,k]
                                pturcoldll[l]+=pturcoldgrid[i,j,k]*volgrid[i,j,k]
                                pturHIll[l]+=pturHIgrid[i,j,k]*volgrid[i,j,k]
                                pthHIll[l]+=pthHIgrid[i,j,k]*volgrid[i,j,k]
                            pkezll[l]+=pkezgrid[i,j,k]*volgrid[i,j,k]
                            pvzll[l]+=vzgrid[i,j,k]*vzgrid[i,j,k]*rhogrid[i,j,k]*volgrid[i,j,k]
                            if havecr>0:
                                pcrll[l]+=pcrgrid[i,j,k]*volgrid[i,j,k]
                            if haveB>0:
                                pBll[l]+=pBgrid[i,j,k]*volgrid[i,j,k]
                            volll[l]+=volgrid[i,j,k]
                            #print 'volll', volll
        gz=gzlist=gzl/volll
        rhol=rholl/volll
        if outHI==1:
            rhohotl=rhohotll/volll
            rhocoldl=rhocoldll/volll
            rhoHIl=rhoHIll/volll
            pturhotl=pturhotll/volll
            pturcoldl=pturcoldll/volll
            pturHIl=pturHIll/volll
            pthHIl=pthHIll/volll
        pthl=pthll/volll
        pturl=pturll/volll
        pkezl=pkezll/volll
        pvzl=pvzll/volll
        pcrl=pcrll/volll
        pBl=pBll/volll
        #print 'pturl, volll', pturl, volll
        rhog = gz*rhol
        rhol[~np.isfinite(rhol)] = 0; pthl[~np.isfinite(pthl)] = 0;
        pturl[~np.isfinite(pturl)] = 0;
        if outHI==1:
            rhohotl[~np.isfinite(rhohotl)] = 0; 
            rhocoldl[~np.isfinite(rhocoldl)] = 0;
            rhoHIl[~np.isfinite(rhoHIl)] = 0;
            pturhotl[~np.isfinite(pturhotl)] = 0;
            pturcoldl[~np.isfinite(pturcoldl)] = 0;
            pturHIl[~np.isfinite(pturHIl)] = 0;
            pthHIl[~np.isfinite(pthHIl)] = 0;
        pkezl[~np.isfinite(pkezl)] = 0;
        pvzl[~np.isfinite(pvzl)] = 0;
        if havecr>0:
            pcrl[~np.isfinite(pcrl)] = 0;
        if haveB>0:
            pBl[~np.isfinite(pBl)] = 0; 
        rhog[~np.isfinite(rhog)] = 0;
        outdict = {'rlist':rlist, 'xlist':xlist,'ylist':ylist,'zlist':zlist,'rhog':rhog,'dPhidz':gz,'rhol':rhol,\
               'pthl':pthl, 'pturl':pturl, 'pkezl':pkezl, 'pcrl':pcrl, 'pBl':pBl,'gzlist':gzlist, 'pvzl':pvzl, 'volll':volll}
        #print 'SSF pturl', pturl
        if outHI==1:
            outdict['rhohotl']=rhohotl
            outdict['rhocoldl']=rhocoldl
            outdict['rhoHIl']=rhoHIl
            outdict['pturhotl']=pturhotl
            outdict['pturcoldl']=pturcoldl
            outdict['pturHIl']=pturHIl
            outdict['pthHIl']=pthHIl
        return outdict