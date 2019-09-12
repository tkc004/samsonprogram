from stdmodandoption import *
import collections
import gastracklib as GTL

def Tztrack_out(ssdict):
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
    Tcut_t = ssdict['Tcut_t'] #K
    highTcut_t = ssdict['highTcut_t']
    vcut_t= ssdict['vcut_t'] #outflow velocity cut
    vhcut_t = ssdict['vhcut_t'] #upper outflow velocity cut
    withinr_t = ssdict['withinr_t']
    zup_t = ssdict['zup_t']
    zdown_t = ssdict['zdown_t']
    zup = ssdict['zup']
    zdown = ssdict['zdown']
    trackstart =ssdict['trackstart']
    print 'runtodo', runtodo
    if wanted=='Tztrack' or wanted=='rhoztrack' or wanted=='vzztrack' or wanted=='Tvztrack':
        rcParams['figure.figsize'] = 5, 4
        from matplotlib.patches import FancyArrowPatch 
        def arrow(x,y,ax,n,color):
            d = len(x)//(n+1)    
            ind = np.arange(d,len(x),d)
            for i in ind:
                ar = FancyArrowPatch ((x[i-1],y[i-1]),(x[i],y[i]), 
                      arrowstyle='->', mutation_scale=20,color=color)
            ax.add_patch(ar)
        userad=0
        snaptrack = Nsnap
        if trackstart==1:
            snaptrack = startno
        if wanted=='Tztrack':
                extent = [np.log10(zdown),np.log10(zup),1.0,6.5]
                #extent = [zdown,zup,1.0,6.5]
        #try:
        info = SSF.outdirname(runtodo, Nsnap)
        runtitle=info['runtitle']
        ptitle=title
        ymin=extent[2]
        ymax = extent[3]
        xpoints=[]; ypoints=[];
        #ax.set_aspect("equal")
        for i in range(startno,Nsnap+1,snapsep):
            info = SSF.outdirname(runtodo, i)
            color=info['color']
            dclabel=info['dclabel']
            haveB=info['haveB']
            newlabel=info['newlabel']
            labelneed=newlabel
            if haveB==1:
                lsn='dashed'
            else:
                lsn='solid'
            data = GTL.gastrack(runtodo,snaptrack,i,Tcut_t=Tcut_t,highTcut_t=highTcut_t,\
                vcut_t=vcut_t,vhcut_t=vhcut_t,withinr_t=withinr_t,zup_t=zup_t,zdown_t=zdown_t,userad=userad)
            Gmass = data['Gmass']*1e10  
            if wanted=='Tztrack':
                TrueTemp = data['TrueTemp']
                partZ = np.absolute(data['partZ'])
                xm = np.average(partZ,weights=Gmass)
                xmed = np.log10(xm)
                ym = np.average(TrueTemp,weights=Gmass)
                ymed = np.log10(ym)
            xpoints = np.append(xpoints,xmed)
            ypoints = np.append(ypoints,ymed)
    if wanted == 'Tztrack':
        plotdict[wanted]['xlab'] = r"$\mathrm{Log}_{10} (z {\rm [kpc]})$"
        plotdict[wanted]['xnl'] = xpoints
        plotdict[wanted]['ylab'] = r"$\mathrm{Log}_{10} (T {\rm [K]})$"
        plotdict[wanted]['ynl'] = ypoints
        plotdict[wanted]['runtodo'] = runtodo
        plotdict[wanted]['labelneed'] = labelneed
        plotdict[wanted]['linelab']= r'$\rho g$'
        plotdict[wanted]['lsn'] = lsn
        plotdict[wanted]['lw'] = 2
        plotdict[wanted]['marker'] = 'o'
        plotdict[wanted]['color'] = color
        plotdict[wanted]['runtitle'] = runtitle
        plotdict[wanted]['ptitle'] = ptitle
        filename=plotloc+'/CRplot/Tztrack/Tztrack_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+fmeat+'.pdf'
        plotdict[wanted]['filename'] = filename
        plotdict[wanted]['legendneed'] = 1                
    return plotdict
