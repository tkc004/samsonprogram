import matplotlib
matplotlib.use('Agg')
from stdmodandoption import *
import plot_setup as PS
import collections
import gastracklib as GTL




def calTz(runtodo,wanted,startno,Nsnap,snapsep,fmeat,trackstart=0):
    title='MW'
    titleneed=title
    needlog=0
    dclabelneed=1
    correctIa=0
    useM1=1
    the_prefix='snapshot'
    the_suffix='.hdf5'
    withinr=15.0
    nogrid=40
    maxlength=10.0
    med = -0.1 
    wholeboxgas=1
    diskgas=1
    Rfrac = 0.5
    nosum=0
    xaxis_snapno=0
    nested_dict = lambda: collections.defaultdict(nested_dict)
    plotdict = nested_dict()
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
        Tcut_t=1e5 #K
        highTcut_t = 1e10
        vcut_t=2e2 #outflow velocity cut
        vhcut_t = 1e10 #upper outflow velocity cut
        withinr_t=10.0
        zup_t=2.0
        zdown_t=1.0
        zup=100.0
        zdown=0.1
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
            #print 'zup_t', zup_t
            #print 'zdown_t', zdown_t
            #print 'withinr_t', withinr_t
            data = GTL.gastrack(runtodo,snaptrack,i,Tcut_t=Tcut_t,highTcut_t=highTcut_t,\
                vcut_t=vcut_t,vhcut_t=vhcut_t,withinr_t=withinr_t,zup_t=zup_t,zdown_t=zdown_t,userad=userad)
            Gmass = data['Gmass']*1e10  
            if wanted=='Tztrack':
                TrueTemp = data['TrueTemp']
                partZ = np.absolute(data['partZ'])
                print 'partZ, np.amax(partZ), np.amin(partZ)', partZ, np.amax(partZ), np.amin(partZ)
                xm = np.average(partZ,weights=Gmass)
                xmed = np.log10(xm)
                ym = np.average(TrueTemp,weights=Gmass)
                ymed = np.log10(ym)
                print 'xm, ym', xm, ym
            xpoints = np.append(xpoints,xmed)
            ypoints = np.append(ypoints,ymed)
    print 'xpoints', xpoints
    print 'ypoints', ypoints
    if wanted == 'Tztrack':
        totalname = plotloc+'/CRplot/Tztrack/Tztrack_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    if wanted=='Tztrack':
        plotdict[wanted]['xlab'] = r"$\mathrm{Log}_{10} (z {\rm [kpc]})$"
        plotdict[wanted]['xnl'] = xpoints
        plotdict[wanted]['ylab'] = r"$\mathrm{Log}_{10}\;(T\;{\rm [K]})$"
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
        filename=plotloc+'/CRplot/Tztrack/Tztrack_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        plotdict[wanted]['filename'] = filename
        plotdict[wanted]['legendneed'] = 0                
        if runtitle=='SMC': plotdict[wanted]['legendneed'] = 1
    return plotdict

                
def main():
    #startno=450
    #Nsnap=500
    startno=550
    Nsnap=600
    snapsep=50
    wanted='Tztrack'
    trackstart=1
    dirdict = collections.defaultdict(dict)
    keylist=['plot1']
    dirdict['plot1']=['m12imhdcvhr']
    #dirdict['plot1']=['m12icr_700hr']
    noplots=len(dirdict.keys())
    fig, ax = PS.setupfig(nrows=noplots, ncols=1,sharex=True,sharey=False)
    fmeat=''
    for i, key in enumerate(keylist):
        items=dirdict[key]
        for j, runtodo in enumerate(items):
            print 'runtodo', runtodo
            plotdict = calTz(runtodo,wanted,startno,Nsnap,snapsep,fmeat,trackstart=trackstart)
            xlab = plotdict[wanted]['xlab'];
            ylab = plotdict[wanted]['ylab'];
            ptitle = plotdict[wanted]['ptitle']
            filename = plotdict[wanted]['filename'];
            xpoints = plotdict[wanted]['xnl'];
            ypoints = plotdict[wanted]['ynl'];
            labelneed = plotdict[wanted]['labelneed'];
            color = plotdict[wanted]['color'];
            lsn = plotdict[wanted]['lsn'];
            lw = plotdict[wanted]['lw'];
            marker = plotdict[wanted]['marker'];
            linelabel = plotdict[wanted]['linelab'];
            legendneed = plotdict[wanted]['legendneed']
            print 'lw', lw
            print 'lsn', lsn
            print 'color', color
            print 'marker', marker
            if noplots==1:
                ax.plot(xpoints,ypoints,label=labelneed,lw=lw,ls=lsn,color=color,marker=marker)
            else:
                ax[i].plot(xpoints,ypoints,label=labelneed,lw=lw,ls=lsn,color=color,marker=marker)
            ax.scatter(xpoints[1:-1],ypoints[1:-1],color=color)
            ax.scatter(xpoints[0],ypoints[0],color=color,marker='s',s=4*rcParams['lines.markersize']**2) 
            ax.scatter(xpoints[-1],ypoints[-1],color=color,marker='>',s=4*rcParams['lines.markersize']**2)            
        if i<noplots-1: xlab=''
        if noplots==1:
            ax.text(0.25, 0.95, ptitle, horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes,fontsize=22)
            PS.miscsetup(ax,logx=1,logy=1,xlab=xlab,ylab=ylab,legendneed=legendneed,labfs=22,legfs=12)
        else:
            ax[i].text(0.25, 0.95, ptitle, horizontalalignment='center',
            verticalalignment='center', transform=ax[i].transAxes,fontsize=22)
            PS.miscsetup(ax[i],logx=1,logy=1,xlab=xlab,ylab=ylab,legendneed=legendneed,labfs=22,legfs=12)
    PS.finishsave(plt,filename)


        
if __name__ == '__main__':
    main()
