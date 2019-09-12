from stdmodandoption import *
import plot_setup as PS
import collections
from kslaw import kslaw

keystore=['plot1','plot2','plot3','plot4','plot5','plot6']

def kslawtestinput(subdict):
    startno=subdict['startno']
    Nsnap=subdict['Nsnap']
    snapsep=subdict['snapsep']
    wanted=subdict['wanted']
    dirneed=subdict['dirneed']
    fmeat=subdict['fmeat']
    dirdict = collections.defaultdict(dict)
    keylist = np.array([])
    noplots=1
    fig, ax = PS.setupfig(nrows=noplots, ncols=1,sharex=True,sharey=False)
    for j, runtodo in enumerate(dirneed):
        ssdict = collections.defaultdict(dict)
        ssdict = subdict
        ssdict['runtodo'] = runtodo
        print 'runtodo', runtodo
        plotdict = kslaw(ssdict)
        xlab = plotdict[wanted]['xlab'];
        ylab = plotdict[wanted]['ylab'];
        ptitle = plotdict[wanted]['ptitle']
        filename = plotdict[wanted]['filename'];
        xnl = plotdict[wanted]['xnl'];
        ynl = plotdict[wanted]['ynl'];
        labelneed = plotdict[wanted]['labelneed'];
        color = plotdict[wanted]['color'];
        lsn = plotdict[wanted]['lsn'];
        lw = plotdict[wanted]['lw'];
        marker = plotdict[wanted]['marker'];
        legendneed = 1
        ax.plot(xnl,ynl,label=labelneed,lw=lw,ls=lsn,color=color,marker=marker, alpha=0.7)
    PS.miscsetup(ax,logx=1,logy=1,xlab=xlab,ylab=ylab,legendneed=0,labfs=22,legfs=12)
    PS.finishsave(plt,filename)
    return None