from stdmodandoption import *
import plot_setup as PS
import collections
from hgashalf_out import hgashalf_out

keystore=['plot1','plot2','plot3','plot4','plot5','plot6']

def hgashalftestinput(subdict):
    startno=subdict['startno']
    Nsnap=subdict['Nsnap']
    snapsep=subdict['snapsep']
    wanted=subdict['wanted']
    dirneed=subdict['dirneed']
    fmeat=subdict['fmeat']
    dirdict = collections.defaultdict(dict)
    keylist = np.array([])
    for ir, runtodo in enumerate(dirneed):
        keylist=np.append(keylist,keystore[ir])
        dirdict[keystore[ir]]=[runtodo]
    #noplots=len(dirdict.keys())
    noplots=1
    print 'noplots', noplots
    fig, ax = PS.setupfig(nrows=noplots, ncols=1,sharex=True,sharey=False)

    for i, key in enumerate(keylist):
        items=dirdict[key]
        for j, runtodo in enumerate(items):
            ssdict = collections.defaultdict(dict)
            ssdict = subdict
            ssdict['runtodo'] = runtodo
            print 'runtodo', runtodo
            plotdict = hgashalf_out(ssdict)
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
            linelabel = plotdict[wanted]['linelab'];
            legendneed = 1
            if i==0 and j==0: 
                label = linelabel
            else:
                label = '_nolegend_' 
            print 'xnl', xnl
            print 'ynl', ynl
            if noplots==1:
                ax.plot(xnl,ynl,label=labelneed,lw=lw,ls=lsn,color=color,marker=marker)
            else:
                ax[i].plot(xnl,ynl,label=labelneed,lw=lw,ls=lsn,color=color,marker=marker)
        if i<noplots-1: xlab='' 
        if i>1: legendneed=0
        if noplots==1:
            ax.set_title(ptitle, fontsize=16)
        else:
            ax[i].set_title(ptitle, fontsize=16)        
        #ax.text(0.25, 0.95, ptitle, horizontalalignment='center',
        #verticalalignment='center', transform=ax.transAxes,fontsize=22)            
        if noplots==1:
            PS.miscsetup(ax,logx=0,logy=1,xlab=xlab,ylab=ylab,legendneed=legendneed,labfs=22,legfs=12)
        else:
            PS.miscsetup(ax[i],logx=0,logy=1,xlab=xlab,ylab=ylab,legendneed=legendneed,labfs=22,legfs=12)        
    PS.finishsave(plt,filename)
    return None