from stdmodandoption import *
import plot_setup as PS
import collections
from crdenfromdata import crdenfromdata

keystore=['plot1','plot2','plot3','plot4','plot5','plot6']

def crdenfromdata_testinput(subdict):
    startno=subdict['startno']
    Nsnap=subdict['Nsnap']
    snapsep=subdict['snapsep']
    wanted=subdict['wanted']
    dirneed=np.array(subdict['dirneed'])
    fmeat=subdict['fmeat']
    dirdict = collections.defaultdict(dict)
    dshape=dirneed.shape
    if dirneed.ndim==1:
        ncols=1
        nrows=dshape[0]
    else:
        ncols=dshape[1]
        nrows=dshape[0]
    keylist = dirneed
    if dirneed.ndim==1:
        enulist = list(enumerate(dirneed))
    else:
        enulist = list(np.ndenumerate(dirneed))
    for (index, runtodo) in enulist:
        dirdict[keylist[index]]=[runtodo]
    noplots=len(dirdict.keys())
    fig, gs = PS.setupgs(nrows=nrows, ncols=ncols)
    for (index, key) in enulist:
        items=dirdict[key]
        for j, runtodo in enumerate(items):
            ssdict = collections.defaultdict(dict)
            ssdict = subdict
            ssdict['runtodo'] = runtodo
            plotdict = crdenfromdata(ssdict)
            xlab = plotdict[wanted]['xlab'];
            ylab = plotdict[wanted]['ylab'];
            ptitle = plotdict[wanted]['ptitle']
            runtitle = plotdict[wanted]['runtitle']
            filename = plotdict[wanted]['filename'];
            linelist = plotdict[wanted]['linelist']
            for k,inkey in enumerate(linelist):
                xnl = plotdict[wanted]['xnl'][inkey];
                ynl = plotdict[wanted]['ynl'][inkey];
                labelneed = plotdict[wanted]['labelneed'];
                color = plotdict[wanted]['color'][inkey];
                lsn = plotdict[wanted]['lsn'][inkey];
                lw = plotdict[wanted]['lw'][inkey];
                marker = plotdict[wanted]['marker'][inkey];
                linelabel = plotdict[wanted]['linelab'][inkey];
                legendneed = 1
                ax = plt.subplot(gs[index])
                if dirneed.ndim==1:
                    if index==nrows-1: 
                        label = linelabel
                    else:
                        label = '_nolegend_'
                else:
                    if index[0]==nrows-1: 
                        label = linelabel
                    else:
                        label = '_nolegend_'
                ax.plot(xnl,ynl,label=label,lw=lw,ls=lsn,color=color,marker=marker)
        if dirneed.ndim==1:
            if index<nrows-1: 
                xlab=''
            if index>0:
                title=''
            else:
                title=runtitle
        else:
            if index[1]>0: 
                ylab=''
                ax.tick_params(labelleft='off')
            if index[0]<nrows-1: 
                xlab=''
                ax.tick_params(labelbottom='off')
            if index[0]>0:
                title=''
            else:
                title=runtitle            
        ax.text(0.8, 0.9, ptitle, horizontalalignment='center',
        verticalalignment='center', transform=ax.transAxes,fontsize=22)
        logx=0; logy=1;
        locneed='lower left'
        if wanted=='vr':
            logy=1; locneed='upper left'
        PS.miscsetup(ax,logx=logx,logy=logy,xlab=xlab,ylab=ylab,legendneed=legendneed,\
                    labfs=22,legfs=12,title=title,locneed=locneed)
        if wanted=='pr':
            ax.set_ylim([1e-14,5e-9])
        if wanted=='pz':
            ax.set_ylim([1e-17,5e-10])
        if wanted=='vr':
            ax.set_ylim([1.0,100.0])
        if wanted=='vz':
            ax.set_ylim([0.1,600.0])
        if wanted=='vturr':
            ax.set_ylim([1.0,500.0])
        
        if dirneed.ndim==1:        
            if not index==nrows-1:
                if not noplots==1:
                    ax.legend().set_visible(False)
        else:       
            if not index[0]==nrows-1:
                if not noplots==1:
                    ax.legend().set_visible(False)
    PS.finishsave(plt,filename,subplotadjust=0,tightbbox=1)
    return None
