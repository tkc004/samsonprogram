from stdmodandoption import *
import plot_setup as PS
import collections
from gasdenThis_out import gasdenThis_out
import matplotlib.lines as mlines
keystore=['plot1','plot2','plot3','plot4','plot5','plot6']

def gasdenThistestinput(subdict):
    startno=subdict['startno']
    Nsnap=subdict['Nsnap']
    snapsep=subdict['snapsep']
    wanted=subdict['wanted']
    dirneed=np.array(subdict['dirneed'])
    fmeat=subdict['fmeat']
    title=subdict['title']
    dirdict = collections.defaultdict(dict)
    keylist = np.array([])
    for ir, subdirneed in enumerate(dirneed):
        keylist=np.append(keylist,keystore[ir])
        if dirneed.ndim>1:
            dirdict[keystore[ir]]=subdirneed
        else:
            dirdict[keystore[ir]]=[subdirneed]
    if dirneed.ndim>1:
        noplots=len(keylist)
    else:
        noplots=1
    print 'noplots', noplots
    fig, ax = PS.setupfig(nrows=noplots, ncols=1,sharex=True,sharey=False)
    linelablist=[]
    lxlist=[]
    for i, key in enumerate(keylist):
        items=dirdict[key]
        for j, runtodo in enumerate(items):
            ssdict = collections.defaultdict(dict)
            ssdict = subdict
            ssdict['runtodo'] = runtodo
            print 'runtodo', runtodo
            if dirneed.ndim>1: 
                ssdict['fmeat']=fmeat[-1]
                ssdict['title']=title[i]
            plotdict = gasdenThis_out(ssdict)
            xlab = plotdict[wanted]['xlab'];
            ylab = plotdict[wanted]['ylab'];
            ptitle = plotdict[wanted]['ptitle']
            filename = plotdict[wanted]['filename'];
            for k,inkey in enumerate(plotdict[wanted]['xnl']):
                xnl = plotdict[wanted]['xnl'][inkey];
                ynl = plotdict[wanted]['ynl'][inkey];
                labelneed = plotdict[wanted]['labelneed'];
                color = plotdict[wanted]['color'];
                lsn = plotdict[wanted]['lsn'][inkey];
                lw = plotdict[wanted]['lw'];
                marker = plotdict[wanted]['marker'];
                linelabel = plotdict[wanted]['linelab'][inkey];
                legendneed = 1
                if i==0 and j==0: 
                    label = linelabel
                else:
                    label = '_nolegend_' 
                print 'xnl', xnl
                print 'ynl', ynl
                if noplots==1:
                    ax.plot(xnl,ynl,label=label,lw=lw,ls=lsn,color=color,marker=marker)
                else:
                    ax[i].plot(xnl,ynl,label=label,lw=lw,ls=lsn,color=color,marker=marker)
                if k==0 and i==noplots-1:
                    lxlist.append(mlines.Line2D([], [], color=color, ls=lsn,label=labelneed))
                    linelablist = np.append(linelablist,labelneed)
        if i<noplots-1: xlab=''
        if i>1: legendneed=0
        if noplots==1:
            ax.set_title(ptitle, fontsize=22)
            legend1 = ax.legend(lxlist, linelablist, loc='lower left',fontsize=12,ncol=1)
            ax.add_artist(legend1)
        else:
            ax[i].set_title(ptitle, fontsize=22)
            if i==noplots-1:
                legend1 = ax[i].legend(lxlist, linelablist, loc='lower left',fontsize=12,ncol=1)
                ax[i].add_artist(legend1)
        if noplots==1:
            PS.miscsetup(ax,logx=0,logy=1,xlab=xlab,ylab=ylab,legendneed=legendneed,labfs=22,legfs=12)
        else:
            PS.miscsetup(ax[i],logx=0,logy=1,xlab=xlab,ylab=ylab,legendneed=legendneed,labfs=22,legfs=12)
        if not (i==0 or i==noplots-1):
            if noplots==1:
                ax.legend().set_visible(False)
            else:
                ax[i].legend().set_visible(False)
    #plt.legend(loc='upper right',fontsize=12)
    PS.finishsave(plt,filename)
    return None