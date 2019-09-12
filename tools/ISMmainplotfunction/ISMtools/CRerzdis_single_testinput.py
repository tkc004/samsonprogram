from stdmodandoption import *
import plot_setup as PS
import collections
from CR_rz_energy_distribution_single import CR_rz_energy_distribution_single
import matplotlib.lines as mlines

keystore=['plot1','plot2','plot3','plot4','plot5','plot6']

def CRerzdis_single_testinput(subdict):
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
    noplots=1
    print 'noplots', noplots
    fig, ax = PS.setupfig(nrows=noplots, ncols=1,sharex=True,sharey=False)
    linelablist=[]
    lxlist=[]
    for i, key in enumerate(keylist): #i for different runs
        items=dirdict[key]
        for j, runtodo in enumerate(items): #j for different CR models
            ssdict = collections.defaultdict(dict)
            ssdict = subdict
            ssdict['runtodo'] = runtodo
            print 'runtodo', runtodo
            plotdict = CR_rz_energy_distribution_single(ssdict)
            xlab = plotdict[wanted]['xlab'];
            ylab = plotdict[wanted]['ylab'];
            ptitle = plotdict[wanted]['ptitle']
            filename = plotdict[wanted]['filename'];
            print 'filename', filename
            for k,inkey in enumerate(plotdict[wanted]['xnl']): #k for different energies
                xnl = plotdict[wanted]['xnl'][inkey];
                ynl = plotdict[wanted]['ynl'][inkey];
                labelneed = plotdict[wanted]['labelneed']; #CR models
                color = plotdict[wanted]['color'];
                lsn = plotdict[wanted]['lsn'];
                lw = plotdict[wanted]['lw'][inkey];
                marker = plotdict[wanted]['marker'];
                linelabel = plotdict[wanted]['linelab']; #runlabel
                legendneed = 1
                if j==1: 
                    label = linelabel
                else:
                    label = '_nolegend_' 
                print 'j,label', j,label
                if noplots==1:
                    ax.plot(xnl,ynl,label=label,lw=lw,ls=lsn,color=color,marker=marker)
                else:
                    ax[i].plot(xnl,ynl,label=label,lw=lw,ls=lsn,color=color,marker=marker)
                if i==0:
                    lxlist.append(mlines.Line2D([], [], color=color, ls='solid',label=labelneed))
                    linelablist = np.append(linelablist,labelneed)
        if i<noplots-1: xlab=''
        legend1 = ax.legend(lxlist, linelablist, loc='lower left',fontsize=12,ncol=1)
        ax.add_artist(legend1)
        if noplots==1:
            PS.miscsetup(ax,logx=0,logy=1,xlab=xlab,ylab=ylab,legendneed=legendneed,labfs=22,legfs=12,locneed='upper right')
        else:
            PS.miscsetup(ax[i],logx=0,logy=1,xlab=xlab,ylab=ylab,legendneed=legendneed,labfs=22,legfs=12,locneed='upper right')        
    PS.finishsave(plt,filename)
    return None