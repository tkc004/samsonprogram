from stdmodandoption import *
import plot_setup as PS
import collections
from CR_rz_energy_distribution_volave import CR_rz_energy_distribution_volave

keystore=['plot1','plot2','plot3','plot4','plot5','plot6']

def CRerzdis_volave_testinput(subdict):
    startno=subdict['startno']
    Nsnap=subdict['Nsnap']
    snapsep=subdict['snapsep']
    wanted=subdict['wanted']
    dirneed=subdict['dirneed']
    fmeat=subdict['fmeat']
    title=subdict['title']
    dirdict = collections.defaultdict(dict)
    keylist = np.array([])
    for ir, runtodo in enumerate(dirneed):
        keylist=np.append(keylist,keystore[ir])
        dirdict[keystore[ir]]=[runtodo]
    noplots=len(dirdict.keys())
    print 'noplots', noplots
    fig, ax = PS.setupfig(nrows=noplots, ncols=1,sharex=True,sharey=False)
    for i, key in enumerate(keylist):
        items=dirdict[key]
        for j, runtodo in enumerate(items):
            ssdict = collections.defaultdict(dict)
            ssdict = subdict
            ssdict['runtodo'] = runtodo
            print 'runtodo', runtodo
            plotdict = CR_rz_energy_distribution_volave(ssdict)
            xlab = plotdict[wanted]['xlab'];
            ylab = plotdict[wanted]['ylab'];
            ptitle = plotdict[wanted]['ptitle']
            filename = plotdict[wanted]['filename'];
            print 'filename', filename
            for k,inkey in enumerate(plotdict[wanted]['xnl']):
                xnl = plotdict[wanted]['xnl'][inkey];
                ynl = plotdict[wanted]['ynl'][inkey];
                labelneed = plotdict[wanted]['labelneed'];
                color = plotdict[wanted]['color'];
                lsn = plotdict[wanted]['lsn'][inkey];
                lw = plotdict[wanted]['lw'][inkey];
                marker = plotdict[wanted]['marker'][inkey];
                linelabel = plotdict[wanted]['linelab'][inkey];
                legendneed = 1
                if i==noplots-1 and j==0: 
                    label = linelabel
                else:
                    label = '_nolegend_' 
                if noplots==1:
                    ax.plot(xnl,ynl,label=label,lw=lw,ls=lsn,color=color,marker=marker)
                else:
                    ax[i].plot(xnl,ynl,label=label,lw=lw,ls=lsn,color=color,marker=marker)
        if i<noplots-1: xlab=''
        if i>0: title=''
        if noplots==1:
            ax.text(0.25, 0.1, ptitle, horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes,fontsize=22)            
        else:
            ax[i].text(0.25, 0.1, ptitle, horizontalalignment='center',
            verticalalignment='center', transform=ax[i].transAxes,fontsize=22)
        if noplots==1:
            PS.miscsetup(ax,logx=0,logy=1,xlab=xlab,ylab=ylab,legendneed=legendneed,labfs=22,legfs=12,title=title)
        else:
            PS.miscsetup(ax[i],logx=0,logy=1,xlab=xlab,ylab=ylab,legendneed=legendneed,labfs=22,legfs=12,title=title)
        if not (i==noplots-1): 
            if noplots>1:
                ax[i].legend().set_visible(False)
    PS.finishsave(plt,filename)
    return None