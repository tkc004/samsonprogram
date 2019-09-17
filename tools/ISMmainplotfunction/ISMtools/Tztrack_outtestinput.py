from stdmodandoption import *
import plot_setup as PS
import collections
from Tztrack_out import Tztrack_out

keystore=['plot1','plot2','plot3','plot4','plot5','plot6']

def Tztrack_outtestinput(subdict):
    dirneed=subdict['dirneed']
    wanted=subdict['wanted']
    Tcut_t = subdict['Tcut_t'] #K
    highTcut_t = subdict['highTcut_t']
    hotwarmmode = subdict['hotwarmmode']
    trackstart =subdict['trackstart']
    if hotwarmmode==1:
        modedict = collections.defaultdict(dict)
        modelist=['hot','warm']
        modedict['warm']['Tcut_t']=1.0
        modedict['warm']['highTcut_t']=1.0e5
        modedict['hot']['Tcut_t']=1.0e5
        modedict['hot']['highTcut_t']=1.0e10
        modedict['warm']['lsn']='dashed'
        modedict['hot']['lsn']='dotted'
    else:
        modedict = collections.defaultdict(dict)
        modelist=['normal']
        modedict['normal']['Tcut_t']=Tcut_t
        modedict['normal']['highTcut_t']=highTcut_t  
    dirdict = collections.defaultdict(dict)
    keylist = np.array([])
    for ir, subdirneed in enumerate(dirneed):
        keylist=np.append(keylist,keystore[ir])
        dirdict[keystore[ir]]=subdirneed
    noplots=1
    fig, ax = PS.setupfig(nrows=noplots, ncols=1,sharex=True,sharey=False)
    for i, key in enumerate(keylist):
        items=dirdict[key]
        ssdict = collections.defaultdict(dict)
        ssdict = subdict
        for j, runtodo in enumerate(items):
            ssdict['runtodo'] = runtodo
            print 'runtodo', runtodo
            for labelcount, modekey in enumerate(modelist):
                ssdict['Tcut_t']=modedict[modekey]['Tcut_t']
                ssdict['highTcut_t']=modedict[modekey]['highTcut_t']
                print 'modekey', modekey
                print 'ssdict[Tcut_t]', ssdict['Tcut_t']
                print 'ssdict[highTcut_t]', ssdict['highTcut_t']
                plotdict = Tztrack_out(ssdict)
                xlab = plotdict[wanted]['xlab'];
                ylab = plotdict[wanted]['ylab'];
                ptitle = plotdict[wanted]['ptitle']
                filename = plotdict[wanted]['filename'];
                xpoints = plotdict[wanted]['xnl'];
                ypoints = plotdict[wanted]['ynl'];
                labelneed = plotdict[wanted]['labelneed'];
                color = plotdict[wanted]['color'];
                lsn = plotdict[wanted]['lsn'];
                if hotwarmmode==1:
                    lsn = modedict[modekey]['lsn']
                lw = plotdict[wanted]['lw'];
                marker = plotdict[wanted]['marker'];
                linelabel = plotdict[wanted]['linelab'];
                legendneed = plotdict[wanted]['legendneed']
                print 'xpoints,ypoints', xpoints,ypoints
                if labelcount==0:
                    ax.plot(xpoints,ypoints,label=labelneed,lw=lw,ls=lsn,color=color)
                else:
                    ax.plot(xpoints,ypoints,lw=lw,ls=lsn,color=color,marker=marker)
                ax.scatter(xpoints[1:-1],ypoints[1:-1],color=color)
                if labelcount==0 and j==0:
                    ax.scatter(xpoints[0],ypoints[0],color=color,marker='s',s=4*rcParams['lines.markersize']**2,label='Start') 
                    ax.scatter(xpoints[-1],ypoints[-1],color=color,marker='>',s=4*rcParams['lines.markersize']**2,label='End')
                else:
                    ax.scatter(xpoints[0],ypoints[0],color=color,marker='s',s=4*rcParams['lines.markersize']**2) 
                    ax.scatter(xpoints[-1],ypoints[-1],color=color,marker='>',s=4*rcParams['lines.markersize']**2) 
        logx=1
        if trackstart==1:
            logx=0
        PS.miscsetup(ax,logx=logx,logy=1,xlab=xlab,ylab=ylab,legendneed=legendneed,labfs=22,legfs=20)
    if hotwarmmode==1:
        plt.xlim([np.log10(5.0),2.0])
    if trackstart==1:
        plt.xlim([np.log10(0.5),2.0])
    PS.finishsave(plt,filename)
    return None