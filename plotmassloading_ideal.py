import matplotlib as mpl
mpl.use('Agg')
from stdmodandoption import * 
import plot_setup as PS
import calmassloading as CML
import collections
#dirneed=['m12icr_b_70', 'm12imhdcv', 'm12icr_700', 'f553']
dirneed=[
'bwsmclr','bwmwmr','bwsbclr',
#'bwsmclrmhd', 'bwmwmrmhd', 'bwsbclrmhd',
#'bwsmclrdc28', 'bwmwmrdc28', 'bwsbclrdc28',
#'bwsmclrdc29', 'bwmwmrdc29', 'bwsbclrdc29',
#'m11dmhdcv',\
#'m11bcr_b_70','m11bmhdcv','m11bcr_700',\
#'476',\
#'bwsmclr','bwsmclrdc0','bwsmclrdc27','bwsmclrdc28','bwsmclrdc29','bwsmclrstr','bwsmclrdc28mhd'
#'m09','m10','m11','m12v','B1','m12qq','383','476',\
#'fm10qmd','fm10vmd','fm11dmd','fm11emd', 'fm11hmd','fm11imd', 'fm11qmd','fm12fmd','fm12imd',\
#'m10qcr_b_70','m10vcr_b_70', 'm11bcr_b_70','m11dcr_b_70','m11vcr_b_70','m11v1cr_b_70',\
#'m12icr_b_70','m11fcr_b_70','m11hcr_b_70','m11gcr_b_70',\
#'m12fmhdcv','m12mmhdcv',\
#'m10qmhdcv','m10vmhdcv', 'm11bmhdcv','m11dmhdcv','m12imhdcv','m11fmhdcv','m11hmhdcv','m11gmhdcv',\
#'m12fcr_700','m12mcr_700',\
#'m10vcr_700', 'm11bcr_700','m11dcr_700','m11fcr_700','m11gcr_700','m11hcr_700','m12icr_700',
#'fm12m','fm12b','fm10q','fm11d','f553','f573','f476','f383','fm11q','fm11v','fm11v1','fm11v2','f1146','f46','f61',
]
galcen=0
hubble=0.702
#Nsnaplist=range(300,441,10)
#Nsnaplist=range(590,601,10)
#Nsnaplist=range(500,601,10)
Nsnaplist=range(400,501,5)
#Nsnaplist=[570,580,590,600]
fmeat='test'
#fmeat='md'
#fmeat='all'
#fmeat='FIRE1'
#fmeat='FIRE2'
#fmeat='dc28mhd'
#fmeat='hydro'
#fmeat='mhdcv'
#fmeat='cr_b_70'
#fmeat='cr_700'
#fmeat='crhydroave'
#fmeat='hydrosn600'
#Nsnaplist=[580,590,600]
neederr=1
avesnap=1
figdir='figures/Massloading/'

plotlist = ['mlmv','mlms','msmv','SFRMs'] #mlmv: massloading Mvir mlms: massloading Mstar; #msmv: mstar mvir #SFRMs: SFR Mstar

mldata= CML.calmassloading(dirneed,Nsnaplist,avesnap=1)
dcl = mldata['dcl']; sfrl=mldata['sfrl']
mvl=mldata['mvl']; msl = mldata['msl']; mll = mldata['mll']; mlld=mldata['mlld']; mllu=mldata['mllu'];
markerl=['o','s','^','d','*']
colorl=['y','b','g','k','brown']
labell=['Hydro','MHD',r'$\kappa$=3e28',r'$\kappa$=3e29','metal diffusion']

'''
build up dictionary for plots
'''

plotdict = collections.defaultdict(dict)
type(plotdict)
for plotneed in plotlist:
    if plotneed == 'mlmv':
        plotdict[plotneed]['xlab'] = r'$M_{\rm vir}\,[{\rm M}_{\odot}]$'
        plotdict[plotneed]['xnl'] = mvl
        plotdict[plotneed]['ylab'] = r'$\eta=\dot{M}_{\rm w}/\dot{M}_{*}$'
        plotdict[plotneed]['ynl'] = mll
        plotdict[plotneed]['selectl'] = dcl
        plotdict[plotneed]['filename'] = figdir+'/Massloading_'+fmeat+'.pdf'
        if neederr==1:
            plotdict[plotneed]['yld'] = mlld
            plotdict[plotneed]['ylu'] = mllu
    if plotneed == 'mlms':
        plotdict[plotneed]['xlab'] = r'$M_{\rm *}\,[{\rm M}_{\odot}]$'
        plotdict[plotneed]['xnl'] = msl
        plotdict[plotneed]['ylab'] = r'$\eta=\dot{M}_{\rm w}/\dot{M}_{*}$'
        plotdict[plotneed]['ynl'] = mll
        plotdict[plotneed]['selectl'] = dcl
        plotdict[plotneed]['filename'] = figdir+'/MlMs_'+fmeat+'.pdf'
        if neederr==1:
            plotdict[plotneed]['yld'] = mlld
            plotdict[plotneed]['ylu'] = mllu
    if plotneed == 'msmv':
        plotdict[plotneed]['xlab'] = r'$M_{\rm vir}\,[{\rm M}_{\odot}]$'
        plotdict[plotneed]['xnl'] = mvl
        plotdict[plotneed]['ylab'] = r'$M_{*}\,[{\rm M}_{\odot}]$'
        plotdict[plotneed]['ynl'] = msl
        plotdict[plotneed]['selectl'] = dcl
        plotdict[plotneed]['filename'] = figdir+'/MsMv_'+fmeat+'.pdf'
    if plotneed == 'SFRMs':
        plotdict[plotneed]['xlab'] = r'$M_{*}\,[{\rm M}_{\odot}]$'
        plotdict[plotneed]['xnl'] = msl
        plotdict[plotneed]['ylab'] = r'SFR $[{\rm M}_{\odot}]$'
        plotdict[plotneed]['ynl'] = sfrl
        plotdict[plotneed]['selectl'] = dcl
        plotdict[plotneed]['filename'] = figdir+'/SFRMs_'+fmeat+'.pdf' 


def doplot(plotdict,plotneed,markerl,colorl,labell,neederr=0):
    fig, ax = PS.setupfig()
    xnl = plotdict[plotneed]['xnl'];
    ynl = plotdict[plotneed]['ynl'];
    xlab = plotdict[plotneed]['xlab'];
    ylab = plotdict[plotneed]['ylab'];
    selectl = plotdict[plotneed]['selectl'];
    filename=plotdict[plotneed]['filename'];
    if neederr==1 and (plotneed=='mlmv' or plotneed=='mlms'):
        yld = plotdict[plotneed]['yld'];
        ylu = plotdict[plotneed]['ylu'];
    for i in range(len(labell)):
        if neederr==1 and (plotneed=='mlmv' or plotneed=='mlms'):
            ax.errorbar(xnl[selectl==i], ynl[selectl==i], yerr=[yld[selectl==i],ylu[selectl==i]]
                        ,color=colorl[i],fmt=markerl[i],markersize=7,label=labell[i])
        else:

            ax.plot(xnl[selectl==i],ynl[selectl==i],marker=markerl[i],ls='none',color=colorl[i],label=labell[i])
    PS.miscsetup(ax,logx=1,logy=1,xlab=xlab,ylab=ylab,legendneed=1,labfs=22,legfs=12)
    PS.finishsave(plt,filename)

    
for plotneed in plotlist:
        doplot(plotdict,plotneed,markerl,colorl,labell,neederr=neederr)
