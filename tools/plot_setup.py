from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
from mpl_toolkits.axes_grid1 import AxesGrid
from textwrap import wrap
import matplotlib.gridspec as gridspec

def setupfig(nrows=1, ncols=1, sharex='col', sharey='row', squeeze=True
    ): 
    figsize= (7.*ncols,5.*nrows)
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex=sharex, sharey=sharey,figsize=figsize, squeeze=squeeze)
    #xinch = 10.*nrows; yinch = 4*ncols
    #fig.set_size_inches(xinch,yinch)
    return fig, ax

def setupgs(nrows=1, ncols=1):
    figsize= (7.*ncols,5.*nrows)
    #fig = plt.subplots(sharex=sharex,sharey=sharey,figsize=figsize)
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(nrows, ncols,
         #width_ratios=[1, 1, 1],
         wspace=0.0, hspace=0.0, top=0.95, bottom=0.15, left=0.17, right=0.845) 
    return fig,gs


def setupaxesgrid(nrows=1, ncols=1):
    fig = plt.figure(figsize=(10.*ncols,6.*nrows))
    gridag = AxesGrid(fig, 111,
                nrows_ncols = (nrows, ncols),
                axes_pad = 0.05,
                aspect = True,
                label_mode = "L",
                share_all = "yaxis",
                cbar_location="right",
                cbar_mode="single",
                cbar_size="7%",
                cbar_pad="0%")
    return fig, gridag

def finishaxesgrid(im,cblabel,gridag,plotit,plotupt,loc='lower left',addtitleinbox=0):
    cbar = gridag.cbar_axes[0].colorbar(im, extend='max')
    cbar.ax.set_ylabel(cblabel, rotation=270,labelpad=25)
    cbar.ax.tick_params()
    for ax, im_title, up_title in zip(gridag, plotit, plotupt):
        t = add_inner_title(ax, "\n".join(wrap(im_title,17)), loc=loc)
        t.patch.set_alpha(0.5)
        if addtitleinbox==1:
            t1 = add_inner_title(ax, up_title, loc='upper left')
            t1.patch.set_alpha(0.5)            
        else:
            ax.set_title(up_title)


def add_inner_title(ax, title, loc, size=None, **kwargs):
    from matplotlib.offsetbox import AnchoredText
    from matplotlib.patheffects import withStroke
    if size is None:
        size = dict(size=plt.rcParams['legend.fontsize'])
    at = AnchoredText(title, loc=loc, prop=size,
                      pad=0., borderpad=0.5,
                      frameon=False, **kwargs)
    ax.add_artist(at)
    at.txt._text.set_path_effects([withStroke(foreground="w", linewidth=3)])
    return at


def miscsetup(ax,logx=0,logy=0,xlab='',ylab='',title='',legendneed=0,titfs=24,labfs=22,legfs=18,
    legloc='best',ncol=1,locneed='best'):
    if logx==1: ax.set_xscale('log')
    if logy==1: ax.set_yscale('log')
    if xlab: ax.set_xlabel(xlab,fontsize=labfs)
    if ylab: ax.set_ylabel(ylab,fontsize=labfs)
    if title: ax.set_title(title,fontsize=titfs)
    if legendneed==1: ax.legend(loc=legloc,ncol=ncol)
    if legfs>0. and legendneed==1: ax.legend(fontsize=legfs,loc=locneed)
    ax.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=labfs)
        
def finishsave(plt,filename,clear=1,ptitle='',subplotadjust=1,nosubplotspace=0,tightbbox=0):
    print 'filename', filename
    if subplotadjust==1: plt.subplots_adjust(left  = 0.18, bottom=0.15, right=0.95, top=0.9)
    if nosubplotspace==1: plt.subplots_adjust(wspace=0, hspace=0)
    #if nosubplotspace==1: plt.subplots_adjust(left  = 0.0, bottom=0.0, right=1.0, top=1.0)
    if tightbbox==1:
        plt.savefig(filename, bbox_inches='tight')
    else:
        plt.savefig(filename)
    if clear==1:
        plt.clf()
        
def dopointplot(plotdict,plotneed,markerl,colorl,labell,neederr=0):
    fig, ax = setupfig()
    xnl = plotdict[plotneed]['xnl'];
    ynl = plotdict[plotneed]['ynl'];
    xlab = plotdict[plotneed]['xlab'];
    ylab = plotdict[plotneed]['ylab'];
    selectl = plotdict[plotneed]['selectl'];
    filename=plotdict[plotneed]['filename'];
    if neederr==1:
        yld = plotdict[plotneed]['yld'];
        ylu = plotdict[plotneed]['ylu'];
    for i in range(len(labell)):
        if not xnl[selectl==i]: continue
        if neederr==1:
            ax.errorbar(xnl[selectl==i], ynl[selectl==i], yerr=[yld[selectl==i],ylu[selectl==i]],
                        color=colorl[i],fmt=markerl[i],markersize=7,label=labell[i])
        else:
            ax.plot(xnl[selectl==i],ynl[selectl==i],marker=markerl[i],ls='none',color=colorl[i],label=labell[i])
    miscsetup(ax,logx=1,logy=1,xlab=xlab,ylab=ylab,legendneed=1,labfs=22,legfs=12)
    finishsave(plt,filename)
