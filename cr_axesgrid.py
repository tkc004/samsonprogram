from samson_const import *
import matplotlib as mpl
mpl.use('Agg')
from readsnap_cr import readsnapcr
import Sasha_functions as SF
import graphics_library as GL
import gas_temperature as GT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from samson_functions import *
from matplotlib import rcParams
import pylab 
from textwrap import wrap
from mpl_toolkits.axes_grid1 import AxesGrid
#rcParams['figure.figsize'] = 5, 5
rcParams['figure.figsize'] = 10, 5
rcParams['legend.fontsize']=16
rcParams['font.size']=18
rcParams['font.family']='serif'
rcParams['text.usetex']=True
rcParams.update({'figure.autolayout': True})
import matplotlib.patches as patches
rcParams['axes.linewidth'] = 2
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True
rcParams['axes.unicode_minus']=False
colortable = [ 'b', 'g', 'r']

dirneed=['bwmwmr','bwmwmrdc27','bwmwmrdc28','bwmwmrdc29']
wanted='Tvz'
startno=450
Nsnap=500
snapsep=10
title=''
nrows=1
ncols=len(dirneed)
#if len(dirneed)>3:
#        nrows = 2
#        ncols = 3
#if len(dirneed)>6:
#        nrows = 3
#        ncols = 3

print 'nrows, ncols', nrows, ncols

fig = plt.figure()
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

plotupt=[]
plotit=[]


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

for idir, runtodo in enumerate(dirneed):
    if wanted=='Tz' or wanted=='rhoz' or wanted=='vzz' or wanted=='Tvz' or wanted=='pcrpth':
        newlabelneed=1
        from crtestfunction import *
        Tcut=-1.0 #K
        highTcut = 1e10
        vcut=0.0 #outflow velocity cut
        vhcut = 1e10 #upper outflow velocity cut
        withinr=20.0
        zup=30.0
        zdown=20.0
        trackgas=0
        tlabel=''
        userad=0
        if trackgas==1:
            snaptrack = Nsnap
            Tcut_t=1.0 #K
            highTcut_t = 1e10
            vcut_t=0.0 #outflow velocity cut
            vhcut_t = 1e10 #upper outflow velocity cut
            withinr_t=20.0
            zup_t=30.0
            zdown_t=20.0
            tlabel='track'
        if userad==1:
            tlabel+='rad'
        if wanted=='Tz':
            extent = [zdown,zup,1,9]
        if wanted=='rhoz':
            extent = [zdown,zup,-7,4]
        if wanted=='vzz':
            extent = [zdown,zup,0,3.2]
        if wanted=='Tvz':
            extent = [0,3.2,1,9]
            #extent = [0,9,0,9]
        if wanted=='pcrpth':
            extent = [-20,-11,-20,-11]
        nobin=51
        needcontour=1
        cbcolor = 'hot'
        info=outdirname(runtodo, Nsnap)
        cosmo=info['cosmo']
        runtitle=info['runtitle']
        havecr=info['havecr']
        dclabel=info['dclabel']
        newlabel=info['newlabel']
        if wanted=='pcrpth' and havecr==0:
            continue
        ptitle=title
        if runtitle=='SMC':
            ptitle='Dwarf'
        elif runtitle=='SBC':
            ptitle='Starburst'
        elif runtitle=='MW':
            ptitle=r'$L\star$ Galaxy'
        ptitle=''
        labelneed=dclabel
        if newlabelneed==1:
            labelneed="\n".join(wrap(newlabel,17))
        plotupt=np.append(plotupt,ptitle)
        plotit =np.append(plotit,labelneed)
        #if cosmo==1:
    #               userad=1
        Hadd = np.zeros((nobin-1,nobin-1))
        if trackgas==1:
            data = gaswindphase(runtodo,snaptrack,Tcut=Tcut_t,highTcut=highTcut_t,\
            vcut=vcut_t,vhcut=vhcut_t,withinr=withinr_t,zup=zup_t,zdown=zdown_t,userad=userad)
            Gid_t = data['Gid']
        nooftimes=0
        for i in range(startno,Nsnap,snapsep):
            if trackgas==1:
                vcut = -1e10 #if we track gas particles, we should consider all velocities
                if nooftimes==1:
                    continue
            data = gaswindphase(runtodo,i,Tcut=Tcut,highTcut=highTcut,\
            vcut=vcut,vhcut=vhcut,withinr=withinr,zup=zup,zdown=zdown,userad=userad)
            partZ = data['partZ']
            partR = data['partR']
            vz = data['vz']
            vr = data['vr']
            TrueTemp = data['TrueTemp']
            Gu = data['Gu']
            rho = data['rho']
            cregy = data['cregy']
            converted_rho = data['convertedrho']
            Gmass = data['Gmass']
            vmax = data['vmax']
            #vmax=8
            vmin = data['vmin']
            #vmin=0
            Gid = data['Gid']
            if trackgas==1:
                idint0=np.in1d(Gid,Gid_t)
                if userad==1:
                    partR=partR[idint0]
                    vr=vr[idint0]
                else:
                    partZ=partZ[idint0]
                    vz=vz[idint0]
                TrueTemp=TrueTemp[idint0]
                converted_rho=converted_rho[idint0]
                Gmass = Gmass[idint0]
            if userad==1:
                x = partR
            else:
                x = partZ
            if wanted=='Tvz':
                if userad==1:
                    x = np.log10(vr)
                else:
                    x = np.log10(vz)
            if wanted=='pcrpth':
                pth = (GAMMA-1.0)*Gu*u_codetocgs*(rho*rho_codetocgs) #cgs
                x = np.log10(pth)
                print 'x', x
            if wanted=='Tz' or wanted=='Tvz':
                y = np.log10(TrueTemp)
            if wanted=='rhoz':
                y = np.log10(converted_rho)
            if wanted=='vzz':
                if userad==1:
                    y = np.log10(vr)
                else:
                    y = np.log10(vz)
            if wanted=='pcrpth':
                pcr = (CRgamma-1.0)*cregy*cregy_codetocgs*(rho*rho_codetocgs)/(Gmass*m_codetocgs) #cgs
                y = np.log10(pcr)
                print 'y', y
            gridxx = np.linspace(extent[0],extent[1],nobin)
            gridyy = np.linspace(extent[2],extent[3],nobin)
            H, xedges, yedges = np.histogram2d(x, y, bins=[gridxx, gridyy],weights=Gmass*1e10)
            H=H.T
            Hadd += H
            nooftimes += 1
        Hadd=Hadd/nooftimes
        if needcontour==1:
            numaspect = (extent[0]-extent[1])/(extent[2]-extent[3])
            gridag[idir].set_aspect(numaspect)
            print 'np.sum(Hadd)', np.sum(Hadd)
            print 'xedges', xedges
            print 'extent', extent
            levels = np.linspace(vmin,vmax,num=9)
            im = gridag[idir].contourf(xedges[:-1], yedges[:-1], np.log10(Hadd),\
            aspect='auto',cmap=cbcolor,levels=levels, extent=extent,origin='lower')
        else:
            im = gridag[idir].imshow(np.log10(Hadd), extent=extent,\
            aspect='auto',interpolation='nearest', cmap=cbcolor, origin='lower',vmax=vmax,vmin=vmin)
        cblabel = r'$\mathrm{Log}_{10} (M {\rm [M_\odot\;per\;grid]})$'
        if userad==1:
            gridag[idir].set_xlabel(r"$r {\rm [kpc]}$")
        else:
            gridag[idir].set_xlabel(r"$z {\rm [kpc]}$")
        if wanted=='Tvz':
            if userad==1:
                gridag[idir].set_xlabel(r"$\mathrm{Log}_{10} (v_r {\rm [km/s]})$")
            else:
                gridag[idir].set_xlabel(r"$\mathrm{Log}_{10} (v_z {\rm [km/s]})$")
        if wanted=='pcrpth':
            gridag[idir].set_xlabel(r"$\mathrm{Log}_{10} (P_{\rm th} {\rm [erg/cm^3]})$")
            gridag[idir].set_ylabel(r"$\mathrm{Log}_{10} (P_{\rm cr} {\rm [erg/cm^3]})$")
        if wanted=='rhoz':
            gridag[idir].set_ylabel(r"$\mathrm{Log}_{10} (n {\rm [cm^{-3}]})$")
        if wanted=='Tz' or wanted=='Tvz':
            gridag[idir].set_ylabel(r"$\mathrm{Log}_{10} (T {\rm [K]})$")
        if wanted=='vzz':
            if userad==1:
                gridag[idir].set_ylabel(r"$\mathrm{Log}_{10} (v_r {\rm [km/s]})$")
            else:
                gridag[idir].set_ylabel(r"$\mathrm{Log}_{10} (v_z {\rm [km/s]})$")
        if wanted == 'Tz':
            if needcontour==1:
                totalname = 'CRplot/Tz/Tz_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+'_contour'+tlabel+'ag.pdf'
            else:
                totalname = 'CRplot/Tz/Tz_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+tlabel+'ag.pdf'
        if wanted == 'rhoz':
            if needcontour==1:
                totalname = 'CRplot/rhoz/rhoz_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+'_contour'+tlabel+'ag.pdf'
            else:
                totalname = 'CRplot/rhoz/rhoz_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+tlabel+'ag.pdf'
        if wanted == 'vzz':
            totalname = 'CRplot/vzz/vzz_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+tlabel+'ag.pdf'
        if wanted == 'Tvz':
            totalname = 'CRplot/Tvz/Tvz_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+tlabel+'ag.pdf'
        if wanted == 'pcrpth':
            totalname = 'CRplot/pcrpth/pcrpth_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+tlabel+'ag.pdf'
        del Hadd


#cbar = grid.cbar_axes[0].colorbar(im,extend='both', norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax)
cbar = gridag.cbar_axes[0].colorbar(im, extend='max')
cbar.ax.set_ylabel(cblabel, rotation=270,fontsize=12, labelpad=15)
cbar.ax.tick_params(labelsize=10)

for ax, im_title, up_title in zip(gridag, plotit, plotupt):
      t = add_inner_title(ax, "\n".join(wrap(im_title,17)), loc=3)
      t.patch.set_alpha(0.5)
      ax.set_title(up_title)

#plt.axes().set_aspect('equal', 'datalim')

print 'totalname', totalname
plt.savefig(totalname,bbox_inches='tight',pad_inches = 0.1)
plt.clf()
