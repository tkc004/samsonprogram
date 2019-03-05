from samson_const import *
from pathloc import *
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
from crtestfunction import *
from matplotlib import rcParams
from pylab import *
from textwrap import wrap
from scipy.optimize import curve_fit
#rcParams['figure.figsize'] = 5, 5
rcParams['figure.figsize'] = 10, 5
rcParams['font.size']=12
rcParams['font.family']='serif'
rcParams['text.usetex']=True
#rcParams.update({'figure.autolayout': True})
import matplotlib.patches as patches
rcParams['axes.linewidth'] = 2
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True
rcParams['axes.unicode_minus']=False
colortable = [ 'b', 'g', 'r']
dirneed=['bwmwmr','bwmwmrstr','bwmwmrdc27','bwmwmrdc28', 'bwmwmrdc29']
#dirneed=['bwmwmr']
#dirneed=['bwsmclrdc28']
#dirneed=['bwmwmrdc29']
#dirneed=['m11gcr_700']
#dirneed=['bwsbclrdc27']
#dirneed=['m11dcr_b_70']
wanted='rhoT'
startno=490
Nsnap=501
snapsep=10
the_prefix='snapshot'
the_suffix='.hdf5'
fmeat=''

if wanted=='rhoT' or wanted=='rhoTwind':
    rcParams['figure.figsize'] = 5, 4
    if wanted=='rhoTwind':
        vcut=0.0 #outflow velocity cut
        withinr=20.0
        zup=25.0
        zdown=15.0
        extent = [-7,-1, 1, 9]
    else:
        extent = [-7,4, 1, 9]
    nooftimes=0
    nobin=51
    needcontour=1
        for runtodo in dirneed:
        Hadd = np.zeros((nobin-1,nobin-1))
        for i in range(startno,Nsnap,snapsep):
            info=outdirname(runtodo, i)
            rundir=info['rundir']
            maindir=info['maindir']
            halostr=info['halostr']
            runtitle=info['runtitle']
            slabel=info['slabel']
            snlabel=info['snlabel']
            dclabel=info['dclabel']
            resolabel=info['resolabel']
            the_snapdir=info['the_snapdir']
            Nsnapstring=info['Nsnapstring']
            havecr=info['havecr']
            Fcal=info['Fcal']
            iavesfr=info['iavesfr']
            timestep=info['timestep']
            color=info['color']
            cosmo=info['cosmo']
            usepep=info['usepep']
            if wanted=='rhoTwind':
                vmax=8.5
                vmin=0.0
                if runtitle=='MW':
                    vmax = 6.0
                    vmin = -1.0
                if runtitle=='SMC':
                    vmax = 5.0
                    vmin = 0.5 
                if runtitle=='SBC':
                    vmax = 6.5
                                        vmin = 0.5
            if wanted=='rhoT':
                                if runtitle=='MW':
                                        vmax = 9.0
                                        vmin = -2.0
                                if runtitle=='SMC':
                                        vmax = 8.0
                                        vmin = 0.5
            if cosmo==1:
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                header=readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1,header_only=1)
                h0 = header['hubble']
                atime = header['time']
                if usepep==1:
                    halosA = read_halo_history_pep(rundir, finalno, beginno=beginno,\
                     singlesnap=0, firever=firever,halonostr=halostr, hubble=h0, comoving=1, maindir=maindir)
                    afactor=atime
                else:
                    halosA = read_halo_history(rundir, maindir=maindir, halonostr=halostr, hubble=h0, comoving=0)
                    afactor=1.0
                redlist = halosA['redshift']
                haloid = halosA['ID']
                a_scale = 1.0/(1.0+redlist)
                xcenl = halosA['x']*afactor
                ycenl = halosA['y']*afactor
                zcenl = halosA['z']*afactor
                xvl = halosA['xv']
                yvl = halosA['yv']
                zvl = halosA['zv']
                Rvirl = halosA['R']*afactor
                xcen = np.interp(atime,a_scale,xcenl)
                ycen = np.interp(atime,a_scale,ycenl)
                zcen = np.interp(atime,a_scale,zcenl)
                                xvcen = np.interp(atime,a_scale,xvl)
                                yvcen = np.interp(atime,a_scale,yvl)
                                zvcen = np.interp(atime,a_scale,zvl)
                Rvir = np.interp(atime,a_scale,Rvirl)
            else:
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                xcen =0; ycen=0; zcen=0; xvcen=0; yvcen=0; zvcen=0;
                h0=1

            Gpos = G['p'][:,:]
            Gvel = G['v'][:,:]
            Gmass = G['m'][:]*1e10
            print 'np.sum(Gmass)', np.sum(Gmass)
            print 'np.amax(Gmass)', np.amax(Gmass)
            print 'np.amin(Gmass)', np.amin(Gmass)
            Tb = G['u'][:]
            rho = G['rho'][:]
            Neb = G['ne'][:]
            partX = Gpos[:,0]-xcen
            partY = Gpos[:,1]-ycen
            partZ = Gpos[:,2]-zcen
            partXV = Gvel[:,0]-xvcen
            partYV = Gvel[:,1]-yvcen
            partZV = Gvel[:,2]-zvcen
            header = G['header'][:]
            redshift = header[3]
            boxsize = header[9]
            if wanted == 'rhoTwind':
                cutz = (np.absolute(partZ)>zdown) & (np.absolute(partZ)<zup) #kpc
                cutxy = partX*partX+partY*partY<withinr*withinr
                if cosmo==1:
                    cutv = partXV*partX+partYV*partY\
                    +partZV*partZ>vcut*np.sqrt(partX*partX+partY*partY+partZ*partZ) #outflowing gas
                else:
                    cutv = partZV*partZ/np.absolute(partZ)>vcut #outflowing gas
                cut = cutz*cutxy*cutv
                Tb = Tb[cut]
                rho = rho[cut]
                Neb = Neb[cut]
                Gmass = Gmass[cut]

            TrueTemp, converted_rho  = SF.convertTemp(Tb, Neb, rho)
            if wanted == 'rhoT':
                if needcontour==1:
                                        totalname = 'CRplot/rhoT/rhoT_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+'_contour.pdf'
                else:
                    totalname = 'CRplot/rhoT/rhoT_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            elif wanted == 'rhoTwind':
                                if needcontour==1:
                    totalname = 'CRplot/rhoTwind/rhoTwind_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+'_contour.pdf'
                else:
                    totalname = 'CRplot/rhoTwind/rhoTwind_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            y = np.log10(TrueTemp)
            x = np.log10(converted_rho)
                        gridxx = np.linspace(extent[0],extent[1],nobin)
                        gridyy = np.linspace(extent[2],extent[3],nobin)
            #gridxx = np.linspace(np.amin(x),np.amax(x),nobin)
            #gridyy = np.linspace(np.amin(y),np.amax(y),nobin)
            H, xedges, yedges = np.histogram2d(x, y, bins=[gridxx, gridyy],weights=Gmass)
            H=H.T
            Hadd += H
            nooftimes += 1 
        Hadd=Hadd/nooftimes
        if needcontour==1:
            levels = np.linspace(vmin,vmax,num=15)
            plt.contourf(xedges[:-1], yedges[:-1], np.log10(Hadd), levels=levels, extent=extent, origin='lower',aspect='auto')
        else:
            plt.imshow(np.log10(Hadd), extent=extent, interpolation='nearest',origin='lower',aspect='auto',vmax=vmax,vmin=vmin)
        cbar = plt.colorbar(extend='both', norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax))
        cbar.set_label(r'$\mathrm{Log}_{10} (M {\rm [M_\odot]})$', rotation=270, labelpad=20)
        plt.xlabel(r"$\mathrm{Log}_{10} (n {\rm [cm^{-3}]})$")
        plt.ylabel(r"$\mathrm{Log}_{10} (T {\rm [K]})$")
        plt.title(titleneed)
        plt.tight_layout()
        print totalname
        plt.savefig(totalname)
        plt.clf()
