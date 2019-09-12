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
import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='inferno', cmap=cmaps.inferno)
plt.register_cmap(name='plasma', cmap=cmaps.plasma)
import collections


def disksurdentestinput(subdict):
        dirneed=subdict['dirneed']
        wanted=subdict['wanted']
        print 'wanted', wanted
        startno=subdict['startno']
        Nsnap=subdict['Nsnap']
        snapsep=subdict['snapsep']
        the_prefix=subdict['the_prefix']
        the_suffix=subdict['the_suffix']
        fmeat=subdict['fmeat']
        usez1kpc = subdict['usez1kpc']
        print 'usez1kpc',usez1kpc
        title = subdict['title']
        maxlength = subdict['maxlength']
        multz = subdict['multz']
        print 'multz', multz
        #title='MW'
        ptitle='LSG'
        nolegend=0
        resoneed=subdict['resoneed']
        newlabelneed=subdict['newlabelneed']
        legendneed = subdict['legendneed']
        plt.figure(figsize=(5,4))
        for runtodo in dirneed:
            print 'runtodo', runtodo
            withinr=15.
            nogrid = 15
            dr = withinr/nogrid
            radlist=np.linspace(0.001,withinr,num=nogrid)
            Slist=0.*radlist
            Glist=0.*radlist
            numoftimes=0
            for i in range(startno,Nsnap,snapsep):
                        snaplist=[]
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
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
                        cosmo=info['cosmo']
                        maindir=info['maindir']
                        color=info['color']
                        haveB=info['haveB']
                        M1speed=info['M1speed']
                        newlabel=info['newlabel']
                        snumadd=info['snumadd']
                        usepep=info['usepep']
                        halostr=info['halostr']
                        ptitle=title
                        if runtitle=='SMC':
                            ptitle='Dwarf'
                        elif runtitle=='SBC':
                            ptitle='Starburst'
                        elif runtitle=='MW':
                            ptitle=r'$L\star$ Galaxy'
                        labelneed=dclabel
                        if newlabelneed==1:
                            labelneed="\n".join(wrap(newlabel,17))
                        rotface=1
                        loccen=1
                        if cosmo==1:
                            h0=1
                        else:
                            h0=0
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        G = readsnapwcen(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix,\
                         havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
                         loccen=loccen,runtodo=runtodo,rundir=rundir,halostr=halostr)
                        Gsurdenlist = calsurden(G,radlist,maxlength)
                        S = readsnapwcen(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix,\
                         havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
                         loccen=loccen,runtodo=runtodo,rundir=rundir,halostr=halostr)
                        Ssurdenlist = calsurden(S,radlist,maxlength)
                        Glist += Gsurdenlist
                        Slist += Ssurdenlist
                        numoftimes+=1
            plt.plot(radlist[:-1],(Slist[:-1]+Glist[:-1])/numoftimes,label=labelneed,lw=2,ls='solid',color=color)
            plt.plot(radlist[:-1],(Glist[:-1])/numoftimes,lw=2,ls='dashed',color=color)
            plt.plot(radlist[:-1],(Slist[:-1])/numoftimes,lw=2,ls='dashdot',color=color)
        import matplotlib.lines as mlines
        linelablist=[]
        lxlist=[]
        lslist=['solid','dashed','dashdot']  
        linelablist = ['tot','gas','star']  
        for icount in range(int(len(lslist))):
            lineobj = mlines.Line2D([], [])
            lxlist.append(lineobj)
        for ils, ls in enumerate(lslist):
            linelabel = linelablist[ils]
            lxlist[ils] = mlines.Line2D([], [], color='k', ls=ls,label=linelabel)
            linelablist = np.append(linelablist,linelabel)
        legend1 = plt.legend(lxlist, linelablist, loc='lower left',fontsize=8,ncol=1)
        plt.gca().add_artist(legend1)
        plt.xlabel(r'${\rm r [kpc]}$', fontsize=18)
        plt.xlim(xmax=np.amax(radlist))
        plt.ylabel(r'$\Sigma\;[{\rm M_\odot/pc^2}]$',fontsize=18)
        plt.yscale('log')
        plt.legend(loc='best',fontsize=10,ncol=2)
        filename=plotloc+'CRplot/gassurdenM_pc2/surdenM_pc2_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if cosmo==1:
            plt.legend(loc='upper right',fontsize=10,ncol=2)
        print 'filename', filename
        plt.title(ptitle,fontsize=18)
        plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
        plt.subplots_adjust(left  = 0.18, bottom=0.15, right=0.95, top=0.9)
        plt.savefig(filename)
        plt.clf()
        return None
