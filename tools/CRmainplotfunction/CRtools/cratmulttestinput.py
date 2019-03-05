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

def cratmulttestinput(subdict):
        dirneed=subdict['dirneed']
        wanted=subdict['wanted']
        startno=subdict['startno']
        Nsnap=subdict['Nsnap']
        snapsep=subdict['snapsep']
        the_prefix=subdict['the_prefix']
        the_suffix=subdict['the_suffix']
        fmeat=subdict['fmeat']
        title='MW'
        ptitle='LSG'
        nolegend=0


        print 'wanted', wanted
        print 'fmeat', fmeat
        print 'runtodo', dirneed
        print 'startno', startno
        print 'Nsnap', Nsnap
        print 'snapsep', snapsep


        if wanted=='cratmult':
            plt.figure(figsize=(5,4))
            normalizedsm=subdict['normalizedsm']
            M1labelneed=subdict['M1labelneed']
            M1runlabelneed=subdict['M1runlabelneed']
            resoneed=subdict['resoneed']
            diffusionsolverneed=subdict['diffusionsolverneed']
            newlabelneed=subdict['newlabelneed']
            strlabelneed=subdict['strlabelneed']
            legendneed=subdict['legendneed']
            runlabelneed=subdict['runlabelneed']
            withincr=subdict['withincr']
            checkcr=0
            twolegend=0
            rateneed=subdict['rateneed']
            #withincr=0
            snon=1
            coolon=1
            adon=1
            numon=0
            stroff=0
            avetime=0 #number of time bins
            outputcrout=0
            ratiocrout_sn=0
            ratiocrout_sou=subdict['ratiocrout_sou']
            if ratiocrout_sou==1:
                snon=coolon=adon=numon=0
            if twolegend==1:
                import matplotlib.lines as mlines
                lxlist=[]
                lslist=[]
                dclablist=[]
                for icount in range(int(len(dirneed))):
                    lineobj = mlines.Line2D([], [])
                    lxlist.append(lineobj)
            plotn=0
            for runtodo in dirneed:
                snaplist=[]
                enclist=[]
                englist=[]
                enllist=[]
                endlist=[]
                enplist=[]
                enalist=[]
                avesfrl=[]
                prel=0
                preg=0
                prec=0
                pred=0
                prep=0
                prea=0
                presm = 0
                presnap = []
                crecuml = []
                crecumg = []
                crecuma = []
                crecumd = []
                crecum  = []
                crecump = []
                crecumlout = []
                crecumgout = []
                crecumaout = []
                crecumdout = []
                crecumout  = []
                crecumpout = []
                crecumain = []
                timel = []
                info=outdirname(runtodo, Nsnap)
                havecr = info['havecr']
                if havecr==0:
                    continue
                for i in range(startno,Nsnap,snapsep):
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
                    exceptcool=info['exceptcool']
                    havemetal=info['havemetal']
                    the_prefix=info['the_prefix']
                    the_suffix=info['the_suffix']
                    stron = info['stron']
                    newlabel = info['newlabel']
                    strlabel = info['strlabel']
                    haveB = info['haveB']
                    color=info['color']
                    ptitle=title
                    labelneed=dclabel
                    if newlabelneed==1:
                        labelneed="\n".join(wrap(newlabel,17))
                        if strlabelneed==1:
                                labelneed="\n ".join(wrap(strlabel,40))
                    if runtitle=='SMC':
                        ptitle='Dwarf'
                        rin = 6 
                    elif runtitle=='SBC':
                        ptitle='Starburst'
                        rin = 10
                    elif runtitle=='MW':
                        ptitle=r'$L\star$ Galaxy'
                        rin = 10
                    print 'the_snapdir', the_snapdir
                    print 'Nsnapstring', Nsnapstring
                    print 'havecr', havecr
                    print 'exceptcool', exceptcool
                    try:
                        if cosmo==1:
                                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr, cosmological=1,h0=1,exceptcool=exceptcool,havemetal=havemetal)
                                header=G['header']
                                timeneed=header[2]
                        else:
                                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,exceptcool=exceptcool, havemetal=havemetal)
                        Gpos = G['p']
                        Gx = Gpos[:,0]
                        Gy = Gpos[:,1]
                        Gz = Gpos[:,2]
                        Gv = G['v']
                        Gvx = Gv[:,0]
                        Gvy = Gv[:,1]
                        Gvz = Gv[:,2]
                        Grho = G['rho']
                        Gm = G['m']*1e10
                        Gu = G['u']
                        Gr = np.sqrt(np.square(Gx)+np.square(Gy)+np.square(Gz))
                        cutr = Gr<rin
                        codetocgs = 1e10*Msun_in_g*km_in_cm*km_in_cm #from code unit to erg
                        cregyl = G['cregyl']*codetocgs
                        cregyg = G['cregyg']*codetocgs
                        cregy  = G['cregy']*codetocgs
                        cregyd = G['cregyd']*codetocgs
                        if withincr==1 or outputcrout==1 or ratiocrout_sou==1:
                            cregylout = cregyl[~cutr]
                            cregygout = cregyg[~cutr]
                            cregyout = cregy[~cutr]
                            cregydout = cregyd[~cutr]
                        if withincr==1:
                            cregyl=cregyl[cutr]
                            cregyg=cregyg[cutr]
                            cregy=cregy[cutr]
                            cregyd=cregyd[cutr]
                        if havecr>4:
                            cregyp = G['cregyp']*codetocgs
                            if withincr==1 or outputcrout==1 or ratiocrout_sou==1:
                                cregypout = cregyp[~cutr]
                            if withincr==1:
                                cregyp=cregyp[cutr]
                        if havecr>5:
                            cregya = G['cregya']*codetocgs
                            if withincr==1 or outputcrout==1 or ratiocrout_sou==1:
                                cregyaout = cregya[~cutr]
                            if withincr==1:
                                cregya=cregya[cutr]
                            if ratiocrout_sou==1:
                                cregyain = cregya[cutr]
                                cregyain = cregyain[cregyain>0.] #take gain only
                    except KeyError:
                            print 'KeyError'
                            continue
                    crecum = np.append(crecum, np.sum(cregy))
                    crecuml = np.append(crecuml,np.sum(cregyl))
                    crecumg = np.append(crecumg, np.sum(cregyg))
                    crecumd = np.append(crecumd, np.sum(cregyd))
                    if havecr>4:
                        crecump = np.append(crecump, np.sum(cregyp))
                    if havecr>5:
                        crecuma = np.append(crecuma, np.sum(cregya))
                    if withincr==1 or  outputcrout==1 or ratiocrout_sou==1:
                        crecumout = np.append(crecumout, np.sum(cregyout))
                        crecumlout = np.append(crecumlout,np.sum(cregylout))
                        crecumgout = np.append(crecumgout, np.sum(cregygout))
                        crecumdout = np.append(crecumdout, np.sum(cregydout))
                        if havecr>4:
                            crecumpout = np.append(crecumpout, np.sum(cregypout))
                        if havecr>5:
                            crecumaout = np.append(crecumaout, np.sum(cregyaout))
                    if ratiocrout_sou==1:
                        if havecr>5:
                            crecumain = np.append(crecumain, np.sum(cregyain))
                    print 'cregy', np.sum(cregy)
                    print 'cregyp', np.sum(cregyp)
                    if cosmo==1:
                            readtimelist=readtime(firever=2)
                            snap2list=readtimelist['snaplist']
                            time2list=readtimelist['tiimelist']
                            a2list=readtimelist['alist']
                            timenow = np.interp(timeneed,a2list,time2list)*1e3
                    else:
                            timenow=i*0.98
                            timel = np.append(timel, timenow)
                            prec=np.sum(cregy)
                            prel=np.sum(cregyl)
                            preg=np.sum(cregyg)
                    if havecr>5:
                        prea=np.sum(cregya)
                        pred=np.sum(cregyd)
                        prep=np.sum(cregyp)
                crecumnum = crecum - crecumg - crecuml - crecuma
                if stron==1:
                    crecumpold = crecump 
                    #crecump = crecumpold - (crecum - crecumg - crecuml - crecumd)
                    crecumnum = crecumnum-crecumpold
                if withincr==1 or outputcrout==1 or ratiocrout_sou==1:
                    #crecumesc = -crecumout + crecumaout + crecumgout + crecumlout + crecumpout
                    crecumesc = -np.array(crecumout)+ np.array(crecumaout)+np.array(crecumlout)
                if plotn ==0:
                    snelabel = 'Supernovae'
                    losslabel = 'Loss'
                    crecumalabel = 'Adiabatic'
                    strlabel = 'Streaming'
                    numlabel = 'Numerical'
                    if withincr==1:
                        esclabel = 'Escape'
                    elif outputcrout==1:
                        esclabel = r'$\dot{E}_{\rm cr,esc}(r\,>$'+str(int(rin))+'kpc$)$'
                else:
                    snelabel = losslabel = crecumalabel=strlabel=esclabel= '_nolegend_' 
                if np.mod(plotn,4)==0:
                    ls='-'
                elif np.mod(plotn,4)==1:
                    ls='dashed'
                elif np.mod(plotn,4)==2:
                    ls='dashdot'
                elif np.mod(plotn,4)==3:
                    ls='dotted'
                if rateneed==1:
                    xcr = -(crecum[1:]-crecum[:-1])/(timel[1:]-timel[:-1])/1e6/yr_in_sec
                    xcrg = (crecumg[1:]-crecumg[:-1])/(timel[1:]-timel[:-1])/1e6/yr_in_sec
                    xcrl = (crecuml[1:]-crecuml[:-1])/(timel[1:]-timel[:-1])/1e6/yr_in_sec
                    xcra = (crecuma[1:]-crecuma[:-1])/(timel[1:]-timel[:-1])/1e6/yr_in_sec
                    if stron ==1:
                            xcrp = (crecump[1:]-crecump[:-1])/(timel[1:]-timel[:-1])/1e6/yr_in_sec
                    if withincr==1 or outputcrout==1 or ratiocrout_sou==1:
                        xcresc = (crecumesc[1:]-crecumesc[:-1])/(timel[1:]-timel[:-1])/1e6/yr_in_sec
                    if ratiocrout_sou==1:
                        xcrain = (crecumain[1:]-crecumain[:-1])/(timel[1:]-timel[:-1])/1e6/yr_in_sec
                    xcrn = (crecumnum[1:]-crecumnum[:-1])/(timel[1:]-timel[:-1])/1e6/yr_in_sec
                    timel = (timel[1:]+timel[:-1])/2.
                else:
                    xcrg = crecumg
                    xcrl = crecuml
                    xcra = crecuma
                    xcr = crecum
                    xcrn = crecumnum
                    if stron ==1:
                        xcrp = crecump
                    if withincr==1 or outputcrout==1 or ratiocrout_sou==1:
                        xcresc = crecumesc
                    if ratiocrout_sou==1:
                        xcrain = crecumain
                if avetime>1:
                    notimebin = len(xcrg)
                    noelement=notimebin/avetime
                    newelementno = avetime*noelement
                    xcrg = xcrg[0:newelementno]
                    xcrl = xcrl[0:newelementno]
                    xcra = xcra[0:newelementno]
                    xcr = xcr[0:newelementno]
                    xcrn = xcrn[0:newelementno]
                    if stron ==1:
                        xcrp = xcrp[0:newelementno]
                    if withincr==1 or outputcrout==1:
                        xcresc = xcresc[0:newelementno]
                    if ratiocrout_sou==1:
                        xcrain = xcrain[0:newelementno]
                        timel = timel[0:newelementno]
                        xcrg = np.average(xcrg.reshape((noelement,avetime)),axis=1)
                        xcrl = np.average(xcrl.reshape((noelement,avetime)),axis=1)
                        xcra = np.average(xcra.reshape((noelement,avetime)),axis=1)
                        xcr = np.average(xcr.reshape((noelement,avetime)),axis=1)
                        xcrn = np.average(xcrn.reshape((noelement,avetime)),axis=1)
                    if stron ==1:
                            xcrp = np.average(xcrp.reshape((noelement,avetime)),axis=1)
                    if withincr==1 or outputcrout==1:
                            xcresc = np.average(xcresc.reshape((noelement,avetime)),axis=1)
                    if ratiocrout_sou==1:
                            xcrain = np.average(xcrain.reshape((noelement,avetime)),axis=1)
                    timel = np.average(timel.reshape((noelement,avetime)),axis=1)
                if runlabelneed==1:
                    if snon==1:
                        plt.plot(timel, xcrg, lw=2,ls=ls, label=labelneed,color=cmaps.inferno(0.1))
                    if coolon==1:
                        plt.plot(timel, xcrl, lw=2,ls=ls, color=cmaps.inferno(0.4))
                    if adon==1:
                        plt.plot(timel, xcra, lw=2, ls=ls, color=cmaps.inferno(0.7))
                    if checkcr==1:
                        plt.plot(timel, xcr, lw=2, ls=ls, color='k')
                    if stron ==1 and stroff==0:
                            plt.plot(timel, xcrp, lw=2, ls=ls,color=cmaps.inferno(0.9))
                    if numon ==1:
                        plt.plot(timel, xcrp, lw=2, ls=ls,color='g')
                    if (withincr==1 or outputcrout==1) and ratiocrout_sn==0 and ratiocrout_sou==0:
                        plt.plot(timel, xcresc, lw=2, ls=ls,color='m')
                    if ratiocrout_sn==1:
                        plt.plot(timel, np.absolute(xcresc/xcrg), lw=2, ls=ls,label=labelneed,color='k')
                    if ratiocrout_sou==1:
                            ls='solid'
                            if haveB>0:
                                    ls='dashed'
                            plt.plot(timel, np.absolute(xcresc/(xcrg+xcrain+xcrn)), lw=2, ls=ls, label=labelneed,color=color)
                else:
                    if snon==1:
                        plt.plot(timel, xcrg, lw=2,ls=ls, label=snelabel,color=cmaps.inferno(0.1))
                    if coolon==1:
                        plt.plot(timel, xcrl, lw=2,ls=ls, label=losslabel,color=cmaps.inferno(0.4))
                    if adon==1:
                        plt.plot(timel, xcra, lw=2, ls=ls, label=crecumalabel,color=cmaps.inferno(0.7))
                    if checkcr==1:
                        plt.plot(timel, xcr, lw=2, ls=ls,color='k')
                    if stron ==1 and stroff==0:
                        plt.plot(timel, xcrp, lw=2, ls=ls, label=strlabel,color=cmaps.inferno(0.9))
                    if numon==1:
                            plt.plot(timel, xcrn, lw=2,ls=ls, label=numlabel,color='g')
                    if (withincr==1 or outputcrout==1) and ratiocrout_sn==0 and ratiocrout_sou==0:
                            plt.plot(timel, xcresc, lw=2, ls=ls,label=esclabel,color='m')
                    if ratiocrout_sn==1:
                            plt.plot(timel, np.absolute(xcresc/xcrg), lw=2, ls=ls,color='k')
                    if ratiocrout_sou==1:
                        ls='solid'
                        if haveB>0:
                            ls='dashed'
                        plt.plot(timel, np.absolute(xcresc/(xcrg+xcrain+xcrn)), lw=2, ls=ls,color=color)
                    if twolegend==1:
                        lxlist[plotn] = mlines.Line2D([], [], color='r', ls=ls,label=dclabel)
                        lslist = np.append(lslist,ls)
                        dclablist = np.append(dclablist,labelneed)
                plotn+=1
        if legendneed==1:
            if runlabelneed==1:
                plt.legend(loc='best',ncol=2,fontsize=10)
            else:
                if twolegend==1:
                    legend1 = plt.legend(lxlist, dclablist, loc=2,fontsize=8,ncol=2)
                    plt.gca().add_artist(legend1)
                    plt.legend(loc=3, fontsize=8)
                else:
                    plt.legend(loc='best', ncol=2, fontsize=10)
        if rateneed==1:
                    plt.ylabel(r'$ \dot{E}_{\rm cr} [ {\rm erg/s}]$',fontsize=18)
        else:
            plt.ylabel(r'${\rm CR\; energy\; [ erg}]$',fontsize=18)
        if ratiocrout_sn==1 and rateneed==1:
            plt.ylabel(r'$\dot{E}_{\rm esc}/\dot{E}_{\rm SN}$',fontsize=18)
        if ratiocrout_sn==1 and rateneed==0:
                plt.ylabel(r'$E_{\rm esc}/E_{\rm SN}$',fontsize=18)
        if ratiocrout_sou==1 and rateneed==1:
                plt.ylabel(r'$\dot{E}_{\rm esc}/\dot{E}_{\rm source}$',fontsize=18)
        if ratiocrout_sou==1 and rateneed==0:
                plt.ylabel(r'$E_{\rm esc}/E_{\rm source}$',fontsize=18)
        if rateneed==1 and ratiocrout_sn==0 and ratiocrout_sou==0:
            if runtitle=='SMC':
                #plt.ylim(ymax=2.3e39,ymin=-2.3e39) 
                plt.ylim(ymax=2.3e39,ymin=-1.0e39) 
            if runtitle=='MW':
                #plt.ylim(ymax=2.5e41,ymin=-2.5e41)
                plt.ylim(ymax=2.7e41,ymin=-1.2e41)
            if runtitle=='SBC':
                #plt.ylim(ymax=4.4e41,ymin=-4.4e41)
                plt.ylim(ymax=4.4e41,ymin=-3.0e41)
        if ratiocrout_sn==1 or ratiocrout_sou==1:
            if rateneed==1:
                plt.ylim(ymin=0.0,ymax=1.5)
            else:
                plt.ylim(ymin=0.0,ymax=1.1)
        plt.xlabel(r'$t [{\rm Myr}]$', fontsize=18)
        plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
        plt.title(ptitle, fontsize=16)
        fnsuffix = ''
        if withincr==1:
            fnsuffix = '_withinr'
        if outputcrout==1:
                    fnsuffix = '_outputcrout'
        if ratiocrout_sn==1:
            fnsuffix = '_ratiocrout_sn'
        if ratiocrout_sou==1:
                fnsuffix = '_ratiocrout_sou'
        if rateneed==1:
            filename=plotloc+'CRplot/cratmult/'+fmeat+'_rate'+fnsuffix+'.pdf'
        else:
            filename=plotloc+'CRplot/cratmult/'+fmeat+'_time'+fnsuffix+'.pdf'
        print 'filename', filename
        plt.subplots_adjust(left  = 0.18, bottom=0.15, right=0.95, top=0.9)
        plt.savefig(filename)
        plt.clf()
        return None
