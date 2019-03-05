from samson_const import *
from pathdir import *
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
import collections

def gammasfrsnaptestinpput(subdict):
        dirneed=subdict['dirneed']
        wanted=subdict['wanted']
        startno=subdict['startno']
        Nsnap=subdict['Nsnap']
        snapsep=subdict['snapsep']
        the_prefix=subdict['the_prefix']
        the_suffix=subdict['the_suffix']
        fmeat=subdict['fmeat']


        if wanted=='dirgammasfr' or wanted=='dirgamma' or wanted=='dirsfr' or wanted=='dirsm':
            Rfrac=0.25
            xaxis_snapno=0
            normalizedsm=subdict['normalizedsm']
            M1labelneed=subdict['M1labelneed']
            M1runlabelneed=subdict['M1runlabelneed']
            resoneed=subdict['resoneed']
            diffusionsolverneed=subdict['diffusionsolverneed']
            newlabelneed=subdict['newlabelneed']
            strlabelneed=subdict['strlabelneed']
            showstarburst=subdict['showstarburst']
            legendneed=subdict['legendneed']
            correctIa=subdict['correctIa']
            refereelabelneed=subdict['refereelabelneed']
            for runtodo in dirneed:
                print 'runtodo', runtodo
                snaplist=[]
                enclist=[]
                englist=[]
                enllist=[]
                sml=[]
                diskml=[]
                bulgeml=[]
                nsml=[]
                timel=[]
                Lgcalist=[]
                pretime=0
                presnap=startno
                havecr=0
                if wanted=='dirgammasfr' or wanted=='dirgamma':     
                    info=outdirname(runtodo, Nsnap)
                    havecr=info['havecr']
                    if havecr==0:
                        continue
                for i in range(startno,Nsnap, snapsep):
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
                    cosmo=info['cosmo']
                    color=info['color']
                    withinRv=info['withinRv']
                    usepep=info['usepep']
                    beginno=info['beginno']
                    finalno=info['finalno']
                    firever=info['firever']     
                    initsnap=info['initsnap']   
                    haveB=info['haveB']
                    M1speed=info['M1speed']
                    Rvirguess=info['Rvir']
                    newlabel=info['newlabel']
                    strlabel=info['strlabel']
                    labelneed=dclabel
                    if newlabelneed==1:
                        labelneed="\n".join(wrap(newlabel,17))
                    if strlabelneed==1:
                        labelneed="\n".join(wrap(strlabel,40))
                    ptitle=title
                    if runtitle=='SMC':
                        ptitle='Dwarf'
                    elif runtitle=='SBC':
                        ptitle='Starburst'
                    elif runtitle=='MW':
                        ptitle=r'$L\star$ Galaxy'
                    if cosmo==0:
                        inittime=initsnap
                    else:
                        inittime=0
                    print 'set initial time'    
                    print 'the_snapdir', the_snapdir
                    print 'Nsnapstring', Nsnapstring
                    print 'havecr', havecr
                    print 'withinRv', withinRv
                    print 'cosmo', cosmo
                    try:
                        if cosmo==1:
                            G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                            S = readsnapcr(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                            if correctIa==1:
                                Disk = readsnapcr(the_snapdir, Nsnapstring, 2, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                                Bulge = readsnapcr(the_snapdir, Nsnapstring, 3, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                        else:
                            G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                            S = readsnapcr(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                            if correctIa==1:
                                Disk = readsnapcr(the_snapdir, Nsnapstring, 2, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                                Bulge = readsnapcr(the_snapdir, Nsnapstring, 3, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        if havecr>4:
                            cregyl = G['cregyl']*1e10*solar_mass_in_g*km_in_cm*km_in_cm #in erg (originally in code unit: 1e10Msun*(km/s)^2)
                            cregyg = G['cregyg']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                        if havecr > 0:
                            cregy_codeunit = G['cregy']
                            cregy  = cregy_codeunit*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                        Grho = G['rho']
                        Neb = G['ne']
                    except KeyError:
                        print 'Keyerror'
                        break
                    try:
                        header=S['header']
                        timeneed=header[2]
                        print 'timeneed', timeneed
                        Smi=S['m']
                        Sage=S['age']
                        if withinRv ==1 and cosmo==1:
                            Sp = S['p']
                            Sx = Sp[:,0]
                            Sy = Sp[:,1]
                            Sz = Sp[:,2]
                            header=readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1,header_only=1)
                            h0 = header['hubble']
                            atime = header['time']
                            if usepep==1:
                                halosA = SF.read_halo_history_pep(rundir, finalno, beginno=beginno,\
                                 singlesnap=0, firever=firever,halonostr=halostr, hubble=h0, comoving=1, maindir=maindir)
                                afactor=atime
                            else:
                                halosA = SF.read_halo_history(rundir, maindir=maindir, halonostr=halostr, hubble=h0, comoving=0)
                                afactor=1.0
                            redlist = halosA['redshift']
                            haloid = halosA['ID']
                            a_scale = 1.0/(1.0+redlist)
                            xcenl = halosA['x']*afactor
                            ycenl = halosA['y']*afactor
                            zcenl = halosA['z']*afactor
                            Rvirl = halosA['R']*afactor
                            xcen = np.interp(atime,a_scale,xcenl)
                            ycen = np.interp(atime,a_scale,ycenl)
                            zcen = np.interp(atime,a_scale,zcenl)
                            Rvir = np.interp(atime,a_scale,Rvirl)
                            Sxrel = Sx-xcen
                            Syrel = Sy-ycen
                            Szrel = Sz-zcen
                            Sr = np.sqrt(Sxrel*Sxrel+Syrel*Syrel+Szrel*Szrel)
                            cutrvs = Sr<Rvir*Rfrac
                            Smi = Smi[cutrvs]
                            Sage = Sage[cutrvs]
                        Sm = np.sum(Smi)*1e10 #in solar mass
                        tcut=Sage>pretime
                        Nsm = np.sum(Smi[tcut])*1e10
                        if correctIa==1:
                            diskm = np.sum(Disk['m'])*1e10
                            bulgem = np.sum(Bulge['m'])*1e10
                    except KeyError:
                        print 'key error'
                        Sm = 0.
                        Nsm = 0.
                        timeneed=0
                        if correctIa==1:
                            diskm = 0.
                            bulgem = 0.
                    if withinRv ==1:
                        Gp = G['p']
                        Gx = Gp[:,0]
                        Gy = Gp[:,1]
                        Gz = Gp[:,2]
                        header=readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1,header_only=1)
                        h0 = header['hubble']
                        atime = header['time']
                        if cosmo==1:
                            if usepep==1:
                                halosA = SF.read_halo_history_pep(rundir, finalno,\
                                beginno=beginno, singlesnap=0, firever=firever,halonostr=halostr,\
                                hubble=h0, comoving=1, maindir=maindir, snapsep=snapsep)
                                afactor=atime
                            else:
                                halosA = SF.read_halo_history(rundir, maindir=maindir, halonostr=halostr, hubble=h0, comoving=0)
                                afactor=1.0
                            redlist = halosA['redshift']
                            haloid = halosA['ID']
                            a_scale = 1.0/(1.0+redlist)
                            xcenl = halosA['x']*afactor
                            ycenl = halosA['y']*afactor
                            zcenl = halosA['z']*afactor
                            Rvirl = halosA['R']*afactor
                            xcen = np.interp(atime,a_scale,xcenl)
                            ycen = np.interp(atime,a_scale,ycenl)
                            zcen = np.interp(atime,a_scale,zcenl)
                            Rvir = np.interp(atime,a_scale,Rvirl)
                        else:
                            xcen=0
                            ycen=0
                            zcen=0
                            Rvir=Rvirguess
                        Gxrel = Gx-xcen
                        Gyrel = Gy-ycen
                        Gzrel = Gz-zcen
                        Gr = np.sqrt(Gxrel*Gxrel+Gyrel*Gyrel+Gzrel*Gzrel)
                        cutrv = Gr<Rvir*Rfrac
                        Grho = Grho[cutrv]
                        Neb = Neb[cutrv]
                        if havecr>0:
                            cregy = cregy[cutrv]
                            cregy_codeunit = cregy_codeunit[cutrv]
                            if havecr>4:
                                    cregyl = cregyl[cutrv]
                                    cregyg = cregyg[cutrv]
                    if cosmo==1:
                            readtimelist=readtime(firever=2)
                            snap2list=readtimelist['snaplist']
                            time2list=readtimelist['timelist']
                            a2list=readtimelist['alist']
                            tnow = np.interp(timeneed,a2list,time2list)*1e9
                            pret = np.interp(pretime,a2list,time2list)*1e9
                    if havecr>4:
                        cregygt=np.sum(cregyg)
                        cregylt=np.sum(cregyl)
                    if havecr>0:
                        cregyt =np.sum(cregy)
                        Lout = outLgamma_nism(Grho,Neb,cregy_codeunit)
                        Lgcal = np.sum(Lout['Lgamma'])
                        print 'Lgcal', Lgcal
                    snaplist.append(float(i))
                    if cosmo==1:
                        timel.append(tnow)
                    else:
                        timel.append(float(i)*0.98*1e6)
                    if havecr>0:
                        enclist.append(cregyt)
                    if havecr>4:
                        englist.append(cregygt)
                        enllist.append(cregylt)
                    if havecr>0:
                        Lgcalist.append(Lgcal)
                    sml.append(Sm)
                    if correctIa==1:
                        diskml.append(diskm)
                        bulgeml.append(bulgem)
                    nsml.append(Nsm)
                    pretime=timeneed
                    del G, S
                sml=np.array(sml)
                if correctIa==1:
                    diskml=np.array(diskml)
                    bulgeml=np.array(bulgeml)
                nsml=np.array(nsml)
                if havecr>0:
                    enclist=np.array(enclist)
                if havecr>4:
                    englist=np.array(englist)
                    enllist=np.array(enllist)
                snaplist=np.array(snaplist)
                if havecr>0:
                    Lgcalist=np.array(Lgcalist)
                timel=np.array(timel) #in yr
                #above is the coefficient for Salpeter only; for Kroupa, the coefficient is 50% larger:
                avesfrl=(nsml[1:])/(timel[1:]-timel[:-1]) #in Msun/yr
                print 'nsml',nsml
                print 'avesfrl', avesfrl
                Lsfr = Kroupa_Lsf*avesfrl*solar_mass_in_g/yr_in_sec*cspeed_in_cm_s*cspeed_in_cm_s #times 1.5? 6.2e-4 for Kroupa 3.8e-4 for Salpeter
                if correctIa==1:
                    IaeffSFR=5.3e-8/3.0e-4*(sml[1:]+diskml[1:]+bulgeml[1:])/0.03753/1e9
                    print 'IaeffSFR', IaeffSFR
                    LsfrnoIa = Kroupa_Lsf*(avesfrl+IaeffSFR)*solar_mass_in_g/yr_in_sec*cspeed_in_cm_s*cspeed_in_cm_s
                print 'Lsfr', Lsfr
                print 'timel', timel
                if havecr>4:
                    Lgamma = (enllist[1:]-enllist[:-1])/((timel[1:]-timel[:-1])*yr_in_sec)/(hadronicdecayrate+coulombdecayrate)*hadronicdecayrate*betapi/nopi_per_gamma
                    Lgamma_sfr = Lgamma/Lsfr
                if havecr>0:
                    Lgcal_sfr = Lgcalist[1:]/Lsfr
                if correctIa==1:
                    Lgamma_sfr_noIa = Lgamma/LsfrnoIa
                if haveB>0:
                    lsn='dashed'
                else:
                    lsn='solid'
                if M1labelneed>1 or M1runlabelneed==1:
                    if M1speed>499:
                        lsn='solid'
                    if M1speed>999:
                        lsn='dashed'
                    if M1speed>1999:
                        lsn='dashdot'
                    if M1speed>3999:
                        lsn='dotted'
                if M1labelneed==1:
                    if M1speed>499:
                            labelneed=r'$\tilde{c}=500$'
                    if M1speed>999:
                            labelneed=r'$\tilde{c}=1000$'
                    if M1speed>1999:
                            labelneed=r'$\tilde{c}=2000$'
                    if M1speed>3999:
                            labelneed=r'$\tilde{c}=4000$'

                if resoneed==1:
                    if resolabel=='llr':
                        labelneed='Lowest res'
                        lsn = 'solid'
                    if resolabel=='lr':
                        labelneed='Low res'
                        lsn = 'dashed'
                    if resolabel=='mr':
                        labelneed='Standard res'
                        lsn = 'dashdot'
                if diffusionsolverneed==1:
                    if runtodo=='bwmwlrdc27ds':
                        labelneed='Zeroth moment'
                        lsn = 'dashed'
                    else:
                        labelneed='Two moment'
                        lsn = 'solid'
                if refereelabelneed==1:
                    if runtodo=='bwsmclrdc28mhdref' or runtodo=='bwmwlrdc28mhdref':
                        labelneed='Modified'
                    else:
                        labelneed='Original'
                    
                if xaxis_snapno==1:
                    xaxisl = snaplist[1:]
                    xlab = 'snapshot number'
                else:
                    xaxisl = timel[1:]/1e6+inittime
                    xlab = 'Myr'
                if showstarburst==1 and runtitle=='SBC':
                    if runtodo=='bwsbclr':
                        pstartno=600
                    if runtodo=='bwsbclrdc0':
                        pstartno=590
                    if runtodo=='bwsbclrdc27':
                        pstartno=440
                    if runtodo=='bwsbclrdc28':
                        pstartno=505
                    if runtodo=='bwsbclrdc29':
                        pstartno=370
                print 'labelneed', labelneed
                if wanted=='dirgammasfr':
                    yaxisl = np.absolute(Lgamma_sfr)
                    if havecr>4:
                        if (M1labelneed==1 and color=='r'):
                            labelneed='_nolegend_'
                    elif havecr>0:
                        yaxisl = np.absolute(Lgcal_sfr)
                    if correctIa==1 and havecr>0:
                        yaxisl = np.absolute(Lgamma_sfr_noIa)
                if wanted=='dirgamma':
                    if havecr>4:
                        yaxisl = np.absolute(Lgamma)
                    elif havecr>0:
                        yaxisl = np.absolute(Lgcalist[1:])
                if wanted=='dirsfr':
                    yaxisl = avesfrl
                    if M1labelneed==1 and color=='r':
                        labelneed='_nolegend_'
                    if showstarburst==1 and runtitle=='SBC':
                        plt.axvline(x=pstartno, color=color,ls=lsn,lw=1)
                if wanted=='dirsm':
                    if M1labelneed==1 and color=='r':
                        if normalizedsm==1:
                            plt.plot(xaxisl, sml[:-1]-sml[0], color=color,ls=lsn,lw=2)
                        else:
                            plt.plot(xaxisl, sml[1:], color=color,ls=lsn,lw=2)
                    else:
                        if normalizedsm==1:
                                plt.plot(xaxisl, sml[:-1]-sml[0], color=color,ls=lsn,lw=2, label=labelneed)
                        else:
                                plt.plot(xaxisl, sml[1:], color=color,ls=lsn,lw=2, label=labelneed)
                    if showstarburst==1 and runtitle=='SBC':
                            plt.axvline(x=pstartno, color=color,ls='dashed',lw=1)
                if wanted=='dirgammasfr':
                    plt.axhline(y=0.0002,ls='--',color='k')
                    if runtitle=='SMC' and legendneed==1:
                        plt.legend(loc='best', fontsize=8)
                    elif legendneed==1 and M1runlabelneed==1:
                        plt.legend(loc='best', fontsize=8,ncol=2)
                    ylab=r'$L_{\gamma}/L_{\rm SF}$'
                    figname=plotloc+'CRplot/dirgammasfr/gammasfrsnap_'+fmeat+'.pdf'
                if wanted=='dirgamma':
                    plt.yscale('log')
                    if runtitle=='SMC' or legendneed==1:
                            plt.legend(loc='best', fontsize=8,ncol=3)
                    if legendneed==1 and M1runlabelneed:
                            plt.legend(loc='best', fontsize=8,ncol=2)
                    plt.xlabel(xlab, fontsize=16)
                    plt.ylabel(r'$L_{\gamma} {\rm erg/s}$', fontsize=16)
                    figname=plotloc+'CRplot/dirgamma/gammasnap_'+fmeat+'.pdf'
                if wanted=='dirsfr':
                    plt.yscale('log')
                    if runtitle=='SMC' or legendneed==1:
                        plt.legend(loc='best', fontsize=8,ncol=3)
                    if M1labelneed==1 and legendneed==1:
                        plt.legend(loc='best', fontsize=8,ncol=2)
                    plt.xlabel(xlab, fontsize=16)
                    plt.ylabel(r'${\rm SFR (M_{\odot}/yr)} $', fontsize=16)
                    figname=plotloc+'CRplot/dirsfr/sfrsnap_'+fmeat+'.pdf'
                if wanted=='dirsm':
                    if runtitle=='SMC' and legendneed==1:
                        if strlabelneed==1 or M1labelneed==1:
                            plt.legend(loc='best', fontsize=10,ncol=2)
                        else:
                            plt.legend(loc='best', fontsize=8,ncol=3)
                    plt.xlabel('Myr', fontsize=16)
                    plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
                    if normalizedsm==1:
                        plt.ylabel(r'${\rm M_{*,new} (M_{\odot})} $', fontsize=16)
                    else:
                        plt.ylabel(r'${\rm M_* (M_{\odot})} $', fontsize=16)
                    if normalizedsm==1:
                        figname=plotloc+'CRplot/dirsm/nsm_'+fmeat+'.pdf'
                    else:
                        figname=plotloc+'CRplot/dirsm/sm_'+fmeat+'.pdf'
        plotdict[wanted]['xlab'] = xlab
        plotdict[wanted]['ylab'] = ylab
        plotdict[wanted]['xnl'] = xaxisl
        plotdict[wanted]['ynl'] = thdenlist[:-1]
        plotdict[wanted]['lw'] = 2    
        plotdict[wanted]['marker'] = 'None'    
        plotdict[wanted]['linelab'] = 'thermal'
        plotdict[wanted]['runtodo'] = runtodo
        plotdict[wanted]['labelneed'] = labelneed
        plotdict[wanted]['lsn'] = lsn
        plotdict[wanted]['color'] = color
        plotdict[wanted]['runtitle'] = runtitle
        plotdict[wanted]['ptitle'] = ptitle
        plotdict[wanted]['filename'] = figname
        return plotdict
