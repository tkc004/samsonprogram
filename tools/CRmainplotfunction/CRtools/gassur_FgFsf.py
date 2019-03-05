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

#dirneed=['bwsmclrdc0','bwsmclrdc27','bwsmclrdc28','bwsmclrdc29']

#dirneed=['bwmwlrdc27','bwmwlrdc28','bwmwlrdc29','bwmwlrdc28mhd']

dirneed=['bwsmclrstr','bwsmclrdc0','bwsmclrdc27','bwsmclrdc28','bwsmclrdc29','bwsmclrdc28mhd','bwsmclrdc28str',\
'bwmwmrdc0','bwmwmrdc27','bwmwmrdc28','bwmwmrdc29','bwmwmrstr','bwmwmrdc28mhd','bwmwmrdc28str',\
'bwsbclrdc0','bwsbclrdc27','bwsbclrdc28','bwsbclrdc29','bwsbclrstr','bwsbclrdc28mhd','bwsbclrdc28str']
#]

fmeat=dirneed[-1]

startno=400

Nsnap=501

snapsep=10

wanted='gassur_FgFsf'

title='MW'
ptitle='LSG'
nolegend=0


print 'wanted', wanted
print 'fmeat', fmeat
print 'runtodo', dirneed


if wanted=='gassur_FgFsf':
    import time
    rcParams['figure.figsize'] = 6,6
    varR = 0
    sepave=0
    useFcal=1
    atstarburst=0
    sbgassurmin=0.08
    newlabelneed=1
    if varR==1:
        radgas=0.0
    else:
        radgas=1.0 #kpc 
    Rfrac=0.03
    for runtodo in dirneed:
        gassurl=[]
        gassurxl=[]
        gassuryl=[]
        gassurzl=[]
        info=outdirname(runtodo, Nsnap)
        havecr=info['havecr']
        dclabel=info['dclabel']
        runtitle=info['runtitle']
        print 'havecr', havecr
        print 'runtitle', runtitle
        haveB=info['haveB']
        cosmo=info['cosmo']
        color=info['color']
        correctIa=info['correctIa']
        if havecr==0:
            continue
        runtitle=info['runtitle']
        newlabel=info['newlabel']
        starburst=0
        shiftz=0
        mks=10
        pstartno = startno
        pNsnap = Nsnap
        psnapsep = snapsep
        crshiftz=1
        marker='^'
        if runtitle=='SMC':
            physicalR=3.0
            shiftz=1
            marker='o'
            radgas=2.0
        elif runtitle=='MW':
            physicalR=10.0
            shiftz=1
            marker='^'
            mks=14
            radgas=4.0
        elif runtitle=='SBC':
            physicalR=10.0 
            starburst=1
            pstartno=240
            pNsnap=650
            psnapsep=5
            radgas=0.25
            if atstarburst==1:
                psnapsep=1
                if runtodo=='bwsbclr':
                    pNsnap=528
                    pstartno=518
                if runtodo=='bwsbclrdc0':
                    pNsnap=594
                    pstartno=584
                if runtodo=='bwsbclrdc27':
                    pNsnap=143
                    pstartno=133
                if runtodo=='bwsbclrdc28':
                    pNsnap=525
                    pstartno=515
                if runtodo=='bwsbclrdc29':
                    pNsnap=438
                    pstartno=428
                if runtodo=='bwsbclrstr':
                    pNsnap=250
                    pstartno=240
                if runtodo=='bwsbclrdc28mhd':
                    pNsnap=563
                    pstartno=553
                if runtodo=='bwsbclrdc28str':
                    pNsnap=375
                    pstartno=365
            shiftz=1
            marker='D'
        if cosmo == 1:
            physicalR=6.0
            if runtodo=='m10qcr_b_70':
                physicalR=3.0
        outdata = dirout(runtodo,pstartno,pNsnap,psnapsep, Rfrac,physicalR=physicalR,shiftz=crshiftz)
        avesfrl = outdata['avesfrl']    
        Lsfr = outdata['Lsfr']
        if havecr>4:
            Lgamma_sfr=outdata['Lgamma_sfr']
            Lgamma = outdata['Lgamma']
        elif havecr>0:
            Lgamma_sfr=outdata['Lgcal_sfr']
            Lgamma = outdata['Lgcal']
        if useFcal==1:
            Lgamma_sfr=outdata['Lgcal_sfr']
            Lgamma = outdata['Lgcal']
        if correctIa==1 and havecr>0:
            Lgamma_sfr=outdata['Lgamma_sfr_noIa']
        print 'Lgamma_sfr', Lgamma_sfr
        Lgamma_sfr = Lgamma_sfr[np.isfinite(Lgamma_sfr)]
                #if starburst==1:
                #        sfrcut=avesfrl>6.0
        #   Lgamma_sfr=Lgamma_sfr[sfrcut]
        if starburst==0:
            prol = ['x','y','z']
            for projection in prol:
                outgas = gassurlout(runtodo,pstartno,pNsnap,psnapsep, radgas,projection,shiftz=shiftz)
                gassurl = np.concatenate((gassurl,outgas['gassurl']))
        else:
            outgas = gassurlout(runtodo,pstartno,pNsnap,psnapsep, radgas,'x',shiftz=shiftz)
            gassurl = outgas['gassurl']
            outgasy = gassurlout(runtodo,pstartno,pNsnap,psnapsep, radgas,'y',shiftz=shiftz)
            gassuryl = outgasy['gassurl']
            outgasz = gassurlout(runtodo,pstartno,pNsnap,psnapsep, radgas,'z',shiftz=shiftz)
            gassurzl = outgasz['gassurl']
            gassurcut =sbgassurmin
            print 'radgas', radgas
            print 'gassurl', gassurl
            print 'Lgamma_sfr', Lgamma_sfr
            print 'gassurl.shape', gassurl.shape
            print 'Lgamma_sfr.shape', Lgamma_sfr.shape
            #gassurlold=np.fmax(gassurl,gassuryl)
            #gassurlold=np.fmax(gassurlold,gassurzl)
            gassurlold = (gassurl+gassuryl+gassurzl)/3.*2.
            sbcut=gassurlold>gassurcut
            gassurl=gassurlold[sbcut]
            gassurzl=gassurzl[sbcut]
            gassuryl=gassuryl[sbcut]
            sbcut=gassurlold[1:]>gassurcut
            Lgamma_sfr=Lgamma_sfr[sbcut]
            #gassurl = np.append(gassurl,outgas['gassurl'])
            gassurl=np.concatenate((gassurl,gassurzl,gassuryl))
        print 'gassurl', gassurl
        print 'Lgamma_sfr', Lgamma_sfr
        if len(Lgamma_sfr)==0:
            continue
        print 'correctIa', correctIa
        print 'Lgamma_sfr',Lgamma_sfr
        gmed = np.median(gassurl)
        g1sigu = np.percentile(gassurl,84)-gmed
        g1sigd = gmed-np.percentile(gassurl,15.8)
        Lmed = np.median(Lgamma_sfr)
        L1sigu = np.percentile(Lgamma_sfr,84)-Lmed
        L1sigd = Lmed-np.percentile(Lgamma_sfr,15.8)
        print 'Lmed', Lmed
        if sepave==1:
            Lmed = np.median(Lgamma)/np.median(Lsfr)
        print 'gmed', gmed, g1sigu, g1sigd
        print 'Lmed', Lmed, L1sigu, L1sigd
        labelneed=dclabel
        if newlabelneed==1:
            labelneed="\n".join(wrap(newlabel,17))
        print 'labelneed', labelneed
        if haveB>0:
        #if runtodo=='bwmwmrdc28mhd' or runtodo=='bwsmclrdc28mhd' or runtodo=='bwsbclrdc28mhd':
            fillstyle='top'
        else:
            fillstyle='full'
        if not ( runtitle=='MW' or runtitle=='COSMO'):
            labelneed=''
        plt.errorbar(gmed, Lmed, xerr=[[g1sigd],[g1sigu]], yerr=[[L1sigd],[L1sigu]]\
,color=color,fmt=marker,markersize=7,label=labelneed, fillstyle=fillstyle)
        #else:
        #   plt.errorbar(gmed, Lmed, xerr=[[g1sigd],[g1sigu]], yerr=[[L1sigd],[L1sigu]],color=color,fmt=marker,markersize=7)
    obsdatal = [[-2.29505,-2.60616,-1.99283,-5.07909,-5.12171,-5.01514]\
    #NGC 4945
    ,[-0.721246,-0.721246,-0.721246,-3.762,-3.93339,-3.63943]\
    #NGC 1068
    ,[-1.69897,-1.69897,-1.69897,-4.41986,-4.56119,-4.3134]\
    ,[-2.99719,-2.99719,-2.99719,-5.05151,-5.15809,-4.96626]\
    ,[-2.51958,-2.572,-2.484, -5.63,-5.85677,-5.41]\
    ,[-2.7049,-2.7049,-2.7049,-5.40359,-5.52615,-5.30768]\
    ,[-0.81742,-0.81742,-0.81742,-4.37769,-4.62812,-4.21782]\
    ,[-0.59429,-0.59429,-0.59429,-4.07425,-4.21278,-3.96768]]
    stardatal =[[-0.81742,-0.81742,-0.81742,-4.07394,-4.18584,-3.99403]]
    for data in obsdatal:
        plt.errorbar(np.power(10,data[0]), np.power(10,data[3]),\
        xerr=[[np.power(10,data[0])-np.power(10,data[1])],[np.power(10,data[2])-np.power(10,data[0])]],\
        yerr=[[np.power(10,data[3])-np.power(10,data[4])],[np.power(10,data[5])-np.power(10,data[3])]],\
        fmt='s', color='0.5',ls='dashed', mfc='0.9',markersize=6)
    for data in stardatal:
        plt.errorbar(np.power(10,data[0]), np.power(10,data[3]),\
        xerr=[[np.power(10,data[0])-np.power(10,data[1])],[np.power(10,data[2])-np.power(10,data[0])]],\
        yerr=[[np.power(10,data[3])-np.power(10,data[4])],[np.power(10,data[5])-np.power(10,data[3])]],\
        fmt='*', color='0.5',ls='dashed',  mfc='0.9',markersize=mks)
    ftxt = open(programdir+'/data/LTQ_line.txt', 'r')
    ftxt.readline()
    dars = ftxt.readlines()
    ftxt.close()
    LTQx=[]
    LTQy=[]
    for line in dars:
        xsd = line.split()
        LTQx=np.append(LTQx, float(xsd[0]))
        LTQy=np.append(LTQy, float(xsd[1]))

        ftxt = open(programdir+'/data/LTQ_up.txt', 'r')
        ftxt.readline()
        dars = ftxt.readlines()
        ftxt.close()
        LTQupx=[]
        LTQupy=[]
        for line in dars:
                xsd = line.split()
                LTQupx=np.append(LTQupx, float(xsd[0]))
                LTQupy=np.append(LTQupy, float(xsd[1]))

        ftxt = open(programdir+'/data/LTQ_down.txt', 'r')
        ftxt.readline()
        dars = ftxt.readlines()
        ftxt.close()
        LTQdownx=[]
        LTQdowny=[]
        for line in dars:
                xsd = line.split()
                LTQdownx=np.append(LTQdownx, float(xsd[0]))
                LTQdowny=np.append(LTQdowny, float(xsd[1]))
    LTQfy = np.interp(LTQupx,LTQdownx, LTQdowny)
    plt.plot(np.power(10,LTQx),np.power(10,LTQy), color='k')
    plt.fill_between(np.power(10,LTQupx), np.power(10,LTQupy), np.power(10,LTQfy),color='0.1',facecolor='0.1',alpha=0.2)
    import matplotlib.lines as mlines
    lackiline = mlines.Line2D([], [], color='k', ls='solid',label='')
    lackigray = mlines.Line2D([], [], lw=9, ls='solid',color='0.1',alpha=0.2)
    obsbar = plt.errorbar([],[], fmt='s', color='0.5', mfc='0.9',markersize=6)
    dwarfpoint = plt.errorbar([],[], fmt='o', color=cmaps.plasma(0.01), mfc=cmaps.plasma(0.01),markersize=6)
    lstarpoint = plt.errorbar([],[], fmt='^', color=cmaps.plasma(0.01), mfc=cmaps.plasma(0.01),markersize=6)
    sbpoint    = plt.errorbar([],[], fmt='D', color=cmaps.plasma(0.01), mfc=cmaps.plasma(0.01),markersize=6)
    lxlist = [(lackiline,lackigray),obsbar,dwarfpoint,lstarpoint,sbpoint]
    dclablist = ['Lacki+11 models','Observations','Dwarf',r'L$\star$ Galaxy', 'Starburst']
    legend1 = plt.legend(lxlist, dclablist, loc=2,fontsize=8,ncol=3)
    plt.gca().add_artist(legend1)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(xmax=1.0)
    plt.axhline(y=0.0002,ls='--',color='k')
    plt.xlabel(r'$\Sigma_{\rm g}\;{\rm [g/cm^2]}$', fontsize=16)
    plt.ylabel(r'${\rm L_{\gamma}/L_{\rm SF}}$', fontsize=18)
    plt.legend(loc=4,fontsize=10,ncol=2, numpoints=1)
    figname=plotloc+'CRplot/gassur_FgFsf/gassur_FgFsf_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    print 'figname', figname
    plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
    plt.savefig(figname,bbox_inches='tight')
    plt.clf()
