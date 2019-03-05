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


dirneed=['bwsmclrdc0','bwsmclrdc27','bwsmclrdc28','bwsmclrdc29','bwsmclrstr','bwsmclrdc28mhd','bwsmclrdc28str',\
 'bwmwmrdc0','bwmwmrdc27','bwmwmrdc28', 'bwmwmrdc29','bwmwmrstr','bwmwmrdc28mhd','bwmwmrdc28str',\
 'bwsbclrdc0','bwsbclrdc27', 'bwsbclrdc28','bwsbclrdc29','bwsbclrstr','bwsbclrdc28mhd','bwsbclrdc28str']

fmeat=dirneed[-1]

startno=400

Nsnap=501

snapsep=10

wanted='SFRLg'

title='MW'
ptitle='LSG'
nolegend=0


print 'wanted', wanted
print 'fmeat', fmeat
print 'runtodo', dirneed

if wanted=='SFRLg':
        import time
        rcParams['figure.figsize'] = 6,6
        varR = 0
        sepave=0
        useFcal=1
        atstarburst=1
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
                        pstartno=300
                        pNsnap=600
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
                #       Lgamma_sfr=Lgamma_sfr[sfrcut]
                print 'Lgamma_sfr', Lgamma_sfr
                if len(Lgamma_sfr)==0:
                        continue
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
                if starburst==1:
                    avesfrold = avesfrl
                    sfrcut=avesfrold>3.0
                    avesfrl = avesfrl[sfrcut]
                    Lgamma=Lgamma[:-1]
                    Lgamma = Lgamma[sfrcut]
                SFRmed = np.median(avesfrl)
                SFR1sigu = np.percentile(avesfrl,84)-SFRmed
                SFR1sigd = SFRmed-np.percentile(avesfrl,15.8)
                Lmed = np.median(Lgamma)
                L1sigu = np.percentile(Lgamma,84)-Lmed
                L1sigd = Lmed-np.percentile(Lgamma,15.8)
                print 'SFRmed, Lmed', SFRmed, Lmed
                plt.errorbar(SFRmed, Lmed, xerr=[[SFR1sigd],[SFR1sigu]], yerr=[[L1sigd],[L1sigu]]\
                ,color=color,fmt=marker,markersize=7,label=labelneed, fillstyle=fillstyle)
        obsdatal = [[0.2,0.2,0.2,1.7e37,1.3e37,2.1e37]\
        #LMC
        ,[0.068,0.036,0.1,4.3e36,4.3e36,4.3e36]\
        #SMC
        ,[2.0,1.2,3.4,6.5e38,6.5e38,6.5e38]\
        #MW
        ,[1.0,1.0,1.0,1.9e38,1.5e38,2.3e38]\
    #M31
        ,[6.18,6.18,6.18,2e40,1.6e40,2.4e40]\
    #M82
        ,[2.82,2.82,2.82,9.5e39,6.0e39,14.0e39]\
    #NGC253
        ,[3.49,3.49,3.49,1.7e40,1.2e40,2.2e40]\
    #NGC4945
        ,[37,37,37,8e40,6e40,10e40]]
    #NGC1068
    
        ftxt = open(programdir+'/data/nondetectGammaray.txt', 'r')
        ftxt.readline()
        dars = ftxt.readlines()
        ftxt.close()
        SFRlnd=[]
        LGRl=[]


        for line in dars:
            xsd = line.split()
            SFRlnd =np.append(SFRlnd, float(xsd[0]))
            LGRl=np.append(LGRl, float(xsd[1]))
        
        plt.plot(np.power(10.,SFRlnd),np.power(10.,LGRl),markersize=4,mfc='None'\
                 ,markeredgecolor='red', marker='s',ls='None')

        for data in obsdatal:
                plt.errorbar(data[0], data[3],\
                xerr=[[data[0]-data[1]],[data[2]-data[0]]],\
                yerr=[[data[3]-data[4]],[data[5]-data[3]]],\
                fmt='s', color='0.5',ls='dashed', mfc='0.9',markersize=6)

        ndpoint = plt.errorbar([],[], fmt='s', color='red', mfc='None',markersize=6)        
        obsbar = plt.errorbar([],[], fmt='s', color='0.5', mfc='0.9',markersize=6)
        dwarfpoint = plt.errorbar([],[], fmt='o', color='g', mfc='g',markersize=6)
        lstarpoint = plt.errorbar([],[], fmt='^', color='g', mfc='g',markersize=6)
        sbpoint    = plt.errorbar([],[], fmt='D', color='g', mfc='g',markersize=6)
        lxlist = [ndpoint,obsbar,dwarfpoint,lstarpoint,sbpoint]
        dclablist = ['Non detection','Observations','Dwarf',r'L$\star$ Galaxy', 'Starburst']
        legend1 = plt.legend(lxlist, dclablist, loc=2,fontsize=10,ncol=3)
        plt.gca().add_artist(legend1)
        #calorimetric line:
        sfrref = np.power(10,np.linspace(-3,1.7))
        Lgcal = sfrref*6e39
        plt.plot(sfrref,Lgcal,color='k',ls='dashed')
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel(r'${\rm SFR}\;{\rm [M_\odot/yr]}$', fontsize=20)
        plt.ylabel(r'$L_{\gamma}$', fontsize=20)
        plt.legend(loc=4,fontsize=10,ncol=2, numpoints=1)
        figname=plotloc+'CRplot/SFRLg/SFRLg_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        print 'figname', figname
        plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=18)
        plt.savefig(figname,bbox_inches='tight')
        plt.clf()

