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
dirneed=['bwmwmr','bwmwmrdc28', 'bwmwmrdc29']
#dirneed=['bwmwmr','bwmwmrdc27','bwmwmrdc28', 'bwmwmrdc29']
#dirneed=['bwmwmr']
#dirneed=['bwsmclrdc28']
#dirneed=['bwsmclrdc0']
#dirneed=['bwsbclrdc27']
#dirneed=['bwmwmrdc28']
#dirneed=['bwmwmrdc29']
#dirneed=['bwmwmrmhd','bwmwmrstr','bwmwmrdc28mhd','bwmwmrdc28str']
#dirneed=['m11bcr_700']
#dirneed=['m11bcr_b_70']
#dirneed=['f46']
#dirneed=['fm12m']
#dirneed=['m12imhdcv']
#dirneed=['m12icr_700']
#dirneed=['m12icr_b_70']
#dirneed=['m11gcr_b_70']
#dirneed=['m11bmhdcv']
#dirneed=['m12fcr_b_70']
#dirneed=['m12fcr_700']
#dirneed=['m12fmhdcv']
#dirneed=['m11gmhdcv']
#dirneed=['m12mcr_b_70']
#dirneed=['m12mcr_700']
#dirneed=['m12mmhdcv']
#dirneed=['bwsbclrdc0']
#dirneed=['bwsmclrdc0','bwsmclrdc27','bwsmclrdc28','bwsmclrdc29']
#dirneed=['bwsmclrstr','bwsmclrdc28mhd','bwsmclrdc28str']
#dirneed=['bwmwmrdc29']
#dirneed=['bwsbclrdc29']
#dirneed=['bwsmclrdc29']
#dirneed=['bwsmclrdc0']
#dirneed=['bwsmclrdc0','bwsmclrdc27','bwsmclrdc28','bwsmclrdc29']
#dirneed=['bwsmclrdc27']
#dirneed=['bwmwmrmhd']
#dirneed=['bwmwmrdc28mhd']
#dirneed=['bwsmclrdc0','bwsmclrdc27','bwsmclrdc28','bwsmclrdc29','bwsmclrstr','bwsmclrdc28mhd','bwsmclrdc28str']
#dirneed=['bwmwmrdc0','bwmwmrdc27','bwmwmrdc28', 'bwmwmrdc29']
#dirneed=['bwmwmrstr','bwmwmrdc28mhd','bwmwmrdc28str']
#dirneed=['bwmwmrdc0','bwmwmrdc27','bwmwmrdc28', 'bwmwmrdc29','bwmwmrstr','bwmwmrdc28mhd','bwmwmrdc28str']
#dirneed=['bwmwmr','bwmwmrdc0','bwmwmrdc27','bwmwmrdc28', 'bwmwmrdc29','bwmwmrmhd','bwmwmrstr','bwmwmrdc28str']
#dirneed=['bwmwmrdc27']
#dirneed=['bwmwmrdc28']
#dirneed=['bwmwmrdc0','bwmwmrdc27','bwmwmrdc28', 'bwmwmrdc29']
#dirneed=['bwsbclrdc0','bwsbclrdc27','bwsbclrdc28','bwsbclrdc29']
#dirneed=['bwsbclrstr','bwsbclrdc28mhd','bwsbclrdc28str']
#dirneed=['bwsbclr','bwsbclrdc28']
#dirneed=['bwsbclrdc29']
#dirneed=['bwsmclrdc29']
#dirneed=['bwsmclr','bwsmclrdc28']
#dirneed=['bwmwmrdc28', 'bwmwmrdc29']
#dirneed=['bwmwmr','bwmwmrdc28', 'bwmwmrdc29']
#dirneed=['bwmwmr','bwmwmrdc0','bwmwmrdc27','bwmwmrdc28', 'bwmwmrdc29']
#dirneed=['mw_mr_2_21','bwmwmrdc0','bwmwmrdc28']
#dirneed=['bwsmclrstr']
#dirneed=['bwmwlrdc27ds']
#dirneed=['bwsmclr','bwsmclrdc0','bwsmclrdc27','bwsmclrdc28','bwsmclrdc29']
#dirneed=['bwsmclrstr','bwsmclrdc0','bwsmclrdc27','bwsmclrdc28','bwsmclrdc29','bwsmclrdc28mhd','bwsmclrdc28str',\
#'bwmwmrdc0','bwmwmrdc27','bwmwmrdc28','bwmwmrdc29','bwmwmrstr','bwmwmrdc28mhd','bwmwmrdc28str',\
#'bwsbclrdc0','bwsbclrdc27','bwsbclrdc28','bwsbclrdc29','bwsbclrstr','bwsbclrdc28mhd','bwsbclrdc28str']
#]
#dirneed=['bwsbclrdc0']
#dirneed=['bwsbclrdc28mhd']
#dirneed=['bwsbclrdc0','bwsbclrdc27','bwsbclrdc28','bwsbclrdc29']
#dirneed=['bwsbclrdc27']
#dirneed=['bwmwlrmhd']
#dirneed=['bwmwmrdc28']
#dirneed=['bwsmcmr']
#dirneed=['bwsmclr','bwsmclrdc0','bwsmclrdc27','bwsmclrdc28','bwsmclrdc29']
#dirneed=['bwsmclr','bwsmclrdc0','bwsmclrdc27','bwsmclrdc28','bwsmclrdc29','bwsmclrmhd','bwsmclrstr','bwsmclrdc28str']
#dirneed=['bwsmclrdc27','bwsmclrdc28','bwsmclrdc29','bwsmclrstr']
#dirneed=['bwsmclrmhd','bwsmclrstr','bwsmclrdc28mhd','bwsmclrdc28str']
#dirneed=['bwsmclrmhd','bwsmclrdc0','bwsmclrstr']
#dirneed=['bwsmclr','bwsmclrdc27','bwsmclrdc28','bwsmclrdc29']
#dirneed=['bwsmclrdc0', 'bwsmclrdc27','bwsmclrdc28','bwsmclrdc29']
#dirneed=['bwmwlrdc27ehh', 'bwmwlrdc28ehh']
#dirneed=['bwmwlrstr']
#dirneed=['bwmwlrmhd','bwmwlrstr','bwmwlrdc28str','bwmwlrdc28mhd']
#dirneed=['bwmwlr','mw_cr_lr_dc0_2_21','bwmwlrdc27','bwmwlrdc28','bwmwlrdc29','bwmwlrstr','bwmwlrdc28str']
#dirneed=['mw_cr_lr_dc0_2_21','bwmwlrdc27','bwmwlrdc28','bwmwlrdc29','bwmwlrmhd','bwmwlrstr','bwmwdc28mhd','bwmwlrdc28str']
#dirneed=['bwmwlrdc28str']
#dirneed=['bwsmclr','bwsmclrdc0','bwsmclrdc27','bwsmclrdc28','bwsmclrdc29','bwsmclrstr','bwsmclrdc28str']
#dirneed=['bwmwlrdc27ehh']
#dirneed=['bwmwlr']
#dirneed=['bwmwlrdc29']
#dirneed=['bwmwlrdc28']
#dirneed=['bwmwlr','bwmwlrdc27','bwmwlrdc28', 'bwmwlrdc29']
#dirneed=['mw_cr_lr_dc0_2_21','bwmwlrdc27','bwmwlrdc28', 'bwmwlrdc29']
#dirneed=['bwmwlr','mw_cr_lr_dc0_2_21','bwmwlrdc27','bwmwlrdc28', 'bwmwlrdc29','bwmwlrmhd', 'bwmwlrstr','bwmwlrdc28mhd','bwmwlrdc28str']
#dirneed=['bwmwlrdc28mhd']
#dirneed=['mw_cr_lr_dc0_2_21']
#dirneed=['bwmwlrmhd','bwmwlrstr','bwmwlrdc28str']
#dirneed=['bwsmclrdc27']
#dirneed=['mw_mr_2_21_IC']
#dirneed=['mw_cr_mr_dc27_2_21_M1']
#dirneed=['bwmwmrdc28']
#dirneed=['bwsbclrdc28']
#dirneed=['bwsbclr','bwsbclrdc0','bwsbclrdc27','bwsbclrdc28','bwsbclrdc29']
#dirneed=['m12icr_b_70','m12imhdcv']
#dirneed=['m12imhdcv']
#dirneed=['m10qmhdcv','m10qcr_b_70']
#dirneed=['m10qcr_b_70','m11bcr_b_70','m12icr_b_70','m11fcr_b_70']
#dirneed=['m10qmhdcv','m11bmhdcv','m12imhdcv','m11fmhdcv']
#dirneed=['m09cr_b_70']
#dirneed=['m11bcr_b_70','m11dcr_b_70','m12icr_b_70']
#dirneed=['m11bcr_b_70']
#dirneed=['m11bcr_700']
#dirneed=['m11bcr_700','m11dcr_700','m12icr_700']
#dirneed=['m09cr_b_70','m10qcr_b_70','m11bcr_b_70','m11dcr_b_70','m11gcr_b_70','m11hcr_b_70','m12icr_b_70','m11fcr_b_70']
#dirneed=['m11bcr_700','m11dcr_700','m12icr_700']
#dirneed=['m11bcr_b_70','m11dcr_b_70','m11fcr_b_70','m11gcr_b_70','m11hcr_b_70','m12icr_b_70']
#dirneed=['m11bmhdcv']
#dirneed=['m11bcr_700']
#dirneed=['bwsmclrstr','bwsmclrstrll','smclrstrll4va']
#dirneed=['m12imhdcv']
#dirneed=['m10qcr_b_70']
#dirneed=['m11bcr_700','m11dcr_700','m12icr_700']
#dirneed=['m12icr_700']
#dirneed=['m10vcr_700', 'm11bcr_700','m11dcr_700','m11fcr_700', 'm11gcr_700', 'm11hcr_700','m12icr_700']
#dirneed=['m10qcr_b_70','m10vcr_b_70','m11bcr_b_70','m11dcr_b_70','m11gcr_b_70','m11hcr_b_70','m12icr_b_70','m11fcr_b_70']
#dirneed=['bwsmclrstr','bwsmclrstrll','smclrstrllva','smclrstrll4va']
fmeat=dirneed[-1]
#fmeat='m12i'
#fmeat='mwmr3kpc'
#fmeat='testleg'
#fmeat='sbc'
#fmeat='sbcmhd'
#fmeat='mw'
#fmeat='mwmhd'
#fmeat='smcmhd'
#fmeat='smc'
#fmeat='5Myr'
#fmeat='bwsbclr'
#fmeat='bwsmclrdc29'
#fmeat='bwsmclrdc28'
#fmeat='bwsbclrdc0'
#fmeat='bwsbclrdc28mhd'
#fmeat='bwmwmrdc29'
#fmeat='bwsmclrdc27'
#fmeat='bwmwmrdc28mhd'
#fmeat='bwmwmrdc27'
#fmeat='bwmwmrdc0'
#fmeat='gasdenxyz'
#fmeat='gasdenxyz2'
#fmeat='sbcgasdenz'
#fmeat='sbcmhd'
#fmeat='smcmhd'
#fmeat='m10vcr_700'
#fmeat = 'm10vcr_b_70'
#fmeat = 'smclrstrll4va'
#fmeat='m12icr700r7_9h1'
#fmeat='m12imhdcvr10h1'
#fmeat='m11bcr_700r5h1'
#fmeat='m12imhdcv'
#fmeat='cosmo'
#fmeat='cosmo1kpc'
#fmeat='cosmo6kpctest'
#fmeat='cosmo6kpcgas1kpczdc28'
#fmeat='cosmo6kpcgas1kpcxdc29'
#fmeat='cosmo6kpcgas1kpczdc29'
#fmeat='smclrstrll4va'
#fmeat='sbcdc29'
#fmeat='cosmo_m09'
#fmeat = 'mwmrdc28'
#fmeat = 'mwmrinner'
#fmeat = 'mwmhd_Fcal'
#fmeat='smcmhd20kpc'
#fmeat='mwmrmhd5_20kpc'
#fmeat='mwmhd79kpc'
#fmeat='test'
#fmeat='mw35kpc'
#fmeat='mw79kpc'
#fmeat='sbc'
#fmeat='smcmhd'
#fmeat='smcsmall'
#fmeat='mwmr'
#fmeat='mwmrsmall'
#fmeat='mwmrmhd'
#fmeat='mwdc27ds'
#fmeat = 'test'
#fmeat='atstarburst'
#fmeat='sbc_totke'
#fmeat='mwmhdnocr'
#fmeat='mwtest'
#fmeat='smcinner'
#fmeat = 'mwinner'
#fmeat = 'mwmhdinner'
#fmeat = 'mwmhdsolarcircle'
#fmeat='mwmhd'
#fmeat='mwlr'
#fmeat='smc'
#fmeat='m12icr'
#fmeat='mwmrmhdtest'
#fmeat='cosmo'
#fmeat='testave'
#fmeat='testLTQ'
#fmeat='idealsmallgamma'
#fmeat='ideal1kpc'
#fmeat='smcmhd'
#fmeat='smcdc28'
#fmeat='mwehh'
#startno=0
#startno=400
#startno=440
#startno=450
#startno=480
#startno=550
#startno=450
#startno=200
#startno=400
#startno=550
#startno=600
#startno=580
#startno=590
startno=500
#Nsnap=300
#Nsnap=601
#Nsnap=651
#Nsnap=501
#Nsnap=651
#Nsnap=190
#Nsnap=600
Nsnap=501
#Nsnap=10
#snapsep=50
#snapsep=20
snapsep=10
#wanted='SFRLg'
#wanted='crdata'
#wanted='dpcrz_rhog'
#wanted='dpcrzmrhog'
#wanted='dGmden'
#wanted='dpcrz_dpth'
#wanted='vcr_vth'
#wanted='dpcrz'
wanted='rhog'
#wanted='pcr_pth'
#wanted='gassurden'
#wanted='eturdentime'
#wanted='pretime'
#wanted='etur_mtime'
#wanted='eturtime'
#wanted='crturmap'
#wanted='gassur_FgFsf'
#wanted='cratmult'
#title='SBG'
#title='Cosmo'
#title=''
#title='Dwarf'
#t='m12imhdcvr10h1'
title='MW'
title='LSG'
nolegend=0


print 'wanted', wanted
print 'fmeat', fmeat
print 'runtodo', dirneed




if wanted=='crdata':
    if title=='Dwarf':
        sfrten = 1e-2
        Smten = 1e6
        tdten = 1e8
        dcgten = 1e8
        Lgten = 1e37
    elif title=='MW':
        sfrten=1.0
        Smten = 1e8
        tdten = 1e9
        dcgten = 1e9
        Lgten = 1e39
        
    runtodol=[]
    Smnewl=[]
    totdiskgml=[]
    diskcoldgml=[]
    Lgammal = []
    z12gasl = []
    r12crl = []
    z12dengasl = []
    gassurdenl = []
    avesfrl = []
    crdenl = []
    for runtodo in dirneed:
        info=outdirname(runtodo, Nsnap)
        dclabel=info['dclabel']
        crdata=outcrdata(runtodo, Nsnap, startno,snapsep)
        Smnew = crdata['Smnew']
        totdiskgm = crdata['totdiskgm']
        diskcoldgm = crdata['diskcoldgm']
        Lgamma = np.absolute(crdata['Lgamma'])
        z12gas = crdata['z12gas']
        r12cr = crdata['r12cr']
        z12dengas = crdata['z12dengas']
        gassurden = crdata['gassurden']
        avesfr = crdata['avesfr']
        crden = crdata['crden']
        #print 'r12cr', r12cr
        runtodol = np.append(runtodol,dclabel)
        Smnewl = np.append(Smnewl, Smnew)
        totdiskgml = np.append(totdiskgml,totdiskgm)
        diskcoldgml = np.append(diskcoldgml,diskcoldgm)
        Lgammal = np.append(Lgammal,Lgamma)
        z12gasl = np.append(z12gasl,z12gas)
        r12crl = np.append(r12crl,r12cr)
        z12dengasl = np.append(z12dengasl,z12dengas)
        gassurdenl = np.append(gassurdenl,gassurden)
        avesfrl = np.append(avesfrl, avesfr)
        crdenl = np.append(crdenl,crden)
    for i in range(len(dirneed)):
        Smnewround = np.around((Smnewl[i])/Smten, decimals=2)
        totdiskgmround = np.around((totdiskgml[i])/tdten, decimals=2)
        diskcoldgmround = np.around((diskcoldgml[i])/dcgten, decimals=2)
        Lgammaround = np.around((Lgammal[i])/Lgten, decimals=2)
        z12gasround = np.around((z12gasl[i]), decimals=2)
        r12crround = np.around((r12crl[i]), decimals=2)
        z12denround = np.around((z12dengasl[i]), decimals=2)
        coldgasfrac =  np.around((diskcoldgml[i])/(totdiskgml[i]), decimals=2)
        gassurdenround = np.around(gassurdenl[i], decimals=3)
        crdenround = np.around(crdenl[i], decimals=3)
        avesfrround = np.around(avesfrl[i]/sfrten, decimals=2)
        strneed = runtodol[i]+'&      $'+str(Smnewround)+'$'\
        '&      $'+str(avesfrround)+'$'\
        +'&      $'+str(totdiskgmround)+'$'\
        +'&      $'+str(diskcoldgmround)+'$'\
                +'&      $'+str(z12gasround)+'$'\
                +'&      $'+str(gassurdenround)+'$'\
        +'&      $'+str(crdenround)+'$'\
        +'&      $'+str(r12crround)+'$'\
                +'&      $'+str(Lgammaround)+'$'\
        +'\\\\'
        print strneed
                #+'&      $'+str(z12denround)+'$'\
                #+'&      $'+str(coldgasfrac)+'$'

                



if wanted=='gassurden':
    projection='y'
    for runtodo in dirneed:
        info=outdirname(runtodo, Nsnap)
        dclabel=info['dclabel']
        color=info['color']
        haveB=info['haveB']
        if haveB>0:
                lsn='dashed'
        else:
                lsn='solid'
        outdata = gassurout(runtodo,Nsnap,projection)
        rneed=outdata['rneed']
        gassurdenl=outdata['gassurdenl']
    plt.plot(rneed, gassurdenl,label=dclabel,color=color,lw=2,ls=lsn)
    plt.xlabel('r [kpc]')
    plt.ylabel(r'gas surface density $(<r) [{\rm g/cm^2}]$')
    plt.legend(loc='best')
    #plt.yscale('log')
    figname = 'CRplot/gassurden/gassurden_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
    print 'figname', figname
    plt.savefig(figname)
    plt.clf()


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
        plt.errorbar(gmed, Lmed, xerr=[[g1sigd],[g1sigu]], yerr=[[L1sigd],[L1sigu]],\
                     color=color,fmt=marker,markersize=7,label=labelneed, fillstyle=fillstyle)
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
    ftxt = open('/home/tkc004/samsonprogram/data/LTQ_line.txt', 'r')
    ftxt.readline()
    dars = ftxt.readlines()
    ftxt.close()
    LTQx=[]
    LTQy=[]
    for line in dars:
        xsd = line.split()
        LTQx=np.append(LTQx, float(xsd[0]))
        LTQy=np.append(LTQy, float(xsd[1]))

        ftxt = open('/home/tkc004/samsonprogram/data/LTQ_up.txt', 'r')
        ftxt.readline()
        dars = ftxt.readlines()
        ftxt.close()
        LTQupx=[]
        LTQupy=[]
        for line in dars:
                xsd = line.split()
                LTQupx=np.append(LTQupx, float(xsd[0]))
                LTQupy=np.append(LTQupy, float(xsd[1]))

        ftxt = open('/home/tkc004/samsonprogram/data/LTQ_down.txt', 'r')
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
    dwarfpoint = plt.errorbar([],[], fmt='o', color='g', mfc='g',markersize=6)
    lstarpoint = plt.errorbar([],[], fmt='^', color='g', mfc='g',markersize=6)
    sbpoint    = plt.errorbar([],[], fmt='D', color='g', mfc='g',markersize=6)
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
    figname='CRplot/gassur_FgFsf/gassur_FgFsf_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    print 'figname', figname
    plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
    plt.savefig(figname,bbox_inches='tight')
    plt.clf()




if wanted=='SFRLg':
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
        for data in obsdatal:
                plt.errorbar(data[0], data[3],\
                xerr=[[data[0]-data[1]],[data[2]-data[0]]],\
                yerr=[[data[3]-data[4]],[data[5]-data[3]]],\
                fmt='s', color='0.5',ls='dashed', mfc='0.9',markersize=6)

        obsbar = plt.errorbar([],[], fmt='s', color='0.5', mfc='0.9',markersize=6)
        dwarfpoint = plt.errorbar([],[], fmt='o', color='g', mfc='g',markersize=6)
        lstarpoint = plt.errorbar([],[], fmt='^', color='g', mfc='g',markersize=6)
        sbpoint    = plt.errorbar([],[], fmt='D', color='g', mfc='g',markersize=6)
        lxlist = [obsbar,dwarfpoint,lstarpoint,sbpoint]
        dclablist = ['Observations','Dwarf',r'L$\star$ Galaxy', 'Starburst']
        legend1 = plt.legend(lxlist, dclablist, loc=2,fontsize=12,ncol=3)
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
        figname='CRplot/SFRLg/SFRLg_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        print 'figname', figname
        plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=18)
        plt.savefig(figname,bbox_inches='tight')
        plt.clf()




if wanted=='crturmap':
    rcParams['figure.figsize'] = 5,4
    nbin=50
    cylr=20
    for runtodo in dirneed:
        info=outdirname(runtodo, Nsnap)
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
        maindir=info['maindir']
        haveB=info['haveB']
        egydata=energyoutput(runtodo,Nsnap)
        Begycutz=egydata['Begyl']
        cregycutz=egydata['cregyl']
        GEintcutz=egydata['therml']
        turl=egydata['turl']
        Gxpl=egydata['Gxpl']
        Gypl=egydata['Gypl']
        Gzpl=egydata['Gzpl'] 
        Gxcl=egydata['xcell']
        Gycl=egydata['ycell']
        Gzcl=egydata['zcell']
        cylz=egydata['cylz']
        Grxycz=np.sqrt(Gxpl*Gxpl+Gypl*Gypl)
        rxycell=np.sqrt(Gxcl*Gxcl+Gycl*Gycl)
        print 'np.amax(rxycell)', np.amax(rxycell)
        ran = np.linspace(0,cylr,num=nbin+1)
        volan = np.pi*(ran[1:]*ran[1:]-ran[:-1]*ran[:-1])*cylz*kpc_in_cm*kpc_in_cm*kpc_in_cm
        if haveB>0:
            histB, bin_edges = np.histogram(Grxycz, bins=nbin, range=[0,cylr], weights=np.absolute(Begycutz))
        if havecr>0:
            histcr, bin_edges = np.histogram(Grxycz, bins=nbin, range=[0,cylr], weights=np.absolute(cregycutz))
        histE, bin_edges = np.histogram(Grxycz, bins=nbin, range=[0,cylr], weights=np.absolute(GEintcutz))
        histtur, bin_edges = np.histogram(rxycell, bins=nbin, range=[0,cylr], weights=np.absolute(turl))
        if haveB>0:
            plt.plot((bin_edges[:-1]+bin_edges[1:])*0.5,np.cumsum(histB),label='B field',color='b',lw=2)
        if havecr>0:
            plt.plot((bin_edges[:-1]+bin_edges[1:])*0.5,np.cumsum(histcr),label='cosmic ray',color='y',lw=2)
        plt.plot((bin_edges[:-1]+bin_edges[1:])*0.5,np.cumsum(histtur),label='turbulence',color='g',lw=2)
        plt.plot((bin_edges[:-1]+bin_edges[1:])*0.5,np.cumsum(histE),label='thermal',color='r',lw=2)
        plt.legend(loc=4)
        plt.yscale('log')
        plt.xlabel(r'${\rm r\;[kpc]}$',fontsize=18)
        plt.ylabel(r'${\rm E(<r)\; [erg]}$',fontsize=18)
        plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
        plt.tight_layout(w_pad=0.5)
        plt.savefig('CRplot/crturmap/'+runtodo+'_rxy_cregy_cell.pdf')
        plt.clf()
        if haveB>0:
            histB=histB*erg_in_eV/volan
        if havecr>0:
            histcr=histcr*erg_in_eV/volan
        histtur=histtur*erg_in_eV/volan
        histE=histE*erg_in_eV/volan
        if haveB>0:
            plt.plot((bin_edges[:-1]+bin_edges[1:])*0.5,histB,label='B field',color='b',lw=2)
        if havecr>0:
            plt.plot((bin_edges[:-1]+bin_edges[1:])*0.5,histcr,label='cosmic ray',color='y',lw=2)
        plt.plot((bin_edges[:-1]+bin_edges[1:])*0.5,histtur,label='turbulence',color='g',lw=2)
        plt.plot((bin_edges[:-1]+bin_edges[1:])*0.5,histE,label='thermal',color='r',lw=2)
        plt.legend(loc=4)
        plt.yscale('log')
        plt.xlabel(r'${\rm r\;[kpc]}$',fontsize=18)
        plt.ylabel(r'${\rm Energy\; density\; [eV/cm^3]}$',fontsize=18)
        plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
        plt.tight_layout(w_pad=0.5)
        plt.savefig('CRplot/crturmap/'+runtodo+'_rxy_egyden_cell.pdf')
        plt.clf()

if wanted=='pretime' or wanted=='etur_mtime' or wanted=='eturtime' or wanted=='eturdentime':
    rcParams['figure.figsize'] = 5,4
    ncount=0
    ncrcount=0
    nbcount=0
    energylabel=1
    usesolarcircle=0
    usecentral=1
    if usesolarcircle==1:
        fmeat+='_usesolarcircle'
    elif usecentral==1:
        fmeat+='_usecentral'
    for runtodo in dirneed:
        info=outdirname(runtodo, Nsnap)
        rundir=info['rundir']
        runtitle=info['runtitle']
        slabel=info['slabel']
        snlabel=info['snlabel']
        dclabel=info['dclabel']
        resolabel=info['resolabel']
        the_snapdir=info['the_snapdir']
        the_prefix = info['the_prefix']
        the_suffix = info['the_suffix']
        Nsnapstring=info['Nsnapstring']
        havecr=info['havecr']
        Fcal=info['Fcal']
        iavesfr=info['iavesfr']
        timestep=info['timestep']
        maindir=info['maindir']
        haveB=info['haveB']
        cosmo=info['cosmo']
        ptitle=title
        if runtitle=='SMC':
            ptitle='Dwarf'
        elif runtitle=='SBC':
            ptitle='Starburst'
        elif runtitle=='MW':
            ptitle=r'$L\star$ Galaxy'
        turt=[]
        thermt=[]
        crt=[]
        bt=[]
        timel=[]
        Gmt=[]
        for i in range(startno,Nsnap,snapsep):
            try:
                egydata=energyoutput(runtodo,i,usesolarcircle=usesolarcircle,usecentral=usecentral)
                Begycutz=egydata['Begyl']
                cregycutz=egydata['cregyl']
                GEintcutz=egydata['therml']
                turl=egydata['turl']
                Gxpl=egydata['Gxpl']
                Gypl=egydata['Gypl']
                Gzpl=egydata['Gzpl']
                Gxcl=egydata['xcell']
                Gycl=egydata['ycell']
                Gzcl=egydata['zcell']
                gminl=egydata['gminl']
                timen=egydata['timen']
                cylrin=egydata['cylrin']
                cylr=egydata['cylr']
                cylz=egydata['cylz']
                nbin=egydata['nbin']
                Gmt = np.append(Gmt, np.sum(gminl))
                turt=np.append(turt, np.sum(turl))
                thermt=np.append(thermt,np.sum(GEintcutz))
                if havecr>0:
                    crt = np.append(crt,np.sum(cregycutz))
                if haveB>0:
                    bt = np.append(bt,np.sum(Begycutz))
                timel = np.append(timel,timen/1e6) #Myr 
            except (ZeroDivisionError,IndexError):
                continue
            turt=np.array(turt)
            thermt=np.array(thermt)
            crt=np.array(crt)
            bt=np.array(bt)
            timel=np.array(timel)
            Gmt=np.array(Gmt)       
        print 'Gmt', np.sum(Gmt)
        if energylabel==1:
            if ncount==0:
                ls='solid'
                tlb = 'turb'
                thlb = 'therm'
                clb = 'CR'
                Blb = 'B'
            else:
                tlb=thlb=clb=Blb='_nolegend_'
            if ncrcount==0: 
                clb = 'CR'
            if nbcount==0:
                Blb = 'B'
            if ncount==1:
                ls='dashed'
            if ncount==2:
                ls='dashdot'
            if ncount==3:
                ls='dotted'
        else:
            if ncount==0:
                    ls='solid'
            if ncount==1:
                    ls='dashed'
            if ncount==2:
                    ls='dashdot'
            if ncount==3:
                    ls='dotted'
        tlb = dclabel
        thlb=clb=Blb='_nolegend_'
        
        if wanted =='etur_mtime': 
            plt.plot(timel,turt/Gmt,ls=ls,lw=2,label=tlb,color='g')
            plt.plot(timel,thermt/Gmt,ls=ls,lw=2,label=thlb,color='r')
            vol = np.pi*(cylr*cylr-cylrin*cylrin)*cylz
            print 'thden', thermt*erg_in_eV/vol/kpc_in_cm/kpc_in_cm/kpc_in_cm
            if havecr>0:
                print 'clb', clb
                plt.plot(timel,crt/Gmt,ls=ls,lw=2,label=clb, color='y')
                print 'crden', crt*erg_in_eV/vol/kpc_in_cm/kpc_in_cm/kpc_in_cm
                ncrcount+=1
            if haveB>0:
                print 'Blb', Blb
                plt.plot(timel,bt/Gmt,ls=ls,lw=2,label=Blb, color='b')
                nbcount+=1
        if wanted =='eturtime':
                        plt.plot(timel,turt,ls=ls,lw=2,label=tlb,color='g')
                        plt.plot(timel,thermt,ls=ls,lw=2,label=thlb,color='r')
                        if havecr>0:
                                plt.plot(timel,crt,ls=ls,lw=2,label=clb, color='y')
                        if haveB>0:
                                plt.plot(timel,bt,lw=2,ls=ls,label=Blb, color='b')
        if wanted =='eturdentime':
            vol = np.pi*(cylr*cylr-cylrin*cylrin)*cylz*kpc_in_cm*kpc_in_cm*kpc_in_cm
                        plt.plot(timel,turt*erg_in_eV/vol,ls=ls,lw=2,label=tlb,color='g')
                        plt.plot(timel,thermt*erg_in_eV/vol,ls=ls,lw=2,label=thlb,color='r')
                        if havecr>0:
                                print 'clb', clb
                                plt.plot(timel,crt*erg_in_eV/vol,ls=ls,lw=2,label=clb, color='y')
                                ncrcount+=1
                        if haveB>0:
                                print 'Blb', Blb
                                plt.plot(timel,bt*erg_in_eV/vol,ls=ls,lw=2,label=Blb, color='b')
                if wanted =='pretime':
                        vol = np.pi*(cylr*cylr-cylrin*cylrin)*cylz*kpc_in_cm*kpc_in_cm*kpc_in_cm
                        plt.plot(timel,(GAMMA-1.0)*thermt/vol,ls=ls,lw=2,label=thlb,color='r')
                        if havecr>0:
                                print 'clb', clb
                                plt.plot(timel,(CRgamma-1.0)*crt/vol,ls=ls,lw=2,label=clb, color='y')
                                ncrcount+=1
                        if haveB>0:
                                print 'Blb', Blb
                                plt.plot(timel,bt/vol,ls=ls,lw=2,label=Blb, color='b')
        ncount+=1
    if wanted =='etur_mtime':
        plt.yscale('log')
        plt.xlabel(r'${\rm Myr}$', fontsize=18)
        plt.ylabel(r'$E/m {\rm [erg/g]}$', fontsize=18)
        if nolegend==0: 
            plt.legend(loc=2,fontsize=12,ncol=3,frameon=False)
        plt.title(ptitle,fontsize=16)
        figname='CRplot/etur_mtime/etur_mtime_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    if wanted =='eturtime':
        plt.yscale('log')
        plt.xlabel(r'${\rm Myr}$', fontsize=18)
        plt.ylabel(r'$E {\rm [erg]}$', fontsize=18)
        if nolegend==0:
            plt.legend(loc='best',fontsize=12,ncol=3,frameon=False)
        figname='CRplot/eturtime/eturtime_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted =='eturdentime':
                plt.yscale('log')
                plt.xlabel(r'${\rm Myr}$', fontsize=18)
                plt.ylabel(r'$e {\rm [eV/cm^{3}]}$', fontsize=18)
                if nolegend==0:
                        plt.legend(loc='best',fontsize=12,ncol=3,frameon=False)
                figname='CRplot/eturdentime/eturdentime_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted =='pretime':
                plt.yscale('log')
                plt.xlabel(r'${\rm Myr}$', fontsize=18)
                plt.ylabel(r'$e {\rm [dyne/cm^{2}]}$', fontsize=18)
                if nolegend==0:
                        plt.legend(loc='best',fontsize=12,ncol=3,frameon=False)
                figname='CRplot/pretime/pretime_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    print figname
    plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
    plt.tight_layout()
    plt.savefig(figname)
    plt.clf()   
        



if wanted=='crez' or wanted=='dpcrz' or wanted=='dpcrz_rhog' or wanted=='rhog'\
     or wanted=='dpcrzmrhog' or wanted=='dpcrz_dpth' or wanted=='pcr_pth'\
     or wanted=='vcr_vth' or  wanted=='dGmden':
        def func(x, a, b, c):
                return a+b*x+c*x*x
        def d10func(r, a, b, c):
                return np.power(10.,func(np.log10(r), a, b, c))*(b+2*c*np.log(r)/np.log(10.))/r
        rcParams['figure.figsize'] = 5,4
        resoneed=0
    rotface=1
    needlogz=0
    thermalneed=1
    needfit=0
    labcount=0
    linestylabel=1
    newlabelneed=1
    twolegend=1
    ptotneed = 1
    if twolegend==1:
        import matplotlib.lines as mlines
        lxlist=[]
        lslist=[]
        dclablist=[]
        for icount in range(int(len(dirneed))):
            lineobj = mlines.Line2D([], [])
            lxlist.append(lineobj)
        for runtodo in dirneed:
        print 'runtodo', runtodo
        maxlength = 50.
        minlength=5.
        nogrid=8
        withinr=5.0
        diskr=10
        diskh = 1.0
        if needlogz==1:
            zl=np.logspace(-1,np.log10(maxlength),num=nogrid)
        else:
            zl=np.linspace(minlength,maxlength,num=nogrid)
            dz = maxlength/nogrid
            crden=zl*0.0
        thden=zl*0.0
                Gmden=zl*0.0
                mtotl=zl*0.0
        mdtot = 0.0
                numoftimes=0
                info=outdirname(runtodo, Nsnap)
                havecr=info['havecr']
        color=info['color']
        runtitle=info['runtitle']
        ptitle=title
        if runtitle=='SMC':
            ptitle='Dwarf'
        elif runtitle=='SBC':
            ptitle='Starburst'
        elif runtitle=='MW':
            ptitle=r'$L\star$ Galaxy'
        if wanted=='crez':
            if havecr==0:
                continue
                for i in range(startno,Nsnap,snapsep):
                        info=outdirname(runtodo, i)
                        M1speed=info['M1speed']
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        the_prefix = info['the_prefix']
                        the_suffix = info['the_suffix']
                        Nsnapstring=info['Nsnapstring']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        haveB=info['haveB']
                        newlabel=info['newlabel']
            cosmo=info['cosmo']
            halostr=info['halostr']
            firever=info['firever']
            maindir=info['maindir']
            usepep=info['usepep']
            snumadd=info['snumadd']
            h0=cosmo
            if cosmo==1:
                datasup=0
            else:
                datasup=1
                        Gextra = readsnapwcen(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix,\
                         havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
                         datasup=datasup,runtodo=runtodo,rundir=rundir,halostr=halostr)
                        Gx = Gextra['x']; Gy = Gextra['y']; Gz = Gextra['z'];
                        Gvx = Gextra['vx']; Gvy = Gextra['vy']; Gvz = Gextra['vz'];
                        Grho = Gextra['rho']; Gu = Gextra['u']; Gm = Gextra['m']
                        if havecr>0:
                                cregy = Gextra['cregy']*1e10*Msun_in_g*km_in_cm*km_in_cm #cosmic ray energy in (originally) 1e10Msun km^2/sec^2

                        if haveB==0:
                                lsn='solid'
                        else:
                                lsn='dashed'
                        #if M1speed>999:
                        #       lsn='dotted'
                        #if M1speed>1999:
                        #       lsn='dashdot'
            if wanted=='dpcrz_rhog' or wanted=='rhog' or wanted=='dpcrzmrhog' or wanted=='vcr_vth' or wanted=='dGmden':
                mdisktot=0.0
                print 'diskr, diskh', diskr, diskh
                for ipd in [0,2,4]:
                    try:
                        Pextra = readsnapwcen(the_snapdir, Nsnapstring, ipd, snapshot_name=the_prefix, extension=the_suffix,\
                         havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
                         datasup=datasup,runtodo=runtodo,rundir=rundir,halostr=halostr)
                        Px = Pextra['x']; Py = Pextra['y']; Pz = Pextra['z'];   
                        Pm = Pextra['m']
                        Prxy = np.sqrt(Px*Px+Py*Py)
                        cutxyz = (Prxy<diskr)*(np.absolute(Pz)<diskh) #kpc
                        mdisktot += np.sum(Pm[cutxyz]*1e10*Msun_in_g) #in g
                        print 'mdisktot', mdisktot
                    except KeyError:
                        continue
                mdtot += mdisktot
                        for iz in range(len(zl)-1):
                                cutxy = (Gx*Gx+Gy*Gy < withinr*withinr)
                                cutzu = Gz<zl[iz+1]
                                cutzd = Gz>zl[iz]
                                cutz = cutzu*cutzd
                                cut = cutxy*cutz
                #print 'withinr', withinr
                #print 'Gm[cut]', Gm[cut]
                #print 'zl', zl
                if havecr>0:
                    crecut = cregy[cut]
                    crden[iz] += np.sum(crecut)/(zl[iz+1]-zl[iz])/kpc_in_cm/(np.pi*withinr*withinr*kpc_in_cm*kpc_in_cm)
                thermalcut = Gm[cut]*Gu[cut]*1e10*Msun_in_g*km_in_cm*km_in_cm
                thden[iz] += np.sum(thermalcut)/(zl[iz+1]-zl[iz])/kpc_in_cm/(np.pi*withinr*withinr*kpc_in_cm*kpc_in_cm)
                if wanted=='dpcrz_rhog' or wanted=='rhog' or wanted=='dpcrzmrhog' or wanted=='vcr_vth' or wanted=='dGmden':
                    Gmcut=Gm[cut]*1e10*Msun_in_g
                    Gmden[iz] += np.sum(Gmcut)/(zl[iz+1]-zl[iz])/kpc_in_cm/(np.pi*withinr*withinr*kpc_in_cm*kpc_in_cm)
                    mintot = 0
                    for ipa in [0,1,2,3,4]:
                        try:
                            Pextra = readsnapwcen(the_snapdir, Nsnapstring, ipd, snapshot_name=the_prefix, extension=the_suffix,\
                             havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
                             datasup=datasup,runtodo=runtodo,rundir=rundir,halostr=halostr)
                            Px = Pextra['x']; Py = Pextra['y']; Pz = Pextra['z'];
                            Pm = Pextra['m']
                            Pr = np.sqrt(Px*Px+Py*Py+Pz*Pz)
                            cutpr = Pr<zl[iz]
                            if ipa==0 or ipa==2 or ipa==4:
                                Prxy = np.sqrt(Px*Px+Py*Py)
                                if zl[iz]<diskr:
                                    drc = zl[iz]
                                else:
                                    drc = diskr
                                cutdisk = (Prxy<drc)*(np.absolute(Pz)<diskh)
                                mdiskcut = np.sum(Pm[cutdisk]*1e10*Msun_in_g)
                            Pmtot = np.sum(Pm[cutpr]*1e10*Msun_in_g)-mdiskcut #in g
                            print 'Pmtot', Pmtot
                        except KeyError:
                            continue
                        mintot += Pmtot
                        mtotl[iz] += mintot
                        numoftimes+=1
        if havecr>0:
            crden = crden/numoftimes
            thden = thden/numoftimes
            Gmden = Gmden/numoftimes
            mtotl = mtotl/numoftimes
            mdtot = mdtot/numoftimes
        if havecr>0:
            print 'crden', crden
            print 'thden', thden
            print 'mtotl', mtotl
            print 'mdtot', mdtot
        labelneed=dclabel
        if newlabelneed==1:
            labelneed="\n".join(wrap(newlabel,17))
                if resoneed==1:
                        if resolabel=='llr':
                                labelneed='low res'
                                lsn = 'solid'
                        if resolabel=='lr':
                                labelneed='std res'
                                lsn = 'dashed'
                        if resolabel=='mr':
                                labelneed='high res'
                                lsn = 'dashdot'
        zlm = (zl[1:]+zl[:-1])/2.
        if needfit==1:
            x0 = [0.1,-1.0,0.0]
            thden1 = thden[:-1]
            xdata=np.log10(zlm[~np.isnan(thden1)])
            ydata=np.log10(thden1[~np.isnan(thden1)])
            outfit=optimize.curve_fit(func, xdata, ydata, x0)
            afit=outfit[0][0]
            bfit=outfit[0][1]
            cfit=outfit[0][2]
            dthcr = -d10func(zlm,afit,bfit,cfit)*(GAMMA-1.0)/kpc_in_cm
        else:
            dthcr = -(thden[:-1]-thden[1:])/(zl[:-1]-zl[1:])*(GAMMA-1.0)/kpc_in_cm
        if havecr>0:
            if needfit==1:
                x0 = [0.1,-1.0,0.0]
                crden1 = crden[:-1]
                xdata=np.log10(zlm)
                ydata=np.log10(crden1)
                xdata = xdata[~np.isinf(ydata)]
                ydata = ydata[~np.isinf(ydata)]
                print 'xdata', xdata
                print 'ydata', ydata
                outfit=optimize.curve_fit(func, xdata, ydata, x0)
                afit=outfit[0][0]
                bfit=outfit[0][1]
                cfit=outfit[0][2]
                dpcr = -d10func(zlm,afit,bfit,cfit)*(CRgamma-1.0)/kpc_in_cm
            else:
                dpcr = -(crden[:-1]-crden[1:])/(zl[:-1]-zl[1:])*(CRgamma-1.0)/kpc_in_cm
                if wanted=='crez':
                        plt.plot(zlm, crden[:-1]*erg_in_eV, label=labelneed,lw=2,ls=lsn)
                if wanted=='dpcrz':
            if havecr>0:
                plt.plot(zlm, dpcr, label=labelneed,lw=2,ls=lsn,color=color)
                        plt.plot(zlm, dthcr,lw=1,color=color)
        if wanted=='dpcrz_dpth':
            if havecr>0:
                plt.plot(zlm, (crden[:-1]-crden[1:])/(thden[:-1]-thden[1:]),ls=lsn,color=color,marker='o')
                plt.plot(zlm, dpcr/dthcr, label=labelneed,lw=2,ls=lsn,color=color)
                if wanted=='pcr_pth':
                        if havecr>0:
                                plt.plot(zlm, (CRgamma-1.0)*crden[:-1]/(GAMMA-1.0)/thden[:-1], label=labelneed,lw=2,ls=lsn,color=color)
        if wanted=='vcr_vth':
            dden = -(Gmden[:-1]-Gmden[1:])/(zl[:-1]-zl[1:])/kpc_in_cm
            gdisk = 2.0*NewtonG_cgs*mdtot/diskr/diskr*(1.0-zl/np.sqrt(diskr*diskr+zl*zl))/kpc_in_cm/kpc_in_cm
            gsphere = mtotl*NewtonG_cgs/(zl*kpc_in_cm)/(zl*kpc_in_cm)
            vg  = np.sqrt((gdisk+gsphere)*zl*kpc_in_cm)/km_in_cm
            vth = np.sqrt(dthcr/dden)/km_in_cm
            if havecr>0:
                vcr = np.sqrt(dpcr/dden)/km_in_cm
                vcom = np.sqrt((dthcr+dpcr)/dden)/km_in_cm
            plt.plot(zl[:-1], vg[:-1], label=labelneed,lw=2,ls=lsn,color=color)
            if havecr>0:
                plt.plot(zlm, vcr,ls=lsn,marker='o',color=color,markeredgewidth=0.0,lw=1)
            plt.plot(zlm, vth,ls=lsn,marker='^',color=color,markeredgewidth=0.0,lw=1)   
                if wanted=='dGmden':
            dden = -(Gmden[:-1]-Gmden[1:])/(zl[:-1]-zl[1:])/kpc_in_cm
            plt.plot(zlm,dden,ls='none',color=color,marker='o')
                        plt.plot(zlm,dGmden, label=labelneed,lw=2,ls=lsn,color=color)
                if wanted=='dpcrz_rhog':
            print 'runtodo', runtodo
            gdisk = 2.0*NewtonG_cgs*mdtot/diskr/diskr*(1.0-zl/np.sqrt(diskr*diskr+zl*zl))/kpc_in_cm/kpc_in_cm
            gsphere = mtotl*NewtonG_cgs/(zl*kpc_in_cm)/(zl*kpc_in_cm)
                        rhog = Gmden*(gdisk+gsphere)
            print 'rhog', rhog
            if havecr>0:
                plt.plot(zlm, dpcr/rhog[:-1], label=labelneed,lw=2,ls=lsn,color=color)
            if thermalneed>0:
                plt.plot(zlm, dthcr/rhog[:-1],lw=1,ls=lsn,color=color)
                if wanted=='dpcrzmrhog':
                        if havecr>0:
                                pcrden = (CRgamma-1.0)*crden
                        pthden = (GAMMA-1.0)*thden
                        print 'Gmden, mtotl', Gmden, mtotl
                        gdisk = 2.0*NewtonG_cgs*mdtot/diskr/diskr*(1.0-zl/np.sqrt(diskr*diskr+zl*zl))/kpc_in_cm/kpc_in_cm
                        gsphere = mtotl*NewtonG_cgs/(zl*kpc_in_cm)/(zl*kpc_in_cm)
                        gacc= gdisk+gsphere
                        if havecr>0:
                                print 'pcrden[:-1]', pcrden[:-1]
                                dpcr_rho = (pcrden[1:]-pcrden[:-1])/(zl[1:]-zl[:-1])/kpc_in_cm/Gmden[:-1]
                        dthcr_rho = (pthden[1:]-pthden[:-1])/(zl[1:]-zl[:-1])/kpc_in_cm/Gmden[:-1]
                        if havecr>0:
                                plt.plot(zlm[:-1], (dpcr_rho[:-1]+gacc[:-2])/km_in_cm*yr_in_sec*1e6, label=labelneed,lw=2,ls=lsn,color=color)
                        plt.plot(zlm[:-1], (dthcr_rho[:-1]+gacc[:-2])/km_in_cm*yr_in_sec*1e6,lw=1,ls=lsn,color=color)
                if wanted=='rhog':
                        gdisk = 2.0*NewtonG_cgs*mdtot/diskr/diskr*(1.0-zl/np.sqrt(diskr*diskr+zl*zl))/kpc_in_cm/kpc_in_cm
                        gsphere = mtotl*NewtonG_cgs/(zl*kpc_in_cm)/(zl*kpc_in_cm)
                        rhog = Gmden*(gdisk+gsphere)
            if thermalneed>0:
                if linestylabel==1 and labcount==1:
                    plt.plot(zlm[:-1], dthcr[:-1],ms=4,ls='dashed',lw=1,color=color,marker='^',label=r'${\rm d} p_{\rm th}/{\rm d} z$')
                else:
                    plt.plot(zlm[:-1], dthcr[:-1],ms=4,ls='dashed',lw=1,color=color,marker='^')
                        if havecr>0:
                                if linestylabel==1 and labcount==1:
                                        plt.plot(zlm[:-1], dpcr[:-1],ms=4,ls='solid',lw=1,color=color,marker='o',label=r'${\rm d} p_{\rm CR}/{\rm d} z$')
                    if ptotneed==1:
                        plt.plot(zlm[:-1], dpcr[:-1]+dthcr[:-1],ms=4,ls='solid',lw=1,color=color,marker='d',label=r'${\rm d} p_{\rm tot}/{\rm d} z$')
                else:
                    plt.plot(zlm[:-1], dpcr[:-1],ms=4,ls='solid',lw=1,color=color,marker='o')
                    if ptotneed==1:
                        plt.plot(zlm[:-1], dpcr[:-1]+dthcr[:-1],ms=4,ls='solid',lw=1,color=color,marker='d')
            if linestylabel==1:
                if labcount==0:
                    plt.plot(zl[:-1], rhog[:-1],color=color,label=r'$\rho g$',ls=lsn,lw=2)
                else:
                    plt.plot(zl[:-1], rhog[:-1],color=color,ls=lsn,lw=2)
            else:
                plt.plot(zl[:-1], rhog[:-1],color=color,label=labelneed,ls=lsn,lw=2)
                        if twolegend==1:
                                lxlist[labcount] = mlines.Line2D([], [], color=color, ls='solid',label=labelneed)
                                dclablist = np.append(dclablist,labelneed)
            labcount+=1
        if wanted=='crez':
                plt.yscale('log')
                plt.xlabel('z [kpc]')
                plt.ylabel(r'$e_{\rm cr} [{\rm eV/cm^3}]$')
                plt.legend(loc='best')
                plt.title('cylinder centered at disk with a radius ' + str(withinr) + ' kpc')
                filename='CRplot/crez/crez_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted=='dpcrz':
                plt.yscale('log')
                plt.xlabel(r'${\rm z\; [kpc]}$',fontsize=18)
                plt.ylabel(r'${\rm d}P/{\rm d}z [{\rm erg/cm^4}]$',fontsize=18)
                plt.legend(loc='best')
                #plt.title('cylinder centered at disk with a radius ' + str(withinr) + ' kpc')
                filename='CRplot/dpcrz/dpcrz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted=='dpcrz_dpth':
                plt.yscale('log')
                #plt.xscale('log')
                plt.xlabel(r'${\rm z\; [kpc]}$',fontsize=18)
                plt.ylabel(r'${\rm d}P_{\rm cr}/{\rm d}P_{\rm th}$',fontsize=18)
                #plt.legend(loc='best')
                #plt.title('cylinder centered at disk with a radius ' + str(withinr) + ' kpc')
                filename='CRplot/dpcrz_dpth/dpcrz_dpth_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted=='pcr_pth':
                plt.yscale('log')
                #plt.xscale('log')
                plt.xlabel(r'${\rm z\; [kpc]}$',fontsize=18)
                plt.ylabel(r'$P_{\rm cr}/P_{\rm th}$',fontsize=18)
                #plt.legend(loc='best')
                #plt.title('cylinder centered at disk with a radius ' + str(withinr) + ' kpc')
                filename='CRplot/pcr_pth/pcr_pth_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
                plt.title(ptitle)
        if wanted=='vcr_vth':
                #plt.yscale('log')
                #plt.xscale('log')
                plt.xlabel(r'${\rm z\; [kpc]}$',fontsize=18)
                plt.ylabel(r'$v\;{\rm [km/s]}$',fontsize=18)
                #plt.legend(loc='best')
                #plt.title('cylinder centered at disk with a radius ' + str(withinr) + ' kpc')
                filename='CRplot/vcr_vth/vcr_vth_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
                plt.title(ptitle,fontsize=16)
        if wanted=='dGmden':
                #plt.yscale('log')
                #plt.xscale('log')
                plt.xlabel(r'${\rm z\; [kpc]}$',fontsize=18)
                plt.ylabel(r'$\;{\rm d} \rho/{\rm d} z{\rm [g/cm^4]}$',fontsize=18)
                #plt.legend(loc='best')
                #plt.title('cylinder centered at disk with a radius ' + str(withinr) + ' kpc')
                filename='CRplot/dGmden/dGmden_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
                plt.title(ptitle,fontsize=16)
        if wanted=='dpcrz_rhog':
                #plt.yscale('log')
                #plt.axhspan(-2, 1, alpha=0.5,color='0.5')
                plt.ylim(ymin=-0.5)
                plt.xlabel(r'${\rm z\; [kpc]}$',fontsize=18)
                plt.ylabel(r'${\rm d}P/{\rm d}z/(\rho g)$',fontsize=18)
                #plt.legend(loc='best')
                #plt.title('cylinder centered at disk with a radius ' + str(withinr) + ' kpc')
                filename='CRplot/dpcrz_rhog/dpcrz_rhog_hd_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted=='dpcrzmrhog':
                plt.xlabel(r'${\rm z\; [kpc]}$',fontsize=18)
                plt.ylabel(r'$g+{\rm d}_zP/\rho\;  [{\rm km/s/Myr}]$',fontsize=18)
                #plt.legend(loc='best')
                #plt.title('cylinder centered at disk with a radius ' + str(withinr) + ' kpc')
                filename='CRplot/dpcrzmrhog/dpcrzmrhog_hd_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'

        if wanted=='rhog':
                if twolegend==1:
                        #print 'lxlist', lxlist
                        #print 'dclablist', dclablist
                        legend1 = plt.legend(lxlist, dclablist, loc=3,fontsize=8,ncol=2)
                        plt.gca().add_artist(legend1)
                        plt.legend(loc=3, fontsize=8)
                plt.yscale('log')
                #plt.xscale('log')
                plt.xlabel(r'${\rm z\; [kpc]}$',fontsize=18)
                plt.ylabel(r'$\rho g [{\rm erg/cm^4}]$',fontsize=18)
        if runtitle=='SMC' or linestylabel==1:
            plt.legend(loc='best',ncol=3,fontsize=10)
            plt.title(ptitle,fontsize=16)
            #plt.title('cylinder centered at disk with a radius ' + str(withinr) + ' kpc')
            filename='CRplot/rhog/rhog_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf' 
    print 'filename', filename
        plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
        plt.tight_layout()
        plt.savefig(filename)
        plt.clf()







if wanted=='cratmult':
    rcParams['figure.figsize'] = 5,4
    strlabelneed=0  
    newlabelneed=0
    runlabelneed=0
    legendneed=1
    checkcr=1
    twolegend=0
    rateneed=0
    withincr=0
    snon=1
    coolon=1
    adon=1
    numon=1
    stroff=0
    avetime=0 #number of time bins
    outputcrout=0
    ratiocrout_sn=0
        ratiocrout_sou=0
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
                    if withincr==1 or outputcrout==1:
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
                                time2list=readtimelist['timelist']
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
            crecumnum = crecumnum-crecump
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
            if withincr==1 or outputcrout==1:
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
            if withincr==1 or outputcrout==1:
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
                plt.plot(timel, xcrg, lw=2,ls=ls, label=labelneed,color='r')
            if coolon==1:
                plt.plot(timel, xcrl, lw=2,ls=ls, color='y')
            if adon==1:
                plt.plot(timel, xcra, lw=2, ls=ls, color='b')
            if checkcr==1:
                plt.plot(timel, xcr, lw=2, ls=ls, color='k')
                        if stron ==1 and stroff==0:
                                plt.plot(timel, xcrp, lw=2, ls=ls,color='c')
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
                plt.plot(timel, xcrg, lw=2,ls=ls, label=snelabel,color='r')
            if coolon==1:
                plt.plot(timel, xcrl, lw=2,ls=ls, label=losslabel,color='y')
            if adon==1:
                plt.plot(timel, xcra, lw=2, ls=ls, label=crecumalabel,color='b')
            if checkcr==1:
                plt.plot(timel, xcr, lw=2, ls=ls,color='k')
            if stron ==1 and stroff==0:
                plt.plot(timel, xcrp, lw=2, ls=ls, label=strlabel,color='c')
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
                #print 'lxlist', lxlist
                #print 'dclablist', dclablist
                legend1 = plt.legend(lxlist, dclablist, loc=2,fontsize=8,ncol=2)
                plt.gca().add_artist(legend1)
                plt.legend(loc=3, fontsize=8)
            else:
                plt.legend(loc='best', ncol=2, fontsize=10)
    if rateneed==1:
                plt.ylabel(r'$ \dot{E}_{\rm mr} [ {\rm erg/s}]$',fontsize=18)
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
            plt.ylim(ymax=2.3e39,ymin=-2.3e39)  
        if runtitle=='MW':
            plt.ylim(ymax=2.5e41,ymin=-2.5e41)
        if runtitle=='SBC':
            plt.ylim(ymax=4.4e41,ymin=-4.4e41)
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
        filename='CRplot/cratmult/'+fmeat+'_rate'+fnsuffix+'.pdf'
    else:
        filename='CRplot/cratmult/'+fmeat+'_time'+fnsuffix+'.pdf'
    print 'filename', filename
    plt.savefig(filename,bbox_inches='tight')
    plt.clf()

