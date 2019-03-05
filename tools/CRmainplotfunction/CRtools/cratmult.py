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
#dirneed=['bwmwmr','bwmwmrdc28', 'bwmwmrdc29']
#dirneed=['bwmwmr','bwmwmrdc27','bwmwmrdc28', 'bwmwmrdc29']
#dirneed=['bwmwmr']
dirneed=['bwmwmrdc28str']
#dirneed=['bwmwmrstreve']
#dirneed=['bwsmclrdc28mhdref']
#dirneed=['bwmwlrdc28']
#dirneed=['bwmwlrdc28ref']
#dirneed=['bwsmclrdc0','bwsmclrdc27','bwsmclrdc28','bwsmclrdc29']
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
startno=0
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
#wanted='rhog'
#wanted='pcr_pth'
wanted='cratmult'
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
ptitle='LSG'
nolegend=0


print 'wanted', wanted
print 'fmeat', fmeat
print 'runtodo', dirneed


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
    numon=0
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
    filename=plotloc+'CRplot/cratmult/'+fmeat+'_rate'+fnsuffix+'.pdf'
else:
    filename=plotloc+'CRplot/cratmult/'+fmeat+'_time'+fnsuffix+'.pdf'
print 'filename', filename
plt.savefig(filename,bbox_inches='tight')
plt.clf()

