import matplotlib as mpl
from readsnap_cr import *
from Sasha_functions import *
from samson_functions import *
from samson_const import *
import sys
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pylab import *
#rcParams['figure.figsize'] = 5, 5
rcParams['figure.figsize'] = 6, 5
rcParams['font.size']=18
rcParams['font.family']='serif'
#rcParams['text.usetex']=True
rcParams.update({'figure.autolayout': True})
import matplotlib.patches as patches
rcParams['axes.linewidth'] = 2
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True
rcParams['axes.unicode_minus']=False
colortable = [ 'b', 'g', 'r']



ic=0
#filename='box_CR_trapezoidal'
#the_prefix = '/home/tkc004/oasis/testICs/'+filename
#the_suffix = '.hdf5'
#the_snapdir='/home/tkc004/oasis/testICs/'
#the_snapdir='/home/tkc004/oasis/testICs/box_CR_trapezoidal'
#filename='box_CR_Gauss_eps_0_5_1D'
#the_prefix = the_snapdir+filename
#the_snapdir ='/home/tkc004/oasis/CRdifftest/dc27g1D_eps_0_5/output'
#filename='dc28g1D_eps_0_5'
#filename='M1dc28g1D_eps_0_5_c1000_testlimiter'
#filename='M1dc1_27g1D_eps_0_5_c100'
#filename='M1dc27g1D_eps_0_5'
#filename='M1dc28g1D_eps_0_5_c500_testlimiter'
#filename='M1dc28g1D_eps_0_5_c2000_testlimiter'
#filename='M1dc29g1D_eps_0_5_c1000_testlimiter'
#filename='M1dc28g1D_eps_0_5_c4000_testlimiter'
#filename='M1strg1D_testflag'
#filename='M1strg1D_eps_0_5_c1000'
#filename='M1strg1D_va'
#filename='M1dc27g1D_eps_0_5_c2000_testlimiter3'
#filename='M1dc27g1D_eps_0_5_c2000_testlimiter4'
#filename='M1dc28g1D_eps_0_5_c500_corHLL'
#filename='M1strg1D_va_c4000_noHLL'
#filename='M1strg1D_va_c4000_corHLL'
#filename='M1strg1D_fva_c4000_corHLL'
#filename='M1strg1D_va10_c2000_useflux_Jianhll'
#filename='M1strg1D_va10_c500_useflux_Jianhll'
#filename='M1strg1D_va10_c2000_useflux_oldhll'
#filename='M1strg1D_va10_c500_useflux_oldhll'
#filename='M1strt1D_va_c500_corHLL'
#filename='M1strt1D_va_c500_noHLL'
#filename='M1strt1D_va_c500_test0vstr'
#filename='M1strg1D_va_c500_noHLL'
#filename='M1dc28g1D_eps_0_5_c100_noHLL'
#filename='M1strtra1D_va_c500_noHLL'
#filename='M1dc27g1D_eps_0_5_c100_testfluxlimit'
#filename='M1dc27g1D_eps_0_5_c100_noHLL'
#filename='M1dc27g1D_eps_0_5_c500_noHLL'
#filename='M1dc27g1D_eps_0_5_c500_Jianhll'
#filename='M1dc27g1D_eps_0_5_c500_oldhll'
#filename='M1dc27g1D_eps_0_5_c500_corHLL'
#filename='M1strg1D_va10_c2000_useflux_Jianhll'
#filename='M1strg1D_va10_c2000_useflux_Jianhll_referee'
#filename='M1dc25g1D_eps_0_5_c500_Jianhll'
#filename='M1dc25g1D_eps_0_5_c500_oldhll'
#filename='M1dc25g1D_eps_0_5_c500_corHLL'
#filename='M1dc28g1D_eps_0_5_c1000_noHLL'
#filename='M1strt1D_va_c1000_noHLL'
#filename='M1dc27g1D_eps_0_5_c1000_testlimiter5'
#filename='M1dc0g1D_eps_0_5_c1000'
#filename='M1dc1_27g1D_eps_0_5_c100'
#filename='M1dc27g1D_eps_0_5_c2000_testlimiter2'
#filename='M1strg1D_va_c4000'
#filename='M1dc27g1D_eps_0_5_c100'
#filename='M1dc27g1D_eps_0_5_c300_testlimiter'
#filename='M1dc27g1D_eps_0_5_c1000_testlimiter'
##filename='M1dc27g1D_eps_0_5_c2000_testlimiter'


filename='M1strt1D_va_c1000_useflux_2048'
#filename='M1strt1D_va_c500_useflux_1024'
#filename='M1strt1D_va_c500_useflux'
#the_snapdir='/home/tkc004/oasis/CRdifftest/M1strt1D_va_c500_useflux/box_CR_triflux'
the_snapdir ='/home/tkc004/oasis/CRdifftest/'+filename+'/output'

#the_snapdir = '/home/tkc004/oasis/bw/CRdifftest/M1dc28g3D_c1000_usefluxy2000_Jianhll_referee/output'

the_prefix = 'snapshot'
the_suffix = '.hdf5'
#the_prefix = the_snapdir
#Nsnapstring = '080'
Nsnapstring = '010'
the_ptype = 0
havecr=1
haveB=1
needana=2
neednum=1
readflux=1
kappa = 0.
#kappa =  10. #diffusion coefficient #in code unit #kpc^2/Gyr
#kappa = 0.1
#kappa =  3.3/0.78 #diffusion coefficient #in code unit #kpc^2/Gyr
sep=20


def streamingsol(va0,timecr,Ec0,xglin):
    xm = np.square(CRgamma*va0*timecr+1e-15)+(CRgamma*va0*timecr)*2.*Ec0
    #print 'xm', xm
    xm = np.sqrt(xm)
    cre_g = xglin*0.0
    cutin=np.logical_and(xglin<xm,xglin>-xm)
    cutout=np.logical_or(xglin>xm,xglin<-xm)
    cre_g[cutin] = Ec0-xm+CRgamma*va0*timecr
    cre_g[cutout] = Ec0-np.absolute(xglin[cutout])+CRgamma*va0*timecr
    cre_g[cre_g<0.] = 0.
    xm=0.
    cutin=np.logical_and(xglin<xm,xglin>-xm)
    cutout=np.logical_or(xglin>xm,xglin<-xm)
    cre_g_0 = xglin*0.0
    cre_g_0[cutin] = Ec0-xm
    cre_g_0[cutout] = Ec0-np.absolute(xglin[cutout])
    cre_g_0[cre_g<0.] = 0.
    return cre_g, cre_g_0


if neednum>0:
        G = readsnapcr(the_snapdir, Nsnapstring, the_ptype, snapshot_name=the_prefix,\
        extension=the_suffix,havecr=havecr,ic=ic,readflux=readflux)
        hsml = G['h']
        pos=G['p']
        vel=G['v']
        if readflux==1:
                cregy=G['cregy']
                crflux=G['crflux']
                crfx = crflux[:,0]
                crfy = crflux[:,1]
                print 'crfy', crfy, np.min(crfy), np.max(crfy)
        xg = pos[:,0]
        #rho=G['rho']
        m =np.array(G['m'])
        u = np.array(G['u'])
        rho = np.array(G['rho'])
        cregy = np.array(G['cregy'])

# for needana==1, i.e. diffusion solution
DIMS=1; # number of dimensions 
Ngas=len(m); # 1D particle number (so total particle number is N_1D^DIMS)
Lbox = 10.0 # box side length #kpc
rho_desired = 0.00245 # box average initial gas density #1e10Msun/kpc^3 #6.805e-22g/cm^3 #407 cm^-3
P_desired = 1.0 # initial gas pressure
gamma_eos = 5./3. # polytropic index of ideal equation of state the run will assume
CRE = 1. # total CR energy #1e10Msun*(km/s)^2 #2e53 erg
eps = 0.5 # small offset in initial Gaussian packet # in code unit #kpc

# for needana==2, i.e. streaming and initial condition = triangle
Ec0=1.0 #initial peak CR energy density
#va0 =22.0/kpc_in_cm*km_in_cm*Gyr_in_sec #alfven velocity #in kpc/Gyr
va0=1.0
#va0=10.0
#va0 =1.8*3/kpc_in_cm*km_in_cm*Gyr_in_sec #alfven velocity #in kpc/Gyr
xm0 =1.0 #box frame considered (count from center) #in kpc

Ec0=2.0 #initial peak CR energy density
#va0 =1.0 #alfven velocity #in kpc/Gyr
#xm0 =1.0 #box frame considered (count from center) #in kpc


timecr = 0.00098*float(Nsnapstring) #by default #in Gyr
Nnum=500
xglin=np.linspace(np.amin(xg),np.amax(xg),num=Nnum)
#xglin=np.linspace(-1.,1.,num=Nnum)
#B=1e-4 #gauss
#nism=rho*1e10*Msun_in_g/protonmass_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm
#if haveB==1:
#       GB = G['B']; Bx=GB[:,0]; By=GB[:,1]; Bz=GB[:,2];
#       B= np.sqrt(Bx*Bx+By*By+Bz*Bz)
#else:
#       B=1e-4
#print 'np.amax(B),np.amin(B)', np.amax(B), np.amin(B)
#print 'B', B
#va=B/np.sqrt(np.pi*4.*nism*protonmass_in_g)/km_in_cm
if needana==1:
        cre_g=diffusionsol(CRE,kappa,timecr,xglin,offset=eps,dim=DIMS)
elif needana==2:
        cre_g, cre_g_0 = streamingsol(va0,timecr,Ec0,xglin)
        #print 'np.absolute(xglin[cutout])', np.absolute(xglin[cutout])
if neednum>0:   
        nism=rho*1e10*Msun_in_g/protonmass_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm
        if haveB==1:
               GB = G['B']; Bx=GB[:,0]; By=GB[:,1]; Bz=GB[:,2];
               B= np.sqrt(Bx*Bx+By*By+Bz*Bz)
        else:
               B=1e-4
        print 'np.amax(B),np.amin(B)', np.amax(B), np.amin(B)
        print 'B', B
        va=B/np.sqrt(np.pi*4.*nism*protonmass_in_g)/km_in_cm
        creden = cregy*Ngas/Lbox
        pre = u*gamma_eos/(gamma_eos-1.)*rho
        vx = vel[:,0]
        print 'h', hsml
        print 'va', va
        print 'vx', vx
        print 'xg', xg
        print 'cregy', cregy
        print 'm', m
        print 'rho',rho
        print 'pre/rho*gamma', u*gamma_eos/(gamma_eos-1.)
        print 'pre', pre
        print 'pcr', cregy*(1./3.)
        print 'm/rho', m/rho, np.amax(m/rho), np.amin(m/rho)
        print 'np.amax(cregy)', np.amax(cregy)
        print 'np.sum(cregy)', np.sum(cregy)
        #print 'pre', pre
        #print 'rholeft', rho[xg<250]
        #print 'rhoright', rho[xg>250]
        argx = np.argsort(xg)
        xg = xg[argx]
        creden = creden[argx]
        xgsep=xg[::sep]
        credensep=creden[::sep]
        crehis,edges=np.histogram(xg, bins=xglin,weights=creden*Lbox/Ngas)
        credenhis = crehis/(edges[1]-edges[0])
        xghis =(edges[:-1]+edges[1:])/2.
if needana>0:
        plt.plot(xglin, cre_g,lw=1, label='Analytic')
        #plt.plot(xglin, cre_g_0,lw=1,color='r')
#plt.plot(xg, crfx,ls='none',marker='o',ms=1, label='Sim')
if neednum>0:
        plt.plot(xgsep, credensep,ls='none',marker='o',ms=5, label='Sim')
#plt.plot(xg,va,ls='none',marker='o',ms=2)
#plt.plot(xgsep, credensep,lw=1,ls='dashed', label='Sim')
plt.xlim([-1,1])
plt.ylim(ymin=1.0)
plt.xlabel('x (kpc)', fontsize=16)
plt.ylabel(r'$e_{\rm cr}\;{\rm [code\;unit]}$', fontsize=16)
plt.legend()
plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
figname='/home/tkc004/samsonprogram/CRplot/CRIC/'+filename+'_'+Nsnapstring+'.pdf'
#figname='/home/tkc004/samsonprogram/CRplot/CRIC/M1dc27g1D_2_21_smallts_'+Nsnapstring+'.pdf'
#figname='/home/tkc004/samsonprogram/CRplot/CRIC/'+filename+'.pdf'
#figname='/home/tkc004/samsonprogram/CRplot/CRIC/box_CR_Gauss_eps_0_1_1D.pdf'
print 'figname', figname
plt.savefig(figname,bbox_inches='tight')
#plt.clf()

if readflux==1:
        crfluxsep = crfx[::sep]*Ngas/Lbox
        plt.plot(xgsep, crfluxsep,ls='none',marker='o',ms=1, label='Sim')
        plt.xlabel('x (kpc)', fontsize=16)
        plt.ylabel(r'$F_{\rm cr}$', fontsize=16)
        #plt.legend()
        #plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
        figname='/home/tkc004/samsonprogram/CRplot/CRIC/Fcr_'+filename+'_'+Nsnapstring+'.pdf'
        #figname='/home/tkc004/samsonprogram/CRplot/CRIC/M1dc27g1D_2_21_smallts_'+Nsnapstring+'.pdf'
        #figname='/home/tkc004/samsonprogram/CRplot/CRIC/'+filename+'.pdf'
        #figname='/home/tkc004/samsonprogram/CRplot/CRIC/box_CR_Gauss_eps_0_1_1D.pdf'
        plt.savefig(figname,bbox_inches='tight')
        plt.clf()


        vst = crfluxsep/credensep/CRgamma
        plt.plot(xgsep, vst,ls='none',marker='o',ms=1, label='Sim')
        plt.xlabel('x (kpc)', fontsize=16)
        plt.ylabel(r'$V_{\rm st}$', fontsize=16)
        plt.legend()
        plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
        figname='/home/tkc004/samsonprogram/CRplot/CRIC/Vst_'+filename+'_'+Nsnapstring+'.pdf'
        #figname='/home/tkc004/samsonprogram/CRplot/CRIC/M1dc27g1D_2_21_smallts_'+Nsnapstring+'.pdf'
        #figname='/home/tkc004/samsonprogram/CRplot/CRIC/'+filename+'.pdf'
        #figname='/home/tkc004/samsonprogram/CRplot/CRIC/box_CR_Gauss_eps_0_1_1D.pdf'
        plt.savefig(figname,bbox_inches='tight')
        plt.clf()

#       vasep=va[::sep]
#       print 'va', vasep
 #       plt.plot(xgsep, vasep,ls='none',marker='o',ms=1, label='Sim')
  #      plt.xlabel('x (kpc)', fontsize=16)
   #     plt.ylabel(r'$V_{\rm st}$', fontsize=16)
    #    plt.legend()
#        plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
#        figname='/home/tkc004/samsonprogram/CRplot/CRIC/VA_'+filename+'_'+Nsnapstring+'.pdf'
        #figname='/home/tkc004/samsonprogram/CRplot/CRIC/M1dc27g1D_2_21_smallts_'+Nsnapstring+'.pdf'
        #figname='/home/tkc004/samsonprogram/CRplot/CRIC/'+filename+'.pdf'
        #figname='/home/tkc004/samsonprogram/CRplot/CRIC/box_CR_Gauss_eps_0_1_1D.pdf'
#        plt.savefig(figname,bbox_inches='tight')
#        plt.clf()

