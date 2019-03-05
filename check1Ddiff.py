import matplotlib as mpl
mpl.use('macOsX')
mpl.use('Agg')
from readsnap_cr import *
from Sasha_functions import *
from samson_functions import *
from samson_const import *
from pathloc import *
import sys
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pylab import *
#rcParams['figure.figsize'] = 5, 5
rcParams['figure.figsize'] = 6, 5
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



ic=0
filenamec100='M1dc27g1D_eps_0_5_c100_noHLL'
the_snapdirc100 ='/home/tkc004/oasis/CRdifftest/'+filenamec100+'/output'
filenamec500='M1dc27g1D_eps_0_5_c500_noHLL'
the_snapdirc500 ='/home/tkc004/oasis/CRdifftest/'+filenamec500+'/output'
the_prefix = 'snapshot'
the_suffix = '.hdf5'
#Nsnapstring = '050'
Nsnapstring = '005'
the_ptype = 0
havecr=1
kappa =  10. #diffusion coefficient #in code unit #kpc^2/Gyr
#kappa =  3.3/0.78 #diffusion coefficient #in code unit #kpc^2/Gyr
sep=10


ftxt = open(programdir+'/data/Ecr_c100_t5Myr.dat', 'r')
ftxt.readline()
dars = ftxt.readlines()
ftxt.close()
xcrc100l=[]
Ecrc100l=[]

for line in dars:
        xsd = line.split()
        xcrc100l=np.append(xcrc100l, float(xsd[0]))
        Ecrc100l=np.append(Ecrc100l, float(xsd[1]))


ftxt = open(programdir+'/data/Ecr_c500_t5Myr.dat', 'r')
ftxt.readline()
dars = ftxt.readlines()
ftxt.close()
xcrc500l=[]
Ecrc500l=[]

for line in dars:
        xsd = line.split()
        xcrc500l=np.append(xcrc500l, float(xsd[0]))
        Ecrc500l=np.append(Ecrc500l, float(xsd[1]))

ftxt = open(programdir+'/data/Ecr_cc_t5Myr.dat', 'r')
ftxt.readline()
dars = ftxt.readlines()
ftxt.close()
xcrccl=[]
Ecrccl=[]

for line in dars:
        xsd = line.split()
        xcrccl=np.append(xcrccl, float(xsd[0]))
        Ecrccl=np.append(Ecrccl, float(xsd[1]))

#calculate pure diffusion solution

DIMS=1; # number of dimensions 
Ngas=2049; # 1D particle number (so total particle number is N_1D^DIMS)
Lbox = 10.0 # box side length #kpc
rho_desired = 0.00245 # box average initial gas density #1e10Msun/kpc^3 #6.805e-22g/cm^3 #407 cm^-3
P_desired = 1.0 # initial gas pressure
gamma_eos = 5./3. # polytropic index of ideal equation of state the run will assume
CRE = 1. # total CR energy #1e10Msun*(km/s)^2 #2e53 erg
eps = 0.5 # small offset in initial Gaussian packet # in code unit #kpc
timecr = 0.001*float(Nsnapstring) #by default
Nnum=300
xglin=np.linspace(0.0,5.0,num=Nnum)
cre_g=diffusionsol(CRE,kappa,timecr,xglin,offset=eps,dim=DIMS)


Gc100 = readsnapcr(the_snapdirc100, Nsnapstring, the_ptype, snapshot_name=the_prefix, extension=the_suffix,havecr=havecr,ic=ic)
posc100=Gc100['p']
xgc100 = posc100[:,0]
cregyc100 = np.array(Gc100['cregy'])
credenc100 = cregyc100*Ngas/Lbox
argxc100 = np.argsort(xgc100)
xgc100 = xgc100[argxc100]
credenc100 = credenc100[argxc100]
xgsepc100=xgc100[::sep]
credensepc100=credenc100[::sep]

Gc500 = readsnapcr(the_snapdirc500, Nsnapstring, the_ptype, snapshot_name=the_prefix, extension=the_suffix,havecr=havecr,ic=ic)
posc500=Gc500['p']
xgc500 = posc500[:,0]
cregyc500 = np.array(Gc500['cregy'])
credenc500 = cregyc500*Ngas/Lbox
argxc500 = np.argsort(xgc500)
xgc500 = xgc500[argxc500]
credenc500 = credenc500[argxc500]
xgsepc500=xgc500[::sep]
credensepc500=credenc500[::sep]

#print 'xcrc100l', xcrc100l
#print 'Ecrc100l', Ecrc100l

plt.plot(xglin, cre_g,lw=1.2, color='k',label='Diffusion')
plt.plot(xcrccl,Ecrccl,lw=1.2,ls='dashed', color='r', label=r'$\tilde{c}=c$')
plt.plot(xcrc500l,Ecrc500l,lw=1.2,ls='dotted', color='g', label=r'$\tilde{c}=$500')
plt.plot(xcrc100l,Ecrc100l,lw=1.2,ls='dashed', color='b', label=r'$\tilde{c}=$100')
plt.plot(xgsepc500, credensepc500,ls='none',color='g',marker='o',ms=5,mfc='none')
plt.plot(xgsepc100, credensepc100,ls='none',color='b',marker='o',ms=5,mfc='none', label='Simulations')


plt.xlim([0.0,2.0])
plt.xlabel('r [kpc]', fontsize=20)
plt.ylabel(r'$e_{\rm cr}$ [code unit]', fontsize=20)
plt.legend()
plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=20)
figname=plotloc+'/CRplot/CRdifftest/'+filenamec100+'_'+Nsnapstring+'.pdf'
plt.savefig(figname,bbox_inches='tight')
plt.clf()

