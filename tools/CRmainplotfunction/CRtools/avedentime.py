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
import collections



dirneed=['bwsbclrmhd','bwsbclrdc0','bwsbclrdc27','bwsbclrdc28']

#dirneed=['bwsbclrdc29','bwsbclrstr','bwsbclrdc28mhd','bwsbclrdc28str']
#dirneed=['bwsbclr','bwsbclrmhd','bwsbclrdc0','bwsbclrdc27','bwsbclrdc28','bwsbclrdc29','bwsbclrstr','bwsbclrdc28mhd','bwsbclrdc28str']
wanted='avedenr'
print 'wanted', wanted
startno=200
Nsnap=651
snapsep=1
the_prefix='snapshot'
the_suffix='.hdf5'
fmeat=''


print 'wanted', wanted
print 'fmeat', fmeat
print 'runtodo', dirneed



if wanted=='gmr' or wanted=='avedenr' or wanted=='avecrdenr'\
or wanted=='avethdenr' or wanted=='aveedenr' or wanted=='FpFgr'\
or wanted=='avekedenr':
    plt.figure(figsize=(5,4))
    def funcab(x, a, b):
            return a+b*x
    def func(x, a, b, c):
            return a+b*x+c*x*x
    def d10func(r, a, b, c):
            return np.power(10.,func(np.log10(r), a, b, c))*(b+2*c*np.log(r)/np.log(10.))/r
    nogrid=20
    withinr=2.0
    atstarburst=1
    linelabelneed=0
    if wanted=='aveedenr':
        linelabelneed=1
    labcount=0
    newlabelneed=0
    legendneed=0
    for runtodo in dirneed:
        snapshot_list=np.array([])
        avedenlist=np.array([])
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
            snumadd=info['snumadd']
            ptitle=title
            labelneed=dclabel
            if newlabelneed==1:
                    labelneed="\n".join(wrap(newlabel,17))
            if runtitle=='SMC':
                ptitle='Dwarf'
            elif runtitle=='SBC':
                ptitle='Starburst'
            elif runtitle=='MW':
                ptitle=r'$L\star$ Galaxy'
            lsn='solid'
            if haveB>0:
                lsn='dashed'
            if havecr==0 and wanted=='avecrdenr':
                continue
            if cosmo==1:
                h0=1
                rotface=1
                datasup=0;
            else:
                h0=0
                rotface=0
                datasup=1;
            datasup=1
            G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
            Gp = G['p']
            Gv = G['v']
            Gu = G['u']
            Gm = G['m']
            Gvx = Gv[:,0]
            Gvy = Gv[:,1]
            Gvz = Gv[:,2]
            GEint = Gu*km_in_cm*km_in_cm*Gm*1e10*Msun_in_g
            cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
            ke = 0.5*(Gvx*Gvx+Gvy*Gvy+Gvz*Gvz)*km_in_cm*km_in_cm*Gm*1e10*Msun_in_g  
            #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
            Gz = Gp[:,2]
            Gx = Gp[:,0]
            Gy = Gp[:,1]
            dr = withinr/nogrid
            from crtestfunction import findcenz
            datasup=1
            xcen,ycen,zcen = findcennew(runtodo,Nsnap,withinr=2.5,dir='x',datasup=datasup,Gx=Gx,Gy=Gy,Gz=Gz,Gm=Gm)
            Gz = Gz-zcen; Gy = Gy-ycen; Gx = Gx-xcen;
            GEint = Gu*km_in_cm*km_in_cm*Gm*1e10*Msun_in_g
            cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
            ke = 0.5*(Gvx*Gvx+Gvy*Gvy+Gvz*Gvz)*km_in_cm*km_in_cm*Gm*1e10*Msun_in_g   
            Gr = np.sqrt(Gx*Gx+Gy*Gy+Gz*Gz)
            dr = withinr/nogrid
            cut = Gr < 0.2 #kpc 
            Gmcut = Gm[cut]
            Gm_in_sun = np.sum(Gmcut*1e10)
            vol = np.power(0.2*kpc_in_cm,3.)*4./3.*np.pi
            den = Gm_in_sun*Msun_in_g/vol/protonmass_in_g
            avedenlist=np.append(avedenlist,den)
            snapshot_list=np.append(snapshot_list,i)
        indexmax = np.argmax(avedenlist)
        avedenmax = avedenlist[indexmax]
        snapshotmax = snapshot_list[indexmax]
        filename = programdir+'/data/'+runtodo+'.txt'
        fh = open(filename, "w")
        fh.write('avedenmax  '+str(avedenmax)+'\n')
        fh.write('snapshotmax  '+str(snapshotmax)+'\n')
        for jno in range(len(snapshot_list)):
            fh.write('snapshot  '+str(snapshot_list[jno])+'  aveden  '+str(avedenlist[jno])+'\n')
        fh.close()
            
