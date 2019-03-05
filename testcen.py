from samson_functions import readsnapwcen, outdirname
import numpy as np


the_prefix='snapshot'
the_suffix='.hdf5'

runtodo = 'f61'
i=600


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
h0=cosmo


Dextra = readsnapwcen(the_snapdir, Nsnapstring, 1, snapshot_name=the_prefix, extension=the_suffix,\
havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,runtodo=runtodo,rundir=rundir,halostr=halostr)
Dx = Dextra['x']; Dy = Dextra['y']; Dz = Dextra['z'];
Dvx = Dextra['vx']; Dvy = Dextra['vy']; Dvz = Dextra['vz'];

Dr = np.sqrt(Dx*Dx+Dy*Dy+Dz*Dz)
Drcut = Dr<200
print 'Dv', np.average(Dvx), np.average(Dvy), np.average(Dvz)


