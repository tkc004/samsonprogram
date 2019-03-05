
# coding: utf-8

### Loading the data

# First we set up our imports:

# In[ ]:
import os
from samson_const import *
import matplotlib
matplotlib.use('Agg')
import yt
import numpy as np
import yt.units as units
import pylab
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib import rcParams
from yt.fields.particle_fields import add_volume_weighted_smoothed_field
from matplotlib.colors import LogNorm
from textwrap import wrap
from readsnap_cr import readsnapcr
import Sasha_functions as SF
from samson_functions import *
from yt.frontends.gizmo.api import GizmoDataset


#rcParams['figure.figsize'] = 15, 5
rcParams['font.size']=18
rcParams['font.family']='serif'
#rcParams.update({'figure.autolayout': True})
rcParams['axes.linewidth'] = 2
rcParams['pdf.fonttype'] = 42
rcParams['text.usetex']=True
rcParams['legend.fontsize']=10
#rcParams['ps.fonttype'] = 42
#rcParams['ps.useafm'] = True
#rcParams['pdf.use14corefonts'] = True
#rcParams['axes.unicode_minus']=False



def add_inner_title(ax, title, loc, size=None, **kwargs):
    from matplotlib.offsetbox import AnchoredText
    from matplotlib.patheffects import withStroke
    if size is None:
        size = dict(size=plt.rcParams['legend.fontsize'])
    at = AnchoredText(title, loc=loc, prop=size,
                      pad=0., borderpad=0.5,
                      frameon=False, **kwargs)
    ax.add_artist(at)
    at.txt._text.set_path_effects([withStroke(foreground="w", linewidth=3)])
    return at


# First we load the data set, specifying both the unit length/mass/velocity, as well as the size of the bounding box (which should encapsulate all the particles in the data set)
# 
# At the end, we flatten the data into "ad" in case we want access to the raw simulation data

# >This dataset is available for download at http://yt-project.org/data/GadgetDiskGalaxy.tar.gz (430 MB).

# In[ ]:




#Nsnap=100
Nsnap=500
#Nsnap=580
#Nsnap=600

fns=[\
#'bwsmclrmhd',\
#'bwmwmr',\
#'bwmwmrmhd',\
#'bwmwmrdc0',\
#'bwmwmrdc27',\
#'bwmwmrdc28',\
'bwmwmrdc29',\
#'bwmwmrstr',\
#'bwmwmrdc28mhd',\
#'bwmwmrdc28str',\
#'bwmwmrdc28mhd',\
#'m11dmhdcv',\
#'m11dcr_b_70',\
#'m11dcr_700',\
#'fm11q',\
#'m11qmhdcv',\
#'f476',\
#'m11bmhdcv',\
#'m11bcr_b_70',\
#'m11bcr_700',\
#'m12imhdcv',\
#'m12icr_b_70',\
#'m12icr_700',\
#'m12mmhdcv',\
#'m12mcr_b_70',\
#'m12mcr_700',\
#'f61'
]


#fns=[\
#'/home/tkc004/oasis/bw/smc/smc_lr_2_21/output/snapshot_490.hdf5',\
#'/home/tkc004/oasis/bw/smc/smc_lr_2_21_mhd/output/snapshot_490.hdf5',\
#'/home/tkc004/oasis/bw/smc/smc_cr_lr_dc0_2_21/output/snapshot_490.hdf5',\
#'/home/tkc004/oasis/bw/smc/smc_cr_lr_dc27_2_21_M1_c500/output/snapshot_490.hdf5',\
#'/home/tkc004/oasis/bw/smc/smc_cr_lr_dc28_2_21_M1_c1000/output/snapshot_490.hdf5',\
#'/home/tkc004/oasis/bw/smc/smc_cr_lr_dc29_2_21_M1_c2000/output/snapshot_490.hdf5',\
#'/home/tkc004/oasis/bw/smc/smc_cr_lr_2_21_purestream/output/snapshot_490.hdf5',\
#'/home/tkc004/oasis/bw/smc/smc_cr_lr_dc28_2_21_M1_mhd_c1000/output/snapshot_490.hdf5',\
#'/home/tkc004/oasis/bw/smc/smc_cr_lr_dc28_2_21_M1_mhd_stream_c1000/output/snapshot_490.hdf5'\
#]


plotit=[]

#plotit = [\
#'Hydro no CR',\
#'MHD no CR',\
#'Advection',\
#r'$\kappa$=3e27',\
#r'$\kappa$=3e28',\
#r'$\kappa$=3e29',\
#r'MHD             Streaming',\
##STRLL',\
#r'MHD $\kappa$=3e28',\
#'MHD stream'
#r'MHD $\kappa$=3e28 Streaming'\
#'Diffusion      Streaming'\
#]

#plotit = [\
#'STREAM',\
#'STREAMLLVA',\
#'STREAMLL4VA'\
#]

#plotit = ['DC=3e27','DC=3e28','DC=3e29']
#plotit = ['DC=0','DC=3e27','DC=3e28']


plotupt=[]

#plotupt=['','SMC','']
#plotupt=['','MHD']
#plotupt=['Milky Way']
#plotupt=['Milky Way', '']
#plotupt=[\
#'',\
#'Starburst',\
#'Dwarf',\
#r'$L\star$ Galaxy',\
#'',\
#'',\
#'',\
#'',\
#'',\
#'',\
#''\
#]
#plotupt=['','Milky Way', '']
#plotupt=['','Dwarf','']


#quivercolor=0
#logquiver=0

#wanted = 'gamma'
#wanted = 'temp'
wanted = 'density'
#wanted = 'cosmic_ray'
#wanted = 'Bfield'

#wanted='ke'

#wanted='zvel'


#remark='MWLRMHDIC_x'


galname=fns[-1]
#galname='m12icr_700'
#galname='m12icr_b_70'
#galname='mw'
#galname='smc'
#galname='sbc'

#remark='mrmhd'
#remark='mrthree'
#remark='mrshort'
remark=str(Nsnap)
#remark='lrmhd'
#remark='streamlr'
#remark='mrz'

#remark = 'CRDC_hole_z'

mode = 'Slice'
#mode = 'projected'

#needquiver=''
needquiver='v'
#needquiver='B'

projectionaxis='x'

rotface=1
faceon= 0 #work only with rotface

newlabelneed=1

if faceon==1:
        remark+='faceon'
elif rotface==1:
        remark+='rotface'

labfs=10
labs=10

if needquiver=='B':
        logquiver=1
        headaxislength=0
        headwidth=0
        headlength=0
else:
        logquiver=0
        headaxislength=4.5
        headwidth=3.
        headlength=5.


fig = plt.figure()

unit_base = {'UnitLength_in_cm'         : 3.08568e+21,
             'UnitMass_in_g'            :   1.989e+43,
             'UnitVelocity_in_cm_per_s' :      100000}


bbox_lim = 1e5 #kpc

bbox = [[-bbox_lim,bbox_lim],
        [-bbox_lim,bbox_lim],
        [-bbox_lim,bbox_lim]]

#newboxsize=75
 
newboxsize=150

projectionaxis='x'

nrows=1
ncols=len(fns)
#if len(fns)>4:
#       nrows = 2
#       ncols = 4
#if len(fns)>8:
#        nrows = 3
#        ncols = 4
if len(fns)>3:
        nrows = 2
        ncols = 3
if len(fns)>6:
        nrows = 3
        ncols = 3

#grid = AxesGrid(fig, (0.07,0.07,0.8,0.8),
grid = AxesGrid(fig, 111,
                nrows_ncols = (nrows, ncols),
                axes_pad = 0.05,
                label_mode = "L",
                share_all = True,
                cbar_location="right",
                cbar_mode="single",
                cbar_size="3%",
                cbar_pad="0%")

def symlog(x):
    """ Returns the symmetric log10 value """
    return np.sign(x) * np.log10(np.abs(x))


for i, runtodo in enumerate(fns):
        #ds = yt.load(fname,unit_base=unit_base)
        info=outdirname(runtodo, Nsnap)
        havecr=info['havecr']
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
        multifile=info['multifile']
        usepep=info['usepep']
        snumadd=info['snumadd']
        labelneed=dclabel
        print 'dclabel', dclabel
        if newlabelneed==1:
                labelneed="\n".join(wrap(newlabel,17))
        plotit.append(labelneed)
        if runtitle=='SMC':
                ptitle='Dwarf'
        elif runtitle=='SBC':
                ptitle='Starburst'
        elif runtitle=='MW':
                ptitle=r'$L\star$ Galaxy'
        if i==1 and cosmo==0:
                plotupt.append(ptitle)
        else:
                plotupt.append('')
        print 'labelneed', labelneed
        if multifile=='y':
                fname=the_snapdir+'/snapdir_'+Nsnapstring+'/snapshot_'+Nsnapstring+'.0.hdf5'
        else:
                fname=the_snapdir+'/snapshot_'+Nsnapstring+'.hdf5'
        print 'fname', fname
        print 'usepep', usepep
        if cosmo==1:
                header=readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1,header_only=1)
                h0 = header['hubble']
                atime = header['time']
                print 'halostr', halostr
                if usepep==1:
                        halosA = SF.read_halo_history_pep(rundir, Nsnap, singlesnap=1, firever=firever,halonostr=halostr, hubble=h0, comoving=1, maindir=maindir)
                        redlist = halosA['redshift']
                        haloid = halosA['ID']
                        xcen = halosA['x']*atime
                        ycen = halosA['y']*atime
                        zcen = halosA['z']*atime
                        vxcen = halosA['xv']
                        vycen = halosA['yv']
                        vzcen = halosA['zv']
                else:
                        halosA = SF.read_halo_history(rundir, maindir=maindir, halonostr=halostr,\
                         hubble=h0, comoving=0, snumadd=snumadd)
                        redlist = halosA['redshift']
                        haloid = halosA['ID']
                        a_scale = 1.0/(1.0+redlist)
                        xcenl = halosA['x']
                        ycenl = halosA['y']
                        zcenl = halosA['z']
                        vxcenl = halosA['xv']
                        vycenl = halosA['yv']
                        vzcenl = halosA['zv']
                        xcen = np.interp(atime,a_scale,xcenl)
                        ycen = np.interp(atime,a_scale,ycenl)
                        zcen = np.interp(atime,a_scale,zcenl)
                        vxcen = np.interp(atime,a_scale,vxcenl)
                        vycen = np.interp(atime,a_scale,vycenl)
                        vzcen = np.interp(atime,a_scale,vzcenl)

                #if rotface==1:
                        #faceon:
                        #projectionaxis=[Lxstar,Lystar,Lzstar]
                        #edgeon
                        #projectionaxis=[1,-Lxstar/Lystar,0]
                        #north_vector = [Lxstar, Lystar, Lzstar]
                        #projectionaxis=[0,0,1]
        else:
                xcen=ycen=zcen=vxcen=vycen=vzcen=0.



        print 'loading', fname
        ds = GizmoDataset(fname,unit_base=unit_base,bounding_box=bbox)
        ds.index
        ad= ds.all_data()
        sorted(ds.field_list)
        #print 'ds.field_list', ds.field_list
        ad = ds.all_data()


        center = ds.arr([xcen,ycen,zcen],'kpc')
        new_box_size = ds.quan(newboxsize,'kpc')
        left_edge = center - new_box_size/2
        right_edge = center + new_box_size/2
        ad2= ds.region(center=center, left_edge=left_edge, right_edge=right_edge)
        sp = ds.sphere(center, (6, 'kpc'))
        #angmom = sp.quantities.angular_momentum_vector(use_gas=True,use_particles=False)
        if rotface==1:
                angmom = sp.quantities.angular_momentum_vector(use_gas=False,use_particles=True,particle_type='PartType4')
                angmom = angmom/np.sqrt(angmom[0]*angmom[0]+angmom[1]*angmom[1]+angmom[2]*angmom[2])
                nx = angmom[2]; ny = 0.0; nz = -angmom[0]; 
                projectionaxis = [nx,ny,nz]
                north_vector = angmom
                if faceon==1:
                        projectionaxis = angmom
                        north_vector = [nx,ny,nz]
                
        


        if (wanted == 'temp'):
                if mode == 'projected':
                        px = yt.ProjectionPlot(ds, projectionaxis, ('gas', 'temperature'), center=center, width=new_box_size, weight_field=('gas','density'))
                        cblabel = r'projected ${\rm T \;[K]}$'
                if mode == 'Slice':
                        if rotface==1:
                                px = yt.OffAxisSlicePlot(ds, projectionaxis, ('gas', 'temperature'), center=center, width=new_box_size, north_vector=north_vector)
                        else:
                                px = yt.SlicePlot(ds, projectionaxis, ('gas', 'temperature'), center=center, width=new_box_size)
                        cblabel = r'${\rm T \;[K]}$'
                px_frb = px.data_source.to_frb((newboxsize, "kpc"), 128)
                px_dens = np.array(px_frb[ ('gas', 'temperature')])
                cbcolor = 'hot'

        if (wanted == 'Bfield'):
                def _Benergydensityfun(field, data):
                    Bf = data['PartType0', 'MagneticField'].in_cgs()
                    print 'Bf.shape', Bf.shape
                    Bdensity = (Bf[:,0]*Bf[:,0]+Bf[:,1]*Bf[:,1]+Bf[:,2]*Bf[:,2])/8./np.pi
                    return Bdensity
                ds.add_field(('PartType0', 'Benergydensity'), function=_Benergydensityfun,
                             units="erg/cm**3", particle_type=True, display_name="Benergydensity")

                add_volume_weighted_smoothed_field('PartType0', 'Coordinates', 'Masses',
                                                   'SmoothingLength', 'Density',
                                                   'Benergydensity', ds.field_info)
                ad2= ds.region(center=center, left_edge=left_edge, right_edge=right_edge)
                if mode == 'projected':
                        px = yt.ProjectionPlot(ds, projectionaxis, ('deposit', 'PartType0_smoothed_Benergydensity'), center=center, width=new_box_size, north_vector=north_vector)
                if mode == 'Slice':
                        if rotface==1:
                                px = yt.OffAxisSlicePlot(ds, projectionaxis, ('deposit', 'PartType0_smoothed_Benergydensity'),\
                                 center=center, width=new_box_size, north_vector=north_vector)
                        else:
                                px = yt.SlicePlot(ds, projectionaxis, ('deposit', 'PartType0_smoothed_Benergydensity'), center=center, width=new_box_size)
                px_frb = px.data_source.to_frb((newboxsize, "kpc"), 128)
                px_dens = np.array(px_frb[('deposit', 'PartType0_smoothed_Benergydensity')])
                cbcolor = 'RdPu'
                cblabel = r'${\rm B\;energy\;[ erg/cm^3}]$'


        if wanted=='density':
                if mode == 'projected':
                        px = yt.Plot(ds, projectionaxis, ('gas', 'density'), center=center, width=new_box_size)
                if mode == 'Slice':
                        if rotface==1:
                                px = yt.OffAxisSlicePlot(ds, projectionaxis, ('gas', 'density'), center, new_box_size,north_vector=north_vector)
                        else:
                                px = yt.SlicePlot(ds, projectionaxis, ('gas', 'density'), center=center, width=new_box_size)    
                px_frb = px.data_source.to_frb((newboxsize, "kpc"), 128)
                px_dens = np.array(px_frb[('gas', 'density')])
                cbcolor='YlOrRd'
                cblabel=r'${\rm \log (\rho [g/cm^3])}$'

        if (wanted == 'cosmic_ray'):
                def _cosmicrayenergydensity(field, data):
                    energy = data['PartType0', 'CosmicRayEnergy']
                    energy = yt.YTArray(1e10*Msun_in_g*km_in_cm*km_in_cm*erg_in_eV*energy,"eV")
                    mass = data['PartType0', 'Masses'].in_cgs()
                    density = data['PartType0', 'Density'].in_cgs()
                    crdensity = energy/mass*density
                    return crdensity
                ds.add_field(('PartType0', 'cosmicrayenergydensity'), function=_cosmicrayenergydensity,
                    units="eV/cm**3", particle_type=True)

                add_volume_weighted_smoothed_field('PartType0', 'Coordinates', 'Masses',
                                   'SmoothingLength', 'Density',
                                   'cosmicrayenergydensity', ds.field_info)
                if mode=='projected':
                        px = yt.ProjectionPlot(ds, projectionaxis, ('deposit', 'PartType0_smoothed_cosmicrayenergydensity'), center=center, width=new_box_size)
                if mode=='Slice':
                        if rotface==1:
                                px = yt.OffAxisSlicePlot(ds, projectionaxis, ('deposit', 'PartType0_smoothed_cosmicrayenergydensity'), center=center, width=new_box_size,north_vector=north_vector)
                        else:
                                px = yt.SlicePlot(ds, projectionaxis, ('deposit', 'PartType0_smoothed_cosmicrayenergydensity'), center=center, width=new_box_size)
                #px.set_unit(('deposit', 'PartType0_smoothed_cosmicrayenergydensity'), 'eV/cm**3')
                px_frb = px.data_source.to_frb((newboxsize, "kpc"), 128)
                px_dens = np.array(px_frb[('deposit', 'PartType0_smoothed_cosmicrayenergydensity')])
                cbcolor = 'YlGnBu'
                cblabel = r'${\rm \log (e_{CR} [eV/cm^3])}$'


        if (wanted=='gamma'):
                def _gammaenergydensity(field, data):
                    #Fcal = 7e-2
                    energy = data['PartType0', 'CosmicRayEnergy']
                    energy = yt.YTArray(energy*2e53, "g*cm**2/s**2")
                    mass = data['PartType0', 'Masses'].in_cgs()
                    density = data['PartType0', 'Density'].in_cgs()
                    crdensity = energy/mass*density
                    tpi = 2e5/density*yt.YTArray(250.0,"g/cm**3")
                    gammadensity = crdensity*0.25/tpi/3.2e7*0.7/3.0
                    return gammadensity
                ds.add_field(('PartType0', 'gammaenergydensity'), function=_gammaenergydensity,
                             units="erg/cm**3", particle_type=True, display_name="Gamma ray Energy Density")
                add_volume_weighted_smoothed_field('PartType0', 'Coordinates', 'Masses',
                                                   'SmoothingLength', 'Density',
                                                   'gammaenergydensity', ds.field_info)
                ds.index
                ad= ds.all_data()
                center = ds.arr([0,0,0],'kpc')
                new_box_size = ds.quan(newboxsize,'kpc')
                left_edge = center - new_box_size/2
                right_edge = center + new_box_size/2
                ad2= ds.region(center=center, left_edge=left_edge, right_edge=right_edge)
                if mode=='projected':
                        px = yt.ProjectionPlot(ds, projectionaxis, ('deposit', 'PartType0_smoothed_gammaenergydensity'), center=center, width=new_box_size)
                #       px.set_zlim('all',1e-28,1e-38)
                        px.set_colorbar_label(('deposit', 'PartType0_smoothed_gammaenergydensity'),r'Projected $\rho_{\gamma}\; \mathrm{(\frac{erg}{cm^2\;s})}$')
                        cblabel = r'Projected $\rho_{\gamma}\; \mathrm{(\frac{erg}{cm^2\; s})}$'
                if mode=='Slice':
                        if rotface==1:
                                px = yt.OffAxisSlicePlot(ds, projectionaxis, ('deposit', 'PartType0_smoothed_gammaenergydensity'), center=center, width=new_box_size,north_vector=north_vector)
                        else:
                                px = yt.SlicePlot(ds, projectionaxis, ('deposit', 'PartType0_smoothed_gammaenergydensity'), center=center, width=new_box_size)
                        px.set_colorbar_label(('deposit', 'PartType0_smoothed_gammaenergydensity'),r'$\rho_{\gamma}\; \mathrm{(\frac{erg}{cm^3})}$')
                        cblabel = r'$\rho_{\gamma}\; \mathrm{(\frac{erg}{cm^3})}$'
                px_frb = px.data_source.to_frb((newboxsize, "kpc"), 128)
                px_dens = np.array(px_frb[('deposit', 'PartType0_smoothed_gammaenergydensity')])
                cbcolor = 'hot'

        x_min, x_max, y_min, y_max = px_frb.bounds
        # Extract X, Y, U, V from the frb
        if needquiver=='B':
                def _Bx(field, data):
                    Bf = data['PartType0', 'MagneticField'].in_cgs()
                    By = Bf[:,0]
                    return By
                ds.add_field(('PartType0', 'Bfield_x'), function=_Bx,
                             units="gauss", particle_type=True, display_name="Bx")
                add_volume_weighted_smoothed_field('PartType0', 'Coordinates', 'Masses',\
                                   'SmoothingLength', 'Density',\
                                   'Bfield_x', ds.field_info)
                def _By(field, data):
                    Bf = data['PartType0', 'MagneticField'].in_cgs()
                    By = Bf[:,1]
                    return By
                ds.add_field(('PartType0', 'Bfield_y'), function=_By,
                             units="gauss", particle_type=True, display_name="By")
                add_volume_weighted_smoothed_field('PartType0', 'Coordinates', 'Masses',\
                                   'SmoothingLength', 'Density',\
                                   'Bfield_y', ds.field_info)
                def _Bz(field, data):
                    Bf = data['PartType0', 'MagneticField'].in_cgs()
                    Bz = Bf[:,2]
                    return Bz
                ds.add_field(('PartType0', 'Bfield_z'), function=_Bz,
                             units="gauss", particle_type=True, display_name="Bz")
                add_volume_weighted_smoothed_field('PartType0', 'Coordinates', 'Masses',\
                                   'SmoothingLength', 'Density',\
                                   'Bfield_z', ds.field_info)

                if projectionaxis=='x':
                        fvx = ('deposit', 'PartType0_smoothed_Bfield_y')
                        fvy = ('deposit', 'PartType0_smoothed_Bfield_z')
                if projectionaxis=='y':
                        fvx = ('deposit', 'PartType0_smoothed_Bfield_x')
                        fvy = ('deposit', 'PartType0_smoothed_Bfield_z')
                if projectionaxis=='z':
                        fvx = ('deposit', 'PartType0_smoothed_Bfield_x')
                        fvy = ('deposit', 'PartType0_smoothed_Bfield_y')
                if rotface==1:
                    def _rotBx(field, data):
                        Bf = data['PartType0', 'MagneticField'].in_cgs()           
                        vpx = -np.cross(projectionaxis,north_vector)
                        Byf = np.dot(Bf,vpx)
                        return Byf
                    ds.add_field(('PartType0', 'Bfx'), function=_rotBx,
                    units="gauss", particle_type=True, display_name="Bfx")
                    add_volume_weighted_smoothed_field('PartType0', 'Coordinates', 'Masses',\
                        'SmoothingLength', 'Density',\
                        'Bfx', ds.field_info)
                    def _rotBy(field, data):
                        Bf = data['PartType0', 'MagneticField'].in_cgs()
                        Bzf = np.dot(Bf,north_vector)
                        return Bzf
                    ds.add_field(('PartType0', 'Bfy'), function=_rotBy,
                        units="gauss", particle_type=True, display_name="Bfy")
                    add_volume_weighted_smoothed_field('PartType0', 'Coordinates', 'Masses',\
                        'SmoothingLength', 'Density',\
                        'Bfy', ds.field_info)
                    fvx = ('deposit', 'PartType0_smoothed_Bfx')
                    fvy = ('deposit', 'PartType0_smoothed_Bfy')
        elif needquiver=='v':
                if projectionaxis=='x':
                        fvx = ('gas', 'velocity_y')
                        fvy = ('gas', 'velocity_z')
                if projectionaxis=='y':
                        fvx = ('gas', 'velocity_x')
                        fvy = ('gas', 'velocity_z')
                if projectionaxis=='z':
                        fvx = ('gas', 'velocity_x')
                        fvy = ('gas', 'velocity_y')
                if rotface==1:
                        vcen=ds.arr([vxcen*1e5,vycen*1e5,vzcen*1e5],'cm/s') 
                        def _rotvx(field, data):
                            vi = data['PartType0', 'Velocities'].in_cgs()           
                            vic = vi-vcen
                            vpx = -np.cross(projectionaxis,north_vector)
                            vyf = np.dot(vic,vpx)
                            vyf[vyf>1e9]=1e9
                            vyf[vyf<-1e9]=-1e9
                            return vyf
                        ds.add_field(('PartType0', 'Vel_x'), function=_rotvx,
                                     units="cm/s", particle_type=True, display_name="Vx")
                        add_volume_weighted_smoothed_field('PartType0', 'Coordinates', 'Masses',\
                                           'SmoothingLength', 'Density',\
                                           'Vel_x', ds.field_info)
                        def _rotvy(field, data):
                            vi = data['PartType0', 'Velocities'].in_cgs()
                            vic = vi-vcen
                            vzf = np.dot(vic,north_vector)
                            vzf[vzf>1e9]=1e9
                            vzf[vzf<-1e9]=-1e9
                            return vzf
                        ds.add_field(('PartType0', 'Vel_y'), function=_rotvy,
                                     units="cm/s", particle_type=True, display_name="Vy")
                        add_volume_weighted_smoothed_field('PartType0', 'Coordinates', 'Masses',\
                                           'SmoothingLength', 'Density',\
                                           'Vel_y', ds.field_info)
                        fvx = ('deposit', 'PartType0_smoothed_Vel_x')
                        fvy = ('deposit', 'PartType0_smoothed_Vel_y')
        if needquiver=='v' or needquiver=='B':
                if rotface==1:
                        p = yt.OffAxisSlicePlot(ds, projectionaxis, [fvx, fvy], center=center,\
                         width=new_box_size,north_vector=north_vector) # Here you can use ProjectionPlot too
                else:
                        p = yt.SlicePlot(ds, projectionaxis, [fvx, fvy], center=center, width=new_box_size) # Here you can use ProjectionPlot too
                p.set_buff_size(128)
                frb = p.data_source.to_frb((newboxsize, "kpc"), 128)
                x_min, x_max, y_min, y_max = frb.bounds
                nx, ny = frb.buff_size
                x_edges = np.linspace(x_min, x_max, nx + 1)
                y_edges = np.linspace(y_min, y_max, ny + 1)
                x = (x_edges[:-1] + x_edges[1:]) / 2
                y = (y_edges[:-1] + y_edges[1:]) / 2
                X, Y = np.meshgrid(x, y)

                # Convert to the units you want
                X = X.to('kpc').v
                Y = Y.to('kpc').v
                if needquiver=='v':
                        U = p.frb[fvx].to('km/s').v
                        V = p.frb[fvy].to('km/s').v
                elif needquiver=='B':
                        U = p.frb[fvx].to('gauss').v
                        V = p.frb[fvy].to('gauss').v
        logpx_dens = np.nan_to_num(np.log10(px_dens))
        if i==0:
                vmax=np.amax(logpx_dens)
                vmin=np.amin(logpx_dens)
                print 'vmax', vmax
                print 'vmin', vmin
                if wanted=='temp':
                        if galname=='smc':
                                vmax=6.0
                                vim=1.0
                        else:
                                vmax=6.5
                                vmin=3.0
                if wanted=='density':
                        vmax=-25.0
                        vmin=vmax-4.0
                if wanted=='cosmic_ray':
                        vmin=vmax-7.0
                if wanted=='gamma':
                        if vmin-vmax<-20:
                                vmin=vmax-20
                if wanted=='Bfield':
                        if vmin-vmax<-20:
                                vmin=vmax-20
        # Use matplotlib from there
        im = grid[i].imshow(logpx_dens, origin='lower', vmax=vmax, vmin=vmin, extent=[x_min,x_max,y_min,y_max], cmap=cbcolor)
        qsep=10
        if len(fns)<4:
                ineed=0
        else:
                ineed=3
        if logquiver==1:
                if needquiver=='B':
                        angles=np.arctan2(V[::qsep, ::qsep],U[::qsep, ::qsep])*180.0/np.pi # calculate angles manually
                        Q=grid[i].quiver(X[::qsep, ::qsep], Y[::qsep, ::qsep], symlog(U[::qsep, ::qsep]/1e-12),\
                         symlog(V[::qsep, ::qsep]/1e-12), units='x', angles=angles,\
                         headwidth=headwidth, headaxislength=headaxislength, headlength=headlength)
        elif (needquiver=='B' or needquiver=='v'):
                Q=grid[i].quiver(X[::qsep, ::qsep], Y[::qsep, ::qsep], U[::qsep, ::qsep], V[::qsep, ::qsep], units='x')
        if i==ineed and (needquiver=='B' or needquiver=='v'):
                if needquiver=='v':
                        grid[i].quiverkey(Q, 0.3, 0.3, 100, r'${\rm 100 km/s}$', labelpos='S',coordinates='inches', fontproperties={'size': 'xx-small'})
                elif needquiver=='B':
                        if logquiver==1:
                                grid[i].quiverkey(Q, 0.4, 0.4, np.log10(100000.0) , r'${\rm log(0.1\mu G)}$', labelpos='S',coordinates='inches', fontproperties={'size': 'xx-small'})
                        elif quivercolor==1:
                                print 'use color not length'
                        else:
                                grid[i].quiverkey(Q, 0.3, 0.3, 1e-8, r'${\rm 10^{-2} \mu G}$', labelpos='S',coordinates='inches', fontproperties={'size': 'xx-small'})
        grid[i].set_xlabel('kpc',fontsize=labfs)
        grid[i].set_ylabel('kpc',fontsize=labfs)
        grid[i].tick_params(axis='x', labelsize=labs)
        grid[i].tick_params(axis='y', labelsize=labs)
        grid[i].set_xticks(np.arange(-newboxsize/2.+10.,newboxsize/2.-10.,20.))

cbar = grid.cbar_axes[0].colorbar(im, extend='max')
cbar.ax.set_ylabel(cblabel, rotation=270,fontsize=12, labelpad=15)
cbar.ax.tick_params(labelsize=10)


print 'plotit, plotupt', plotit, plotupt
for ax, im_title, up_title in zip(grid, plotit, plotupt):
      t = add_inner_title(ax, "\n".join(wrap(im_title,17)), loc=3)
      t.patch.set_alpha(0.5)
      ax.set_title(up_title,fontsize=labfs)
                
if needquiver=='v':
        mode+='_vquiver'
if needquiver=='B':
        mode+='_Bquiver'
if headaxislength==0:
        mode+='_nohead'

#px.save('gas'+wanted+'_'+mode+'_'+remark+'.pdf')
directory = 'YTplot/'+galname
if not os.path.exists(directory):
    os.makedirs(directory)
directory = 'YTplot/'+galname+'/'+wanted
if not os.path.exists(directory):
    os.makedirs(directory)
filename='YTplot/'+galname+'/'+wanted+'/'+'gas'+wanted+'_'+mode+'_'+galname+remark+'.pdf'
print 'filename', filename
plt.savefig(filename,bbox_inches='tight',pad_inches = 0.1)
plt.clf()
