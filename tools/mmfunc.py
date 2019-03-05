import utilities as util
import os.path
import visualization.image_maker as imaker
import visualization.movie_maker as mmaker
import gadget
import numpy as np
import utilities as util
import os.path
import visualization.image_maker as imaker
import gadget
import sys
import samson_functions as SSF
import os, errno
from TKutilities import TKmkdir


def runmmfromrun(runtodo,i_snap_min,i_snap_max,theta0=90.,phi0=0.):   
    info=SSF.outdirname(runtodo, i_snap_max)
    snapdir=info['the_snapdir']
    outputdir_master = '/home/tkc004/scratch/philfig/'+runtodo
    TKmkdir(outputdir_master)
    print 'snapdir', snapdir
    mmfunc(snapdir,outputdir_master,i_snap_min,i_snap_max,theta0=theta0,phi0=phi0)
    
def mmfunc(snapdir,outputdir_master,i_snap_min,i_snap_max,theta0=90.,phi0=0.):
        op= mmaker.movie_maker(\
            xmax_of_box=30.0, #"""scale (in code units) of movie box size"""
            snapdir=snapdir,
            outputdir_master=outputdir_master,
            show_gasstarxray='star', #"""determines if image is 'stars','gas','xr' etc"""
            frames_per_gyr=50., #"""sets spacing of frames (even in time)"""
            cosmological=1, #"""is this is a cosmological (comoving) simulation?"""
            time_min=0, #"""initial time of movie"""
            time_max=0, #"""final time (if<=0, just set to time of maximum snapshot)"""
            frame_min=0, #"""initial frame number to process"""
            frame_max=1.0e10, #"""final frame number to process"""
            i_snap_min=i_snap_min, #"""minimum snapshot number in list to scan"""
            i_snap_max=i_snap_max, #"""maximum snapshot number to scan (if =0, will just use largest in snapdir)"""
            pixels=1440, #"""images will be pixels*pixels"""
            show_time_label=1, #"""do or don't place time label on each frame"""
            show_scale_label=1, #"""do or don't place physical scale label on each frame"""
            temp_max=1.0e6, #"""maximum gas temperature for color-weighted gas images (as 'temp_cuts' in 3-color)"""
            temp_min=3.0e2, #"""minimum gas temperature for color-weighted gas images (as 'temp_cuts' in 3-color)"""
            center_on_bh=0, #"""if have a bh particle, uses it for the centering"""
            center_on_com=0, #"""simple centering on center of mass"""
            set_fixed_center=0, #"""option to set a fixed center for the whole movie"""
            use_h0=1, #"""correct snapshots to physical units"""
            four_char=0, #"""snapshots have 4-character filenames"""
            skip_bh=0, #"""skips reading bh information in snapshots"""
            add_extension='', #"""extension for movie snapshot images"""
            theta_0=theta0, #viewing angles
            phi_0=phi0)

        
def runimfromrun(runtodo,snapnum,theta0=90.,phi0=0.):   
    sdir=''
    info=SSF.outdirname(runtodo, snapnum)
    snapdir_master=info['the_snapdir']
    print 'snapdir_master', snapdir_master
    outdir_master = '/home/tkc004/scratch/philfig/'+runtodo
    TKmkdir(outdir_master)
    imfunc(sdir,snapnum,outdir_master=outdir_master,theta=theta0,phi=phi0)
        
def imfunc(sdir,snapnum,snapdir_master='',outdir_master='',snapshot_subdirectory='',theta=90.,phi=0.):
    imaker.image_maker( sdir, snapnum, \
    snapdir_master=snapdir_master,
    outdir_master=outdir_master,
    snapshot_subdirectory=snapshot_subdirectory,
    theta=theta, phi=phi,show_gasstarxray = 'star')

