################################################################################
###### This is an example script to generate HDF5-format ICs for GIZMO
######  The specific example below is obviously arbitrary, but could be generalized
######  to whatever IC you need. 
################################################################################
################################################################################

## load libraries we will use 
import numpy as np
import h5py as h5py
from samson_functions import *
# the main routine. this specific example builds an N-dimensional box of gas plus 
#   a collisionless particle species, with a specified mass ratio. the initial 
#   gas particles are distributed in a uniform lattice; the initial collisionless 
#   particles laid down randomly according to a uniform probability distribution 
#   with a specified random velocity dispersion
#
def makeIC_gaussian():
    datatype = 'float64'
    DIMS=3; # number of dimensions 
    N_1D=257; # 1D particle number (so total particle number is N_1D^DIMS)
    fname='/home/tkc004/oasis/testICs/box_CR_Gauss_offset0_15.hdf5'; # output filename 

    Lbox = 10.0 # box side length #kpc
    rho_desired = 0.00245 # box average initial gas density #1e10Msun/kpc^3 #6.805e-22g/cm^3 #407 cm^-3
    P_desired = 1.0 # initial gas pressure
    gamma_eos = 5./3. # polytropic index of ideal equation of state the run will assume
    CRE = 1. # total CR energy #1e10Msun*(km/s)^2 #2e53 erg
    eps = 0.15 # small offset in initial Gaussian packet # in code unit #kpc
    kappa =  10#diffusion coefficient #in code unit #kpc^2/Gyr
    timecr = 0 #by default
    # first we set up the gas properties (particle type 0)
    
    # make a regular 1D grid for particle locations (with N_1D elements and unit length)
    x0=np.arange(-0.5,0.5,1./N_1D); x0+=0.5*(0.5-x0[-1]);
    # now extend that to a full lattice in DIMS dimensions
    if(DIMS==3):
        xv_g, yv_g, zv_g = np.meshgrid(x0,x0,x0, sparse=False, indexing='xy')
    if(DIMS==2):
        xv_g, yv_g = np.meshgrid(x0,x0, sparse=False, indexing='xy'); zv_g = 0.0*xv_g
    if(DIMS==1):
        xv_g=x0; yv_g = 0.0*xv_g; zv_g = 0.0*xv_g; 
    # the gas particle number is the lattice size: this should be the gas particle number
    Ngas = xv_g.size
    # flatten the vectors (since our ICs should be in vector, not matrix format): just want a
    #  simple list of the x,y,z positions here. Here we multiply the desired box size in
    xv_g=xv_g.flatten()*Lbox; yv_g=yv_g.flatten()*Lbox; zv_g=zv_g.flatten()*Lbox; 
    # set the initial velocity in x/y/z directions (here zero)
    vx_g=0.*xv_g; vy_g=0.*xv_g; vz_g=0.*xv_g;
    # set the initial magnetic field in x/y/z directions (here zero). 
    #  these can be overridden (if constant field values are desired) by BiniX/Y/Z in the parameterfile
    bx_g=0.*xv_g; by_g=0.*xv_g; bz_g=0.*xv_g;
    # set the particle masses. Here we set it to be a list the same length, with all the same mass
    #   since their space-density is uniform this gives a uniform density, at the desired value
    mv_g=rho_desired/((1.*Ngas)/(Lbox*Lbox*Lbox)) + 0.*xv_g
    # set the initial internal energy per unit mass. recall gizmo uses this as the initial 'temperature' variable
    #  this can be overridden with the InitGasTemp variable (which takes an actual temperature)
    uv_g=P_desired/((gamma_eos-1.)*rho_desired)
    rho_g=rho_desired+0.*xv_g
    rho_g=np.array(rho_g)
    hsml_g=np.power(mv_g/rho_g,1./3.)
    # set the gas IDs: here a simple integer list
    id_g=np.arange(1,Ngas+1)
    cre_g=0.*xv_g
    rg = np.sqrt(xv_g*xv_g+yv_g*yv_g+zv_g*zv_g)
    cre_g=diffusionsol(CRE,kappa,timecr,rg,offset=eps)/((1.*Ngas)/(Lbox*Lbox*Lbox))
    cre_g=np.array(cre_g)
    

    # now we get ready to actually write this out
    #  first - open the hdf5 ics file, with the desired filename
    file = h5py.File(fname,'w') 

    # set particle number of each type into the 'npart' vector
    #  NOTE: this MUST MATCH the actual particle numbers assigned to each type, i.e.
    #   npart = np.array([number_of_PartType0_particles,number_of_PartType1_particles,number_of_PartType2_particles,
    #                     number_of_PartType3_particles,number_of_PartType4_particles,number_of_PartType5_particles])
    #   or else the code simply cannot read the IC file correctly!
    #
    npart = np.array([Ngas,0,0,0,0,0]) # we have gas and particles we will set for type 3 here, zero for all others

    # now we make the Header - the formatting here is peculiar, for historical (GADGET-compatibility) reasons
    h = file.create_group("Header");
    # here we set all the basic numbers that go into the header
    # (most of these will be written over anyways if it's an IC file; the only thing we actually *need* to be 'correct' is "npart")
    h.attrs['NumPart_ThisFile'] = npart; # npart set as above - this in general should be the same as NumPart_Total, it only differs 
                                         #  if we make a multi-part IC file. with this simple script, we aren't equipped to do that.
    h.attrs['NumPart_Total'] = npart; # npart set as above
    h.attrs['NumPart_Total_HighWord'] = 0*npart; # this will be set automatically in-code (for GIZMO, at least)
    h.attrs['MassTable'] = np.zeros(6); # these can be set if all particles will have constant masses for the entire run. however since 
                                        # we set masses explicitly by-particle this should be zero. that is more flexible anyways, as it 
                                        # allows for physics which can change particle masses 
    ## all of the parameters below will be overwritten by whatever is set in the run-time parameterfile if
    ##   this file is read in as an IC file, so their values are irrelevant. they are only important if you treat this as a snapshot
    ##   for restarting. Which you shouldn't - it requires many more fields be set. But we still need to set some values for the code to read
    h.attrs['Time'] = 0.0;  # initial time
    h.attrs['Redshift'] = 0.0; # initial redshift
    h.attrs['BoxSize'] = Lbox; # box size
    h.attrs['NumFilesPerSnapshot'] = 1; # number of files for multi-part snapshots
    h.attrs['Omega0'] = 1.0; # z=0 Omega_matter
    h.attrs['OmegaLambda'] = 0.0; # z=0 Omega_Lambda
    h.attrs['HubbleParam'] = 1.0; # z=0 hubble parameter (small 'h'=H/100 km/s/Mpc)
    h.attrs['Flag_Sfr'] = 0; # flag indicating whether star formation is on or off
    h.attrs['Flag_Cooling'] = 0; # flag indicating whether cooling is on or off
    h.attrs['Flag_StellarAge'] = 0; # flag indicating whether stellar ages are to be saved
    h.attrs['Flag_Metals'] = 0; # flag indicating whether metallicity are to be saved
    h.attrs['Flag_Feedback'] = 0; # flag indicating whether some parts of springel-hernquist model are active
    h.attrs['Flag_DoublePrecision'] = 0; # flag indicating whether ICs are in single/double precision
    h.attrs['Flag_IC_Info'] = 0; # flag indicating extra options for ICs
    ## ok, that ends the block of 'useless' parameters
    
    # Now, the actual data!
    #   These blocks should all be written in the order of their particle type (0,1,2,3,4,5)
    #   If there are no particles of a given type, nothing is needed (no block at all)
    #   PartType0 is 'special' as gas. All other PartTypes take the same, more limited set of information in their ICs
    
    # start with particle type zero. first (assuming we have any gas particles) create the group 
    p = file.create_group("PartType0")
    # now combine the xyz positions into a matrix with the correct format
    q=np.zeros((Ngas,3)); q[:,0]=xv_g; q[:,1]=yv_g; q[:,2]=zv_g;
    # write it to the 'Coordinates' block
    p.create_dataset("Coordinates",data=q,dtype=datatype)
    # similarly, combine the xyz velocities into a matrix with the correct format
    q=np.zeros((Ngas,3)); q[:,0]=vx_g; q[:,1]=vy_g; q[:,2]=vz_g;
    # write it to the 'Velocities' block
    p.create_dataset("Velocities",data=q,dtype=datatype)
    # write particle ids to the ParticleIDs block
    p.create_dataset("ParticleIDs",data=id_g)
    # write particle masses to the Masses block
    p.create_dataset("Masses",data=mv_g)
    # write internal energies to the InternalEnergy block
    p.create_dataset("InternalEnergy",data=uv_g)
    # combine the xyz magnetic fields into a matrix with the correct format
    p.create_dataset("Density",data=rho_g, dtype=datatype)
    p.create_dataset("SmoothingLength",data=hsml_g, dtype=datatype)
    p.create_dataset("CosmicRayEnergy",data=cre_g, dtype=datatype)
    q=np.zeros((Ngas,3)); q[:,0]=bx_g; q[:,1]=by_g; q[:,2]=bz_g;
    # write magnetic fields to the MagneticField block. note that this is unnecessary if the code is compiled with 
    #   MAGNETIC off. however, it is not a problem to have the field there, even if MAGNETIC is off, so you can 
    #   always include it with some dummy values and then use the IC for either case
    p.create_dataset("MagneticField",data=q)

    # no PartType1 for this IC
    # no PartType2 for this IC


    # no PartType4 for this IC
    # no PartType5 for this IC

    # close the HDF5 file, which saves these outputs
    file.close()
    # all done!
makeIC_gaussian()
