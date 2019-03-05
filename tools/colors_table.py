def colors_table( age_in_Gyr, metallicity_in_solar_units, 
    BAND_ID=0, SALPETER_IMF=0, CHABRIER_IMF=0, NEW_CHABRIER_IMF=0, KROUPA_IMF=1, QUIET=0, CRUDE=0, 
    RETURN_NU_EFF=0, RETURN_LAMBDA_EFF=0):

    import utilities as util
    import numpy as np
    import math
    import scipy.ndimage.interpolation as interpolate
    import struct


#    print 'age_in_Gyr', age_in_Gyr
#    print 'metallicity_in_solar_units', metallicity_in_solar_units

    age_in_Gyr=np.array(age_in_Gyr,ndmin=1);
    metallicity_in_solar_units=np.array(metallicity_in_solar_units,ndmin=1);

    band=BAND_ID; # default=bolometric
    j = [  0,  6,  7,  8,  9, 10, 11, 12, 13,  1,   2,   3,   4,   5] # ordering I'm used to
    i = [  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10,  11,  12,  13] # ordering of this
    band_standardordering = band
    band = j[band]
    if (band > 13): 
        print 'BAND_ID must be < 13'; 
        return 0;
    
    b=['Bolometric', \
    'Sloan u','Sloan g','Sloan r','Sloan i','Sloan z', \
    'Johnsons U','Johnsons B', 'Johnsons V','Johnsons R','Johnsons I', \
    'Cousins J','Cousins H','Cousins K']
    if (QUIET==0): print 'Calculating M/L in '+str(b[band])+' ('+str(band)+','+str(band_standardordering)+')'

    print 'Calculating M/L in '+str(b[band])+' ('+str(band)+','+str(band_standardordering)+')'   
 
    if (RETURN_NU_EFF==1) or (RETURN_LAMBDA_EFF==1):
        lam_eff=np.array([1.e-5, 3541., 4653., 6147., 7461., 8904., 3600., 4400., \
   		    5556., 6940., 8700., 12150., 16540., 21790.]);
        nu_eff = 2.998e18 / lam_eff;
        if (RETURN_NU_EFF==1): return nu_eff[band];
        if (RETURN_LAMBDA_EFF==1): return lam_eff[band];
    
    froot = 'tools/'; # directory in which the data binaries are stored
 #   print 'KROUPA_IMF', KROUPA_IMF
    if (CHABRIER_IMF==1): fname=froot+'colors.chabrier.dat'
    if (SALPETER_IMF==1): fname=froot+'colors.salpeter.dat'
    if (NEW_CHABRIER_IMF==1): fname=froot+'colors.chabrier_fsps.dat'
    if (KROUPA_IMF==1): fname=froot+'colors.kroupa_fsps.dat'


    # be careful! use 'r' but not 'rb' for the new EOS table 
    lut = open(fname,'r');
    lut_dat = lut.read();
    Nl,Na,Nz = struct.unpack('3i',lut_dat[0:12])
    z_grid = np.array(struct.unpack(str(Nz)+'d',lut_dat[12:12+8*Nz]))
    age_grid = np.array(struct.unpack(str(Na)+'d',lut_dat[12+8*Nz:12+8*Nz+8*Na]))
    l_all_l = np.array(struct.unpack(str(Nl*Na*Nz)+'d',lut_dat[12+8*Nz+8*Na:12+8*Nz+8*Na+8*Nl*Na*Nz]))
    l_all = np.transpose(l_all_l.reshape(Nz,Na,Nl))
    lut.close()
    
    l_band = np.zeros((Na,Nz),dtype=np.float64);
    for iz in range(Nz): l_band[:,iz]=l_all[band,:,iz]
    
    # allow for extreme metallicities (extrapolate linearly past table)
    push_metals = 1;
    if (push_metals==1):
        Nz = Nz + 1;
        z_ext = [1000.0];
        z_grid = np.concatenate([z_grid,z_ext])
        lb1 = l_band[:,Nz-3]
        lb2 = l_band[:,Nz-2]
        lbx = np.zeros((Na,Nz),dtype=np.float64)
        lbx[:,0:Nz-1] = l_band
        lbx[:,Nz-1] = (lb2 - lb1) / (np.log10(z_grid[Nz-2]/z_grid[Nz-3])) * \
            np.log10(z_grid[Nz-1]/z_grid[Nz-2])
        l_band = lbx;

    # get the x-axis (age) locations of input points
    ia_pts=np.interp(np.log10(age_in_Gyr)+9.0,age_grid,np.arange(0,Na,1)); 
    # this returns the boundary values for points outside of them (no extrapolation)
    #f=interp.interp1d(age_grid,np.arange(0,Na,1),kind='linear'); 
    #ia_pts=f(np.log10(age_in_Gyr)+9.0);

    # get the y-axis (metallicity) locations of input points
    zsun = 0.019;
    iz_pts=np.interp(np.log10(metallicity_in_solar_units*zsun),np.log10(z_grid),np.arange(0,Nz,1)); 
    #f=interp.interp1d(np.log10(z_grid),np.arange(0,Nz,1),kind='linear'); 
    #iz_pts=f(np.log10(metallicity_in_solar_units*zsun));

    #all new FSPS tables are in unit log10(L(cgs))-log10(Mstar(Msun))

  #  print 'l_band', l_band

    if (CRUDE==1):
        ia_pts=np.around(ia_pts).astype(int);
        iz_pts=np.around(iz_pts).astype(int);
        l_b=l_band[ia_pts,iz_pts];
    else:
        l_b = interpolate.map_coordinates(l_band, (ia_pts,iz_pts), order=1);
    l_b = 10.**l_b
    	
    return l_b;

