betapi = 0.7 ##should be between 0.5-0.9 from Lacki 2011 #fraction of pi has enough energy (>1GeV)
au_in_cm = 1.496e+13
degree_in_arcsec = 3600. 
yr_in_sec = 3.2e7
Myr_in_sec = yr_in_sec*1.0e6
Gyr_in_sec = yr_in_sec*1.0e9
nopi_per_gamma = 3.0
solar_mass_in_g = 1.989e33
Msun_in_g = solar_mass_in_g
km_in_cm = 1e5
m_in_cm = 1e2
kg_in_g = 1e3
kpc_in_pc = 1e3
kpc_in_cm = 3.085678e21
kpc_in_m = kpc_in_cm/m_in_cm
pc_in_cm = kpc_in_cm/1000.
Mpc_in_cm = kpc_in_cm*1000.
erg_in_eV = 6.242e11
erg_in_GeV = erg_in_eV/1.0e9
eV_in_erg = 1.0/erg_in_eV
cspeed_in_cm_s = 3e10
cspeed_in_cgs = cspeed_in_cm_s
cspeed_in_m_s = cspeed_in_cm_s/m_in_cm
cspeed_in_km_s = cspeed_in_cm_s/km_in_cm
proton_mass_in_g = 1.6726e-24
protonmass_in_g = proton_mass_in_g
electronmass_in_g = 9.1095e-28
electronmass_in_kg = electronmass_in_g/kg_in_g
Kroupa_Lsf=6.2e-4
Salpeter_Lsf=3.8e-4
pidecay_fac = 2.0e5*250.0*yr_in_sec #over neff in cm_3 for pi decay time in s
hadronicdecayrate = 5.86e-16 #from Guo 2008. The hadronic decay rate of CR
coulombdecayrate = 1.65e-16 #Coulomb loss of CR
NewtonG_cgs = 6.67259e-8
NewtonG_in_cgs = NewtonG_cgs
G_in_cgs = NewtonG_cgs
H_MASSFRAC = 0.76
PROTONMASS  = proton_mass_in_g
GAMMA = 5.0/3.0
BOLTZMANN  = 1.3806e-16
kB = BOLTZMANN
Avogadro_const = 6.022140857e23
CRgamma = 4.0/3.0
rho_codetocgs=1e10*Msun_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm
m_codetocgs=1e10*Msun_in_g
u_codetocgs=km_in_cm*km_in_cm
x_codetocgs=kpc_in_cm
area_codetocgs=x_codetocgs*x_codetocgs
vol_codetocgs=x_codetocgs*x_codetocgs*x_codetocgs
cregy_codetocgs=1e10*Msun_in_g*km_in_cm*km_in_cm
t_codetocgs=0.98*1e9*yr_in_sec
echarge_in_statC = 4.80320425e-10
echarge_in_C = 1.60e-19
mu_H_in_g = 2.3e-24 #'average' mass of hydrogen nucleus from Krumholz&Gnedin 2011
m_codetoMsun = 1e10
Myr_in_yr = 1e6
Gyr_in_yr = 1e9
kpc_in_pc = 1e3
Gyr_in_Myr = 1e3
Tesla_in_Gauss = 1e4
degree_in_turn = 360.
arcmin_in_degree = 60.
arcsec_in_degree = 3600.
microarcsec_in_degree = arcsec_in_degree*1e6
microarcsec_in_turn = microarcsec_in_degree*degree_in_turn
arcsec_in_turn = arcsec_in_degree*degree_in_turn
Joule_in_eV = 6.24152272506e18
Joule_in_GeV = 6.24152272506e18/1.0e9
Jy_in_SI = 1e-26
Jy_in_cgs = 1e-23
muJy_in_SI = 1e-32
muJy_in_cgs = 1e-29
npionmass_in_g = 2.406e-25 
propconstfromLackitoFermi = 1.369e-16 #in cgs: multiply by nuclei density and CR energy density to get luminosity per volume
Lackiprop = 1.35e-16
