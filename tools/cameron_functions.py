from samson_const import *
import numpy as np

def calnH(Gmas,Gnh,KernalLengths,density,Gmetal):
    #Inputs = Mass, NeutralHydrogenAbundance, Kernal Length, density, metallicity
    #Unit: Mass (1e10 Msun), .., kpc, 1e10Musn/kpc^3, ...
    #Convert to cgs
    Gmas = Gmas*1e10*Msun_in_g
    KernalLengths = KernalLengths*kpc_in_cm
    density = density*1e10*Msun_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm
    Z = Gmetal[:,0] #metal mass fraction (everything not H, He)
    ZHe = Gmetal[:,1] # He mass fraction
    M_H = proton_mass_in_g
    #M_H = 1.67353269159582103*10**(-27)*(1000.0/unit_M) ##kg to g to unit mass
    #mu_H = 2.3*10**(-27)*(1000.0/unit_M) #'average' mass of hydrogen nucleus from Krumholz&Gnedin 2011
    mu_H = mu_H_in_g
    #sigma_d21 = 1 #dust cross section per H nucleus to 1000 A radiation normalized to MW val
    #R_16 = 1 #rate coefficient for H2 formation on dust grains, normalized to MW val. Dividing these two cancels out metallicity dependence, so I'm keeping as 1

    #chi = 71*(sigma_d21/R_16)*G_o/Gnh #Scaled Radiation Field


    Z_MW = 0.02 #Assuming MW metallicity is ~solar on average (Check!!)
    Z_solar = 0.02 #From Gizmo documentation

    sobColDens =np.multiply(KernalLengths,density) #Cheesy approximation of Column density
    sigma_d = 1e-21*Z/Z_MW #cross section in cm^2
    tau = np.multiply(sobColDens,sigma_d)*1/(mu_H)

    chi = 3.1 * (1+3.1*np.power(Z/Z_solar,0.365)) / 4.1 #Approximation

    s = np.divide( np.log(1.0+0.6*chi+0.01*np.square(chi)) , (0.6 *tau) )
    fH2 = np.divide((1.0 - 0.5*s) , (1+0.25*s)) #Fraction of Molecular Hydrogen (over total complex mass) from Krumholz & Knedin
    fH2[fH2<0] = 0 #Nonphysical negative molecular fractions set to 0
    fHI = 1.0-fH2  #Gnh doesn't account for He and metals
    #nHI =  np.multiply(Gmas,fHI) / M_H #Number of HI atoms in particles
    mHI =  Gmas*fHI*(1.0-(Z+ZHe))*Gnh #Mass of HI atoms in particles
    mH2 =  Gmas*fH2*(1.0-(Z+ZHe))*Gnh #Mass of H2 in particles
    mHII = (1.0-(Z+ZHe))*Gmas*(1.0-Gnh)
    return {'NHI':mHI/M_H, 'NH2':mH2/M_H, 'NHII':mHII/M_H,
        'mHII_in_g':mHII, 'mH2_in_g':mH2, 'mHI_in_g':mHI,'fH2':fH2,'fHI':fHI}
