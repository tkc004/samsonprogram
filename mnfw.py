import numpy as np
reffp=1

zr=0
rhocrit=127*(0.27*(1+zr)*(1+zr)*(1+zr)+0.73)
xLam=-0.73/(0.27*(1+zr)*(1+zr)*(1+zr)+0.73)
delta=18.0*np.pi*np.pi+82.0*(xLam)-39*(xLam)*(xLam)
samlogMv=np.arange(8.0,12.0,0.05)
logc=0.537+(1.025-0.537)*np.exp(-0.718*np.power(zr,1.08))+(0.024*zr-0.097)*(samlogMv-np.log10(0.7)-np.log10(1e12))
cont=np.power(10,logc)
Rvirn=np.power(np.power(10,samlogMv)*3.0/4.0/np.pi/rhocrit/delta,1.0/3.0)

#finding the NFW profile (at z=0) that fits rho at reff and plot it
#x=reffp/Rvirn
#rho = rhocrit*delta/cont/cont/cont/x/(1/cont+x)/(1/cont+x)
Rs = Rvirn/cont
cdelta = cont*cont*cont*delta/3./(np.log(1+cont)-cont/(1+cont))
menc = 4*np.pi*rhocrit*cdelta*Rs*Rs*Rs*(np.log((Rs+reffp)/Rs)-reffp/(Rs+reffp))
