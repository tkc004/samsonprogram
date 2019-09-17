import collections
from inputplotdataISM import inputplotdict
plotlist=['figagfrompardist_m12ihr10kpc']
#plotlist=['figvzzgrid_m12_z200_cr100']
#plotlist=['figvzzgrid_m12_r200']
#plotlist=['figTvznologgrid_m12_z50_100kpc_510']
#plotlist=['figTvznologgrid_m12_z50_100kpc']
#plotlist=['figvzzgrid_m12_z100']
#plotlist=['figTztrackstart_m12i']
#plotlist=['figTztrack_m12i_inflow']
#plotlist=['figintrhogfrompardist_m12hr0_25kpcr0_9_usekez']
#plotlist=['figrhoTgrid_m12iCGM']
#plotlist=['figanHz_m12i']
#plotlist=['figvz_m12hr0_25kpc_usekezHI']
#plotlist=['outssm12mcr_700580']
#plotlist=['gpm12icr_700']
#plotlist=['figintrhogfrompardist_m12fhrtest']
#plotlist=['figintrhogfrompardist_m12fhr10kpc']
for plotneed in plotlist:
    inputplotdict[plotneed]['_plotfunction'](inputplotdict[plotneed])
