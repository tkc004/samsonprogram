import collections
from inputplotdataISM import inputplotdict
plotlist=['figintrhogfrompardist_m12hr0_25kpcr0_9_usekez']
#plotlist=['figrhoTgrid_m12iCGM']
#plotlist=['figanHz_m12i']
#plotlist=['figvz_m12hr0_25kpc_usekezHI']
#plotlist=['outssm12mcr_700580']
#plotlist=['gpm12icr_700']
#plotlist=['figintrhogfrompardist_m12fhrtest']
#plotlist=['figintrhogfrompardist_m12fhr10kpc']
for plotneed in plotlist:
    inputplotdict[plotneed]['_plotfunction'](inputplotdict[plotneed])
