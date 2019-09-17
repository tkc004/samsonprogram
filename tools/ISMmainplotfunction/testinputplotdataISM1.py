import collections
from inputplotdataISM import inputplotdict
#plotlist = ['figsfrmocktur_m12imhdcv',
#                'figsfrmocktur_m12i']
#]
plotlist=['figpcrpthgrid_m12']
#plotlist=['figrhoTgrid_m12']
#plotlist=['figTztrackstart_m12i_warmssout']
#plotlist=['figTztrackstart_m12i_580']
#plotlist=['figintrhogfrompardist_m12hr0_25kpcr0_9_usekez']
#plotlist = ['figsfrmocktur_m11gcr_700']
#plotlist=['figsfrmocktur_m11bcr_700']
#plotlist = ['figsfrmocktur_m12i']
#plotlist=['outssm12f']
#plotlist=['figrhogfrompardist_m12ihr0kpc']
#plotlist=['figsfrtur_m11_m12']
#plotlist=['figCRpr_volave_m12f']
#plotlist = ['figCRpr_volave_m12i']
for plotneed in plotlist:
    inputplotdict[plotneed]['_plotfunction'](inputplotdict[plotneed])
