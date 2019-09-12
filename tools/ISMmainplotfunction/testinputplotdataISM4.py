import collections
from inputplotdataISM import inputplotdict
plotlist=['figgasdenTz_m12']
#plotlist=['figanHz_m12m']
#plotlist=['figsfr_m12']
#plotlist = ['outssm12mcr_700590']
#plotlist=['figrhoTgrid_m12out']
#plotlist=['figintrhogfrompardist_m12ihrtest']
#plotlist=['figrhoTgrid_m12out']
#plotlist=['figintrhogfrompardist_m12imhdcvhrtest']
#plotlist=['figintrhogfrompardist_m12mhr10kpcz10kpc']
#plotlist=['figintrhogfrompardist_m12fhr10kpc']
#plotlist=['figrhoTgrid_m12']
#plotlist=['figCRpr_volave_m11f','figCRpr_volave_m11g']
#plotlist=['figCRvtur_single_m12']
#plotlist=['figrhogfrompardist_m12mhr0kpc']
#plotlist=['figsfrtur_m11_m12']
#plotlist=['figsfrtur_m12i']
#plotlist = ['figetur_m12mhr']
#plotlist = ['figmgasz_m12f']
#plotlist = ['figSFRsur_m12ihr']
#plotlist = ['figSFRsur_m12ihr']
for plotneed in plotlist:
    inputplotdict[plotneed]['_plotfunction'](inputplotdict[plotneed])
