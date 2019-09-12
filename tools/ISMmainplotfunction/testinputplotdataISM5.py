import collections
from inputplotdataISM import inputplotdict
plotlist=['figgasdenTz_m12_atr8kpc']
#plotlist=['figkslaw_m11_m12_out']
#plotlist=['figkslaw_m11_m12']
#plotlist=['figkslaw_m12mhr']
#plotlist=['figCRpr_volave_m11f','figCRpr_volave_m11g']
#plotlist=['figCRpcr_single_m12']
#plotlist=['figrhogfrompardist_m12mhr0kpc']
#plotlist=['figsfrtur_m11_m12']
#plotlist=['figsfrtur_m12i']
#plotlist = ['figetur_m12mhr']
#plotlist = ['figmgasz_m12f']
#plotlist = ['figSFRsur_m12ihr']
#plotlist = ['figSFRsur_m12ihr']
for plotneed in plotlist:
    inputplotdict[plotneed]['_plotfunction'](inputplotdict[plotneed])
