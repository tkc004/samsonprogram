import collections
from inputplotdataISM import inputplotdict
plotlist=[
'figTztrack_m12i_580',
#'figanHz_m12f',
#'figgasdenTz_m12f_atr8kpc',
#'figgasdenTz_m12',
#'figrhoTvoloutz_m12'
#'figTvznologgrid_m12'
#'figgasdenTr_m12'
#'figgasdenTz_m12',
#'figvz_m12hr0_25kpc_usekezHI'
#'figsfrmockturtot_m11_m12'
#'gpm12imhdcvgrid0_125kpc6580HI',\
#'figtdsvz_m12hr0_5kpc_usekez',\
#'figtdspz_m12hr0_5kpc_usekez',\
]
#plotlist=['figintrhogfrompardist_m12hr0_25kpcr0_9_usekez']
#plotlist=['figvr_m12hr0_5kpc_usekez']
#plotlist=['figpr_m12hr0_5kpc_usekez']
#plotlist=['figintrhogfrompardist_m12hr0_25kpcr0_9_usekez']
#plotlist=['outssm12fcr_700590']
#plotlist=['outssm12f590']
#plotlist=['gf12imhdcv']
#plotlist=['outssm12m']
#plotlist=['figintrhogfrompardist_m12ihr10kpc']
#plotlist=['figrhoTgrid_m12out']
for plotneed in plotlist:
    inputplotdict[plotneed]['_plotfunction'](inputplotdict[plotneed])
