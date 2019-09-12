import collections
from importplotfunctionsCR import *
from inputplotdataCR import inputplotdict
plotlist = ['figcrden_m12mhrlog']
#plotlist = ['figsurden_m12fhr2']
for plotneed in plotlist:
    inputplotdict[plotneed]['_plotfunction'](inputplotdict[plotneed])
