import collections
from importplotfunctionsCR import *
from inputplotdataCR import inputplotdict
plotlist = ['figgasden_m12fhr2']
for plotneed in plotlist:
    inputplotdict[plotneed]['_plotfunction'](inputplotdict[plotneed])
