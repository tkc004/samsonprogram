import collections
from importplotfunctionsCR import *
from inputplotdataCR import inputplotdict
plotlist = ['figBrms_m12fhr','figBrms_m12mhr']
#plotlist = ['figBrms_m12bhr','figBrms_m12chr']
for plotneed in plotlist:
    inputplotdict[plotneed]['_plotfunction'](inputplotdict[plotneed])
