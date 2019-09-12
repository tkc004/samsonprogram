import collections
from importplotfunctionsCR import *
from inputplotdataCR import inputplotdict
plotlist = ['figBrms_m12ihr2']
for plotneed in plotlist:
    inputplotdict[plotneed]['_plotfunction'](inputplotdict[plotneed])
