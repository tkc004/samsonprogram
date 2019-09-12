import collections
from importplotfunctionsCR import *
from inputplotdataCR import inputplotdict
plotlist = ['figCRgasmid_m12mhr', 'figCRgasmid_m12bhr']
#plotlist = ['figCRgasmid_m12bhr', 'figCRgasmid_m12chr']
#plotlist = ['figBrms_m12bhr','figBrms_m12chr']
for plotneed in plotlist:
    inputplotdict[plotneed]['_plotfunction'](inputplotdict[plotneed])
