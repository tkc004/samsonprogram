import collections
from inputplotdataISM import inputplotdict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--plotitem', default='gf12imhdcv')
args = parser.parse_args()

print args.plotitem
plotlist=[args.plotitem]
for plotneed in plotlist:
    inputplotdict[plotneed]['_plotfunction'](inputplotdict[plotneed])
