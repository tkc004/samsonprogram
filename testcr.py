from samson_const import *
import matplotlib as mpl
mpl.use('Agg')
from readsnap_cr import readsnapcr
import Sasha_functions as SF
import graphics_library as GL
import numpy as np
import matplotlib.pyplot as plt
from samson_functions import *
from matplotlib import rcParams
from pylab import *
from Sasha_functions import *
#rcParams['figure.figsize'] = 5, 5
rcParams['figure.figsize'] = 10, 5
rcParams['font.size']=12
rcParams['font.family']='serif'
#rcParams.update({'figure.autolayout': True})
import matplotlib.patches as patches
rcParams['axes.linewidth'] = 2
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True
rcParams['axes.unicode_minus']=False
colortable = [ 'b', 'g', 'r']
factor=1.0/proton_mass_in_g*1e10*solar_mass_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm
factor = factor/pidecay_fac
factor = factor*1e10*solar_mass_in_g*km_in_cm*km_in_cm*betapi/nopi_per_gamma
print 'factor',factor

factorSF=Kroupa_Lsf*solar_mass_in_g/yr_in_sec*cspeed_in_cm_s*cspeed_in_cm_s
print 'factorSF', factorSF
data = outdirname('fm11qmd')
print 'data[rundir]',data['rundir'], data['maindir']
