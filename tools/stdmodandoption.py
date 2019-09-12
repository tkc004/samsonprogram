import os
import numpy as np
import matplotlib
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm
import Sasha_functions as SF
import graphics_library as GL
import gas_temperature as GT
from readsnap_cr import readsnapcr
from pylab import *
import matplotlib.colors
import matplotlib.cm
import math
import h5py
import re
import sys
import glob
from textwrap import wrap
from numpy.linalg import inv
rcParams.update({'figure.autolayout': True})
rcParams.update({'font.size': 24})
rcParams['axes.unicode_minus'] = False
import matplotlib.patches as patches
rcParams['axes.linewidth'] = 2
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True
rcParams['xtick.major.width'] = 2
rcParams['ytick.major.width'] = 2
rcParams['xtick.minor.width'] = 2
rcParams['ytick.minor.width'] = 2
rcParams['text.usetex'] = True
colortable = [ 'b', 'g', 'r']
import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='inferno', cmap=cmaps.inferno)
plt.register_cmap(name='plasma', cmap=cmaps.plasma)
from gadget_lib.cosmo import *
import samson_functions as SSF
import crtestfunction as CRTF
from samson_const import *
from pathloc import *
import readsnipshot as RSS