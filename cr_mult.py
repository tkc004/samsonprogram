from samson_const import *
import matplotlib as mpl
mpl.use('Agg')
from readsnap_cr import readsnapcr
import Sasha_functions as SF
import graphics_library as GL
import gas_temperature as GT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from samson_functions import *
from matplotlib import rcParams
from pylab import *
from textwrap import wrap
#rcParams['figure.figsize'] = 5, 5
rcParams['figure.figsize'] = 10, 5
rcParams['font.size']=18
rcParams['font.family']='serif'
rcParams['text.usetex']=True
#rcParams.update({'figure.autolayout': True})
import matplotlib.patches as patches
rcParams['axes.linewidth'] = 2
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True
rcParams['axes.unicode_minus']=False
colortable = [ 'b', 'g', 'r']
#dirneed=['bwsbclrdc28']
#dirneed=['m12imhdcv','m12icr_b_70','m12icr_700','m12icr_70va']
#dirneed=['m12mmhdcv','m12mcr_b_70','m12mcr_700','m12mcr_70va']
#dirneed=['m11dmhdcv','m11dcr_b_70','m11dcr_700','m11dcr_70va']
#dirneed = ['bwmwlrstr','bwmwlrstrc1000']
#dirneed=['bwmwlrdc27ds','bwmwlrdc27']
#dirneed=['mw_mr_2_21','bwmwmrdc0','bwmwmrdc28']
#dirneed=['m09mhdcv','m09cr_b_70','m10qmhdcv','m10qcr_b_70','m11bmhdcv','m11bcr_b_70','m11fmhdcv','m11fcr_b_70','m12imhdcv','m12icr_b_70']
#dirneed=['m10vmhdcv','m10vcr_b_70', 'm10vcr_700']
#dirneed=['m11gcr_b_70','m11gcr_700','m11gmhdcv','f46']
#dirneed=['m12icr_b_70','m12icr_700','m12imhdcv','fm12i']
#dirneed=['m12icr_b_70']
#dirneed=['m12mcr_b_70']
#dirneed=['m11bcr_700']
#dirneed=['bwsmclrdc28']
#dirneed=['bwmwmrdc28']
dirneed=['m12mcr_700']
#dirneed=['m11bcr_b_70','m11bcr_700','m11bmhdcv','f476']
#dirneed=['bwsmclrdc29']
#dirneed=['bwmwmr']
#dirneed=['fm12m']
#dirneed=['m10vmhdcv','m11bmhdcv','m11dmhdcv','m11fmhdcv','m11gmhdcv','m11hmhdcv','m12imhdcv',\
#'m12mmhdcv','m12fmhdcv']
#dirneed=['m10vcr_b_70','m11bcr_b_70','m11dcr_b_70','m11fcr_b_70','m11gcr_b_70','m11hcr_b_70','m12icr_b_70',\
#'m12mcr_b_70','m12fcr_b_70']
#dirneed=['m10vcr_700','m11bcr_700','m11dcr_700','m11fcr_700','m11gcr_700','m11hcr_700','m12icr_700',\
#'m12mcr_700','m12fcr_700']
#dirneed=['m10vcr_b_70','m11bcr_b_70','m11dcr_b_70','m11fcr_b_70','m11gcr_b_70','m11hcr_b_70','m12icr_b_70']
#dirneed=['m10vmhdcv','m11bmhdcv','m11dmhdcv','m11fmhdcv','m11gmhdcv','m11hmhdcv','m12imhdcv']
#dirneed=['m10qmhdcv']
#dirneed=['m10qcr_b_70']
#dirneed=['m10vmhdcv']
#dirneed=['bwmwlrstr','bwmwlrstrll','mwlrstrllva','mwlrstrll4va']
#dirneed=['bwsmclrmhd']
#dirneed=['bwsmclrmhd']
#dirneed=['bwsmclrstr','bwsmclrstrll','smclrstrllva','smclrstrll4va']
#dirneed=['bwmwlrdc27','bwmwlrdc27nc']
#dirneed=['bwmwlrdc27ehh']
#dirneed=['dc27g','M1dc27g_2_21_c100_smallts', 'M1dc27g_2_21_c1000_smallts']
#dirneed=['bwsmclrdc27str']
#dirneed=['m12icr_b_70','m12imhdcv']
#dirneed=['m10qcr_b_70','m10qmhdcv']
#dirneed=['FIRE_2_0_or_h573_criden1000_noaddm_sggs']
#dirneed=['mw_lr_2_21','mw_lr_2_21_mhd']
#dirneed=['bwmwlrmhd']
#dirneed=['mw_lr_2_21','mw_cr_lr_dc27_2_21_M1_smallts','mw_cr_lr_dc28_2_21_M1_smallts']
#dirneed=['mw_cr_lr_dc27_2_21_M1_fhfcr_ets_c100']
#dirneed=['mw_cr_lr_dc27_2_21_M1']
#dirneed=['mw_cr_lr_dc27_2_21_M1_smallts','mw_cr_lr_dc28_2_21_M1_smallts','mw_cr_lr_dc28_2_21_M1_mhd']
#dirneed=['mw_lr_2_21','mw_cr_lr_dc27_2_21_M1_smallts','mw_cr_lr_dc28_2_21_M1_smallts','mw_cr_lr_dc28_2_21_M1_mhd']
#dirneed=['mw_cr_lr_dc27_2_21_M1_smallts','mw_cr_lr_dc28_2_21_M1_smallts','mw_cr_lr_dc28_2_21_M1_mhd','mw_cr_lr_dc28_2_21_M1_mhd_stream','mw_cr_lr_2_21_purestream']
#dirneed=['smc_cr_lr_dc28_2_21_M1_smallts','smc_cr_lr_dc28_2_21_M1_c1000_smallts']
#dirneed=['bwmwlrdc27ehh', 'bwmwlrdc28ehh']
#dirneed=['smc_lr_2_21_mhd_su','smc_lr_2_21_weakB_p10G_su']
#dirneed=['smc_mr_2_21_vtd_t1_mhd']
#dirneed=['mw_cr_lr_dc0_2_21']
#dirneed=['bwmwlrdc27ds','bwmwlrdc27']
#dirneed=['m10qcr_b_70','m11bcr_b_70','m12icr_b_70']
#dirneed=['bwsmclrdc28']
#dirneed=['FIRE_2_0_h573_cr_2_21_M1_z0','FIRE_2_0_h573_2_21']
#dirneed=['FIRE_2_0_h573_cr_2_21_M1_z0']
#dirneed=['fm12f']
#dirneed=['fm12c']
#dirneed=['fm12b']
#dirneed=['m11gcr_b_70']
#dirneed=['fm12q']
#dirneed=['fm12i']
#dirneed=['f46']
#dirneed=['m11bmhdcv']
#dirneed=['f476']
#dirneed=['fm12m']
#dirneed = ['m12mcr_b_70']
#dirneed = ['m12mmhdcv']
#dirneed = ['m12mcr_700']
#dirneed = ['m12imhdcv']
#dirneed = ['m10vmhdcv', 'm11bmhdcv','m11dmhdcv','m12imhdcv','m11fmhdcv','m11hmhdcv',\
#'m11gmhdcv','m10vcr_b_70', 'm11bcr_b_70','m11dcr_b_70','m12icr_b_70','m11fcr_b_70',\
#'m11hcr_b_70','m11gcr_b_70','m10vcr_700', 'm11bcr_700','m11dcr_700',\
#'m11fcr_700','m11gcr_700','m11hcr_700','m12icr_700',\
#'fm11v','fm11d','f383','fm11q','f1146','f61']
#dirneed=['m11bcr_700','m11bcr_b_70','m11bmhdcv']
#dirneed=['m09cr_b_70', 'm09mhdcv','m11dcr_b_70','m11dmhdcv','m11hcr_b_70','m11hmhdcv','m11gmhdcv','m11gcr_b_70']
#dirneed=['m10qmhdcv','m11bmhdcv','m12imhdcv','m10qcr_b_70','m11bcr_b_70','m12icr_b_70']
#dirneed=['m10qmhdcv','m11bmhdcv','m12imhdcv']
#dirneed=['m10qmhdcv','m10qcr_b_70','m11bmhdcv','m11bcr_b_70','m12imhdcv','m12icr_b_70']
#dirneed = ['bwmwllrdc28', 'bwmwlrdc28','bwmwmrdc28']
#dirneed=['bwsmclrdc28', 'bwsmclrdc28cutcr','bwsmclrdc29', 'bwsmclrdc29cutcr']
#dirneed=['bwmwlrdc28', 'bwmwlrdc28cutcr','bwmwlrdc29', 'bwmwlrdc29cutcr']
#dirneed=['bwsbclr','bwsbclrmhd']
#dirneed=['bwsbclrdc28']
#dirneed=['bwsbclr']
#dirneed=['bwsbclrdc28mhd']
#dirneed=['bwsbclrdc29']
#dirneed=['bwsbclr','bwsbclrdc0','bwsbclrdc27','bwsbclrdc28','bwsbclrdc29']
#dirneed=['bwsbclrdc27','bwsbclrdc28','bwsbclrdc29']
#dirneed=['bwsmclrdc0','bwsmclrdc27','bwsmclrdc28','bwsmclrdc29']
#dirneed=['bwsmclr','bwsmclrmhd','bwsmclrdc0','bwsmclrdc27','bwsmclrdc28','bwsmclrdc29','bwsmclrstr','bwsmclrdc28mhd','bwsmclrdc28str']
#dirneed=['bwmwmr','bwmwmrmhd','bwmwmrdc0','bwmwmrdc27','bwmwmrdc28', 'bwmwmrdc29','bwmwmrstr','bwmwmrdc28mhd','bwmwmrdc28str']
#dirneed=['bwmwmr','bwmwmrdc0','bwmwmrdc27','bwmwmrdc28','bwmwmrdc29']
#dirneed=['bwsmclrdc0','bwsmclrdc27','bwsmclrdc28','bwsmclrdc29']
#dirneed=['bwsmclr','bwsmclrdc0','bwsmclrdc27','bwsmclrdc28','bwsmclrdc29']
#dirneed=['bwsmclrdc0','bwsmclrdc27','bwsmclrdc28','bwsmclrdc29','bwsmclrstr','bwsmclrdc28mhd','bwsmclrdc28str']
#dirneed=['bwsmclr','bwsmclrdc0','bwsmclrdc27','bwsmclrdc28','bwsmclrdc29','bwsmclrstr'] 
#dirneed=['bwmwlrmhd']
#dirneed=['bwmwlrstr', 'bwmwlrstrll']
#dirneed=['bwsmclrstr', 'bwsmclrstrll']
#dirneed=['M1strg_2_21_c500_smallts']
#dirneed=['bwsbclr']
#dirneed=['bwsbclrmhd']
#dirneed=['bwsbclrdc0']
#dirneed=['bwsbclrdc27']
#dirneed=['bwsbclrdc28']
#dirneed=['bwsbclrdc29']
#dirneed=['bwsbclrstr']
#dirneed=['bwsbclrdc28mhd']
#dirneed=['bwsbclrdc28str']
#dirneed=['bwmwmrstrll4va']
#dirneed=['bwsbclr','bwsbclrmhd','bwsbclrdc0','bwsbclrdc27', 'bwsbclrdc28','bwsbclrdc29','bwsbclrstr','bwsbclrdc28mhd','bwsbclrdc28str']
#dirneed=['sbc_lr_2_21_extend4']
#dirneed=['bwsmclrmhd','bwsmclrdc28','bwsmclrdc28str']
#dirneed=['bwsmclrdc29'] 
#dirneed=['bwmwlrdc28str']
#dirneed=['mw_cr_lr_dc28_2_21_M1_smallts','bwmwlrdc28str']
#dirneed=['mw_lr_2_21','bwmwlrdc28str']
#dirneed=['sbc_lr_2_21_extend3']
#dirneed=['smc_cr_lr_dc28_2_21_M1_smallts','smc_cr_lr_dc28_2_21_M1_smallts_cutcr']
#dirneed=['bwsmclrdc27','bwsmclrdc28'] 
#dirneed=['mw_cr_lr_dc27_2_21_M1_smallts','mw_cr_lr_dc28_2_21_M1_c1000_smallts','mw_cr_lr_dc29_2_21_M1_c2000_smallts']
#dirneed=['mw_lr_2_21','mw_cr_lr_dc0_2_21','mw_cr_lr_dc27_2_21_M1_smallts','mw_cr_lr_dc28_2_21_M1_c1000_smallts','mw_cr_lr_dc29_2_21_M1_c2000_smallts']
#dirneed=['mw_lr_2_21_mhd_IC']
#dirneed=['mw_lr_2_21_mhd','mw_cr_lr_dc28_2_21_M1_mhd','mw_cr_lr_dc28_2_21_M1_mhd_stream','mw_cr_lr_2_21_purestream']
#dirneed=['mw_cr_lr_dc28_2_21_M1_mhd']
#dirneed=['mw_lr_2_21_mhd','mw_cr_lr_dc28_2_21_M1_mhd','mw_cr_lr_dc28_2_21_M1_mhd_stream','mw_cr_lr_2_21_purestream']
#dirneed=['mw_cr_lr_dc28_2_21_M1_c1000_smallts']
#dirneed=['mw_cr_lr_dc28_2_21_M1_smallts','mw_cr_lr_dc28_2_21_M1_c1000_smallts','mw_cr_lr_dc28_2_21_M1_c2000_smallts']
#dirneed=['mw_lr_2_21','mw_lr_2_21_mhd','mw_cr_lr_dc27_2_21_M1_smallts','mw_cr_lr_dc28_2_21_M1_smallts','mw_cr_lr_dc28_2_21_M1_mhd','mw_cr_lr_dc28_2_21_M1_mhd_stream','mw_cr_lr_2_21_purestream']
#dirneed=['mw_cr_lr_dc28_2_21_M1_mhd']
#dirneed=['mw_cr_lr_dc28_2_21_M1_mhd','mw_cr_lr_dc28_2_21_M1_mhd_stream','mw_cr_lr_2_21_purestream']
#dirneed=['smc_lr_2_21_thickdisk']
#dirneed=['mw_lr_2_21_strongB_IC']
#dirneed=['mw_cr_lr_dc28_2_21_M1_mhd','mw_cr_lr_dc28_2_21_M1_mhd_c1000']
#dirneed=['smc_lr_2_21_vtd_t1','smc_lr_2_21_vtd_t1_mhd','smc_mr_2_21_vtd_t1','smc_mr_2_21_vtd_t1_mhd']
#dirneed=['smc_lr_2_21_mhd_su']
#dirneed=['mw_lr_2_21_mhd_IC']
#dirneed=['smc_lr_2_21_vtd','smc_mr_2_21_vtd']
#dirneed=['smc_mr_2_21_vtd_t1']
#dirneed=['smc_lr_2_21_vtd']
#dirneed=['smc_lr_2_21_thickdisk','smc_mr_2_21_thickdisk']
#dirneed=['smc_lr_2_21_thickdisk_mhd']
#dirneed=['mw_lr_2_21_mhd','mw_lr_2_21_mhd_IC']
#dirneed=['mw_lr_2_21_mhd','mw_cr_lr_dc28_2_21_M1_mhd','mw_cr_lr_dc28_2_21_M1_mhd_stream','mw_cr_lr_2_21_purestream']
#dirneed=['mw_cr_lr_dc27_2_21_M1_smallts','mw_cr_lr_dc28_2_21_M1_c1000_smallts','mw_cr_lr_dc29_2_21_M1_c2000_smallts','mw_cr_lr_dc29_2_21_M1_c4000_smallts']
#dirneed=['mw_cr_lr_dc28_2_21_M1_smallts','mw_cr_lr_dc28_2_21_M1_mhd']
#dirneed=['mw_cr_lr_dc28_2_21_M1_smallts','mw_cr_lr_dc28_2_21_M1_c1000_smallts','mw_cr_lr_dc28_2_21_M1_c2000_smallts']
#dirneed=['mw_cr_lr_dc27_2_21_M1_smallts','mw_cr_lr_dc28_2_21_M1_smallts','mw_lr_5_27']
#dirneed=['mw_cr_lr_dc28_2_21_M1_smallts','mw_cr_lr_dc28_2_21_M1_c1000_smallts']
#dirneed=['mw_cr_lr_dc29_2_21_M1_smallts','mw_cr_lr_dc29_2_21_M1_c1000_smallts', 'mw_cr_lr_dc29_2_21_M1_c2000_smallts','mw_cr_lr_dc29_2_21_M1_c4000_smallts']
#dirneed=['mw_cr_lr_dc29_2_21_M1_c1000_smallts', 'mw_cr_lr_dc29_2_21_M1_c2000_smallts']
#dirneed=['smc_cr_lr_2_21_purestream']
#dirneed=['smc_cr_lr_dc28_2_21_M1_mhd_stream_c1000']
#dirneed=['mw_cr_lr_dc28_2_21_M1_mhd_stream_c1000']
#dirneed=['mw_cr_lr_dc29_2_21_M1_c1000_smallts', 'mw_cr_lr_dc29_2_21_M1_c2000_smallts','mw_cr_lr_dc29_2_21_M1_c4000_smallts']
#dirneed=['mw_cr_lr_dc28_2_21_M1_smallts','mw_cr_lr_dc28_2_21_M1_c1000_smallts','mw_cr_lr_dc29_2_21_M1_smallts','mw_cr_lr_dc29_2_21_M1_c1000_smallts']
#dirneed=['mw_cr_lr_dc28_2_21_M1_smallts']
#dirneed=['mw_lr_2_21']
#dirneed=['mw_cr_lr_dc29_2_21_M1_smallts','mw_cr_lr_dc29_2_21_M1_c1000_smallts']
#dirneed=['sbc_lr_2_21_tk','sbc_cr_lr_dc28_2_21_tk_initCR1']
#dirneed=['sbc_cr_lr_dc28_2_21_tk','sbc_cr_lr_dc28_2_21_tk_initCR1','sbc_cr_lr_dc28_2_21_tk_initCR10']
#dirneed=['sbc_lr_2_21_tk','sbc_cr_lr_dc28_2_21_tk']
#dirneed=['mw_cr_lr_dc28_2_21_M1_nomaxCR']
#dirneed=['mw_cr_lr_dc28_2_21_M1_smallts','mw_cr_lr_dc28_2_21_M1_nomaxCR']
#dirneed=['mw_cr_lr_dc28_2_21_M1_smallts','mw_cr_mr_dc28_2_21_M1']
#dirneed=['mw_lr_2_21_IC','mw_mr_2_21_IC']
#dirneed=['mw_lr_2_21','mw_lr_2_21_mhd']
#dirneed=['smc_cr_lr_dc28_2_21_M1_smallts']
#dirneed=['smc_cr_lr_dc28_2_21_M1_smallts','smc_cr_lr_dc27_2_21_M1_smallts']
#dirneed=['smc_cr_lr_dc28_2_21_M1_smallts','smc_cr_lr_dc28_2_21_M1_c1000_smallts','smc_cr_lr_dc28_2_21_M1_c2000_smallts']
#dirneed=['m82_lr_2_21']
#dirneed=['smc_cr_lr_dc28_2_21_M1_smallts','smc_cr_lr_dc28_2_21_M1_mhd']
#dirneed=['smc_lr_2_21_su','smc_cr_lr_dc27_2_21_M1_smallts','smc_cr_lr_dc28_2_21_M1_smallts']
#dirneed=['smc_lr_2_21_mhd','smc_lr_2_21_weakB','smc_lr_2_21_mhd_su']
#dirneed=['mw_cr_lr_dc28_2_21_M1_smallts','mw_cr_lr_dc28_2_21_M1_c1000_smallts']
#dirneed=['bwmwllrdc28', 'bwmwlrdc28', 'bwmwmrdc28']
#dirneed=['bwmwlr','mw_mr_2_21','bwmwllrdc28', 'bwmwlrdc28', 'bwmwmrdc28']
#dirneed=['mw_cr_lr_dc28_2_21_M1_smallts','mw_cr_lr_dc28_2_21_M1_c1000_smallts','mw_cr_lr_dc28_2_21_M1_c2000_smallts','mw_cr_lr_dc29_2_21_M1_smallts','mw_cr_lr_dc29_2_21_M1_c1000_smallts','mw_cr_lr_dc29_2_21_M1_c2000_smallts','mw_cr_lr_dc29_2_21_M1_c4000_smallts']
#dirneed=['mw_cr_lr_dc29_2_21_M1_smallts','mw_cr_lr_dc29_2_21_M1_nomaxCR']
#dirneed=['smc_lr_2_21_thickdisk','smc_lr_2_21_thickdisk_mhd','smc_mr_2_21_thickdisk']
#dirneed=['smc_lr_2_21_su','smc_mr_2_21_su']
#dirneed=['smc_lr_2_21_mhd','smc_mr_2_21_mhd']
#dirneed=['smc_lr_2_21_mhd','smc_cr_lr_dc28_2_21_M1_mhd']
#dirneed=['smc_lr_2_21_mhd','smc_mr_2_21_mhd']
#dirneed=['smc_lr_2_21_su','mw_lr_2_21']
#dirneed=['smc_cr_lr_dc27_2_21_M1_smallts','smc_cr_lr_dc28_2_21_M1_smallts','smc_cr_lr_dc28_2_21_M1_mhd']
#dirneed=['smc_lr_2_21','smc_lr_2_21_mhd','smc_lr_2_21_weakB']
#dirneed=['smc_lr_2_21_su','smc_lr_2_21_mhd_su']
#dirneed=['smc_lr_2_21','smc_lr_2_21_mhd','smc_cr_lr_dc28_2_21_M1_mhd_stream','smc_cr_lr_2_21_purestream']
#dirneed=['smc_cr_lr_dc27_2_21_M1_smallts','smc_cr_lr_dc28_2_21_M1_c1000_smallts','smc_cr_lr_dc29_2_21_M1_c2000_smallts']
#dirneed=['smc_lr_2_21','smc_lr_2_21_mhd','smc_cr_lr_dc27_2_21_M1_smallts','smc_cr_lr_dc28_2_21_M1_smallts','smc_cr_lr_dc28_2_21_M1_mhd','smc_cr_lr_dc28_2_21_M1_mhd_stream','smc_cr_lr_2_21_purestream']
#dirneed=['smc_lr_2_21_mhd', 'mw_lr_2_21_mhd']
#dirneed=['smc_lr_2_21_su','smc_lr_2_21_mhd']
#dirneed=['smc_lr_2_21_mhd']
#dirneed=['smc_lr_2_21_mhd_su']
#dirneed=['smc_lr_2_21_mhd','smc_lr_2_21_weakB']
#dirneed=['smc_lr_2_21_weakB_p10G_su']
#dirneed=['smc_lr_2_21_weakB_p12G']
#dirneed=['smc_cr_lr_dc28_2_21_M1_mhd']
#dirneed=['smc_cr_lr_dc27_2_21_M1_smallts','smc_cr_lr_dc28_2_21_M1_smallts','smc_cr_lr_dc28_2_21_M1_mhd']
#dirneed=['smc_lr_2_21','smc_cr_lr_dc27_2_21_M1_smallts','smc_cr_lr_dc28_2_21_M1_smallts']
#dirneed=['mw_lr_2_21', 'mw_cr_lr_dc27_2_21_M1_smallts','mw_cr_lr_dc28_2_21_M1_smallts']
#dirneed=['smc_cr_lr_dc28_2_21_M1_smallts','smc_cr_lr_dc28_2_21_M1_smallts_cutcr']
#dirneed=['smc_cr_lr_dc28_2_21_M1_smallts_cutcr']
#dirneed=['FIRE_2_0_h573_2_21']
#dirneed=['mw_cr_lr_dc28_6_15_ets']
#dirneed=['smc_lr_2_21_su']
#dirneed=['sbc_cr_lr_dc28_2_21_tk']
#dirneed=['sbc_lr_2_21_tk','sbc_cr_lr_dc28_2_21_tk']
#dirneed=['sbc_lr_2_21_tk','sbc_cr_lr_dc28_2_21_tk']
#dirneed=['sbc_lr_2_21_su','sbc_mr_2_21_su']
#dirneed=['smc_lr_2_21_ehalo_IC']
#dirneed=['M1dc28g_2_21_smallts_cutcr']
#dirneed=['M1dc29g_2_21_smallts_cutcr']
#dirneed=['FIRE_2_0_h573_cr_2_21_M1_z0']
#dirneed=['FIRE_2_0_h573_2_21','FIRE_2_0_h573_cr_2_21_M1_z0','FIRE_2_0_h573_cr_2_21_M1_z0_dc27']
#dirneed=['mw_cr_lr_dc28_2_21_M1_smallts_testsn']
#dirneed=['mw_cr_lr_dc28_6_15_ets_fhfcr']
#dirneed=['mw_cr_lr_dc27_2_21_M1_fhfcr_ets_ssts']
#dirneed=['mw_cr_lr_dc29_2_21_M1_c1000_smallts']
#dirneed=['FIRE_2_0_h573_cr_2_21_M1_z0_ts2']
#dirneed=['FIRE_2_0_h573_cr_2_21_M1_z0_defaultts']
#dirneed=['M1dc28g_2_21_c100_smallts_gassphere']
#dirneed=['M1dc28g_2_21_smallts_gassphere']
#dirneed=['M1dc28g_2_21_c100_smallts_gassphere','M1dc28g_2_21_smallts_gassphere', 'M1dc28g_2_21_c1000_smallts_gassphere']
#dirneed=['mw_cr_lr_dc28_2_21_M1_smallts']
#dirneed=['smc_cr_lr_dc27_2_21_M1_smallts','smc_cr_lr_dc28_2_21_M1_smallts']
#dirneed=['M1dc29g_2_21_smallts_gassphere','M1dc29g_2_21_c1000_smallts_gassphere']
#dirneed=['FIRE_2_0_h573_cr_2_21_M1_z0']
#dirneed=['FIRE_2_0_h573_cr_2_21_M1_z0','FIRE_2_0_h573_cr_2_21_M1_z0_dc27']
#dirneed=['FIRE_2_0_h573_cr_2_21_M1_z0']
#dirneed=['FIRE_2_0_h573_cr_2_21_M1_z0','FIRE_2_0_h573_cr_2_21_M1_z0_cd100']
#dirneed=['mw_cr_lr_dc27_2_21_M1_smallts_ets']
#dirneed=['mw_cr_lr_dc27_2_21_M1_allowimbalance']
#dirneed=['mw_cr_lr_dc28_2_21_M1_allowimbalance']
#dirneed=['mw_cr_lr_dc27_2_21_M1_smallts']
#dirneed=['mw_cr_lr_dc27_2_21_M1_smallts_ets']
#dirneed=['mw_cr_lr_dc28_2_21_M1_fhfcr_ets']
#dirneed=['mw_cr_lr_dc27_2_21_M1_fhfcr_ets']
#dirneed=['M1dc28_2_21_c2000_smallts']
#dirneed=['dc28','M1dc28_2_21_c2000_smallts']
#dirneed=['dc27','M1dc27_2_21','M1dc27_2_21_c100','M1dc27_r128_2_21']
#dirneed=['dc27','M1dc27g_2_21_c100']
#dirneed=['dc27','dc27_slowc']
#dirneed=['M1dc28g_2_21_smallts', 'M1dc28g_2_21_c100_smallts']
#dirneed=['M1dc29g_2_21_smallts','M1dc29g_2_21_c1000_smallts','M1dc29g_2_21_c2000_smallts']
#dirneed=['M1dc27g_2_21_smallts']
#dirneed=['dc27g']
#dirneed=['M1dc29g_2_21_c2000_smallts']
#dirneed=['M1dc27g_2_21_c100_smallts','M1dc27g_2_21_smallts','dc27g']
#dirneed=['M1dc27_2_21_c100_smallts','M1dc27_2_21_smallts']
#dirneed=['M1dc27g_2_21_smallts','M1dc27g_2_21_c100_smallts']
#dirneed=['dc28','M1dc28_2_21','M1dc28_2_21_c1000','M1dc28_2_21_c1000_smallts','M1dc28_2_21_c2000_smallts']
#dirneed=['M1dc28_2_21_c1000','M1dc28_2_21','dc28_2_21']
#dirneed=['M1dc27_2_21','dc27']
#dirneed=['dc28','M1dc28']
#dirneed=['dc26','M1dc26']
#dirneed=['dc27','M1dc27']
#dirneed=['M1dc27']
#dirneed=['FIRE_2_0_h573_CRtest2']
#dirneed=['smc_lr_6_15_sn300']
#dirneed=['sbc_lr_6_15']
#dirneed=['smc_cr_lr_dc28_6_15']
#dirneed=['smc_cr_lr_dc27_6_15']
#dirneed=['smc_cr_lr_dc26_6_15']
#dirneed=['smc_lr_6_15_sn300','smc_cr_lr_dc27_6_15']
#dirneed=['smc_lr_6_15_sn300','smc_cr_lr_dc27_6_15','smc_cr_lr_dc28_6_15']
#dirneed=['FIRE_2_0_h573_cr_6_15_ets_dc27','FIRE_2_0_h573_cr_6_15_ets']
#dirneed=['FIRE_2_0_h573_6_15','FIRE_2_0_h573_cr_6_15_ets_dc0']
#dirneed=['FIRE_2_0_h573_6_15','FIRE_2_0_h573_cr_6_15_ets_dc0','FIRE_2_0_h573_cr_6_15_ets_dc27','FIRE_2_0_h573_cr_6_15_ets']
#dirneed=['FIRE_2_0_or_h573_criden1000_noaddm_sggs']
#dirneed=['FIRE_2_0_h573_6_15']
#dirneed=['FIRE_2_0_h573_cr_6_15_ets_dc0','FIRE_2_0_h573_cr_6_15_ets_dc27']
#dirneed=['FIRE_2_0_h573_cr_6_15_ets_dc0']
#dirneed=['FIRE_2_0_h573_cr_6_15_ets_dc27']
#dirneed=['FIRE_2_0_h573_cr_6_15_ets']
#dirneed=['FIRE_2_0_h573_cr_6_15_M1']
#dirneed=['FIRE_2_0_h573_cr_6_15_M1_z0']
#dirneed=['FIRE_2_0_h573_cr_6_15_M1_c2000']
#dirneed=['FIRE_2_0_h573_cr_6_15_ets','FIRE_2_0_h573_cr_6_15_ets_dc27#']
#dirneed=['mw_cr_lr_dc28_hole_6_15_M1']
#dirneed=['mw_cr_lr_dc28_hole_6_15_ets']
#dirneed=['mw_cr_lr_dc29_6_15_ets']
#dirneed=['mw_cr_lr_dc29_6_15_ets','mw_cr_lr_dc29_6_15_ets_nolimitspeed','mw_cr_lr_dc29_6_15_ets_limitspeed','mw_cr_lr_dc29_6_15_ets_extremelimit']
#dirneed=['mw_cr_lr_dc28_6_15_ets_freeze_nocooling']
#dirneed=['mw_cr_lr_dc27_6_15_ets','mw_cr_lr_dc27_6_15_M1_nofluxterm']
#dirneed=['mw_cr_lr_dc28_6_15_ets','mw_cr_lr_dc28_6_15_M1_nofluxterm','mw_cr_lr_dc28_6_15_M1_c5000_nofluxterm']
#dirneed=['mw_cr_lr_dc28_hole_6_15_ets','mw_cr_lr_dc28_hole_6_15_M1']
#dirneed=['mw_cr_lr_dc28_6_15_ets','mw_cr_lr_dc28_6_15_M1','mw_cr_lr_dc28_6_15_M1_nofluxterm']
#dirneed=['mw_cr_lr_dc28_6_15_ets','mw_cr_lr_dc28_6_15_M1','mw_cr_lr_dc28_6_15_M1_ca']
#dirneed=['mw_cr_lr_dc29_6_15_M1_nofluxterm']
#dirneed=['mw_cr_lr_dc28_6_15_M1_nofluxterm']
#dirneed=['mw_cr_lr_dc27_6_15_M1_nofluxterm']
#dirneed=['mw_cr_lr_dc29_6_15_M1_nofluxterm']
#dirneed=['mw_cr_lr_dc29_6_15_c5000_M1_nofluxterm']
#dirneed=['mw_cr_lr_dc29_6_15_ets','mw_cr_lr_dc29_6_15_M1_nofluxterm','mw_cr_lr_dc29_6_15_c2000_M1_nofluxterm','mw_cr_lr_dc29_6_15_c5000_M1_nofluxterm']
#dirneed=['mw_cr_lr_dc29_6_15_M1_nofluxterm','mw_cr_lr_dc29_6_15_c2000_M1_nofluxterm','mw_cr_lr_dc29_6_15_c5000_M1_nofluxterm']
#dirneed=['mw_cr_lr_dc29_6_15_ets','mw_cr_lr_dc29_6_15_M1_nofluxterm','mw_cr_lr_dc29_6_15_c2000_M1_nofluxterm','mw_cr_lr_dc29_6_15_c5000_M1_nofluxterm','mw_cr_lr_dc29_6_15_c2000_M1_nofluxterm_relaxHll','mw_cr_lr_dc29_6_15_c5000_M1_nofluxterm_relaxHll']
#dirneed=['mw_cr_lr_dc29_6_15_ets','mw_cr_lr_dc29_6_15_ets_limitspeed']
#dirneed=['mw_cr_lr_dc29_6_15_ets','mw_cr_lr_dc29_6_15_M1_nofluxterm','mw_cr_lr_dc29_6_15_c2000_M1_nofluxterm','mw_cr_lr_dc29_6_15_c5000_M1_nofluxterm','mw_cr_lr_dc29_6_15_c2000_M1_nofluxterm_relaxHll']
#dirneed=['mw_cr_lr_dc27_6_15_ets','mw_cr_lr_dc27_6_15_M1_nofluxterm']
#dirneed=['mw_cr_lr_dc29_6_15_ets','mw_cr_lr_dc29_6_15_M1_nofluxterm','mw_cr_lr_dc29_6_15_c2000_M1_nofluxterm','mw_cr_lr_dc29_6_15_c2000_M1_nofluxterm_relaxHll']
#dirneed=['mw_cr_lr_dc29_6_15_c5000_M1_nofluxterm_relaxHll']
#dirneed=['mw_cr_lr_dc29_6_15_ets']
#dirneed=['mw_cr_lr_dc29_6_15_M1_nofluxterm#']
#dirneed=['mw_cr_lr_dc28_6_15_M1_fhfcr_ets']
#dirneed=['mw_cr_lr_dc28_6_15_ets_freeze']
#dirneed=['m12_tmp']
#dirneed=['mw_cr_lr_dc28_6_15_M1_fhfcr_ets_relaxHll']
#dirneed=['mw_cr_lr_dc27_2_10_M1']
#dirneed=['mw_cr_lr_dc28_2_10_M1_fhfcr_ets']
#dirneed=['mw_cr_lr_dc27_6_15_M1_fhfcr_ets']
#dirneed=['mw_cr_lr_dc27_2_21_M1_fhfcr_ets','mw_cr_lr_dc27_6_15_ets_fhfcr']
#dirneed=['mw_cr_lr_dc27_2_10_M1_fhfcr_ets']
#dirneed=['mw_cr_lr_dc29_6_15_ets_fhfcr','mw_cr_lr_dc29_6_15_M1_fhfcr_ets','mw_cr_lr_dc29_2_10_M1_fhfcr_ets']
#dirneed=['mw_cr_lr_dc28_6_15_ets_fhfcr','mw_cr_lr_dc28_6_15_M1_fhfcr_ets','mw_cr_lr_dc28_6_15_M1_fhfcr_ets_relaxHll','mw_cr_lr_dc28_2_10_M1_fhfcr_ets']
#dirneed=['mw_cr_lr_dc28_6_15_ets_fhfcr','mw_cr_lr_dc28_6_15_M1_fhfcr_ets','mw_cr_lr_dc28_6_15_M1_fhfcr_ets_relaxHll']
#dirneed=['mw_cr_lr_dc28_6_15_c2000_M1_freezehydro_freecr']
#dirneed=['mw_cr_lr_dc28_6_15_ets_freeze_nocooling']
#dirneed=['m12_tmp','m12i_res7000']
#dirneed=['m12_tmp','m12i_res7000_output']
#dirneed=['m10q_M1CR']
#dirneed=['mw_effeos_cr']
#dirneed=['mw_lr_5_27','mw_lr_5_27_weaksn']
#dirneed=['m10q_mass250_MHDCR_K100_M1_tkFIX','m10q_mass250_MHDCR_K1000_M1_tkFIX']
#dirneed=['m11q_mass7000_MHDCR_K100_M1_tkFIX']
#dirneed=['m11q_mass7000_MHDCR_K100_M1_tkFIX','m11q_mass7000_MHDCR_K1000_M1_tkFIX']
#dirneed=['mw_effeos','mw_effeos_cr']
#dirneed=['mw_cr_lr_nodiff_6_15']
#dirneed=['mw_cr_lr_dc28_6_15_M1_c5000']
#dirneed=['mw_cr_lr_dc28_6_15_ets','mw_cr_lr_dc28_6_15_M1']
#dirneed=['mw_cr_lr_dc27_6_15_ets']
#dirneed=['mw_cr_lr_dc27_6_15_ets','mw_cr_lr_dc28_6_15_ets']
#dirneed=['mw_cr_lr_dc28_hole_6_15_ets']
#dirneed=['mw_cr_lr_dc28_6_15_ets','mw_cr_lr_dc28_hole_6_15_ets']
#dirneed=['mw_cr_lr_dc27_5_27']
#dirneed=['mw_lr_5_27']
#dirneed=['smc_lr_6_15']
#dirneed=['smc_lr_6_15','smc_lr_Feb']
#dirneed=['smc_lr_6_15_sn300']
#dirneed=['smc_lr_6_15', 'smc_lr_Feb']
#dirneed=['smc_lr_Feb']
#dirneed=['mw_lr_hole_5_27']
#dirneed=['mw_cr_lr_dc28_6_6_ets']
#dirneed=['smc_booth_5_27']
#dirneed=['FIRE_2_0_or_h573_criden1000_noaddm_sggs']
#dirneed=['FIRE_2_0_h573_CRtest2']
#dirneed=['FIRE_2_0_h573_CRtest2_equaltimestep']
#dirneed=['FIRE_2_0_or_h573_criden1000_noaddm_sggs']
#dirneed=['FIRE_2_0_or_h573_criden1000_noaddm_sggs']
#dirneed=['FIRE_2_0_h573_CRtest2']
fmeat=dirneed[-1]
#fmeat='m11d'
#fmeat='mwIC'
#fmeat='573'
#fmeat='m82'
#fmeat='mwstrll'
#fmeat='smccrmhd'
#fmeat='smcdc28M1'
#fmeat='mwdc29M1'
#fmeat='smccrmhd_testsn'
#fmeat='smcmhdsu'
#fmeat='smcmhdp10'
#fmeat='smcweakB'
#fmeat='mwmr'
#fmeat='mwstrc1000'
#fmeat='mwcrmrstrll4va'
#fmeat='sbc'
#fmeat='m12imhdcv'
#fmeat='bwsbclrdc28'
#fmeat='bwsbclr'
#fmeat='bwsbclrdc28mhd'
#fmeat='bwsmclrdc29'
#fmeat='bwmwmrdc29'
#fmeat='bwsbclrdc28'
#fmeat='smc'
#fmeat='mw'
#fmeat='mw1000kpc'
#fmeat='cosmocr70'
#fmeat='sbcmhd'
#fmeat='mwlrstrll4va'
#fmeat='mwdc29'
#fmeat='smclrstrnl4va'
#fmeat='m11g'
#fmeat='bwsmclrmhd'
#fmeat='smcv30'
#fmeat='m10vcr_b_70'
#fmeat='m10qcr_b_70'
#fmeat='m10vmhdcv'
#fmeat='cr_b_70'
#fmeat='smccutcr'
#fmeat='smccrlrmhd'
#fmeat='smctd'
#fmeat='mwstrdc28'
#fmeat='mwmhd'
#fmeat='sbcextend4'
#fmeat='bwmwmhd'
#fmeat='M1str'
#fmeat='mwM1'
#fmeat='mwreso'
#fmeat='smc'
#fmeat='sbc'
#fmeat='mwlr'
#fmeat='smcmr'
#fmeat='mwiso'
#fmeat='bwmwlrehh'
#fmeat='dc28gsn0'
#fmeat='m12i'
#fmeat='m11g'
#fmeat='m10v'
#fmeat='m11b'
#fmeat='m10vcrb70'
#fmeat='m10vcr700'
#fmeat='mw'
#fmeat='mwdc27ds'
#fmeat='m12icr_700'
#fmeat='sbc'
#fmeat='mwcr'
#fmeat='mwcutcr'
#fmeat='smccrlr'
#fmeat='bwsmccrlrdc28str'
#fmeat='mwcrlr'
#fmeat='smctestmhdIC'
#fmeat='smctestmhdmr'
#fmeat='mwcrdc28strc1000'
#fmeat='smctestM1'
#fmeat='smcmhdweakBIC'
#fmeat='smcmrtd'
#fmeat='mwstrongB'
#fmeat='smcM1dc28str'
#fmeat='smcvtdt1'
#fmeat='smctd'
#fmeat='mwdc28M1'
#fmeat='mwdc29M1'
#fmeat='mwdc29nomaxCR'
#fmeat='sbcinitcr'
#fmeat='mwtestM1'
#fmeat='mwcrmhd'
#fmeat='mwcrmhd_testsn'
#fmeat='mwcrmr'
#fmeat='mwmhd'
#fmeat='sbccr'
#fmeat='mwlrmhd'
#fmeat='smcmhdsu'
#fmeat='smcweakB'
#fmeat='smcdc2728'
#fmeat='mwdc2728'
#fmeat='m82_sn20'
#fmeat='mwM1_2_21'
#fmeat='mwdc29M1_2_21_smallts_sn10_cutcr'
#fmeat='mwdc29M1_2_21_smallts_sn1'
#fmeat='mwdc29M1_2_21_c1000'
#fmeat='M1dc29g_2_21_smallts_cutcr_sn15'
#fmeat='M1dc29g_2_21_smallts_sn1'
#fmeat='M1str_'+str(Nsnap)
#fmeat='M1dc27g_2_21_smallts_sn'+str(Nsnap)
#fmeat='M1dc29g_2_21_smallts_wc_sn8'
#fmeat='M1dc28g_2_21_smallts_wc_sn5'
#fmeat='M1dc28g_2_21_smallts_cutcr_sn40'
#fmeat='M1dc27g_2_21_smallts_sn50'
#fmeat = 'dc28_2_21_sn20'
#fmeat='dc27_gaussian_sn5'
#fmeat='dc28M1_2_21_sn30'
#fmeat='dc27M1_2_21_smallts_sn25'
#fmeat='mwdc27M1nf'
#fmeat='mwdc28M1nf'
#fmeat='mwdc28M1nfvsf'
#fmeat='573M1z0_sn548'
#fmeat='mwdc29limit'
#fmeat='mwdc28fh'
#fmeat='mwdc29elimit'
#fmeat='mwdc29fM1_2_10'
#fmeat='CRtestM1dc27_2_21_c100'
#fmeat='CRtestdc28_2_21'
#fmeat='CRtestM1dc27_2_21_smallts_sn25'
#fmeat='CRtestM1dc28_2_21_c2000_smallts'
#fmeat='CRtestdc26_sn100'
#fmeat='mwdc28fM1_2_10'
#fmeat='mwdc29M1nfrelax'
#fmeat='mwdc29M1onlynf'
#fmeat='mwdc28M1ca'
#fmeat='smccr'
#fmeat = 'mwcr'
#fmeat = 'smcmr'
#fmeat='mwcrlr'
#fmeat='mw100kpc'
#fmeat='sbc'
#fmeat='sbc_starburst'
#fmeat='smclr'
#fmeat='mwmr'
#fmeat='mw1cm_3'
#fmeat='smc100kpc'
#fmeat='mwnc'
#fmeat='cosmo'
#fmeat='mwtest'
#fmeat='smc0_03Rv'
#fmeat='sbc0_25kpc'
#fmeat='smc'
#fmeat='smc1kpc'
#fmeat='mw1kpc'
#fmeat='sbc1kpc'
#fmeat='mwreso'
#fmeat='mwreso_z50'
#startno=499
#startno=590
#startno=10
#startno=30
#startno=400
#startno=480
#startno=300
#startno=350
#startno=100
#startno=200
#startno=490
#startno=400
#startno=570
#startno=490
#startno=500
#startno=11
#startno=300
#startno=440
#startno=600
#startno=10
startno=490
#startno=0
#startno=480
#startno=450
#startno=506
#startno=200
#startno=350
#startno=499
#startno=600
#startno=550
#startno=590
#Nsnap=10
#Nsnap=50
#Nsnap=15
#Nsnap=20
#fmeat='mwmr'
#fmeat='mwdc27ds'
#fmeat='smc10kpc'
#fmeat='smc'
#fmeat='sbc'
#fmeat='M1str_'+str(Nsnap)
#fmeat='dc27g_sn'+str(Nsnap)
#fmeat='M1dc29g_2_21_smallts_sn'+str(Nsnap)
#fmeat='M1dc27g_2_21_smallts_sn'+str(Nsnap)
#fmeat='M1dc27_2_21_smallts_sn'+str(Nsnap)
#Nsnap=100
#Nsnap=601
#Nsnap=999
#Nsnap=200
#Nsnap=0
#Nsnap=210
#Nsnap=250
#Nsnap=300
#Nsnap=460
Nsnap=500
#Nsnap=651
#Nsnap=10
#Nsnap=801
#Nsnap=430
#Nsnap=600
#Nsnap=650
#Nsnap=590
#Nsnap=10
#Nsnap=560
#Nsnap=506
#Nsnap=200
#Nsnap=15
#Nsnap=10
#Nsnap=30
#Nsnap=440
#Nsnap=1000
#Nsnap=250
#Nsnap=400
#Nsnap=601
#Nsnap=300
#Nsnap=700
#Nsnap=500
#Nsnap=506
#Nsnap=518
#Nsnap=8
#Nsnap=555
#Nsnap=600
#Nsnap=1
#Nsnap=161
#Nsnap=160
#Nsnap=580
#Nsnap=30
#Nsnap=10
snapsep=10
#snapsep=1
#snapsep=5
#snapsep=490
#snapsep=20
#fmeat='dc27gsn'+str(Nsnap)
#wanted='phaseTwindz'
#wanted='Twindbar'
#wanted='Tbar'
#wanted='phaseTz'
#wanted='rhoT'
#wanted='rhoTwind'
#wanted='Tzmwed'
#wanted='vzz'
#wanted='vzztrack'
#wanted='Tvz'
#wanted='Tvztrack'
#wanted='testCRcum'
#wanted='testCRdis'
#wanted='Tz'
#wanted='rhoztrack'
#wanted='Tztrack'
#wanted='rhoz'
#wanted = 'nismB'
#wanted ='dpcrz_rhog'
#wanted = 'dpcrz'
#wanted = 'massloadingSasha'
#wanted = 'outflowSasha'
#wanted = 'massloadingBooth'
#wanted = 'outflowBooth'
#wanted = 'gassurden'
#wanted = 'gassurdentime'
#wanted = 'outflowall'
#wanted = 'gasdenmidplane'
#wanted = 'outflow'
#wanted = 'sfrv'
#wanted = 'cumnsmrad'
#wanted = 'nsrad'
#wanted = 'sfrrad'
#wanted = 'sfrarad'
#wanted='sncretime'
#wanted='crdenplanes'
#wanted='crdenv'
#wanted='crdenmidplane'
#wanted='gammadensph'
#wanted='credensph'
#wanted='decayratiosph'
#wanted='decaytimesph'
#wanted='gammasph'
#wanted='Gmencsph'
#wanted='gasdensph'
#wanted='gasden'
#wanted='avecrdenr'
#wanted='avekedenr'
#wanted='avethdenr'
#wanted='aveedenr'
#wanted='denrtime'
#wanted='avedenr'
#wanted = 'FpFgr'
#wanted='gmr'
#wanted = 'rhovr'
#wanted='denr'
#wanted = 'rhovr_denr'
#wanted='crdenr'
#wanted='crer'
#wanted='crdpr'
#wanted='crecumz'
#wanted='crecumr'
#wanted='Begyz'
#wanted='kez'
#wanted='ethz'
#wanted='vz'
#wanted='pz'
#wanted='dpz'
#wanted='gammar'
#wanted='gammaz'
#wanted='cresph'
#wanted='gmz'
#wanted='gasdenz'
#wanted='gasdeny'
#wanted='gasdenx'
#wanted='nismcrg'
#wanted = 'rcrad'
#wanted='nismcrl'
#wanted='nismcrad'
#wanted='nismcre'
#wanted='nismgamma'
#wanted='nismcumgamma'
#wanted='Begytime'
#wanted='Brmstime'
#wanted='Begy_mtime'
#wanted='Bcumgm'
#wanted='nismcumgm'
#wanted='gammadecaytime'
#wanted='gammacumdecaytime'
#wanted='credecaytime'
#wanted='crecumdecaytime'
#wanted='nismcumcre'
#wanted='pvolcumcre'
#wanted='pdxcumcre'
#wanted='pdxcumN'
#wanted='pvolcumN'
#wanted = 'crasnap'
#wanted='testCRcum'
#wanted='rhoT'
#wanted='cratime'
#wanted='gasdenv'
#wanted='crdrad'
#wanted='cramap'
#wanted='cramapv'
wanted='Vquiverxy'
#wanted='Vquiverxz'
#wanted='Bquiverxy'
#wanted='Bquiverxz'
#wanted='Smassxy'
#wanted='Smassxz'
#wanted='testquiver'
#wanted='Twindm'
#wanted='vzwindm'
#wanted='Tz'
#wanted='pcrpth'
#wanted='Alfven_sound'
#wanted='cramapv'
#wanted='crturmap'
#wanted='outflowwrite'
#wanted = 'gassurden'
#wanted='gasfrac'
#wanted='printparmass'
#wanted='enchange'
#wanted='dirssfr10_200Myr'
#wanted='dirage10_200Myr'
#wanted='diragehis'
#wanted='dirage200Myr'
#wanted='diragestd'
#wanted='diragestd10_200Myr'
#wanted='diragemax'
#wanted='dirsm'
#wanted='dirsfr'
#wanted='testcr'
#wanted='dirheader'
#wanted='dirgammasfr'
#wanted='dirgamma'
#wanted='cosmic_raywrite'
#wanted='crdensol' #cosmic ray energy density at solar circle (observation ~ 1 eV/cm^3)
#wanted='smchange'
#wanted='cumgamma'
#wanted='smass'
#wanted='crdtime'
#wanted='gaskechange'
#wanted='gaske' #evolution of total gas kinetic energy 
#wanted='crchange'
#wanted='crtime'
#wanted='crcooling'
#wanted='lgratio'
#wanted='dg_sfr'
#wanted='gamma_partit'
#wanted='gamma_aveCR'
#wanted='dlgratio'
#wanted='gsmratio'
#title='LSG'
title='MW'
#title='Cosmo'
#title='Dwarf'
#titleneed='Dwarf'
#title='SBG'
titleneed=title
needlog=0
dclabelneed=1
correctIa=0
useM1=1
the_prefix='snapshot'
the_suffix='.hdf5'
withinr=15.0
nogrid=40
maxlength=10.0
med = -0.1 
wholeboxgas=1
diskgas=1
Rfrac = 0.5
nosum=0
xaxis_snapno=0
#betapi = 0.7 ##should be between 0.5-0.9 from Lacki 2011 #fraction of pi has enough energy (>1GeV)
#yr_in_sec = 3.2e7
#nopi_per_gamma = 3.0
#solar_mass_in_g = 2e33
#km_in_cm = 1e5
#kpc_in_cm = 3.086e21
#erg_in_eV = 6.242e11
#cspeed_in_cm_s = 3e10
#proton_mass_in_g = 1.67e-24
#Kroupa_Lsf=6.2e-4
#Salpeter_Lsf=3.8e-4
#pidecay_fac = 2.0e5*250.0*yr_in_sec #over neff in cm_3 for pi decay time in s
#hadronicdecayrate = 5.86e-16 #from Guo 2008. The hadronic decay rate of CR
#coulombdecayrate = 1.65e-16 #Coulomb loss of CR


print 'wanted', wanted
print 'fmeat', fmeat
print 'runtodo', dirneed



if wanted=='Alfven_sound':
    rcParams['figure.figsize'] = 5,4
        withinr=20.0
        maxlength=30.0
    nbin=20
    needvA=1
    userad=1
        for runtodo in dirneed:
                info=outdirname(runtodo, Nsnap)
            
                haveB=info['haveB']
        rundir=info['rundir']
        runtitle=info['runtitle']
        slabel=info['slabel']
        snlabel=info['snlabel']
        dclabel=info['dclabel']
        resolabel=info['resolabel']
        the_snapdir=info['the_snapdir']
        Nsnapstring=info['Nsnapstring']
        havecr=info['havecr']
        Fcal=info['Fcal']
        iavesfr=info['iavesfr']
        timestep=info['timestep']
        maindir=info['maindir']
        color=info['color']
        haveB=info['haveB']
                ptitle=title
                if runtitle=='SMC':
                        ptitle='Dwarf'
                elif runtitle=='SBC':
                        ptitle='Starburst'
                elif runtitle=='MW':
                        ptitle=r'$L\star$ Galaxy'
        print 'the_snapdir', the_snapdir
        print 'Nsnapstring', Nsnapstring
        print 'havecr', havecr
        cosmo=info['cosmo']
        if cosmo==1:
            h0=1
        else:
            h0=0
        lsn='solid'
        if haveB>0:
            lsn='dashed'
        header = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo, header_only=1)
        ascale = header['time']
        print 'this time', ascale
        thisred = 1./ascale-1.
        hubble = header['hubble']
        print 'hubble', hubble
        G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo)
        if cosmo==1:
            halosingle = read_halo_history(rundir, halonostr='00',hubble=hubble, comoving=0, maindir=maindir, singlesnap=1, atime=ascale)
            xcen = halosingle['x']
            ycen = halosingle['y']
            zcen = halosingle['z']
            xvcen = halosingle['xv']
            yvcen = halosingle['yv']
            zvcen = halosingle['zv']
            Rvirnow = halosingle['R']
            MgAHF = halosingle['Mg']
        else:
            xcen=0
            ycen=0
            zcen=0
            xvcen=0
            yvcen=0
            zvcen=0
        Gpos = G['p']
        Gvel = G['v']
        Tb = G['u'] #internal energy per unit mass # in km^2/s^2
        Gx = Gpos[:,0]-xcen
        Gy = Gpos[:,1]-ycen
        Gz = Gpos[:,2]-zcen
        Gvx = Gvel[:,0]-xvcen   #km/s
        Gvy = Gvel[:,1]-yvcen
        Gvz = Gvel[:,2]-zvcen
        rho = G['rho']
        Grho = rho*1e10*Msun_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm #g/cm^3
        Gm = G['m']*1e10*Msun_in_g #g
        Neb = G['ne']
        if haveB>0 and needvA>0:
            Bfield = G['B']
            Bx = Bfield[:,0]
            By = Bfield[:,1]
            Bz = Bfield[:,2]
        if havecr>0:
            cregy = G['cregy']*1e10*Msun_in_g*km_in_cm*km_in_cm #erg
        cut = (Gx*Gx+Gy*Gy)<withinr*withinr
        Gx=Gx[cut]
        Gy=Gy[cut]
        Gz=Gz[cut]
        Gr = np.sqrt(Gx*Gx+Gy*Gy+Gz*Gz)
        Gvx=Gvx[cut]
        Gvy=Gvy[cut]
        Gvz=Gvz[cut]
        Gm=Gm[cut]
        Grho=Grho[cut]
        rho=rho[cut]
        Tb = Tb[cut]
        Neb= Neb[cut]
        if havecr>0:
            cregy = cregy[cut]
        if haveB>0 and needvA>0:
            Bx = Bx[cut]
            By = By[cut]
            Bz = Bz[cut]
            B2 = Bx*Bx+By*By+Bz*Bz
            B = np.sqrt(B2)
            nism=rho*1e10*Msun_in_g/protonmass_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm
            #print 'np.median(np.sqrt(B2))', np.median(np.sqrt(B2))
            #print 'np.median(Grho)', np.median(Grho)
            valfven=B/np.sqrt(np.pi*4.*nism*protonmass_in_g)/km_in_cm
                        print 'np.amax(B)', np.amax(B)
                        print 'np.amax(nism)', np.amax(nism)
            print 'np.amax(valfven)', np.amax(valfven)
            print 'np.average(valfven)', np.average(valfven)
            #valfven = np.sqrt(B2/4.0/np.pi/Grho)/km_in_cm #km/s
            #print 'np.median(np.sqrt(B2/Grho))', np.median(np.sqrt(B2/Grho))/km_in_cm
            #print 'np.sqrt(np.median(B2)/np.median(Grho))', np.sqrt(np.median(B2)/np.median(Grho))/km_in_cm
            #print 'np.median(valfven)', np.median(valfven)
        TrueTemp, converted_rho  = SF.convertTemp(Tb, Neb, rho)
        Tm = TrueTemp*Gm
        
        zver = np.linspace(0.,maxlength,num=nbin+1)
        vthl = []
        val = []
        vthul = []
        vthdl = []
        vaul = []
        vadl = []
        vcrl = []
        for j in range(nbin):
            if userad==0:
                cutz = (Gz>zver[j]) & (Gz<zver[j+1])
            else:
                cutz = (Gr>zver[j]) & (Gr<zver[j+1])
            Temp = np.sum(Tm[cutz])/np.sum(Gm[cutz])
            MeanWeights_con = np.average(4.0/(3*H_MASSFRAC+1+4*H_MASSFRAC* Neb[cutz]), weights=Gm[cutz])
            #print 'Temp', Temp
            cs = np.sqrt(GAMMA*BOLTZMANN*Temp/MeanWeights_con/protonmass_in_g)/km_in_cm
            vthl = np.append(vthl, cs)
        #   vthl = np.append(vthl, np.median(vthermal[cutz]))
        #   vthul = np.append(vthul, np.percentile(vthermal[cutz],82))
        #   vthdl = np.append(vthdl, np.percentile(vthermal[cutz],15))
            if havecr>0:
                Pcr = (CRgamma-1.0)*cregy*Grho/Gm
                vcr = np.sqrt(CRgamma*Pcr[cutz]/Grho[cutz])/km_in_cm
                #print 'CRgamma', CRgamma
                #print 'vcrest', np.sqrt(cregy/Gm)/km_in_cm
                vcrave = np.average(vcr)
                vcrl = np.append(vcrl, vcrave)
            if haveB>0 and needvA>0:
                vaave = np.average(valfven[cutz])
                print 'vaave', vaave
                val = np.append(val, vaave)
        #       val = np.append(val, np.median(valfven[cutz]))
        #       vaul = np.append(vaul, np.percentile(valfven[cutz],82))
        #       vadl = np.append(vadl, np.percentile(valfven[cutz],15))
        x = 0.5*(zver[1:]+zver[:-1])
        plt.plot(x, vthl, label=dclabel, color=color,ls=lsn,lw=1)
        #if havecr>0:
        #   plt.plot(x, vcrl, color=color,ls=lsn,lw=2)
        #plt.fill_between(x, vthul, vthdl,facecolor='r', alpha=0.3)
        if haveB>0 and needvA>0:
            plt.plot(x, val, color=color,ls='dotted',lw=1)
        #   plt.fill_between(x, vaul, vadl,facecolor='b', alpha=0.3)
        #plt.yscale('log')
        if userad==0:
            plt.xlabel(r'$z\;[{\rm kpc}]$',fontsize=18)
        else:
            plt.xlabel(r'$r\;[{\rm kpc}]$',fontsize=18)
                plt.ylabel(r'$c\;[{\rm km/s}]$',fontsize=18)
        if runtitle=='SMC':
            plt.legend(loc='best',fontsize=8,ncol=3)
                #plt.yscale('log')
    if userad==0:
        figname = 'CRplot/Alfven_sound/Alfven_sound_'+runtodo+'_sn'+str(Nsnap)+'.pdf'
    else:
        figname = 'CRplot/Alfven_sound/Alfven_sound_'+runtodo+'_sn'+str(Nsnap)+'_rad.pdf'
    plt.title(ptitle,fontsize=20)
    print 'figname', figname
    plt.tight_layout()
    plt.savefig(figname)
    plt.clf()










if wanted=='Smassxz' or wanted=='Smassxy':
        withinr=30
        maxlength=0.5
        needcontour=1
        rotface=1
    newstarneed=1
        for runtodo in dirneed:
                snaplist=[]
                presm = 0
                presnap = 0
                info=outdirname(runtodo, Nsnap)
                for i in [Nsnap]:
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        maindir=info['maindir']
                        haveB=info['haveB']
            snumadd=info['snumadd']
            usepep=info['usepep']
            halostr=info['halostr']
            firever = info['firever']
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
            print'usepep', usepep
                        cosmo=info['cosmo']
                        if cosmo==1:
                                h0=1
                        else:
                                h0=0
                        header = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo, header_only=1)
                        ascale = header['time']
                        print 'this time', ascale
                        thisred = 1./ascale-1.
                        hubble = header['hubble']
                        print 'hubble', hubble
                        #test
                        S = readsnapcr(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo)

                        if cosmo==1:
                if usepep==1:
                    halosingle = SF.read_halo_history_pep(rundir, Nsnap, singlesnap=1, firever=firever,halonostr=halostr, hubble=hubble, comoving=0, maindir=maindir)
                else:
                    halosingle = read_halo_history(rundir, halonostr=halostr,hubble=hubble, comoving=0, maindir=maindir, singlesnap=1, atime=ascale,snumadd=snumadd)
                                xcen = halosingle['x']
                                ycen = halosingle['y']
                                zcen = halosingle['z']
                                xvcen = halosingle['xv']
                                yvcen = halosingle['yv']
                                zvcen = halosingle['zv']
                                Rvirnow = halosingle['R']
                                MgAHF = halosingle['Mg']
                        else:
                                xcen=0
                                ycen=0
                                zcen=0
                                xvcen=0
                                yvcen=0
                                zvcen=0
                        Spos = S['p']
                        Svel = S['v']
                        Sx = Spos[:,0]-xcen
                        Sy = Spos[:,1]-ycen
                        Sz = Spos[:,2]-zcen
                        Sr = np.sqrt(Sx*Sx+Sy*Sy+Sz*Sz)
                        Svx = Svel[:,0]-xvcen   #km/s
                        Svy = Svel[:,1]-yvcen
                        Svz = Svel[:,2]-zvcen
                        Sm = S['m']*1e10 #Msun
            Sage = S['age']
            if cosmo==1 and newstarneed==1:
                readtimelist=readtime(firever=2)
                snap2list=readtimelist['snaplist']
                time2list=readtimelist['timelist']
                a2list=readtimelist['alist']
                Sage = np.interp(Sage,a2list,time2list)*1e9
                tnow = np.interp(ascale,a2list,time2list)*1e9
                print 'tnow', tnow
                print 'Sage', Sage
                tcut=Sage>tnow-1e7
                Sage = Sage[tcut]
                Sx=Sx[tcut]; Sy=Sy[tcut]; Sz=Sz[tcut];Svx=Svx[tcut];Svy=Svy[tcut];Svz=Svz[tcut];Sm=Sm[tcut];Sr=Sr[tcut];
                        if correctIa==1:
                                diskm = np.sum(Disk['m'])*1e10
                                bulgem = np.sum(Bulge['m'])*1e10
                        if rotface==1:
                                cutr = Sr < 5. #kpc
                                Smcutr = Sm[cutr]; #to avoid overflow
                                Sxcutr = Sx[cutr]; Sycutr = Sy[cutr]; Szcutr = Sz[cutr];
                                Svxcutr = Svx[cutr]; Svycutr = Svy[cutr]; Svzcutr = Svz[cutr];
                                Lang = [0.,0.,0.]
                                for i in range(len(Sxcutr)):
                                        Lang += Smcutr[i]*np.cross([Sxcutr[i],Sycutr[i],Szcutr[i]],[Svxcutr[i],Svycutr[i],Svzcutr[i]])
                                #test
                                Sx, Sy, Sz = SF.rotateL_to_z(Sx,Sy,Sz,Lang[0],Lang[1],Lang[2])
                                Svx, Svy, Svz = SF.rotateL_to_z(Svx,Svy,Svz,Lang[0],Lang[1],Lang[2])

                        Vxm = Svx*Sm
                        Vym = Svy*Sm
                        Vzm = Svz*Sm
                        V2m = np.sqrt(Svx*Svx+Svy*Svy+Svz*Svz)*Sm
                        cellvol = 2.0*withinr*2.0*withinr*2.0*maxlength/100./100.*kpc_in_cm*kpc_in_cm*kpc_in_cm

                if wanted=='Smassxy':
                        HmSm, xedges, yedges = np.histogram2d(Sy, Sx, bins=100,range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(Sm))
                        HmV2, xedges, yedges = np.histogram2d(Sy, Sx, bins=100,range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(V2m))
                        im=plt.imshow(np.log10(HmV2/HmSm)/2., interpolation='nearest', origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
                if wanted=='Smassxz':
                        HmSm, xedges, yedges = np.histogram2d(Sz, Sx, bins=100,range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(Sm))
                        HmV2, xedges, yedges = np.histogram2d(Sz, Sx, bins=100,range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(V2m))
                        HV2=np.nan_to_num(HmV2/HmSm)+0.1
                        if needcontour==1:
                                im=plt.contourf(xedges[:-1],yedges[:-1], np.log10(HV2)/2., 30, cmap='GnBu')
                        else:
                                im=plt.imshow(np.log10(HmV2/HmGm+1e-2)/2., interpolation='bicubic', cmap='Oranges', origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
                plt.colorbar(im,fraction=0.046, pad=0.04)
                if wanted=='Smassxy':
                        HmVx, xedges, yedges = np.histogram2d(Sy, Sx, bins=20,range=[[-withinr,withinr],[-withinr,withinr]], weights=Vxm)
                        HmVy, xedges, yedges = np.histogram2d(Sy, Sx, bins=20,range=[[-withinr,withinr],[-withinr,withinr]], weights=Vym)
                        HmVxy = np.sqrt(HmVx*HmVx+HmVy*HmVy)
                        U = HmVx/HmVxy
                        V = HmVy/HmVxy
                if wanted=='Smassxz':
                        HmVx, xedges, yedges = np.histogram2d(Sz, Sx, bins=20,range=[[-withinr,withinr],[-withinr,withinr]], weights=Vxm)
                        HmVy, xedges, yedges = np.histogram2d(Sz, Sx, bins=20,range=[[-withinr,withinr],[-withinr,withinr]], weights=Vzm)
                        HmSm, xedges, yedges = np.histogram2d(Sz, Sx, bins=20,range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(Sm))
                        HmVxy = np.sqrt(HmVx*HmVx+HmVy*HmVy)
                        U = HmVx/HmSm
                        V = HmVy/HmSm
                X = (xedges[:-1]+xedges[1:])/2.
                Y = (yedges[:-1]+yedges[1:])/2.
                Q = plt.quiver(X, Y, U, V, scale_units='x',scale=25)
        suffix=''
        if newstarneed==1:
            suffix = '_newstarneed'
                if wanted=='Smassxy' or wanted=='Smassxz':
                        qk = plt.quiverkey(Q, 0.1, 0.95, 100, '100 km/s', labelpos='E',
                                   coordinates='figure')
                if wanted=='Smassxy':
                        fname = 'CRplot/Smassxy/Smassxy'+runtodo+'_sn'+str(Nsnap)+suffix+'.pdf'
                        print 'fname', fname
                        plt.savefig(fname)
                        plt.clf()
                if wanted=='Smassxz':
                        fname='CRplot/Smassxz/Smassxz'+runtodo+'_sn'+str(Nsnap)+suffix+'.pdf'
                        print 'fname', fname
                        plt.savefig(fname)
                        plt.clf()









if wanted=='Bquiverxz' or wanted=='Bquiverxy' or wanted=='Vquiverxz' or wanted=='Vquiverxy':
    withinr=30
    maxlength=0.5
    needcontour=1
    rotface=1
        for runtodo in dirneed:
                snaplist=[]
                enclist=[]
                englist=[]
                enllist=[]
                endlist=[]
                enplist=[]
                enalist=[]
                avesfrl=[]
                prel=0
                preg=0
                prec=0
                pred=0
                prep=0
                presm = 0
                presnap = 0
        info=outdirname(runtodo, Nsnap)
        haveB=info['haveB']
        if haveB<1 and (wanted=='Bquiverxz' or wanted=='Bquiverxy'):
            continue
                for i in [Nsnap]:
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        maindir=info['maindir']
                        haveB=info['haveB']
                        snumadd=info['snumadd']
                        usepep=info['usepep']
                        halostr=info['halostr']
                        firever = info['firever']

                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        cosmo=info['cosmo']
                        if cosmo==1:
                                h0=1
                        else:
                                h0=0
                        header = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo, header_only=1)
                        ascale = header['time']
                        print 'this time', ascale
                        thisred = 1./ascale-1.
                        hubble = header['hubble']
                        print 'hubble', hubble
            if cosmo==0:
                datasup = 1
            else:
                datasup = 0
            #test
            #G = readsnapcr(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo)

                        Gextra = readsnapwcen(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix,\
                         havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
                         datasup=datasup,runtodo=runtodo,rundir=rundir,halostr=halostr)
                        Gx = Gextra['x']; Gy = Gextra['y']; Gz = Gextra['z'];
                        Gvx = Gextra['vx']; Gvy = Gextra['vy']; Gvz = Gextra['vz'];
                        Grho = Gextra['rho']; Gu = Gextra['u']; Gm = Gextra['m']



                        if wanted == 'Bquiverxy' or wanted == 'Vquiverxy' :
                                cutz = np.absolute(Gz)<maxlength
                        elif wanted == 'Bquiverxz' or wanted == 'Vquiverxz':
                                cutz = np.absolute(Gy)<maxlength
                        if haveB>0 and (wanted == 'Bquiverxy' or wanted == 'Bquiverxz'):
                                Bfield = Gextra['B']
                                Bx = Bfield[:,0]
                                By = Bfield[:,1]
                                Bz = Bfield[:,2]
                if rotface==1:
                    Bx, By, Bz = SF.rotateL_to_z(Bx,By,Bz,Lang[0],Lang[1],Lang[2])
                        Gx=Gx[cutz]
                        Gy=Gy[cutz]
                        Gz=Gz[cutz]
                        Gvx=Gvx[cutz]
                        Gvy=Gvy[cutz]
                        Gvz=Gvz[cutz]
                        Gm=Gm[cutz]
                        Grho=Grho[cutz]
                        if haveB>0 and (wanted == 'Bquiverxy' or wanted == 'Bquiverxz'):
                                Bx = Bx[cutz]
                                By = By[cutz]
                                Bz = Bz[cutz]

                        if haveB>0 and (wanted == 'Bquiverxy' or wanted == 'Bquiverxz'):
                                B2 = Bx*Bx+By*By+Bz*Bz
                                Begy = B2/8./np.pi*Gm/Grho*kpc_in_cm*kpc_in_cm*kpc_in_cm
                Bxvol = Bx*Gm/Grho*kpc_in_cm*kpc_in_cm*kpc_in_cm
                Byvol = By*Gm/Grho*kpc_in_cm*kpc_in_cm*kpc_in_cm
                                Bzvol = Bz*Gm/Grho*kpc_in_cm*kpc_in_cm*kpc_in_cm
            Vxm = Gvx*Gm
            Vym = Gvy*Gm
            Vzm = Gvz*Gm
            V2m = np.sqrt(Gvx*Gvx+Gvy*Gvy+Gvz*Gvz)*Gm
            cellvol = 2.0*withinr*2.0*withinr*2.0*maxlength/100./100.*kpc_in_cm*kpc_in_cm*kpc_in_cm

        if wanted=='Vquiverxy':
            HmGm, xedges, yedges = np.histogram2d(Gy, Gx, bins=100,range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(Gm))
            HmV2, xedges, yedges = np.histogram2d(Gy, Gx, bins=100,range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(V2m))
            im=plt.imshow(np.log10(HmV2/HmGm)/2., interpolation='nearest', origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
        if wanted=='Vquiverxz':
            HmGm, xedges, yedges = np.histogram2d(Gz, Gx, bins=100,range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(Gm))
            HmV2, xedges, yedges = np.histogram2d(Gz, Gx, bins=100,range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(V2m))
            
            HV2=np.nan_to_num(HmV2/HmGm)+0.1
                        if needcontour==1:
                                im=plt.contourf(xedges[:-1],yedges[:-1], np.log10(HV2)/2., 30, cmap='GnBu')
                        else:
                im=plt.imshow(np.log10(HmV2/HmGm+1e-2)/2., interpolation='bicubic', cmap='Oranges', origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
        if wanted=='Bquiverxy':
            HmB2, xedges, yedges = np.histogram2d(Gy, Gx, bins=100,range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(Begy/cellvol))
            im=plt.imshow(np.log10(HmB2)/2., interpolation='nearest', origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
        if wanted=='Bquiverxz':
            HmB2, xedges, yedges = np.histogram2d(Gz, Gx, bins=100,range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(Begy/cellvol))
            HmB2+=1e-18
            if needcontour==1:
                im=plt.contourf(xedges[:-1],yedges[:-1], np.log10(HmB2)/2., 20, cmap='GnBu',vmin=-10,vmax=-5.5)
            else:
                im=plt.imshow(np.log10(HmB2)/2., interpolation='bicubic', cmap='GnBu', origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
            
        plt.colorbar(im,fraction=0.046, pad=0.04)

                if wanted=='Vquiverxy':
                        HmVx, xedges, yedges = np.histogram2d(Gy, Gx, bins=20,range=[[-withinr,withinr],[-withinr,withinr]], weights=Vxm)
                        HmVy, xedges, yedges = np.histogram2d(Gy, Gx, bins=20,range=[[-withinr,withinr],[-withinr,withinr]], weights=Vym)
                        HmVxy = np.sqrt(HmVx*HmVx+HmVy*HmVy)
                        U = HmVx/HmVxy
                        V = HmVy/HmVxy
                if wanted=='Vquiverxz':
                        HmVx, xedges, yedges = np.histogram2d(Gz, Gx, bins=20,range=[[-withinr,withinr],[-withinr,withinr]], weights=Vxm)
                        HmVy, xedges, yedges = np.histogram2d(Gz, Gx, bins=20,range=[[-withinr,withinr],[-withinr,withinr]], weights=Vzm)
            HmGm, xedges, yedges = np.histogram2d(Gz, Gx, bins=20,range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(Gm))
                        HmVxy = np.sqrt(HmVx*HmVx+HmVy*HmVy)
                        U = HmVx/HmGm
                        V = HmVy/HmGm
                if wanted=='Bquiverxy':
                        HmBx, xedges, yedges = np.histogram2d(Gy, Gx, bins=20,range=[[-withinr,withinr],[-withinr,withinr]], weights=Bxvol/cellvol)
                        HmBy, xedges, yedges = np.histogram2d(Gy, Gx, bins=20,range=[[-withinr,withinr],[-withinr,withinr]], weights=Byvol/cellvol)
                        HmBxy = np.sqrt(HmBx*HmBx+HmBy*HmBy)
                        U = HmBx/1e6
                        V = HmBy/1e6
                if wanted=='Bquiverxz':
                        HmBx, xedges, yedges = np.histogram2d(Gz, Gx, bins=20,range=[[-withinr,withinr],[-withinr,withinr]], weights=Bxvol/cellvol)
                        HmBy, xedges, yedges = np.histogram2d(Gz, Gx, bins=20,range=[[-withinr,withinr],[-withinr,withinr]], weights=Bzvol/cellvol)
                        HmBxy = np.sqrt(HmBx*HmBx+HmBy*HmBy)
                        U = HmBx/1e6
                        V = HmBy/1e6
                X = (xedges[:-1]+xedges[1:])/2.
                Y = (yedges[:-1]+yedges[1:])/2.
                Q = plt.quiver(X, Y, U, V, scale_units='x',scale=25)
                if wanted=='Bquiverxy' or wanted=='Bquiverxz':
                        qk = plt.quiverkey(Q, 0.1, 0.95, 1e-6, r'$\mu {\rm G}$', labelpos='E',
                                   coordinates='figure')

                if wanted=='Vquiverxy' or wanted=='Vquiverxz':
                        qk = plt.quiverkey(Q, 0.1, 0.95, 100, '100 km/s', labelpos='E',
                                   coordinates='figure')


        if wanted=='Vquiverxy':
            fname = 'CRplot/Vquiverxy/Vquiverxy'+runtodo+'_sn'+str(Nsnap)+'.pdf'
            print 'fname', fname
            plt.savefig(fname)
            plt.clf()
        if wanted=='Vquiverxz':
            fname='CRplot/Vquiverxz/Vquiverxz'+runtodo+'_sn'+str(Nsnap)+'.pdf'
                        print 'fname', fname
                        plt.savefig(fname)
            plt.clf()
        if wanted=='Bquiverxy':
            fname='CRplot/Bquiverxy/Bquiverxy'+runtodo+'_sn'+str(Nsnap)+'.pdf'
                        print 'fname', fname
                        plt.savefig(fname)
            plt.clf()
        if wanted=='Bquiverxz':
            fname='CRplot/Bquiverxz/Bquiverxz'+runtodo+'_sn'+str(Nsnap)+'.pdf'
                        print 'fname', fname
                        plt.savefig(fname)
            plt.clf()




if wanted=='rhoT' or wanted=='rhoTwind':
    rcParams['figure.figsize'] = 5, 4
    if wanted=='rhoTwind':
        vcut=0.0 #outflow velocity cut
        withinr=20.0
        zup=25.0
        zdown=15.0
        extent = [-7,-1, 1, 9]
    else:
        extent = [-7,4, 1, 9]
    nooftimes=0
    nobin=51
    needcontour=1
        for runtodo in dirneed:
        Hadd = np.zeros((nobin-1,nobin-1))
        for i in range(startno,Nsnap,snapsep):
            info=outdirname(runtodo, i)
            rundir=info['rundir']
            maindir=info['maindir']
            halostr=info['halostr']
            runtitle=info['runtitle']
            slabel=info['slabel']
            snlabel=info['snlabel']
            dclabel=info['dclabel']
            resolabel=info['resolabel']
            the_snapdir=info['the_snapdir']
            Nsnapstring=info['Nsnapstring']
            havecr=info['havecr']
            Fcal=info['Fcal']
            iavesfr=info['iavesfr']
            timestep=info['timestep']
            color=info['color']
            cosmo=info['cosmo']
            usepep=info['usepep']
            if wanted=='rhoTwind':
                vmax=8.5
                vmin=0.0
                if runtitle=='MW':
                    vmax = 6.0
                    vmin = -1.0
                if runtitle=='SMC':
                    vmax = 5.0
                    vmin = 0.5 
                if runtitle=='SBC':
                    vmax = 6.5
                                        vmin = 0.5
            if wanted=='rhoT':
                                if runtitle=='MW':
                                        vmax = 9.0
                                        vmin = -2.0
                                if runtitle=='SMC':
                                        vmax = 8.0
                                        vmin = 0.5
            if cosmo==1:
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                header=readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1,header_only=1)
                h0 = header['hubble']
                atime = header['time']
                if usepep==1:
                    halosA = read_halo_history_pep(rundir, finalno, beginno=beginno,\
                     singlesnap=0, firever=firever,halonostr=halostr, hubble=h0, comoving=1, maindir=maindir)
                    afactor=atime
                else:
                    halosA = read_halo_history(rundir, maindir=maindir, halonostr=halostr, hubble=h0, comoving=0)
                    afactor=1.0
                redlist = halosA['redshift']
                haloid = halosA['ID']
                a_scale = 1.0/(1.0+redlist)
                xcenl = halosA['x']*afactor
                ycenl = halosA['y']*afactor
                zcenl = halosA['z']*afactor
                xvl = halosA['xv']
                yvl = halosA['yv']
                zvl = halosA['zv']
                Rvirl = halosA['R']*afactor
                xcen = np.interp(atime,a_scale,xcenl)
                ycen = np.interp(atime,a_scale,ycenl)
                zcen = np.interp(atime,a_scale,zcenl)
                                xvcen = np.interp(atime,a_scale,xvl)
                                yvcen = np.interp(atime,a_scale,yvl)
                                zvcen = np.interp(atime,a_scale,zvl)
                Rvir = np.interp(atime,a_scale,Rvirl)
            else:
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                xcen =0; ycen=0; zcen=0; xvcen=0; yvcen=0; zvcen=0;
                h0=1

            Gpos = G['p'][:,:]
            Gvel = G['v'][:,:]
            Gmass = G['m'][:]*1e10
            print 'np.sum(Gmass)', np.sum(Gmass)
            print 'np.amax(Gmass)', np.amax(Gmass)
            print 'np.amin(Gmass)', np.amin(Gmass)
            Tb = G['u'][:]
            rho = G['rho'][:]
            Neb = G['ne'][:]
            partX = Gpos[:,0]-xcen
            partY = Gpos[:,1]-ycen
            partZ = Gpos[:,2]-zcen
            partXV = Gvel[:,0]-xvcen
            partYV = Gvel[:,1]-yvcen
            partZV = Gvel[:,2]-zvcen
            header = G['header'][:]
            redshift = header[3]
            boxsize = header[9]
            if wanted == 'rhoTwind':
                cutz = (np.absolute(partZ)>zdown) & (np.absolute(partZ)<zup) #kpc
                cutxy = partX*partX+partY*partY<withinr*withinr
                if cosmo==1:
                    cutv = partXV*partX+partYV*partY\
                    +partZV*partZ>vcut*np.sqrt(partX*partX+partY*partY+partZ*partZ) #outflowing gas
                else:
                    cutv = partZV*partZ/np.absolute(partZ)>vcut #outflowing gas
                cut = cutz*cutxy*cutv
                Tb = Tb[cut]
                rho = rho[cut]
                Neb = Neb[cut]
                Gmass = Gmass[cut]

            TrueTemp, converted_rho  = SF.convertTemp(Tb, Neb, rho)
            if wanted == 'rhoT':
                if needcontour==1:
                                        totalname = 'CRplot/rhoT/rhoT_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+'_contour.pdf'
                else:
                    totalname = 'CRplot/rhoT/rhoT_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            elif wanted == 'rhoTwind':
                                if needcontour==1:
                    totalname = 'CRplot/rhoTwind/rhoTwind_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+'_contour.pdf'
                else:
                    totalname = 'CRplot/rhoTwind/rhoTwind_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
            y = np.log10(TrueTemp)
            x = np.log10(converted_rho)
                        gridxx = np.linspace(extent[0],extent[1],nobin)
                        gridyy = np.linspace(extent[2],extent[3],nobin)
            #gridxx = np.linspace(np.amin(x),np.amax(x),nobin)
            #gridyy = np.linspace(np.amin(y),np.amax(y),nobin)
            H, xedges, yedges = np.histogram2d(x, y, bins=[gridxx, gridyy],weights=Gmass)
            H=H.T
            Hadd += H
            nooftimes += 1 
        Hadd=Hadd/nooftimes
        if needcontour==1:
            levels = np.linspace(vmin,vmax,num=15)
            plt.contourf(xedges[:-1], yedges[:-1], np.log10(Hadd), levels=levels, extent=extent, origin='lower',aspect='auto')
        else:
            plt.imshow(np.log10(Hadd), extent=extent, interpolation='nearest',origin='lower',aspect='auto',vmax=vmax,vmin=vmin)
        cbar = plt.colorbar(extend='both', norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax))
        cbar.set_label(r'$\mathrm{Log}_{10} (M {\rm [M_\odot]})$', rotation=270, labelpad=20)
        plt.xlabel(r"$\mathrm{Log}_{10} (n {\rm [cm^{-3}]})$")
        plt.ylabel(r"$\mathrm{Log}_{10} (T {\rm [K]})$")
        plt.title(titleneed)
        plt.tight_layout()
        print totalname
        plt.savefig(totalname)
        plt.clf()



if wanted=='gassurden':
    withinr=3.0
    for runtodo in dirneed:
        info=outdirname(runtodo, Nsnap)
        rundir=info['rundir']
        maindir=info['maindir']
        halostr=info['halostr']
        runtitle=info['runtitle']
        slabel=info['slabel']
        snlabel=info['snlabel']
        dclabel=info['dclabel']
        resolabel=info['resolabel']
        the_snapdir=info['the_snapdir']
        Nsnapstring=info['Nsnapstring']
        havecr=info['havecr']
        Fcal=info['Fcal']
        iavesfr=info['iavesfr']
        timestep=info['timestep']
        color=info['color']
        cosmo=info['cosmo']
        usepep=info['usepep']
        haveB=info['haveB']
        if cosmo==1:
            G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
            header=readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1,header_only=1)
            h0 = header['hubble']
            atime = header['time']
            if usepep==1:
                halosA = read_halo_history_pep(rundir, finalno,\
                 beginno=beginno, singlesnap=0, firever=firever,halonostr=halostr, hubble=h0, comoving=1, maindir=maindir)
                afactor=atime
            else:
                halosA = read_halo_history(rundir, maindir=maindir, halonostr=halostr, hubble=h0, comoving=0)
                                afactor=1.0
            redlist = halosA['redshift']
            haloid = halosA['ID']
            a_scale = 1.0/(1.0+redlist)
            xcenl = halosA['x']*afactor
            ycenl = halosA['y']*afactor
            zcenl = halosA['z']*afactor
            Rvirl = halosA['R']*afactor
            xcen = np.interp(atime,a_scale,xcenl)
            ycen = np.interp(atime,a_scale,ycenl)
            zcen = np.interp(atime,a_scale,zcenl)
            Rvir = np.interp(atime,a_scale,Rvirl)
        else:
            G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
            xcen =0; ycen=0; zcen=0;
        Gpos = G['p']
                Gx = Gpos[:,0]-xcen
                Gy = Gpos[:,1]-ycen
                Gz = Gpos[:,2]-zcen
        print 'np.amax(Gx),np.amin(Gx)', np.amax(Gx),np.amin(Gx)
                print 'np.amax(Gy),np.amin(Gy)', np.amax(Gy),np.amin(Gy)
                print 'np.amax(Gz),np.amin(Gz)', np.amax(Gz),np.amin(Gz)
                Gm =np.array(G['m'])
                Gu = np.array(G['u'])
                rho = np.array(G['rho'])
                Gr = np.sqrt(Gx*Gx+Gy*Gy)
        #Gr = np.sqrt(Gx*Gx+Gz*Gz)
        dr = withinr/nogrid
        gassurdenl=[]
        rneed=[]
        if haveB>0:
            lsn='dashed'
        else:
            lsn='solid'
        for i in range(nogrid):
            rcut = Gr<dr*(i+1)
            Gmwithin = np.sum(Gm[rcut])*1e10*Msun_in_g
            gassurden = Gmwithin/(np.pi*dr*(i+1)*dr*(i+1)*kpc_in_cm*kpc_in_cm)
            gassurdenl = np.append(gassurdenl, gassurden)
            rneed=np.append(rneed,dr*(i+1))
        plt.plot(rneed, gassurdenl,label=dclabel,color=color,lw=2,ls=lsn)
                plt.xlabel(r'${\rm r [kpc]}$',fontsize=18)
                plt.ylabel(r'${\rm gas\; surface\; density\;(<r) [ g/cm^2}]$',fontsize=18)
                plt.legend(loc='best')
                plt.yscale('log')
    plt.tight_layout()
    plt.savefig('CRplot/gassurden/gassurden_'+fmeat+'_sn'+str(Nsnap)+'.pdf')
    plt.clf()



if wanted=='gassurdentime':
    for runtodo in dirneed:
        tneed=[]
        gassurdenl=[]
        withinr = 1.0
        for i in range(startno, Nsnap, snapsep):
            info=outdirname(runtodo, i)
            rundir=info['rundir']
            maindir=info['maindir']
            halostr=info['halostr']
            runtitle=info['runtitle']
            slabel=info['slabel']
            snlabel=info['snlabel']
            dclabel=info['dclabel']
            resolabel=info['resolabel']
            the_snapdir=info['the_snapdir']
            Nsnapstring=info['Nsnapstring']
            havecr=info['havecr']
            Fcal=info['Fcal']
            iavesfr=info['iavesfr']
            timestep=info['timestep']
            color=info['color']
            cosmo=info['cosmo']
            usepep=info['usepep']
            haveB=info['haveB']
            Rvir=info['Rvir']
            if cosmo==1:
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                header=readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1,header_only=1)
                h0 = header['hubble']
                atime = header['time']
                if usepep==1:
                    halosA = read_halo_history_pep(rundir, finalno, beginno=beginno,\
                     singlesnap=0, firever=firever,halonostr=halostr, hubble=h0, comoving=1, maindir=maindir)
                    afactor=atime
                else:
                    halosA = read_halo_history(rundir, maindir=maindir, halonostr=halostr, hubble=h0, comoving=0)
                    afactor=1.0
                redlist = halosA['redshift']
                haloid = halosA['ID']
                a_scale = 1.0/(1.0+redlist)
                xcenl = halosA['x']*afactor
                ycenl = halosA['y']*afactor
                zcenl = halosA['z']*afactor
                Rvirl = halosA['R']*afactor
                xcen = np.interp(atime,a_scale,xcenl)
                ycen = np.interp(atime,a_scale,ycenl)
                zcen = np.interp(atime,a_scale,zcenl)
                Rvir = np.interp(atime,a_scale,Rvirl)
            else:
                from crtestfunction import * 
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                xcen =0; ycen=0; zcen=0;
                Gpos = G['p']
                Gx = Gpos[:,0]
                Gy = Gpos[:,1]
                Gz = Gpos[:,2]
                Gm = G['m']
                                datasup=1
                                xcen = findcenz(runtodo,Nsnap,withinr=5.,dir='x',datasup=datasup,Gx=Gx,Gy=Gy,Gz=Gz,Gm=Gm)
                                ycen = findcenz(runtodo,Nsnap,withinr=5.,dir='y',datasup=datasup,Gx=Gx,Gy=Gy,Gz=Gz,Gm=Gm)
                                zcen = findcenz(runtodo,Nsnap,withinr=5.,dir='z',datasup=datasup,Gx=Gx,Gy=Gy,Gz=Gz,Gm=Gm)
            Gpos = G['p']
            Gx = Gpos[:,0]-xcen
            Gy = Gpos[:,1]-ycen
            Gz = Gpos[:,2]-zcen
            print 'np.amax(Gx),np.amin(Gx)', np.amax(Gx),np.amin(Gx)
            print 'np.amax(Gy),np.amin(Gy)', np.amax(Gy),np.amin(Gy)
            print 'np.amax(Gz),np.amin(Gz)', np.amax(Gz),np.amin(Gz)
            Gm =np.array(G['m'])
            #Gr = np.sqrt(Gx*Gx+Gy*Gy)
            #Rvir cut
            Grr = np.sqrt(Gx*Gx+Gy*Gy+Gz*Gz)
            rvcut = Grr<Rvir*0.25
            Gm = Gm[rvcut]
            Gx = Gx[rvcut]
            Gy = Gy[rvcut]
            Gz = Gz[rvcut]
            Gr = np.sqrt(Gx*Gx+Gz*Gz)
            if haveB>0:
                lsn='dashed'
            else:
                lsn='solid'
            rcut = Gr<withinr
            Gmwithin = np.sum(Gm[rcut])*1e10*Msun_in_g
            gassurden = Gmwithin/(np.pi*withinr*withinr*kpc_in_cm*kpc_in_cm)
            gassurdenl = np.append(gassurdenl, gassurden)
            tneed = np.append(tneed,i*0.98)
        #plt.plot(tneed, gassurdenl,label=dclabel,color=color,lw=2,ls=lsn)
        plt.plot(tneed, gassurdenl,label=dclabel,lw=2,ls=lsn)
        plt.xlabel('time [Myr]')
        plt.ylabel(r'gas surface density $(< 1 {\rm kpc}) [{\rm g/cm^2}]$')
        plt.legend(loc='best')
        #plt.yscale('log')
    plt.savefig('CRplot/gassurdentime/gassurdentime_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf')
    plt.clf()


            




if wanted =='testCRcum' or wanted=='testCRdis':
    rcParams['figure.figsize'] = 5,4    
    needcuty=0
    ycut = 0.001
    for runtodo in dirneed:
        info=outdirname(runtodo, Nsnap)
        rundir=info['rundir']
        maindir=info['maindir']
        halostr=info['halostr']
        runtitle=info['runtitle']
        slabel=info['slabel']
        snlabel=info['snlabel']
        dclabel=info['dclabel']
        resolabel=info['resolabel']
        the_snapdir=info['the_snapdir']
        Nsnapstring=info['Nsnapstring']
        havecr=info['havecr']
        Fcal=info['Fcal']
        iavesfr=info['iavesfr']
        timestep=info['timestep']
        cosmo=info['cosmo']
        color=info['color']
        withinRv=info['withinRv']
        usepep=info['usepep']
        beginno=info['beginno']
        finalno=info['finalno']
        firever=info['firever']
        initsnap=info['initsnap']
        kappa=info['kappa']
        havemetal=info['havemetal']
        print 'the_snapdir', the_snapdir
        G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,havemetal=havemetal)
        pos=G['p']
        vel=G['v']
        header =  readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr, header_only=1)
        timeneed = header['time']
        xg = pos[:,0]
        #rho=G['rho']
        Gx = pos[:,0]
        Gy = pos[:,1]
        Gz = pos[:,2]
        m =np.array(G['m'])
        Gu = np.array(G['u'])
        rho = np.array(G['rho'])
        cregy = np.array(G['cregy'])
        Gr = np.sqrt(Gx*Gx+Gy*Gy+Gz*Gz)
        cutr = Gr<0.01  #kpc 
        if needcuty==1:
            cuty = np.absolute(Gy)<ycut
            Gx = Gx[cuty]
            Gy = Gy[cuty]
            Gz = Gz[cuty]
            m = m[cuty]
            Gu = Gu[cuty]
            rho = rho[cuty]
            cregy = cregy[cuty]
            Gr = Gr[cuty]
        #derive an analytic solution of 3D diffusion equation with a point source:
        #ecr = Ecr0/(4pi*kappa*t)^(3/2)*Exp(-r^2/(4*kappa*t))
        Ecr0 = 1. #code unit, initial total CR energy (1e10Msun*(km/s)^2)
        #kappa = 1./2. #code unit diffusion coefficient (kpc^2/Gyr) 
        kappa = kappa
        print 'kappa', kappa
        #time = 0.1 #code unit time (Gyr)
        time=timeneed
        print 'timeneed', timeneed
        #rneed = 1. #code unit radius measured (kpc)
        #offset = 0.0389
        #offset = 0.05
        offset = 0.15
        kappahalf = kappa/2.
        kappatenth = kappa/10.
        fit = 1 #add an analytic solution to graph



        #calculate the numerical solution
        nrad=20
        arad=100
        rmax=3.
        dr =rmax/nrad
        dra = rmax/arad
        radl=[]
        if wanted=='testCRcum':
            cregyinl=[]
            for irad in range(nrad):
                cut = Gr < dr*(irad+0.5)
                cregycut = cregy[cut]
                #print 'cregycut', cregycut
                cregyin = np.sum(cregycut)
                cregyinl = np.append(cregyinl, cregyin)
                #print 'unit conversion constant', 1e10*Msun_in_g*km_in_cm*km_in_cm/kpc_in_cm/kpc_in_cm/kpc_in_cm*erg_in_eV
                radl = np.append(radl, dr*(irad+0.5))
                        #plt.plot(radl, cregyinl,ls='none',marker='o',label=dclabel)
            plt.plot(radl, cregyinl,ls='none',marker='o',label=dclabel)
        elif wanted== 'testCRdis':
            cregydenl=[]
            for irad in range(nrad):
                cut = (Gr > dr*irad) & (Gr < dr*(irad+1))
                cregycut = cregy[cut]
                #print 'cregycut', cregycut
                shellvol_in_kpc3 = 4.0/3.0*np.pi*(-np.power(dr*irad,3)+np.power(dr*(irad+1),3))
                cregyden = np.sum(cregycut)/shellvol_in_kpc3
                print 'cregyden', cregyden
                print 'np.sum(cregycut)', np.sum(cregycut)
                cregydenl = np.append(cregydenl, cregyden)
                #print 'unit conversion constant', 1e10*Msun_in_g*km_in_cm*km_in_cm/kpc_in_cm/kpc_in_cm/kpc_in_cm*erg_in_eV
                radl = np.append(radl, dr*(irad+0.5))
            plt.plot(radl, cregydenl,ls='none',marker='o',label=dclabel)

    if wanted=='testCRcum':
        if fit==1:
            rneed=radl
            #rneed=np.linspace(0.5*dra, 5.,num=arad)
            Ecrin = diffCRin(Ecr0,kappa,time,rneed,offset=offset)
            #plt.plot(rneed, Ecrin,label='analytic; kappa(code unit)='+str(kappa))
            plt.plot(rneed, Ecrin)
                        #Ecrin = diffCRin(Ecr0,kappa*0.9,time,rneed,offset=offset)
            #plt.plot(rneed, Ecrin,ls='dashed')
            #Ecrin = diffCRin(Ecr0,kappa*1.1,time,rneed,offset=offset)
                        #plt.plot(rneed, Ecrin,ls='dotted')
            #Ecrin = diffCRin(Ecr0,kappahalf,time,rneed,offset=offset)
            #plt.plot(rneed, Ecrin,label='analytic; kappa(code unit)='+str(kappahalf))
                        #Ecrin = diffCRin(Ecr0,kappatenth,time,rneed,offset=offset)
                        #plt.plot(rneed, Ecrin,label='analytic; kappa(code unit)='+str(kappatenth))
        #plt.yscale('log')
        plt.ylim([0.,1.])
        plt.xlim([0.,rmax])
        plt.xlabel(r'${\rm r\; [kpc]}$',fontsize=18)
        plt.ylabel(r'$E_{\rm CR}(<r)\; [{\rm code unit}]$',fontsize=18)
        plt.legend(loc='best')
        #plt.yscale('log')
        #plt.title('CR energy within sphere '+'time = '+str(time*1e3)+' Myr')
        plt.tight_layout()
        plt.savefig('CRplot/testCRcum/cregycum_sph_'+fmeat+'.pdf')
        plt.clf()
    elif wanted=='testCRdis':
                if fit==1:
                        rneed = np.linspace(0.5*dra, 5.,num=arad)
                        ecr_an = diffusionsol(Ecr0,kappa,time,rneed,offset=offset)
                        plt.plot(rneed, ecr_an,label='analytic; kappa(code unit)='+str(kappa))
                        print 'radl, cregydenl', radl, cregydenl
                        print 'rneed, ecr_an', rneed, ecr_an
        plt.xlim([0,rmax])
                plt.xlabel('r [kpc]')
                plt.ylabel(r'$e_{\rm CR} [{\rm code unit}]$')
                plt.legend(loc='best')
                #plt.yscale('log')
                #plt.title('CR density in each spherical shell '+'time ='+str(time*1e3)+' Myr')
                plt.savefig('CRplot/testCRdis/credensity_sph_'+fmeat+'.pdf')
                plt.clf()




if wanted=='printparmass':
    for runtodo in dirneed:
        for i in [0]:
            info=outdirname(runtodo, Nsnap)
            rundir=info['rundir']
            runtitle=info['runtitle']
            slabel=info['slabel']
            snlabel=info['snlabel']
            dclabel=info['dclabel']
            resolabel=info['resolabel']
            the_snapdir=info['the_snapdir']
            Nsnapstring=info['Nsnapstring']
            havecr=info['havecr']
            Fcal=info['Fcal']
            iavesfr=info['iavesfr']
            timestep=info['timestep']
            DM = readsnapcr(the_snapdir, Nsnapstring, 1, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
            D = readsnapcr(the_snapdir, Nsnapstring, 3, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
            G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
            print 'dmmass', np.average(DM['m']) 
            print 'dmass', np.average(D['m'])
            print 'gasmass', np.average(G['m']) 
            print 'dmsummass', np.sum(DM['m'])
            print 'dsummass', np.sum(D['m'])
            print 'gassummass', np.sum(G['m'])

if wanted=='gasfrac':
    for runtodo in dirneed:
        gbflb = []
        gbfl = []
        Nsnapl = []
                for i in range(0,Nsnap,snapsep):
                        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr=outdirname(runtodo, i)
            havecr=0
            G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
            DM = readsnap(the_snapdir, Nsnapstring, 1, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
            S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
            disk = readsnap(the_snapdir, Nsnapstring, 2, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
            bulge = readsnap(the_snapdir, Nsnapstring, 3, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
            Gp = G['p']
            Grho = G['rho']
            Gu = G['u']
            Gm = G['m']
            diskm = disk['m']
            bulgem = bulge['m']
            try:
                Sm = S['m']
            except KeyError:
                Sm=[]
            DMm = DM['m']
            Gz = Gp[:,2]
            Gx = Gp[:,0]
            Gy = Gp[:,1]
            cutxy = Gx*Gx+Gy*Gy < withinr*withinr
            cutz = np.absolute(Gz) < maxlength
            cut = cutxy*cutz
            gasbarf = np.sum(Gm[cut])/(np.sum(Gm[cut])+np.sum(Sm)+np.sum(diskm)+np.sum(bulgem))
                        gasbarfbox = np.sum(Gm)/(np.sum(Gm)+np.sum(Sm)+np.sum(diskm)+np.sum(bulgem))
            Nsnapl = np.append(Nsnapl, i)
            gbfl = np.append(gbfl,gasbarf)
            gbflb = np.append(gbflb, gasbarfbox)
            #cut = Gx*Gx+Gy*Gy < withinr*withinr
        #plt.plot(Nsnapl*0.001, gbfl)
        if wholeboxgas==1:
            plt.plot(Nsnapl, gbflb, label=runtodo+'_box')
        if diskgas==1:
            plt.plot(Nsnapl, gbfl, label=runtodo+'_disk')
    plt.xlabel('Myr')
    plt.ylabel('Gas fraction')
    plt.legend(loc='best')
    if wholeboxgas==1:
        plt.savefig('gasfraction_'+fmeat+'_box.pdf')
    else:
        plt.savefig('gasfraction_'+fmeat+'.pdf')
    plt.clf()   


if wanted=='outflow':
    for signal in [20,8,16,5]:
        for runtodo in dirneed:
                        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr=outdirname(runtodo, Nsnap)
            spmname='outflow_'+runtodo+'.txt'
            spmfile=open(spmname,"r")
            spmfile.readline()
            dars = spmfile.readlines()

            snapno=[]
            outflow=[]
            outflow16=[]
            outflow8=[]
            outflow0_5=[]
            for line in dars:
                xsd = line.split()
                snapno = np.append(snapno, int(xsd[0]))
                outflow = np.append(outflow, float(xsd[1]))
                outflow16 = np.append(outflow16, float(xsd[2]))
                outflow8 = np.append(outflow8, float(xsd[3]))
                outflow0_5 = np.append(outflow0_5, float(xsd[4]))       
            spmfile.close()

            if havecr==0:
                needlabel = slabel
            else:
                needlabel = r' $\kappa_{\rm di}=$'+dclabel

            if signal==20:
                plt.plot(snapno*1e-3, outflow, label=needlabel)
                plt.title('0-20kpc')
            if signal==8:
                plt.plot(snapno*1e-3, outflow8, label=needlabel)
                plt.title('8-12kpc')
            if signal==16:
                plt.plot(snapno*1e-3, outflow16, label=needlabel)
                plt.title('16-24kpc')
            if signal==5:
                plt.plot(snapno*1e-3, outflow0_5, label=needlabel)
                plt.title('>0.5kpc and v > 1000 km/s')
            plt.xlabel('Gyr')
            plt.ylabel('outflow')
            plt.xlim(xmax=Nsnap*1e-3)
        plt.yscale('log')
        plt.legend(loc='best')
        if signal==20:
            plt.savefig('outflow_'+runtodo+'.pdf')
        if signal==8:
            plt.savefig('outflow8_'+runtodo+'.pdf')
        if signal==16:
            plt.savefig('outflow16_'+runtodo+'.pdf')
        if signal==5:
            plt.savefig('outflow0_5_'+runtodo+'.pdf')
        plt.clf()

if wanted=='sfr':
    for runtodo in dirneed:
        try:
            rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr=outdirname(runtodo, Nsnap)
            S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix)
            Sage = S['age']
        except KeyError:
            Nsnap-=1
            rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr=outdirname(runtodo, Nsnap)
                        S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix)
                        Sage = S['age']
        #aveSmass = np.average(S['m'])
        hist, bin_edges = np.histogram(Sage, bins=10, weights=S['m'])
        print 'bin_edges', bin_edges
        print 'hist', hist
        if havecr==0:
            needlabel = slabel
        else:
            needlabel = r' $\kappa_{\rm di}=$'+dclabel
        plt.plot(bin_edges[1:]*0.98, hist/(bin_edges[1]-bin_edges[0])/0.98*10, label=needlabel)
        plt.xlabel('Gyr')
        plt.ylabel(r'SFR ($M_\odot$/yr)')
    plt.yscale('log')
    plt.legend(loc='best')
    plt.title(runtitle)
    plt.savefig('starage_'+fmeat+'.pdf')
    plt.clf()

if wanted=='pregrad':
    for runtodo in dirneed:
        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr , Fcal, iavesfr =outdirname(runtodo, Nsnap)
        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
        Gp = G['p']
        Grho = G['rho']
        Gu = G['u']
        Gm = G['m']
        if (havecr==1):
            Gcregy=G['cregy']
        Gz = Gp[:,2]
        Gx = Gp[:,0]
        Gy = Gp[:,1]
        cut = Gx*Gx+Gy*Gy < withinr*withinr
        if (havecr==1):
            Gcregyc=Gcregy[cut]*2e53
        Gzcut=Gz[cut]
        Grhoc=Grho[cut]*7e-22
        Guc=Gu[cut]*1e10
        Gmc=Gm[cut]*1e10*2e33
        heightlist=[]
        pgradlist=[]
        Gprelist=[]
        gamma=1.5
        preGpre=0.0
        for i in range(nogrid):
            height = i*maxlength/nogrid-maxlength/2.0
            big=Gzcut>height
            small=Gzcut<height+maxlength/nogrid
            within=big*small
            Gpre0=(gamma-1.0)*Guc[within]*Grhoc[within]
            if (havecr==1):
                Gpre=Gpre0+(0.3333)*Gcregyc[within]*Grhoc[within]/Gmc[within]
                #print 'ratio', Gpre0/Gpre
            else:
                Gpre=Gpre0
            heightlist=np.append(heightlist,height)
            Gprelist=np.append(Gprelist,np.average(Gpre))
            #pgradlist=np.append(pgradlist,np.absolute(np.average(Gpre)-preGpre)*nogrid/maxlength/3e21)
            #preGpre=np.average(Gpre)
        #print 'Gprelist', len(Gprelist)
        #print 'heightlist', len(heightlist)
        #plt.plot(Gprelist,heightlist,label=slabel)
                if havecr==0:
                        needlabel = slabel
                else:
                        needlabel = r' $\kappa_{\rm di}=$'+dclabel
        plt.plot(np.absolute((Gprelist[:-1]-Gprelist[1:])/(heightlist[:-1]-heightlist[1:]))/3.09e21,(heightlist[1:]+heightlist[:-1])*0.5,label=needlabel)
    plt.legend(loc='best',fontsize=14)
    plt.title(runtitle)
    #plt.title('pressure gradient in a cylinder of radius '+str(withinr)+'kpc centered on the '+runtitle+' disk')
    #plt.xlim([-300,300])
    plt.ylim([-maxlength/2, maxlength/2])
    plt.ylabel('z (kpc)')
    plt.xlabel(r'${\rm dP/dz \; (g/cm^2/s^2)}$')
    plt.xscale('log')
    plt.savefig('pgrad_disk_'+runtodo+'_'+str(withinr)+'kpc.pdf')

if wanted=='pre':
        for runtodo in dirneed:
                rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr=outdirname(runtodo, Nsnap)
                G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                Gp = G['p']
                Grho = G['rho']
                Gu = G['u']
                Gm = G['m']
                if (havecr==1):
                        Gcregy=G['cregy']
                Gz = Gp[:,2]
                Gx = Gp[:,0]
                Gy = Gp[:,1]
                cut = Gx*Gx+Gy*Gy < withinr*withinr
                if (havecr==1):
                        Gcregyc=Gcregy[cut]*2e53
                Gzcut=Gz[cut]
                Grhoc=Grho[cut]*7e-22
                Guc=Gu[cut]*1e10
                Gmc=Gm[cut]*1e10*2e33
                heightlist=[]
                pgradlist=[]
                Gprelist=[]
                gamma=1.5
                preGpre=0.0
                for i in range(nogrid):
                        height = i*maxlength/nogrid-maxlength/2.0
                        big=Gzcut>height
                        small=Gzcut<height+maxlength/nogrid
                        within=big*small
                        Gpre0=(gamma-1.0)*Guc[within]*Grhoc[within]
                        if (havecr==1):
                                Gpre=Gpre0+(0.3333)*Gcregyc[within]*Grhoc[within]/Gmc[within]
                        else:
                                Gpre=Gpre0
                        heightlist=np.append(heightlist,height)
                        Gprelist=np.append(Gprelist,np.average(Gpre))
                        #pgradlist=np.append(pgradlist,np.absolute(np.average(Gpre)-preGpre)*nogrid/maxlength/3e21)
                        #preGpre=np.average(Gpre)
                #print 'Gprelist', len(Gprelist)
                #print 'heightlist', len(heightlist)
                #plt.plot(Gprelist,heightlist,label=slabel)
                if havecr==0:
                        needlabel = slabel
                else:
                        needlabel = r' $\kappa_{\rm di}=$'+dclabel
                plt.plot(Gprelist,heightlist,label=needlabel)
        plt.legend(loc='best')
        plt.title(runtitle)
        #plt.title('pressure gradient in a cylinder of radius '+str(withinr)+'kpc centered on the '+runtitle+' disk')
        #plt.xlim([-300,300])
        plt.ylim([-maxlength/2, maxlength/2])
        plt.ylabel('z (kpc)')
        plt.xlabel(r'${\rm P \; (g/cm/s^2)}$')
        plt.xscale('log')
        plt.savefig('pre_disk_'+runtodo+'_'+str(withinr)+'kpc.pdf')


if wanted=='zv':
    for runtodo in dirneed:
        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr =outdirname(runtodo, Nsnap)
        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix)
        Gp = G['p']
        Gv = G['v']
        Gm = G['m']
        Gvz = Gv[:,2]
        Gz = Gp[:,2]
        Gx = Gp[:,0]
        Gy = Gp[:,1]
        h = 1

        #find particles within 3kpc cylinder in z plane:
        cut = Gx*Gx+Gy*Gy < withinr*withinr
        Gzcut=Gz[cut]
        Gmc=Gm[cut]
        Gvzc=Gvz[cut]
        heightlist=[]
        medvlist=[]

        for i in range(nogrid):
            height = i*maxlength/nogrid-maxlength/2.0
            big=Gzcut>height
            small=Gzcut<height+maxlength/nogrid
            within=big*small
            Gmcwi=Gmc[within]
            Gvzcwi=Gvzc[within]
            #print 'Gvzcwi', Gvzcwi
            heightlist=np.append(heightlist,height)
            medvlist=np.append(medvlist,np.median(Gvzcwi))
                if havecr==0:
                        needlabel = slabel
                else:
                        needlabel = r' $\kappa_{\rm di}=$'+dclabel
        plt.plot(medvlist,heightlist,label=needlabel)
    plt.legend(loc='best')
    #plt.title('Median Vz in a cylinder of radius 10kpc centered on the galactic disk')
    plt.title(runtitle)
    plt.xlim([-400,400])
    plt.ylim([-maxlength/2, maxlength/2])
    plt.ylabel('z (kpc)')
    plt.xlabel('Vz (km/s)')
    plt.savefig('vzdisk_'+runtodo+'_'+str(withinr)+'kpc.pdf')



if wanted=='massloading':
    for signal in [20,8,16, 5]:
        for runtodo in dirneed:
            rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr=outdirname(runtodo, Nsnap)
            spmname='outflow_'+runtodo+'.txt'
            spmfile=open(spmname,"r")
            spmfile.readline()
            dars = spmfile.readlines()

            snapno=[]
            outflow=[]
            outflow16=[]
            outflow8=[]
            outflow0_5=[]
            smlist=[]

            for line in dars:
                xsd = line.split()
                snapno = np.append(snapno, int(xsd[0]))
                outflow = np.append(outflow, float(xsd[1]))
                outflow16 = np.append(outflow16, float(xsd[2]))
                outflow8 = np.append(outflow8, float(xsd[3]))
                outflow0_5 = np.append(outflow0_5, float(xsd[4]))
                smlist = np.append(smlist, float(xsd[5]))
            spmfile.close()
            sfr = (smlist[1:]-smlist[:-1])/(snapno[1:]-snapno[:-1])/0.98*10
            print 'sfr', sfr
            massloading=outflow[:-1]/sfr
            massloading8 = outflow8[:-1]/sfr
            massloading16 = outflow16[:-1]/sfr
            massloading0_5 = outflow0_5[:-1]/sfr
            if havecr==0:
                needlabel = slabel
            else:
                needlabel = r' $\kappa_{\rm di}=$'+dclabel

            if signal==20:
                plt.plot((snapno[:-1])*0.98, massloading, label=needlabel)
                plt.title(runtitle+' 0-20kpc')
            if signal==8:
                plt.plot((snapno[:-1])*0.98, massloading8, label=needlabel)
                plt.title(runtitle+' 8-12kpc')
            if signal==16:
                plt.plot((snapno[:-1])*0.98, massloading16, label=needlabel)
                plt.title(runtitle+' 16-24kpc')
            if signal==5:
                plt.plot((snapno[:-1])*0.98, massloading0_5, label=needlabel)
                plt.title('>0.5kpc and v > 100 km/s')
            plt.xlabel('Gyr')
            plt.ylabel(r'$\eta$', fontsize=24)
        plt.yscale('log')
        plt.legend(loc='best')
        if signal==20:
            plt.savefig('massloading_'+runtodo+'.pdf')
        if signal==8:
            plt.savefig('massloading8_'+runtodo+'.pdf')
        if signal==16:
            plt.savefig('massloading16_'+runtodo+'.pdf')
        if signal==5:
            plt.savefig('massloading0_5_'+runtodo+'.pdf')
        plt.clf()



if wanted=='smchange':
    for runtodo in dirneed:
        info=outdirname(runtodo, Nsnap)
        rundir=info['rundir']
        runtitle=info['runtitle']
        slabel=info['slabel']
        snlabel=info['snlabel']
        dclabel=info['dclabel']
        resolabel=info['resolabel']
        the_snapdir=info['the_snapdir']
        Nsnapstring=info['Nsnapstring']
        havecr=info['havecr']
        Fcal=info['Fcal']
        iavesfr=info['iavesfr']
        timestep=info['timestep']
        spmname='outflow_'+runtodo+'.txt'
        spmfile=open(spmname,"r")
        spmfile.readline()
        dars = spmfile.readlines()

        snapno=[]
        smlist=[]
        if havecr==0:
            needlabel = slabel
        else:
            needlabel = r' $\kappa_{\rm di}=$'+dclabel

        for line in dars:
            xsd = line.split()
            snapno = np.append(snapno, int(xsd[0]))
            smlist = np.append(smlist, float(xsd[5]))
        spmfile.close()
                smlists=smlist[::snapsep]
                snapnos=snapno[::snapsep]
        sfr = (smlists[1:]-smlists[:-1])/(snapnos[1:]-snapnos[:-1])/0.98*1e4
        print 'sfr', sfr
        #plt.plot((snapno[1:]+snapno[:-1])/2.*0.98, sfr, label=needlabel)
                #plt.plot((snapnos[1:]+snapnos[:-1])/2.*0.98, sfr, label='kappa='+dclabel+'; '+'time step= '+timestep)
        #plt.plot((snapnos[1:]+snapnos[:-1])/2.*0.98, sfr, label='kappa='+dclabel)
        plt.plot((snapnos[1:]+snapnos[:-1])/2.*0.98, sfr, label=runtodo)
    plt.legend(loc='best')
    plt.title(runtitle)
    plt.xlim(xmin=float(startno),xmax=float(Nsnap))
    plt.yscale('log')
    plt.ylabel(r'SFR $({\rm M}_{\odot}/{\rm yr})$')
    plt.xlabel('Myr')
    plt.savefig('CRplot/sfr_'+fmeat+'.pdf')
    plt.clf()

if wanted=='smass':
        for runtodo in dirneed:
                rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr=outdirname(runtodo, Nsnap)
                spmname='outflow_'+runtodo+'.txt'
                spmfile=open(spmname,"r")
                spmfile.readline()
                dars = spmfile.readlines()

                snapno=[]
                smlist=[]
                if havecr==0:
                        needlabel = slabel
                else:
                        needlabel = r' $\kappa_{\rm di}=$'+dclabel

                for line in dars:
                        xsd = line.split()
                        snapno = np.append(snapno, int(xsd[0]))
                        smlist = np.append(smlist, float(xsd[5]))
                spmfile.close()
                #plt.plot((snapno[1:]+snapno[:-1])/2.*0.98, sfr, label=needlabel)
                plt.plot(snapno*0.98, smlist*1e10, label=runtodo)
        plt.legend(loc='best')
        plt.title(fmeat)
        plt.xlim(xmax=float(Nsnap))
        plt.yscale('log')
        plt.ylabel(r' $M_{*}({\rm M}_{\odot})$')
        plt.xlabel('Myr')
        plt.savefig('smass_'+fmeat+'.pdf')
        plt.clf()
    

if wanted=='KSLaw':
        for runtodo in dirneed:
                rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr=outdirname(runtodo, Nsnap)
                spmname='outflow_'+runtodo+'.txt'
                spmfile=open(spmname,"r")
                spmfile.readline()
                dars = spmfile.readlines()

                snapno=[]
                smlist=[]
        sm1kpc=[]
        gm1kpc=[]
                if havecr==0:
                        needlabel = slabel
                else:
                        needlabel = r' $\kappa_{\rm di}=$'+dclabel

                for line in dars:
                        xsd = line.split()
                        snapno = np.append(snapno, int(xsd[0]))
                        smlist = np.append(smlist, float(xsd[5]))
            sm1kpc = np.append(sm1kpc, float(xsd[6]))
            gm1kpc = np.append(gm1kpc, float(xsd[7]))
                spmfile.close()
                #sfr = (smlist[1:]-smlist[:-1])/(snapno[1:]-snapno[:-1])/0.98*1e4
        sfr1kpc = (sm1kpc[1:]-sm1kpc[:-1])/(snapno[1:]-snapno[:-1])/0.98*1e4
        #plt.plot(np.log10((gm1kpc[1:]+gm1kpc[:-1])/2.*1e10/np.pi/1e6),np.log10(sfr/np.pi), ls='none', marker='o')
        plt.plot(np.log10((gm1kpc[1:]+gm1kpc[:-1])/2.*1e10/np.pi/1e6),np.log10(sfr1kpc/np.pi), ls='none', marker='o', label=runtodo)
                #plt.plot((snapno[1:]+snapno[:-1])/2.*0.98, sfr, label=needlabel)
        #plt.plot((snapno[1:]+snapno[:-1])/2.*0.98, sfr1kpc, label='< 1kpc')
        plt.legend(loc='best')
        plt.title(fmeat)
        plt.ylabel(r'$\log (\sum_{\rm SFR} [{\rm M}_{\rm sun}/{\rm yr}/{\rm kpc^2}])$')
        plt.xlabel(r'$\log (\sum_{\rm gas} [{\rm M}_{\rm sun}/{\rm pc^2}])$')
        plt.savefig('KSlaw_1kpc_'+fmeat+'.pdf')
        plt.clf()


if wanted=='outflowwrite':
    for runtodo in dirneed:
        snaplist=[]
        outlist=[]
        outlist16=[]
        outlist8=[]
        outlist0_5=[]
        outlistpb=[]
        mslist=[]
        gm1kpc=[]
        sm1kpc=[]
            for i in range(0,Nsnap,snapsep):
            veccut=100
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
            print 'the_snapdir', the_snapdir
            print 'Nsnapstring', Nsnapstring
            G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix)
            Gp = G['p']
            Gv = G['v']
            Gm = G['m']
            Gvx = Gv[:,0]
            Gvy = Gv[:,1]
            Gvz = Gv[:,2]
            Gx = Gp[:,0]
            Gy = Gp[:,1]
            Gz = Gp[:,2]
            #find particles within 1kpc
            in1kpc = np.square(Gx)+np.square(Gy)+np.square(Gz)<1
            gasmass1kpc = np.sum(Gm[in1kpc])
            

            #find particles within 0-20, 16-24, 8-12kpc in z direction:
            up = Gz>0
            down = Gz<0
            within = np.absolute(Gz) < 20
            withinu = up*within
            withind = down*within

            up16 = Gz>16
            down16 = Gz<-16
            within16 = np.absolute(Gz) < 24
            withinu16 = up16*within16
            withind16 = down16*within16

            up8 = Gz>8
            down8 = Gz<-8
            within8 = np.absolute(Gz) < 12
            withinu8 = up8*within8
            withind8 = down8*within8
            #find particles above 500pc from the disk
            up0_5 = Gz> 0.5
            down0_5 = Gz <-0.5
            Gviu = Gvz[withinu]
            Gvid = Gvz[withind]
            Gmiu = Gm[withinu]
            Gmid = Gm[withind]
            print 'Gviu', Gviu
            upout = Gviu>veccut
            downout = Gvid<veccut

            Gviu16 = Gvz[withinu16]
            Gvid16 = Gvz[withind16]
            Gmiu16 = Gm[withinu16]
            Gmid16 = Gm[withind16]
            upout16 = Gviu16>veccut
            downout16 = Gvid16<veccut



            Gviu8 = Gvz[withinu8]
            Gvid8 = Gvz[withind8]
            Gmiu8 = Gm[withinu8]
            Gmid8 = Gm[withind8]
            upout8 = Gviu8>veccut
            downout8 = Gvid8<veccut

            Gvu0_5 = Gvz[up0_5]
            Gvd0_5 = Gvz[down0_5]
            Gmu0_5 = Gm[up0_5]
            Gmd0_5 = Gm[down0_5]
            upout0_5 = Gvu0_5>100
            downout0_5 = Gvd0_5<-100
            totmom = np.absolute(np.sum(Gviu[upout]*Gmiu[upout]))+np.absolute(np.sum(Gvid[downout]*Gmid[downout]))
            totmom *= 1e10 # convert to solar mass *km/s

            totmom16 = np.absolute(np.sum(Gviu16[upout16]*Gmiu16[upout16]))+np.absolute(np.sum(Gvid16[downout16]*Gmid16[downout16]))
            totmom16 *= 1e10 # convert to solar mass *km/s

            totmom8 = np.absolute(np.sum(Gviu8[upout8]*Gmiu8[upout8]))+np.absolute(np.sum(Gvid8[downout8]*Gmid8[downout8]))
            totmom8 *= 1e10 # convert to solar mass *km/s

            totmom0_5 = np.absolute(np.sum(Gvu0_5[upout0_5]*Gmu0_5[upout0_5]))+np.absolute(np.sum(Gvd0_5[downout0_5]*Gmd0_5[downout0_5]))
            totmom0_5 *= 1e10 # convert to solar mass *km/s

            #totmomb = np.absolute(np.sum(Gvpb*Gmpb))*1e10




            #divide it by dL to get dM/dt:
            outflowrate = totmom/40./3.08567758e16*3.15569e7 #convert kpc to km, and s to yr so the unit is solar mass/yr
            outflowrate16 = totmom16/16./3.08567758e16*3.15569e7
            outflowrate8 = totmom8/8./3.08567758e16*3.15569e7
            outflowrate0_5 = totmom0_5/3.08567758e16*3.15569e7
            #outflowratepb = totmompb/3.08567758e16*3.15569e7
            snaplist = np.append(snaplist,Nsnapstring)
            outlist = np.append(outlist,outflowrate)
            outlist16 = np.append(outlist16,outflowrate16)
            outlist8 = np.append(outlist8, outflowrate8)
            outlist0_5 = np.append(outlist0_5, outflowrate0_5)
            
            #outlistpb = np.append(outlistpb, outflowratepb)
                        S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix)
            try:
                Sm = S['m']
                ms = np.sum(Sm)
                Sp = S['p']
                Sx = Sp[:,0]
                Sy = Sp[:,1]
                Sz = Sp[:,2]
                in1kpc = np.square(Sx)+np.square(Sy)+np.square(Sz)<1
                starmass1kpc = np.sum(Sm[in1kpc])
                del Sm, Sp, Sx, Sy, Sz
            except KeyError:
                ms = 0.0
                starmass1kpc = 0
            mslist=np.append(mslist,ms)
            sm1kpc=np.append(sm1kpc, starmass1kpc)
            gm1kpc=np.append(gm1kpc, gasmass1kpc)
            del S, Gvu0_5, Gmu0_5, upout0_5,Gvd0_5,Gmd0_5,downout0_5,Gviu8,Gmiu8,Gvid8,Gmid8,upout8,downout8,Gviu16,Gmiu16,Gvid16,Gmid16,upout16,downout16,Gviu,Gmiu,upout,Gmid,downout,G,Gp,Gv,Gm,Gvx,Gvy,Gvz,Gx,Gy,Gz
        spmname='outflow_'+runtodo+'.txt'
        spmfile=open(spmname,"w")
        spmfile.write('snapshot no.   ' + 'outflow (<20kpc)   '+'outflow (16-24kpc)    '+'outflow (8-12kpc)   '+'outflow (obs)      '+'stellar mass      '+'stellar mass (1kpc)       ' + 'gas mass (1kpc)       '+'\n')
        for ncount in range(len(snaplist)):
            spmfile.write(str(snaplist[ncount])+'    '+str(outlist[ncount])+'    '+str(outlist16[ncount])+'    '+str(outlist8[ncount])+'      '+str(outlist0_5[ncount])+'        '+str(mslist[ncount])+'              '+str(sm1kpc[ncount])+'         '+str(gm1kpc[ncount])+'\n')
        spmfile.close()


if wanted=='onlygamma':
        for runtodo in dirneed:
                print 'runtodo', runtodo
                rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr,  Fcal, iavesfr=outdirname(runtodo, Nsnap)
                G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix,havecr=1)
                Gegy = G['u']
                Gm = G['m']
                Gp = G['p']
                Gr = np.sqrt(Gp[:,0]*Gp[:,0]+Gp[:,1]*Gp[:,1]+Gp[:,2]*Gp[:,2])
                Grho = G['rho'] # in 10^10 Msun/kpc^3
                Gcregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                Gnism = Grho*1e10*2e33/2.9e64/1.1*6e23 #in cm^-3 assume: mean molecular weight=1.1 consider only nuclei
                tpi = 2e5/Gnism*250.0 #pi decay time in yr
                betapi=0.7 ##should be between 0.5-0.9 from Lacki 2011
                Lgammagev = Gcregy*0.25/tpi*2e53/3.2e7*betapi/3.0 #in erg/s   in proton calorimetric limit
                totLgammagev =np.sum(Lgammagev)
                print 'np.sum(Gcregy) in erg', np.sum(Gcregy)*2e53
                print 'np.sum(Gegy*Gm)', np.sum(Gegy*Gm)
                print 'np.amax(Gnism)', np.amax(Gnism)
                hist, bin_edges = np.histogram(Gnism, density=True, weights=G['m'])
                print 'np.average(Gnism)', np.average(Gnism)
                print 'In proton calorimetric limit'
                print ' Total gamma luminosity in solar (in log10)', np.log10(totLgammagev/3.85e33)





if wanted=='gamma':
        for runtodo in dirneed:
        print 'runtodo', runtodo
                rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr,  Fcal, iavesfr=outdirname(runtodo, Nsnap)
        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix,havecr=1)
        S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix,havecr=1)
        Sage = S['age']
        Sm = S['m']
        Gegy = G['u']
        Gm = G['m']
        Gp = G['p']
        Gr = np.sqrt(Gp[:,0]*Gp[:,0]+Gp[:,1]*Gp[:,1]+Gp[:,2]*Gp[:,2])
        Grho = G['rho'] # in 10^10 Msun/kpc^3
        Gcregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
        Gnism = Grho*1e10*2e33/2.9e64/1.1*6e23 #in cm^-3 assume: mean molecular weight=1.1 consider only nuclei
        recentstar = Sage < 0.2
                aveSmass = np.average(S['m'])
                hist, bin_edges = np.histogram(Sage, bins=20)
                avesfr=hist[-1]*aveSmass/(bin_edges[1]-bin_edges[0])/0.98*10
        print 'avesfr', avesfr
        closer = Gr < 3.0
        print 'np.sum(Sm)', np.sum(Sm)
        tpi = 2e5/Gnism*250.0 #pi decay time in yr
        betapi=0.7 ##should be between 0.5-0.9 from Lacki 2011
        Lgammagev = Gcregy*0.25/tpi*2e53/3.2e7*betapi/3.0 #in erg/s   in proton calorimetric limit
        Lsfr =  3.8e-4 *avesfr*2e33/3.2e7*9e20
        totLgammagev =np.sum(Lgammagev)
        print 'estimated total supernova energy', np.sum(Sm)/0.4*1e10*0.0037*1e51
        print 'np.sum(Gcregy) in erg', np.sum(Gcregy)*2e53
        print 'np.sum(Gegy*Gm)', np.sum(Gegy*Gm)
        print 'np.amax(Gnism)', np.amax(Gnism)
        hist, bin_edges = np.histogram(Gnism, density=True, weights=G['m'])
        print 'np.average(Gnism)', np.average(Gnism)
        Lsfrideal = 3.8e-4 *iavesfr*2e33/3.2e7*9e20
        print 'Fsfr (ideal)', Lsfrideal/np.pi/4.0/(0.06*3.08e24)/(0.06*3.08e24)
        print 'In proton calorimetric limit'
        print ' Total gamma luminosity in solar (in log10)', np.log10(totLgammagev/3.85e33)
        print 'Fgamma/Fsfr', totLgammagev/Lsfr
        print 'with Fcal from LTQ fiducial model'
        print ' Total gamma luminosity in solar (in log10)', np.log10(totLgammagev*Fcal/3.85e33)
        print 'Fgamma/Fsfr', totLgammagev*Fcal/Lsfr
        print 'Fgamma/Fsfr (ideal)', totLgammagev*Fcal/Lsfrideal
        rneed=10
        fraction=0
        TotGcregy=np.sum(Gcregy)
        while (fraction>0.92) or (fraction<0.88):
            within = Gr<rneed
            fraction = np.sum(Gcregy[within])/TotGcregy
            rneed *= np.exp(1.0*(0.9-fraction))
            #print 'fraction, rneed', fraction, rneed
        hist, bin_edges = np.histogram(Gr,weights=Gcregy)
        plt.plot(bin_edges[:-1],hist, label=runtodo)
        plt.axvline(x=rneed)
    plt.ylabel(r'CR energy')
    plt.xlabel('kpc')
    plt.yscale('log')
    plt.xscale('log')
    plt.legend(loc='best')
    plt.savefig('histcr.pdf')
    plt.clf()


if wanted=='cosmic_raywrite':
        for runtodo in dirneed:
                snaplist=[]
                crlist=[]
        galist=[]
        smlist=[]
        aveGnisml=[]
        #for i in [1, 2]:
        #for i in [500]:
        print 'runtodo', runtodo
                for i in range(startno,Nsnap):
            info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        cosmo=info['cosmo']
                        color=info['color']
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
            if cosmo==1:
                                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,cosmological=1, h0=1)
                                S = readsnapcr(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, cosmological=1, h0=1)
            else:
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                S = readsnapcr(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix) 
            Grho = G['rho']
                Tb = G['u']
            cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
            Neb = G['ne']
            try:
                Sm = S['m']
            except KeyError:
                Sm = []
            #TrueTemp, converted_rho  = convertTemp(Tb, Neb, Grho, 1)
            #Gnism = Grho*1e10*2e33/2.9e64/1.1*6e23 #in cm^-3 assume: mean molecular weight=1.1 consider only nuclei
            Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21
            tpi = 2e5/Gnism*250.0 #pi decay time in yr
            betapi = 0.7 ##should be between 0.5-0.9 from Lacki 2011
            #energylost=cregy/tpi*2e53*1e6
            #out=tpi/2<1e6
            #print 'len(out)/len(tot)', len(tpi[out])/float(len(tpi))
            #print 'crenergy in out', np.sum(cregy[out])
            Lgammagev = cregy/tpi*2e53/3.2e7*betapi/3.0 #in erg/s  
#           print 'gamma ray wo cut', np.sum(Lgammagev) 
#           for j in range(len(tpi)):
#               if tpi[j]/2<1e6:
#                   Lgammagev[j]=cregy[j]/2e6*2e53/3.2e7*betapi/3.0
            if np.sum(cregy) >0:
                aveGnism=np.average(Gnism, weights=cregy)
            else:
                aveGnism=0
            snaplist=np.append(snaplist,Nsnapstring)
            crlist=np.append(crlist,np.sum(cregy))
            galist=np.append(galist,np.sum(Lgammagev))
            smlist=np.append(smlist,np.sum(Sm))
            aveGnisml=np.append(aveGnisml, aveGnism)
            print 'Gnism', aveGnism
            print 'cosmic ray energy', np.sum(cregy)
            print 'gamma ray', np.sum(Lgammagev)
            print 'cr energy loss', np.sum(cregy)*2e53/2e5/3.2e7/250*aveGnism
            print 'cr energy inject', 0.8*2.9e40/17*24.3
            del G, S, Grho, Tb, cregy, Neb, Sm, Lgammagev, Gnism, tpi
        spmname='cregy_'+runtodo+'.txt'
        spmfile=open(spmname,"w")
        spmfile.write('snapshot no.   ' + 'cosmic ray energy     '+'gamma ray luminosity     '+' stellar mass  '+'Gnism     '+'\n')
        for ncount in range(len(snaplist)):
            spmfile.write(snaplist[ncount]+'    '+str(crlist[ncount])+'       '+str(galist[ncount])+'       '+str(smlist[ncount])+'       '+str(aveGnisml[ncount])+'\n')
        print 'file name', spmname
        spmfile.close()

if wanted=='gammatime':
    for runtodo in dirneed:
        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr=outdirname(runtodo, Nsnap)
        spmname='cregy_'+runtodo+'.txt'
        spmfile=open(spmname,"r")
        spmfile.readline()
        dars = spmfile.readlines()

        snapno=[]
        cregy=[]
        galum=[]

        for line in dars:
            xsd = line.split()
            snapno = np.append(snapno, int(xsd[0]))
            cregy = np.append(cregy, float(xsd[1]))
            galum = np.append(galum, float(xsd[2]))
        spmfile.close()
        plt.plot(snapno*1e-3,galum,label=r'$\kappa_{\rm di} =$'+dclabel)
    plt.xlabel('Gyr')
    plt.title(runtitle)
    plt.ylabel('Gamma ray luminosity (erg/s)')
    plt.yscale('log')
    plt.legend(loc='best')
    plt.savefig('galum.pdf')


if wanted=='crtime':
        for runtodo in dirneed:
                rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr, timestep=outdirname(runtodo, Nsnap)
                spmname='cregy_'+runtodo+'.txt'
                spmfile=open(spmname,"r")
                spmfile.readline()
                dars = spmfile.readlines()

                snapno=[]
                cregy=[]
                galum=[]

                for line in dars:
                        xsd = line.split()
                        snapno = np.append(snapno, int(xsd[0]))
                        cregy = np.append(cregy, float(xsd[1]))
                        galum = np.append(galum, float(xsd[2]))
                spmfile.close()
                plt.plot(snapno,cregy*2e53,label='kappa='+dclabel+'; '+'time step= '+timestep)
        plt.xlabel('Myr')
        plt.title(runtitle)
        plt.ylabel('CR energy (erg)')
    plt.yscale('log')
        plt.legend(loc='best')
        plt.savefig('CRplot/cregy.pdf')




if wanted=='gammasfr':
        for runtodo in dirneed:
                rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr, timestep=outdirname(runtodo, Nsnap)
                spmname='cregy_'+runtodo+'.txt'
                spmfile=open(spmname,"r")
                spmfile.readline()
                dars = spmfile.readlines()
        print 'file name', spmname
                snapno=[]
                cregy=[]
                galum=[]
        smlist=[]
                for line in dars:
                        xsd = line.split()
                        snapno = np.append(snapno, int(xsd[0]))
                        cregy = np.append(cregy, float(xsd[1]))
                        galum = np.append(galum, float(xsd[2]))
            smlist = np.append(smlist, float(xsd[3]))
                spmfile.close()
        smlists=smlist[::snapsep]
        snapnos=snapno[::snapsep]
        cumgalum = np.cumsum(galum)
        cumgalums=cumgalum[::snapsep]
        sfr = (smlists[1:]-smlists[:-1])/(snapnos[1:]-snapnos[:-1])*1e4 
        avegalums = (cumgalums[1:]-cumgalums[:-1])/(snapnos[1:]-snapnos[:-1])
        #6.2e-4 comes from the different estimate of SNe rate in Lacki 2011 and in Chan 2015 DM paper
        plt.plot(snapnos[:-1],avegalums/(6.2e-4*sfr*2e33/3.2e7*9e20),label='kappa='+dclabel+'; '+'time step= '+timestep)
        print 'SFR', sfr
        print 'snapno', snapnos[:-1]
        print 'FgammaFsf', avegalums/(6.2e-4*sfr*2e33/3.2e7*9e20)
    plt.axhline(y=0.00023,ls='--',color='k')
        plt.xlabel('Myr')
    plt.xlim(xmax=Nsnap)
        plt.title(runtitle)
        plt.ylabel(r'$F_{\gamma}/F_{\rm SF}$')
    plt.yscale('log')
        plt.legend(loc='best')
        plt.savefig('CRplot/galum_sfr_'+fmeat+'.pdf')
    plt.clf()

if wanted=='cumgamma':
        for runtodo in dirneed:
                rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr, timestep=outdirname(runtodo, Nsnap)
                spmname='cregy_'+runtodo+'.txt'
                spmfile=open(spmname,"r")
                spmfile.readline()
                dars = spmfile.readlines()
                print 'file name', spmname
                snapno=[]
                cregy=[]
                galum=[]
                smlist=[]
                for line in dars:
                        xsd = line.split()
                        snapno = np.append(snapno, int(xsd[0]))
                        cregy = np.append(cregy, float(xsd[1]))
                        galum = np.append(galum, float(xsd[2]))
                        smlist = np.append(smlist, float(xsd[3]))
                spmfile.close()
                smlists=smlist[::snapsep]
                snapnos=snapno[::snapsep]
                cumgalum = np.cumsum(galum*3.15e7) #in erg
                cumgalums=cumgalum[::snapsep]/float(snapsep)
                plt.plot(snapnos,cumgalums,label='kappa='+dclabel+'; '+'time step= '+timestep)
        plt.xlabel('Myr')
    plt.ylabel('Total gamma ray energy in erg')
        plt.xlim(xmax=Nsnap)
        plt.title(runtitle)
        plt.yscale('log')
        plt.legend(loc='best')
        plt.savefig('CRplot/cumgalum_'+fmeat+'.pdf')
        plt.clf()


if wanted=='crcooling':
        for runtodo in dirneed:
                rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr, timestep=outdirname(runtodo, Nsnap)
                spmname='cregy_'+runtodo+'.txt'
                spmfile=open(spmname,"r")
                spmfile.readline()
                dars = spmfile.readlines()

                snapno=[]
                cregy=[]
                galum=[]

                for line in dars:
                        xsd = line.split()
                        snapno = np.append(snapno, int(xsd[0]))
                        cregy = np.append(cregy, float(xsd[1]))
                        galum = np.append(galum, float(xsd[2]))
                spmfile.close()
        crcooling = galum*3.0/0.7*2e5*250.0*3.2e7*7.51e-16
        snapnos=snapno[::snapsep]
        crcoolings=crcooling[::snapsep]  
                #plt.plot(snapnos,crcoolings,label=r'$\kappa_{\rm di} =$'+dclabel)
                plt.plot(snapnos,crcoolings,label='kappa='+dclabel+'; '+'time step= '+timestep)
        plt.xlabel('Myr')
        plt.title(runtitle)
        plt.ylabel('CR cooling rate (erg/s)')
        plt.yscale('log')
        plt.legend(loc='best')
        plt.savefig('CRplot/crcooling_'+runtodo+'.pdf')


if wanted=='sncretime':
        for runtodo in dirneed:
                rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr, timestep=outdirname(runtodo, Nsnap)
                S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix)
                Sage = S['age']
                aveSmass = np.average(S['m'])
                hist, bin_edges = np.histogram(Sage, bins=20)
                if havecr==0:
                        needlabel = slabel
                else:
                        needlabel = r' $\kappa_{\rm di}=$'+dclabel
        sncrerate=hist*aveSmass/(bin_edges[1]-bin_edges[0])/0.98*10/0.4*0.0037*1e51/3.16e7*0.1
        print 'sfr', hist*aveSmass/(bin_edges[1]-bin_edges[0])/0.98*10
                plt.plot(Nsnap*1e-3-bin_edges[1:]*0.98, sncrerate, label=needlabel)
                plt.xlabel('Gyr')
                plt.ylabel('CR power from SNe (erg/s)')
        plt.yscale('log')
        plt.legend(loc='best')
        plt.title(runtitle)
        plt.savefig('sncre_'+runtodo+'.pdf')
        plt.clf()



if wanted=='crchange':
        for runtodo in dirneed:
                rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr, timestep=outdirname(runtodo, Nsnap)
                spmname='cregy_'+runtodo+'.txt'
                spmfile=open(spmname,"r")
                spmfile.readline()
                dars = spmfile.readlines()

                snapno=[]
                cregy=[]
                galum=[]

                for line in dars:
                        xsd = line.split()
                        snapno = np.append(snapno, int(xsd[0]))
                        cregy = np.append(cregy, float(xsd[1]))
                        galum = np.append(galum, float(xsd[2]))
                spmfile.close()
        snapno=np.array(snapno[::50])
        cregy=np.array(cregy[::50])
                plt.plot(snapno[:-1],(cregy[1:]-cregy[:-1])/(snapno[1:]-snapno[:-1])*2e53/1e6/3e7,label='kappa='+dclabel+'; '+'time step= '+timestep)
        plt.xlabel('Myr')
        plt.title(runtitle)
        plt.ylabel('rate of change in CR energy (erg/s)')
        plt.yscale('log')
        plt.legend(loc='best')
        plt.savefig('CRplot/crchange_'+runtodo+'.pdf')


if wanted=='enchange':
        for runtodo in dirneed:
                rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr, timestep=outdirname(runtodo, Nsnap)
                spmname='cregy_'+runtodo+'.txt'
                spmfile=open(spmname,"r")
                spmfile.readline()
                dars = spmfile.readlines()
                snapno=[]
                cregy=[]
                galum=[]
        smlist=[]
                for line in dars:
                        xsd = line.split()
                        snapno = np.append(snapno, int(xsd[0]))
                        cregy = np.append(cregy, float(xsd[1]))
                        galum = np.append(galum, float(xsd[2]))
            smlist = np.append(smlist, float(xsd[3]))
                spmfile.close()
                #S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix)
                #Sage = S['age']
                #aveSmass = np.average(S['m'])
                #hist, bin_edges = np.histogram(Sage, bins=20,weights=S['m'])
        smlist5=smlist[::5]
        snapno5=snapno[::5]
        #experiment 24.3/17 is different between Lacki 2011 and my DM paper in SNe rate.
        sncrerate = (smlist5[1:]-smlist5[:-1])/(snapno5[1:]-snapno5[:-1])*1e4/0.98/0.4*0.0037*1e51*0.1/3.2e7*24.3/17
                #sncrerate=hist/(bin_edges[1]-bin_edges[0])/0.98*10/0.4*0.0037*1e51*0.1/3.2e7
        #print 'bin_edges', bin_edges
        plt.plot(snapno5[:-1]/1e3, sncrerate, label=r'0.1$\dot{E}_{SN}$')
        snapno=snapno[::5]
        cregy=cregy[::5]
        galum=galum[::5]
        dcregy=(cregy[1:]-cregy[:-1])/(snapno[1:]-snapno[:-1])*2e53/1e6/3.2e7
        plt.plot(snapno[:-1]/1e3,dcregy,label=r'$\dot{E}_{\rm CR}$')
        if np.amin(dcregy)<0:
            plt.plot(snapno[:-1]/1e3,-dcregy,label=r'$-\dot{E}_{\rm CR}$')
                crcooling = galum*3.0/0.7*2e5*250.0*3.2e7*7.51e-16
                plt.plot(snapno/1e3,crcooling,label=r'$\dot{E}_{cooling}$')
    plt.xlabel('Gyr')
        plt.title(runtitle)
        plt.ylabel(r'$\dot{E}$(erg/s)')        
    plt.yscale('log')
    plt.xlim(xmax=Nsnap*1e-3)
        plt.legend(loc='best', fontsize=14)
        plt.savefig('enchange_'+runtodo+'.pdf')


if wanted=='checksnIItxt':
    for runtodo in dirneed:
        Nsnap=500 #does not have a role; arbitrary number
        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr=outdirname(runtodo, Nsnap)
        spmname='/home/tkc004/scratch/'+rundir+'/output/SNeIIheating.txt'
        spmfile=open(spmname,"r")
        spmfile.readline()
        dars = spmfile.readlines()

        timeinGyr=[]
        ntotal=[]

        for line in dars:
            xsd = line.split()
            timeinGyr = np.append(timeinGyr, float(xsd[0]))
            ntotal = np.append(ntotal, float(xsd[3]))
        spmfile.close()
        noofbins=100
        hist, bin_edges = np.histogram(timeinGyr, bins=noofbins, weights=ntotal)
        print 'sum', np.sum(ntotal)
        barwidth=(np.amax(timeinGyr)-np.amin(timeinGyr))/noofbins
        print 'barwidth', barwidth
        plt.plot(bin_edges[:-1], hist/barwidth*1e51/3.2e16, label=r'$\kappa_{\rm di}=$'+ dclabel)
        plt.xlabel('Gyr')
        plt.ylabel('Sne power (erg/s')
        plt.yscale('log')
        plt.legend(loc='best')
        plt.savefig('sn2egy.pdf')


if wanted=='crcsnap':
        for runtodo in dirneed:
                snaplist=[]
                crlist=[]
                galist=[]
                smlist=[]
                aveGnisml=[]
                #for i in [1, 2]:
                #for i in [500]:
                for i in range(Nsnap):
                        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr=outdirname(runtodo, i)
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
            cregyl = G['cregyl']
            print 'CR energy loss rate (erg/s)', np.sum(cregyl*2e53/1e9/3.2e7)

if wanted=='crcgsnap':
    for runtodo in dirneed:
        snaplist=[]
        enclist=[]
        englist=[]
        enllist=[]
        endlist=[]
        enplist=[]
        prel=0
        preg=0
        prec=0
        pred=0
        prep=0
        for i in range(Nsnap):
            rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr=outdirname(runtodo, i)
            print 'the_snapdir', the_snapdir
            print 'Nsnapstring', Nsnapstring
            print 'havecr', havecr
            G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
            cregyl = G['cregyl']
            cregyg = G['cregyg']
            cregy  = G['cregy']
            cregyd = G['cregyd']
            if havecr>4:
                cregyp = G['cregyp']
            cregygt=np.sum(cregyg)
            cregyt =np.sum(cregy)
            cregylt=np.sum(cregyl)
            cregydt=np.sum(cregyd)
            if havecr>4:
                cregydtp=np.sum(cregyp)
            eng = (cregygt-preg)/1e6/3.2e7*2e53
            enc = (cregyt-prec)/1e6/3.2e7*2e53
            enl = (cregylt-prel)*2e53/1e6/3.2e7
            end = (cregydt-pred)*2e53/1e6/3.2e7
            if havecr>4:
                enp = (cregydtp-prep)*2e53/1e6/3.2e7
            print 'CR energy loss rate (erg/s)', enl
            print 'CR energy gain rate (erg/s)', eng
            print 'CR energy change rate (erg/s)', enc
            print 'CR energy dt rate (erg/s)', end  #including adiabatic heating and streaming
            snaplist.append(i)
            enclist.append(enc)
            englist.append(eng)
            enllist.append(enl)
            endlist.append(end)
            if havecr>4:
                enplist.append(enp)
            preg = cregygt
            prec = cregyt
            prel = cregylt
            pred = cregydt
            if havecr>4:
                prep = cregydtp
        enclist=np.array(enclist)
        englist=np.array(englist)
        endlist=np.array(endlist)
        enllist=np.array(enllist)
        enplist=np.array(enplist)
        plt.plot(snaplist, enclist, label='CR energy change')
        plt.plot(snaplist, englist, label='CR energy gain')
        plt.plot(snaplist, endlist, label='CR energy dt')
        plt.plot(snaplist, enllist, label='CR energy loss')
        plt.plot(snaplist, englist+endlist+enllist, label='CR energy estimate')
        if havecr>4:
            plt.plot(snaplist, enplist, label='pre CR energy dt')
        #plt.yscale('log')
        plt.legend(loc='best')
        plt.xlabel('Myr', fontsize=25)
        plt.ylabel('dE/dt (erg/s)', fontsize=25)
        plt.savefig('cregsnap_'+fmeat+'.pdf')
        plt.clf()
        avesfr=1
        #Lsfr = 3.8e-4 *avesfr*2e33/3.2e7*9e20
        #above is the coefficient for Salpeter only; for Kroupa, the coefficient is 50% larger:
        Lsfr = 3.8e-4 *1.5*avesfr*2e33/3.2e7*9e20
        Lgammae = (enllist+endlist)/7.51e-16/2e5/3.2e7/250.0*0.7/3.0
        #From Guo 2008: 7.51e-16 is rate of total energy loss, including both hadronic and Coulomb
        # hadronic only: 5.86e-16; Coulumb only: 1.65e-16
        # decay rate = 7.51e-16*ne(cm^3)*rhocr(erg/cm^3) in s^-1
        # decay rate from Lacki et al. with the same unit: 6.25e-16 
                Lgamma = (enllist)/7.51e-16/2e5/3.2e7/250.0*0.7/3.0
        Lgammae_sfr = Lgammae/Lsfr
        Lgamma_sfr = Lgamma/Lsfr
        plt.plot(snaplist, np.absolute(Lgamma_sfr), label=r'Ours')
        plt.plot(snaplist, np.absolute(Lgammae_sfr), label='cooling-adiabatic')
        plt.yscale('log')
        plt.legend(loc='best')
                plt.xlabel('Myr', fontsize=25)
                plt.ylabel(r'$\frac{L_{\gamma}}{L_{\rm SF}}$', fontsize=30)
                plt.savefig('gammasfrsnapg_'+fmeat+'.pdf')
                plt.clf()
                
if wanted=='dirgammasfr' or wanted=='dirgamma' or wanted=='dirsfr' or wanted=='dirsm':
    rcParams['figure.figsize'] = 5,4
    Rfrac=0.25
    #Rfrac=0.03
    normalizedsm=1
    M1labelneed=0
    M1runlabelneed=0
    resoneed=0
    diffusionsolverneed=0
    newlabelneed=1
    strlabelneed=1
    showstarburst=0
    legendneed=1
        for runtodo in dirneed:
                snaplist=[]
                enclist=[]
                englist=[]
                enllist=[]
        sml=[]
        diskml=[]
        bulgeml=[]
        nsml=[]
        timel=[]
        Lgcalist=[]
        pretime=0
        presnap=startno
        havecr=0
        if wanted=='dirgammasfr' or wanted=='dirgamma':     
            info=outdirname(runtodo, Nsnap)
            havecr=info['havecr']
            if havecr==0:
                continue
                for i in range(startno,Nsnap, snapsep):
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
            maindir=info['maindir']
            halostr=info['halostr']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
            cosmo=info['cosmo']
            color=info['color']
            withinRv=info['withinRv']
            usepep=info['usepep']
            beginno=info['beginno']
            finalno=info['finalno']
            firever=info['firever']     
            initsnap=info['initsnap']   
            haveB=info['haveB']
            M1speed=info['M1speed']
            Rvirguess=info['Rvir']
            newlabel=info['newlabel']
            strlabel=info['strlabel']
                        labelneed=dclabel
                        if newlabelneed==1:
                                labelneed="\n".join(wrap(newlabel,17))
            if strlabelneed==1:
                labelneed="\n".join(wrap(strlabel,40))
            #snapsep=info['snapsep']
            ptitle=title
            if runtitle=='SMC':
                ptitle='Dwarf'
            elif runtitle=='SBC':
                ptitle='Starburst'
            elif runtitle=='MW':
                ptitle=r'$L\star$ Galaxy'
            if cosmo==0:
                inittime=initsnap
            else:
                inittime=0
                print 'set initial time'    
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
            print 'withinRv', withinRv
            print 'cosmo', cosmo
            try:
                if cosmo==1:
                    G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                    S = readsnapcr(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                    if correctIa==1:
                        Disk = readsnapcr(the_snapdir, Nsnapstring, 2, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                        Bulge = readsnapcr(the_snapdir, Nsnapstring, 3, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                else:
                    G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                    S = readsnapcr(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                    if correctIa==1:
                        Disk = readsnapcr(the_snapdir, Nsnapstring, 2, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Bulge = readsnapcr(the_snapdir, Nsnapstring, 3, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                if havecr>4:
                    cregyl = G['cregyl']*1e10*solar_mass_in_g*km_in_cm*km_in_cm #in erg (originally in code unit: 1e10Msun*(km/s)^2)
                    cregyg = G['cregyg']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                if havecr > 0:
                    cregy_codeunit = G['cregy']
                    cregy  = cregy_codeunit*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                Grho = G['rho']
                Neb = G['ne']
                        except KeyError:
                print 'Keyerror'
                break
            try:
                header=S['header']
                timeneed=header[2]
                print 'timeneed', timeneed
                Smi=S['m']
                Sage=S['age']
                if withinRv ==1 and cosmo==1:
                    Sp = S['p']
                    Sx = Sp[:,0]
                    Sy = Sp[:,1]
                    Sz = Sp[:,2]
                    header=readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1,header_only=1)
                    h0 = header['hubble']
                    atime = header['time']
                    if usepep==1:
                        halosA = SF.read_halo_history_pep(rundir, finalno, beginno=beginno,\
                         singlesnap=0, firever=firever,halonostr=halostr, hubble=h0, comoving=1, maindir=maindir)
                                                afactor=atime
                    else:
                        halosA = SF.read_halo_history(rundir, maindir=maindir, halonostr=halostr, hubble=h0, comoving=0)
                                                afactor=1.0
                    redlist = halosA['redshift']
                    haloid = halosA['ID']
                    a_scale = 1.0/(1.0+redlist)
                    xcenl = halosA['x']*afactor
                    ycenl = halosA['y']*afactor
                    zcenl = halosA['z']*afactor
                    Rvirl = halosA['R']*afactor
                    xcen = np.interp(atime,a_scale,xcenl)
                    ycen = np.interp(atime,a_scale,ycenl)
                    zcen = np.interp(atime,a_scale,zcenl)
                    Rvir = np.interp(atime,a_scale,Rvirl)
                    Sxrel = Sx-xcen
                    Syrel = Sy-ycen
                    Szrel = Sz-zcen
                    Sr = np.sqrt(Sxrel*Sxrel+Syrel*Syrel+Szrel*Szrel)
                    cutrvs = Sr<Rvir*Rfrac
                    Smi = Smi[cutrvs]
                    Sage = Sage[cutrvs]
                Sm = np.sum(Smi)*1e10 #in solar mass
                tcut=Sage>pretime
                Nsm = np.sum(Smi[tcut])*1e10
                if correctIa==1:
                    diskm = np.sum(Disk['m'])*1e10
                    bulgem = np.sum(Bulge['m'])*1e10
            except KeyError:
                print 'key error'
                Sm = 0.
                Nsm = 0.
                timeneed=0
                if correctIa==1:
                    diskm = 0.
                    bulgem = 0.
                        if withinRv ==1:
                                Gp = G['p']
                                Gx = Gp[:,0]
                                Gy = Gp[:,1]
                                Gz = Gp[:,2]
                                header=readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1,header_only=1)
                                h0 = header['hubble']
                                atime = header['time']
                if cosmo==1:
                    if usepep==1:
                        halosA = SF.read_halo_history_pep(rundir, finalno,\
                        beginno=beginno, singlesnap=0, firever=firever,halonostr=halostr,\
                        hubble=h0, comoving=1, maindir=maindir, snapsep=snapsep)
                                                afactor=atime
                    else:
                        halosA = SF.read_halo_history(rundir, maindir=maindir, halonostr=halostr, hubble=h0, comoving=0)
                                                afactor=1.0
                    redlist = halosA['redshift']
                    haloid = halosA['ID']
                    a_scale = 1.0/(1.0+redlist)
                    xcenl = halosA['x']*afactor
                    ycenl = halosA['y']*afactor
                    zcenl = halosA['z']*afactor
                    Rvirl = halosA['R']*afactor
                    xcen = np.interp(atime,a_scale,xcenl)
                    ycen = np.interp(atime,a_scale,ycenl)
                    zcen = np.interp(atime,a_scale,zcenl)
                    Rvir = np.interp(atime,a_scale,Rvirl)
                else:
                    xcen=0
                    ycen=0
                    zcen=0
                    Rvir=Rvirguess
                                Gxrel = Gx-xcen
                                Gyrel = Gy-ycen
                                Gzrel = Gz-zcen
                                Gr = np.sqrt(Gxrel*Gxrel+Gyrel*Gyrel+Gzrel*Gzrel)
                                cutrv = Gr<Rvir*Rfrac
                                Grho = Grho[cutrv]
                                Neb = Neb[cutrv]
                if havecr>0:
                    cregy = cregy[cutrv]
                    cregy_codeunit = cregy_codeunit[cutrv]
                                if havecr>4:
                                        cregyl = cregyl[cutrv]
                                        cregyg = cregyg[cutrv]
                        if cosmo==1:
                                readtimelist=readtime(firever=2)
                                snap2list=readtimelist['snaplist']
                                time2list=readtimelist['timelist']
                                a2list=readtimelist['alist']
                                tnow = np.interp(timeneed,a2list,time2list)*1e9
                                pret = np.interp(pretime,a2list,time2list)*1e9
            if havecr>4:
                cregygt=np.sum(cregyg)
                cregylt=np.sum(cregyl)
            if havecr>0:
                cregyt =np.sum(cregy)
                Lout = outLgamma_nism(Grho,Neb,cregy_codeunit)
                Lgcal = np.sum(Lout['Lgamma'])
                print 'Lgcal', Lgcal
                        snaplist.append(float(i))
            if cosmo==1:
                timel.append(tnow)
            else:
                timel.append(float(i)*0.98*1e6)
            if havecr>0:
                enclist.append(cregyt)
            if havecr>4:
                englist.append(cregygt)
                enllist.append(cregylt)
            if havecr>0:
                Lgcalist.append(Lgcal)
            sml.append(Sm)
            if correctIa==1:
                diskml.append(diskm)
                bulgeml.append(bulgem)
            nsml.append(Nsm)
            pretime=timeneed
            del G, S
        sml=np.array(sml)
        if correctIa==1:
            diskml=np.array(diskml)
            bulgeml=np.array(bulgeml)
        nsml=np.array(nsml)
        if havecr>0:
            enclist=np.array(enclist)
        if havecr>4:
            englist=np.array(englist)
            enllist=np.array(enllist)
        snaplist=np.array(snaplist)
        if havecr>0:
            Lgcalist=np.array(Lgcalist)
        timel=np.array(timel) #in yr
        #above is the coefficient for Salpeter only; for Kroupa, the coefficient is 50% larger:
        avesfrl=(nsml[1:])/(timel[1:]-timel[:-1]) #in Msun/yr
        print 'nsml',nsml
        print 'avesfrl', avesfrl
        Lsfr = Kroupa_Lsf*avesfrl*solar_mass_in_g/yr_in_sec*cspeed_in_cm_s*cspeed_in_cm_s #times 1.5? 6.2e-4 for Kroupa 3.8e-4 for Salpeter
        if correctIa==1:
            IaeffSFR=5.3e-8/3.0e-4*(sml[1:]+diskml[1:]+bulgeml[1:])/0.03753/1e9
            print 'IaeffSFR', IaeffSFR
            LsfrnoIa = Kroupa_Lsf*(avesfrl+IaeffSFR)*solar_mass_in_g/yr_in_sec*cspeed_in_cm_s*cspeed_in_cm_s
        print 'Lsfr', Lsfr
        print 'timel', timel
        if havecr>4:
            Lgamma = (enllist[1:]-enllist[:-1])/((timel[1:]-timel[:-1])*yr_in_sec)/(hadronicdecayrate+coulombdecayrate)*hadronicdecayrate*betapi/nopi_per_gamma
            Lgamma_sfr = Lgamma/Lsfr
        if havecr>0:
            Lgcal_sfr = Lgcalist[1:]/Lsfr
        if correctIa==1:
            Lgamma_sfr_noIa = Lgamma/LsfrnoIa
        if haveB>0:
            lsn='dashed'
        else:
            lsn='solid'
        if M1labelneed>1 or M1runlabelneed==1:
            if M1speed>499:
                                lsn='solid'
            if M1speed>999:
                lsn='dashed'
            if M1speed>1999:
                lsn='dashdot'
            if M1speed>3999:
                lsn='dotted'
        if M1labelneed==1:
                        if M1speed>499:
                                labelneed=r'$\tilde{c}=500$'
                        if M1speed>999:
                                labelneed=r'$\tilde{c}=1000$'
                        if M1speed>1999:
                                labelneed=r'$\tilde{c}=2000$'
                        if M1speed>3999:
                                labelneed=r'$\tilde{c}=4000$'

        if resoneed==1:
            if resolabel=='llr':
                labelneed='Lowest res'
                lsn = 'solid'
            if resolabel=='lr':
                labelneed='Low res'
                lsn = 'dashed'
            if resolabel=='mr':
                labelneed='Standard res'
                lsn = 'dashdot'
        if diffusionsolverneed==1:
            if runtodo=='bwmwlrdc27ds':
                labelneed='Zeroth moment'
                lsn = 'dashed'
            else:
                labelneed='Two moment'
                lsn = 'solid'
        if xaxis_snapno==1:
            xaxisl = snaplist[1:]
            xlab = 'snapshot number'
        else:
            xaxisl = timel[1:]/1e6+inittime
            xlab = 'Myr'
        if showstarburst==1 and runtitle=='SBC':
            if runtodo=='bwsbclr':
                pstartno=600
            if runtodo=='bwsbclrdc0':
                pstartno=590
            if runtodo=='bwsbclrdc27':
                pstartno=440
            if runtodo=='bwsbclrdc28':
                pstartno=505
            if runtodo=='bwsbclrdc29':
                pstartno=370
        print 'labelneed', labelneed
        if wanted=='dirgammasfr':
            if havecr>4:
                #plt.plot(xaxisl, np.absolute(Lgamma_sfr),  color=color,ls=lsn,lw=2, label=labelneed)
                if (M1labelneed==1 and color=='r'):
                                        plt.plot(xaxisl, np.absolute(Lgamma_sfr),  color=color,ls=lsn,lw=2)
                elif (M1labelneed==2 and lsn=='solid'):
                    plt.plot(xaxisl, np.absolute(Lgamma_sfr),  color=color,ls=lsn,lw=2, label=labelneed)
                else:
                    plt.plot(xaxisl, np.absolute(Lgamma_sfr),  color=color,ls=lsn,lw=2, label=labelneed)
                #print 'M1speed', M1speed
                print 'Lgamma_sfr', Lgamma_sfr
            elif havecr>0:
                plt.plot(xaxisl, np.absolute(Lgcal_sfr), color=color,ls=lsn,lw=2, label=labelneed)
                print 'Lgcal_sfr', Lgcal_sfr
            if correctIa==1 and havecr>0:
                plt.plot(xaxisl, np.absolute(Lgamma_sfr_noIa),ls=lsn,lw=2, color=color)
                print 'Lgamma_sfr_noIa', Lgamma_sfr_noIa
                if wanted=='dirgamma':
                        if havecr>4:
                                plt.plot(xaxisl, np.absolute(Lgamma),  color=color,ls=lsn,lw=2, label=labelneed)
            elif havecr>0:
                plt.plot(xaxisl, np.absolute(Lgcalist[1:]), color=color,ls=lsn,lw=2, label=labelneed)
                if wanted=='dirsfr':
            if M1labelneed==1 and color=='r':
                                plt.plot(xaxisl, np.absolute(avesfrl), color=color,ls=lsn,lw=2)
            elif M1speed<501 and (M1runlabelneed==1 or M1labelneed==1):
                plt.plot(xaxisl, np.absolute(avesfrl), color=color,ls=lsn,lw=2, label=labelneed)
            else:
                plt.plot(xaxisl, np.absolute(avesfrl), color=color,ls=lsn,lw=2, label=labelneed)
            if showstarburst==1 and runtitle=='SBC':
                plt.axvline(x=pstartno, color=color,ls=lsn,lw=1)
        if wanted=='dirsm':
            if M1labelneed==1 and color=='r':
                if normalizedsm==1:
                    plt.plot(xaxisl, sml[:-1]-sml[0], color=color,ls=lsn,lw=2)
                else:
                    plt.plot(xaxisl, sml[1:], color=color,ls=lsn,lw=2)
            else:
                if normalizedsm==1:
                                        plt.plot(xaxisl, sml[:-1]-sml[0], color=color,ls=lsn,lw=2, label=labelneed)
                                else:
                                        plt.plot(xaxisl, sml[1:], color=color,ls=lsn,lw=2, label=labelneed)
                        if showstarburst==1 and runtitle=='SBC':
                                plt.axvline(x=pstartno, color=color,ls='dashed',lw=1)
    if wanted=='dirgammasfr':
        plt.axhline(y=0.0002,ls='--',color='k')
        plt.yscale('log')
        if runtitle=='SMC' and legendneed==1:
            plt.legend(loc='best', fontsize=8)
                elif legendneed==1 and M1runlabelneed==1:
                        plt.legend(loc='best', fontsize=8,ncol=2)
        plt.xlabel(xlab, fontsize=16)
        plt.ylabel(r'$L_{\gamma}/L_{\rm SF}$', fontsize=16)
        figname='CRplot/dirgammasfr/gammasfrsnap_'+fmeat+'.pdf'
        if wanted=='dirgamma':
                plt.yscale('log')
                if runtitle=='SMC' or legendneed==1:
                        plt.legend(loc='best', fontsize=8,ncol=3)
        if legendneed==1 and M1runlabelneed:
                plt.legend(loc='best', fontsize=8,ncol=2)
                plt.xlabel(xlab, fontsize=16)
                plt.ylabel(r'$L_{\gamma} {\rm erg/s}$', fontsize=16)
        figname='CRplot/dirgamma/gammasnap_'+fmeat+'.pdf'
        if wanted=='dirsfr':
                plt.yscale('log')
                if runtitle=='SMC' or legendneed==1:
                        plt.legend(loc='best', fontsize=8,ncol=3)
        if M1labelneed==1 and legendneed==1:
                plt.legend(loc='best', fontsize=8,ncol=2)
                plt.xlabel(xlab, fontsize=16)
        plt.ylabel(r'${\rm SFR (M_{\odot}/yr)} $', fontsize=16)
        figname='CRplot/dirsfr/sfrsnap_'+fmeat+'.pdf'
    if wanted=='dirsm':
        if runtitle=='SMC' and legendneed==1:
            if strlabelneed==1 or M1labelneed==1:
                plt.legend(loc='best', fontsize=10,ncol=2)
            else:
                plt.legend(loc='best', fontsize=8,ncol=3)
        plt.xlabel('Myr', fontsize=16)
        plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
        if normalizedsm==1:
            plt.ylabel(r'${\rm M_{*,new} (M_{\odot})} $', fontsize=16)
        else:
            plt.ylabel(r'${\rm M_* (M_{\odot})} $', fontsize=16)
        if normalizedsm==1:
            figname='CRplot/dirsm/nsm_'+fmeat+'.pdf'
        else:
            figname='CRplot/dirsm/sm_'+fmeat+'.pdf'
    print 'Saving ', figname
    plt.title(ptitle,fontsize=16)
    plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
    plt.savefig(figname,bbox_inches='tight')
    plt.clf()






if wanted=='crasnap':
        for runtodo in dirneed:
                snaplist=[]
                enclist=[]
                englist=[]
                enllist=[]
                endlist=[]
                enplist=[]
        enalist=[]
                avesfrl=[]
                prel=0
                preg=0
                prec=0
                pred=0
                prep=0
                presm = 0
                presnap = 0
                for i in range(0,Nsnap, snapsep):
                        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr,timestep=outdirname(runtodo, i)
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        cregyl = G['cregyl']
                        cregyg = G['cregyg']
                        cregy  = G['cregy']
                        cregyd = G['cregyd']
                        if havecr>4:
                                cregyp = G['cregyp']
            if havecr>5:
                cregya = G['cregya']
                        try:
                                Sm = np.sum(S['m'])
                        except KeyError:
                                Sm = 0.
                        cregygt=np.sum(cregyg)
                        cregyt =np.sum(cregy)
                        cregylt=np.sum(cregyl)
                        cregydt=np.sum(cregyd)
                        if havecr>4:
                                cregydtp=np.sum(cregyp)
            if havecr>5:
                cregyat = np.sum(cregya)
                        eng = (cregygt)*2e53
                        enc = (cregyt)*2e53
                        enl = (cregylt)*2e53
                        end = (cregydt)*2e53
                        if havecr>4:
                                enp = (cregydtp)*2e53
                                print 'CR energy dtp', enp
            if havecr>5:
                ena = cregyat*2e53
                print 'CR energy dta', ena
                        snaplist.append(i)
                        enclist.append(enc)
                        englist.append(eng)
                        enllist.append(enl)
                        endlist.append(end)
                        if havecr>4:
                                enplist.append(enp)
            if havecr>5:
                enalist.append(ena)
                        if i>presnap:
                                avesfr=(Sm-presm)*1e4/0.98/(i-presnap)
                        else:
                                avesfr=0
                        print 'avesfr', avesfr
                        avesfrl.append(avesfr)
                avesfrl=np.array(avesfrl)
                enclist=np.array(enclist)
                englist=np.array(englist)
                endlist=np.array(endlist)
                enllist=np.array(enllist)
                enplist=np.array(enplist)
                enalist=np.array(enalist)
                plt.plot(snaplist, enclist, label='CR energy')
                plt.plot(snaplist, englist, label='SNe')
                if havecr > 4:
                        plt.plot(snaplist, endlist, label='Ad')
        #if havecr > 6:
        #   plt.plot(snaplist, enalist, label='Pure Ad')
        #if havecr > 7:
        #   plt.plot(snaplist, enplist, label='Other Ad')
                #plt.plot(snaplist, endlist, label='CR energy dt')
                plt.plot(snaplist, enllist, label='Loss')
                #plt.plot(snaplist, englist+endlist+enllist, label='CR energy estimate')
                if havecr>4:
                        plt.plot(snaplist, enclist-(englist+endlist+enllist), label='Extra')
                #plt.yscale('log')
        print 'havecr', havecr
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                plt.xlabel('Myr', fontsize=25)
                plt.subplots_adjust(left=0.2,bottom=0.2, right=0.75)
                plt.ylabel(r'$E_{\rm tot}$ (erg)', fontsize=25)
                plt.savefig('CRplot/crasnap_'+fmeat+'.pdf')
                plt.clf()


if wanted=='crdelv':
        for runtodo in dirneed:
                snaplist=[]
                enplist=[]
        enalist=[]
                prep=0
        prea=0
                for i in range(Nsnap):
                        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr=outdirname(runtodo, i)
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        cregyp = G['cregyp']
            cregya = G['cregya']
            cregydtp=np.sum(cregyp)
            cregydta=np.sum(cregya)
            enp = (cregydtp-prep)*2e53/1e6/3.2e7
            ena = (cregydta-prea)*2e53/1e6/3.2e7    
                        snaplist.append(i)
            enplist.append(enp)
            enalist.append(ena)
                        prep = cregydtp
            prea = cregydta
                enplist=np.array(enplist)
        enalist=np.array(enalist)
                plt.plot(snaplist, enplist, label='CR energy change')
                plt.plot(snaplist, enalist, label='CR energy delv')
                plt.legend(loc='best')
                plt.xlabel('Myr', fontsize=25)
                plt.ylabel('dE/dt (erg/s)', fontsize=25)
                plt.savefig('credelv.pdf')
                plt.clf()


if wanted=='gasden':
        for runtodo in dirneed:
                for i in [Nsnap]:
                        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr=outdirname(runtodo, i)
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gp = G['p']
                        #Grho = G['rho']
                        Gu = G['u']
                        Gm = G['m']
                        cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                        Neb = G['ne']
                        #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                        Gz = Gp[:,2]
                        Gx = Gp[:,0]
                        Gy = Gp[:,1]
            dr = withinr/nogrid
            Gnism_in_cm_3l=[]
            radl =[]
            for irad in range(nogrid):
                cutxy = (Gx*Gx+Gy*Gy > dr*irad*dr*irad) & (Gx*Gx+Gy*Gy < dr*(irad+1)*dr*(irad+1))
                cutz = np.absolute(Gz-med)< maxlength/2.
                cut = cutxy*cutz
                Nebcut = Neb[cut]
                Gmcut = Gm[cut] 
                Gm_in_g = Gmcut*1e10*2e33
                shellvol_in_cm3 = np.pi*maxlength*(-np.power(dr*irad,2)+np.power(dr*(irad+1),2))*3.086e21*3.086e21*3.086e21
                Grho_in_g_cm_3 = Gm_in_g/shellvol_in_cm3
                Gnism_in_cm_3 = np.sum((0.78+0.22*Nebcut*0.76)/protonmass_in_g*Grho_in_g_cm_3)
                Gnism_in_cm_3l = np.append(Gnism_in_cm_3l, Gnism_in_cm_3)
                radl = np.append(radl, dr*(irad+0.5))
                        plt.plot(radl, Gnism_in_cm_3l, label=runtodo)
        plt.xlabel('r [kpc]')
        plt.ylabel(r'$n_{\rm ISM} [{\rm cm^{-3}}]$')
        plt.legend(loc='best')
    plt.title('cylinder centered at disk '+str(med)+' kpc above z plane with a height ' + str(maxlength) + ' kpc')
        plt.savefig('gasdensity_'+fmeat+'_r.pdf')
        plt.clf()

if wanted=='gmr' or wanted=='avedenr' or wanted=='avecrdenr'\
or wanted=='avethdenr' or wanted=='aveedenr' or wanted=='FpFgr'\
or wanted=='avekedenr':
        def funcab(x, a, b):
                return a+b*x
        def func(x, a, b, c):
                return a+b*x+c*x*x
        def d10func(r, a, b, c):
                return np.power(10.,func(np.log10(r), a, b, c))*(b+2*c*np.log(r)/np.log(10.))/r
        rcParams['figure.figsize'] = 5,4
    withinr=20.0
    atstarburst=0
    linelabelneed=0
    if wanted=='aveedenr':
        linelabelneed=1
    labcount=0
    newlabelneed=1
        for runtodo in dirneed:
        if atstarburst==1:
            if runtodo=='bwsbclr':
                Nsnap=523
                        if runtodo=='bwsbclrmhd':
                                Nsnap=589
                        if runtodo=='bwsbclrdc0':
                                Nsnap=590
                        if runtodo=='bwsbclrdc27':
                                Nsnap=138
                        if runtodo=='bwsbclrdc28':
                                Nsnap=520
                        if runtodo=='bwsbclrdc29':
                                Nsnap=433
                        if runtodo=='bwsbclrstr':
                                Nsnap=245
            if runtodo=='bwsbclrdc28mhd':
                                Nsnap=558
            if runtodo=='bwsbclrdc28str':
                                Nsnap=370
                for i in [Nsnap]:
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        maindir=info['maindir']
                        halostr=info['halostr']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        cosmo=info['cosmo']
                        color=info['color']
                        withinRv=info['withinRv']
                        usepep=info['usepep']
                        beginno=info['beginno']
                        finalno=info['finalno']
                        firever=info['firever']
                        initsnap=info['initsnap']
                        haveB=info['haveB']
                        M1speed=info['M1speed']
                        Rvirguess=info['Rvir']
            newlabel=info['newlabel']
            ptitle=title
                        labelneed=dclabel
                        if newlabelneed==1:
                                labelneed="\n".join(wrap(newlabel,17))
            if runtitle=='SMC':
                ptitle='Dwarf'
            elif runtitle=='SBC':
                ptitle='Starburst'
            elif runtitle=='MW':
                ptitle=r'$L\star$ Galaxy'
            lsn='solid'
            if haveB>0:
                lsn='dashed'
                        if havecr==0 and wanted=='avecrdenr':
                                continue
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gp = G['p']
            Gv = G['v']
                        Gu = G['u']
                        Gm = G['m']
            Gvx = Gv[:,0]
            Gvy = Gv[:,1]
            Gvz = Gv[:,2]
            GEint = Gu*km_in_cm*km_in_cm*Gm*1e10*Msun_in_g
                        cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
            ke = 0.5*(Gvx*Gvx+Gvy*Gvy+Gvz*Gvz)*km_in_cm*km_in_cm*Gm*1e10*Msun_in_g  
                        #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                        Gz = Gp[:,2]
                        Gx = Gp[:,0]
                        Gy = Gp[:,1]
                        dr = withinr/nogrid
            from crtestfunction import findcenz
            datasup=1
                        xcen = findcenz(runtodo,Nsnap,withinr=withinr,dir='x',datasup=datasup,Gx=Gx,Gy=Gy,Gz=Gz,Gm=Gm)
                        ycen = findcenz(runtodo,Nsnap,withinr=withinr,dir='y',datasup=datasup,Gx=Gx,Gy=Gy,Gz=Gz,Gm=Gm)
            zcen = findcenz(runtodo,Nsnap,withinr=withinr,dir='z',datasup=datasup,Gx=Gx,Gy=Gy,Gz=Gz,Gm=Gm)
            Gz = Gz-zcen;Gy = Gy-ycen; Gx = Gx-xcen
                        Gm_in_sunl=[]
            GEl = []
            crel = []
                        radl =[]
            kel = []
                        for irad in range(nogrid):
                                cutxy = (Gx*Gx+Gy*Gy < dr*irad*dr*irad) 
                                cutz = np.absolute(Gz-med)< maxlength/2.
                cut = cutxy*cutz
                crcut = cregy[cut]*1e10*Msun_in_g*km_in_cm*km_in_cm #erg
                crel = np.append(crel, np.sum(crcut))
                Gmcut = Gm[cut]
                Gm_in_sun = Gmcut*1e10
                Gm_in_sunl = np.append(Gm_in_sunl, np.sum(Gm_in_sun))
                GEcut = GEint[cut]
                kecut = ke[cut]
                GEl = np.append(GEl,np.sum(GEcut)) #erg
                kel = np.append(kel,np.sum(kecut))
                                radl = np.append(radl, dr*irad)
            if wanted == 'gmr':
                plt.plot(radl, Gm_in_sunl, label=runtodo)
            if wanted == 'avedenr':
                vol = np.power(radl*kpc_in_cm,3.)*4./3.*np.pi
                den = Gm_in_sunl*Msun_in_g/vol/protonmass_in_g
                plt.plot(radl,den,label=labelneed,lw=2,ls=lsn,color=color)
                print 'radl, den', radl, den
                                x0 = [1.0,-0.1]
                xdata = np.log10(radl[1:]); ydata = np.log10(den[1:])
                outfit=optimize.curve_fit(funcab, xdata, ydata,x0)
                afit=outfit[0][0]
                bfit=outfit[0][1]
                print 'afit, bfit', afit, bfit
                        if wanted == 'avecrdenr' or wanted == 'aveedenr':
                if havecr>0:
                    vol = np.power(radl*kpc_in_cm,3.)*4./3.*np.pi
                    crden = crel*erg_in_eV/vol
                    if linelabelneed==1:
                        if labcount==2:
                            plt.plot(radl,crden,label='CR',lw=2,ls=lsn,color=color)
                        else:
                            plt.plot(radl,crden,lw=2,ls=lsn,color=color)
                    else:
                        plt.plot(radl,crden,label=dclabel,lw=2,ls=lsn,color=color)
            if wanted == 'avethdenr' or wanted == 'aveedenr':
                vol = np.power(radl*kpc_in_cm,3.)*4./3.*np.pi
                                thden = GEl*erg_in_eV/vol
                if linelabelneed==1:
                    if labcount==2:
                        plt.plot(radl,thden,label='Thermal',lw=1,ls=lsn,color=color)
                    else:
                        plt.plot(radl,thden,lw=1,ls=lsn,color=color)
                else:
                    plt.plot(radl,thden,label=dclabel,lw=2,ls=lsn,color=color)
            if wanted == 'avekedenr':
                vol = np.power(radl*kpc_in_cm,3.)*4./3.*np.pi
                keden = kel*erg_in_eV/vol
                plt.plot(radl,keden,label=dclabel,lw=2,ls=lsn,color=color)
            if wanted == 'FpFgr':
                vol = np.power(radl*kpc_in_cm,3.)*4./3.*np.pi
                dv = vol[1:]-vol[:-1]
                drad = (radl[1:]-radl[:-1])*kpc_in_cm
                dpth = (GAMMA-1.0)*(GEl[1:]-GEl[:-1])/dv/drad #cgs
                Fg_m = Gm_in_sunl*Msun_in_g*NewtonG_in_cgs/(radl*kpc_in_cm)/(radl*kpc_in_cm)
                dm = (Gm_in_sunl[1:]-Gm_in_sunl[:-1])*Msun_in_g
                fg = Fg_m[:-1]*dm/dv #force per unit volume
                plt.plot(radl[:-1],dpth/fg,label=dclabel,lw=2,ls='dashdot',color=color)
                if havecr>0:
                    dpcr = (CRgamma-1.0)*(crel[1:]-crel[:-1])/dv/drad #cgs
                    plt.plot(radl[:-1],dpcr/fg,label=dclabel,lw=2,ls='solid',color=color)
            labcount+=1
    plt.yscale('log')
        plt.xlabel(r'$ r\; [{\rm kpc}]$',fontsize=14)
    if wanted =='gmr':
        plt.ylabel(r'enclosed $M_{\rm g} [M_{\odot}]$')
        filename='CRplot/gmr/gasmass_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
    if wanted =='avedenr':
        plt.ylabel(r'$\bar{n}_{\rm ISM}(<r)\; [{\rm cm^{-3}}]$',fontsize=14)
        filename='CRplot/avedenr/avedenr_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
        plt.legend(loc='best',fontsize=8,ncol=2)
        if wanted =='avecrdenr':
                plt.ylabel(r'$\bar{e}_{\rm cr}\; [{\rm eV/cm^{3}}]$',fontsize=14)
                filename='CRplot/avecrdenr/avecrdenr_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
    if wanted =='avethdenr':
                plt.ylabel(r'$\bar{e}_{\rm th}\; [{\rm eV/cm^{3}}]$',fontsize=14)
                filename='CRplot/avethdenr/avethdenr_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
        if wanted =='aveedenr':
        plt.ylabel(r'$\bar{e}(<r)\; [{\rm eV/cm^{3}}]$',fontsize=14)
                filename='CRplot/aveedenr/aveedenr_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
        plt.legend(loc='best',fontsize=14)
    if wanted=='avekedenr':
        plt.ylabel(r'$\bar{e}_{\rm KE}\; [{\rm eV/cm^{3}}]$',fontsize=14)
                filename='CRplot/avekedenr/avekedenr_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
        if wanted =='FpFgr':
                plt.ylabel(r'$\nabla P/\rho g$',fontsize=14)
                filename='CRplot/FpFgr/FpFgr_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
    plt.tight_layout()
    print 'filename', filename
    plt.title(ptitle,fontsize=16)
        plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
        plt.savefig(filename)
        plt.clf()

if wanted=='crer' or wanted=='crdenr' or wanted=='crdpr' or wanted=='denr' or wanted=='rhovr' or wanted=='rhovr_denr':
    startr=5.
    withinr = 200. #kpc
    thermalneed = 1
    nogrid = 10
    from crtestfunction import *
    def funcab(x, a, b):
        return a+b*x
        def func(x, a, b, c):
                return a+b*x+c*x*x
        def d10func(r, a, b, c):
                return np.power(10.,func(np.log10(r), a, b, c))*(b+2*c*np.log(r)/np.log(10.))/r

        for runtodo in dirneed:
        info=outdirname(runtodo, Nsnap)
        havecr=info['havecr']
                for i in [Nsnap]:
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        maindir=info['maindir']
                        halostr=info['halostr']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        cosmo=info['cosmo']
                        color=info['color']
                        withinRv=info['withinRv']
                        usepep=info['usepep']
                        beginno=info['beginno']
                        finalno=info['finalno']
                        firever=info['firever']
                        initsnap=info['initsnap']
                        haveB=info['haveB']
                        M1speed=info['M1speed']
                        Rvirguess=info['Rvir']
            if haveB>0:
                lsn='dashed'
            else:
                lsn='solid'
                        if cosmo==1:
                                header=readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1,header_only=1)
                                hubble = header['hubble']
                                atime = header['time']
                                print 'halostr', halostr
                                halosA = SF.read_halo_history_pep(rundir, i, singlesnap=1, firever=firever,halonostr=halostr, hubble=hubble, comoving=0, maindir=maindir)
                                redlist = halosA['redshift']
                                haloid = halosA['ID']
                                a_scale = 1.0/(1.0+redlist)
                                xcen = halosA['x']
                                ycen = halosA['y']
                                zcen = halosA['z']
                                xvcen = halosA['xv']
                                yvcen = halosA['yv']
                                zvcen = halosA['zv']
                                Rvir = halosA['R']
                        else:
                                xcen=ycen=zcen=0.
                                xvcen=yvcen=zvcen=0.
            if cosmo==1:
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix,\
                 extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
            else:   
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gp = G['p']
                        #Grho = G['rho']
                        Gu = G['u']
                        Gm = G['m']
            Gv = G['v']
                        cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                        Neb = G['ne']
                        #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                        Gz = Gp[:,2]-zcen
                        Gx = Gp[:,0]-xcen
                        Gy = Gp[:,1]-ycen
            Gvz = Gv[:,2]-zvcen
            Gvx = Gv[:,0]-xvcen
            Gvy = Gv[:,1]-yvcen
            Gr = np.sqrt(Gx*Gx+Gy*Gy+Gz*Gz)
            vr = (Gvx*Gz+Gvy*Gy+Gvz*Gz)/Gr*km_in_cm
                        dr = withinr/nogrid
                        crel= []
            thl = [] 
            gml = []
            gmvl=[]
                        radl =np.linspace(startr,withinr,nogrid)
                        for irad in range(nogrid):
                                cut = (np.sqrt(Gx*Gx+Gy*Gy+Gz*Gz) < radl[irad])
                                crecut = cregy[cut]*1e10*Msun_in_g*km_in_cm*km_in_cm
                                crel = np.append(crel, np.sum(crecut))
                Gmcut = Gm[cut]*1e10*Msun_in_g
                gml = np.append(gml,np.sum(Gmcut))
                vrcut = vr[cut]
                gmvl = np.append(gmvl,np.sum(Gmcut*vr[cut]))
                if thermalneed>0:
                    thcut = Gu[cut]*Gmcut*km_in_cm*km_in_cm
                    thl = np.append(thl, np.sum(thcut))
                        if wanted=='denr':
                                radlm = (radl[1:]+radl[:-1])/2.
                                vol = 4./3.*np.pi*(radl[1:]*radl[1:]*radl[1:]-radl[:-1]*radl[:-1]*radl[:-1])*kpc_in_cm*kpc_in_cm*kpc_in_cm
                denl = (gml[1:]-gml[:-1])/vol
                                plt.plot(radlm, denl, label=dclabel,lw=2,color=color,ls=lsn)
                x0 = [1.0,-0.1]
                xdata=np.log10(radlm)
                ydata=np.log10(denl)
                outfit=optimize.curve_fit(funcab, xdata, ydata,x0)
                afit=outfit[0][0]
                bfit=outfit[0][1]
                print 'afit,bfit', afit,bfit
                                filename='CRplot/denr/denr_'+fmeat+'.pdf'
                        if wanted=='rhovr':
                                radlm = (radl[1:]+radl[:-1])/2.
                                vol = 4./3.*np.pi*(radl[1:]*radl[1:]*radl[1:]-radl[:-1]*radl[:-1]*radl[:-1])*kpc_in_cm*kpc_in_cm*kpc_in_cm
                                rhovrl = (gmvl[1:]-gmvl[:-1])/vol
                                plt.plot(radlm, rhovrl, label=dclabel,lw=2,color=color,ls=lsn)
                                x0 = [1.0,-0.1]
                                xdata=np.log10(radlm)
                                ydata=np.log10(np.absolute(rhovrl))
                                outfit=optimize.curve_fit(funcab, xdata, ydata,x0)
                                afit=outfit[0][0]
                                bfit=outfit[0][1]
                                print 'afit,bfit', afit,bfit
                                filename='CRplot/rhovr/rhovr_'+fmeat+'.pdf'
            if wanted=='rhovr_denr':
                radlm = (radl[1:]+radl[:-1])/2.
                rhovr_denrl = (gmvl[1:]-gmvl[:-1])/(gml[1:]-gml[:-1])
                                plt.plot(radlm, rhovr_denrl, label=dclabel,lw=2,color=color,ls=lsn)
                                x0 = [1.0,-0.1]
                                xdata=np.log10(radlm)
                                ydata=np.log10(np.absolute(rhovr_denrl))
                                outfit=optimize.curve_fit(funcab, xdata, ydata,x0)
                                afit=outfit[0][0]
                                bfit=outfit[0][1]
                                print 'afit,bfit', afit,bfit
                                filename='CRplot/rhovr_denr/rhovr_denr_'+fmeat+'.pdf'
            if wanted=='crer':
                plt.plot(radl, crel, label=dclabel,lw=2,color=color,ls=lsn)
                if thermalneed>1:
                    plt.plot(radl, thl,lw=1,color=color,ls=lsn)
                filename='CRplot/crer/crer_'+fmeat+'.pdf'
            if wanted=='crdenr':
                print 'runtodo', runtodo
                radlm = (radl[1:]+radl[:-1])/2.
                vol = 4./3.*np.pi*(radl[1:]*radl[1:]*radl[1:]-radl[:-1]*radl[:-1]*radl[:-1])*kpc_in_cm*kpc_in_cm*kpc_in_cm
                if havecr>0:
                    crden = (crel[1:]-crel[:-1])*erg_in_eV/vol
                    x0 = [1.0,-0.1]
                    xdata=np.log10(radlm)
                    ydata=np.log10(crden)
                    outfit=optimize.curve_fit(funcab, xdata, ydata,x0)
                    afit=outfit[0][0]
                    bfit=outfit[0][1]
                    print 'afit,bfit', afit,bfit
                    crden_ana = np.power(10.,funcab(np.log10(radlm),afit,bfit))
                    plt.plot(radlm, crden_ana, label=dclabel,lw=2,color=color,ls=lsn)
                    plt.plot(radlm, crden, lw=2,color=color,ls='none',marker='o')
                if thermalneed>0:
                    thden = (thl[1:]-thl[:-1])*erg_in_eV/vol
                    x0 = [0.1,-1.0,0.0]
                    xdata=np.log10(radlm)
                    ydata=np.log10(thden)
                    outfit=optimize.curve_fit(func, xdata, ydata, x0)
                    afit=outfit[0][0]
                    bfit=outfit[0][1]
                    cfit=outfit[0][2]
                    thden = (thl[1:]-thl[:-1])*erg_in_eV/vol
                    plt.plot(radlm, thden, lw=1,color=color,ls=lsn)
                    #plt.plot(radlm, thden, label=dclabel,lw=2,color=color)
                filename='CRplot/crdenr/crdenr_'+fmeat+'.pdf'
            if wanted=='crdpr':
                radlm = (radl[1:]+radl[:-1])/2.
                radlmm = (radlm[1:]+radlm[:-1])/2.
                vol = 4./3.*np.pi*(radl[1:]*radl[1:]*radl[1:]-radl[:-1]*radl[:-1]*radl[:-1])*kpc_in_cm*kpc_in_cm*kpc_in_cm
                if havecr>0:
                    crden = (crel[1:]-crel[:-1])*erg_in_eV/vol
                    x0 = [0.1,-1.0,0.0]
                    xdata=np.log10(radlm)
                    ydata=np.log10(crden)
                    outfit=optimize.curve_fit(func, xdata, ydata, x0)
                    afit=outfit[0][0]
                    bfit=outfit[0][1]
                    cfit=outfit[0][2]
                    print 'afit,bfit,cfit', afit,bfit,cfit
                    dpcr = -d10func(radlmm,afit,bfit,cfit)*(CRgamma-1.0)/erg_in_eV/kpc_in_cm
                    plt.plot(radlmm, dpcr, label=dclabel,lw=2,color=color,ls=lsn)
                    #dpcr = (crden[1:]-crden[:-1])*CRgamma/(radlm[1:]-radlm[:-1])
                                if thermalneed>0:
                                        thden = (thl[1:]-thl[:-1])*erg_in_eV/vol
                    x0 = [1.0,-1.0,0.0]
                    xdata=np.log10(radlm)
                    ydata=np.log10(thden)
                    outfit=optimize.curve_fit(func, xdata, ydata, x0)
                    afit=outfit[0][0]
                    bfit=outfit[0][1]
                    cfit=outfit[0][2]
                    print 'afit,bfit,cfit', afit,bfit,cfit
                    #dpth = (thden[1:]-thden[:-1])*GAMMA/(radlm[1:]-radlm[:-1])
                    dpth = -d10func(radlmm,afit,bfit,cfit)*(GAMMA-1.0)/erg_in_eV/kpc_in_cm
                                        plt.plot(radlmm, dpth, lw=1,color=color,ls=lsn)
                                filename='CRplot/crdpr/crdpr_'+fmeat+'.pdf'
        plt.yscale('log')
    plt.xscale('log')
        plt.xlabel(r'$r\; [{\rm kpc}]$')
    if wanted=='crer':
        plt.ylabel(r'enclosed $E_{\rm cr} [{\rm erg}]$')
    if wanted=='crdenr':
        plt.ylabel(r'$e_{\rm cr} [{\rm eV/cm^3}]$')
        if wanted=='crdpr':
                plt.ylabel(r'${\rm d}P_{\rm cr}/{\rm d} r [{\rm erg/cm^2}]$')
        if wanted=='denr':
                plt.ylabel(r'$\rho [{\rm g/cm^3}]$')
        if wanted=='rhovr':
                plt.ylabel(r'$\rho v [{\rm g/cm^2/s}]$')
        if wanted=='rhovr_denr':
                plt.ylabel(r'$v [{\rm cm/s}]$')
        plt.legend(loc='best')
    print 'filename', filename
        plt.savefig(filename)
        plt.clf()


if wanted=='denrtime':
        rcParams['figure.figsize'] = 5,4
        withinr=2.0
        linelabelneed=0
        labcount=0
        newlabelneed=1
        for runtodo in dirneed:
        timel=[]
        denl=[]
                for i in range(startno,Nsnap,snapsep):
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        maindir=info['maindir']
                        halostr=info['halostr']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        cosmo=info['cosmo']
                        color=info['color']
                        withinRv=info['withinRv']
                        usepep=info['usepep']
                        beginno=info['beginno']
                        finalno=info['finalno']
                        firever=info['firever']
                        initsnap=info['initsnap']
                        haveB=info['haveB']
                        M1speed=info['M1speed']
                        Rvirguess=info['Rvir']
                        newlabel=info['newlabel']
                        ptitle=title
                        labelneed=dclabel
                        if newlabelneed==1:
                                labelneed="\n".join(wrap(newlabel,17))
                        if runtitle=='SMC':
                                ptitle='Dwarf'
                        elif runtitle=='SBC':
                                ptitle='Starburst'
                        elif runtitle=='MW':
                                ptitle=r'$L\star$ Galaxy'
                        lsn='solid'
                        if haveB>0:
                                lsn='dashed'
                        if havecr==0 and wanted=='avecrdenr':
                                continue
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gp = G['p']
                        Gv = G['v']
                        Gu = G['u']
                        Gm = G['m']
                        Gvx = Gv[:,0]
                        Gvy = Gv[:,1]
                        Gvz = Gv[:,2]
                        GEint = Gu*km_in_cm*km_in_cm*Gm*1e10*Msun_in_g
                        #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                        Gz = Gp[:,2]
                        Gx = Gp[:,0]
                        Gy = Gp[:,1]
                        dr = withinr/nogrid
                        from crtestfunction import findcenz
                        datasup=1
                        xcen = findcenz(runtodo,Nsnap,withinr=withinr,dir='x',datasup=datasup,Gx=Gx,Gy=Gy,Gz=Gz,Gm=Gm)
                        ycen = findcenz(runtodo,Nsnap,withinr=withinr,dir='y',datasup=datasup,Gx=Gx,Gy=Gy,Gz=Gz,Gm=Gm)
                        zcen = findcenz(runtodo,Nsnap,withinr=withinr,dir='z',datasup=datasup,Gx=Gx,Gy=Gy,Gz=Gz,Gm=Gm)
                        Gz = Gz-zcen;Gy = Gy-ycen; Gx = Gx-xcen
            irad=1
            cutxy = (Gx*Gx+Gy*Gy < dr*irad*dr*irad)
            cutz = np.absolute(Gz-med)< maxlength/2.
            cut = cutxy*cutz
            Gmcut = Gm[cut]
            Gm_in_sun = np.sum(Gmcut)*1e10
            vol = np.power(dr*irad*kpc_in_cm,3.)*4./3.*np.pi
            den = Gm_in_sun*Msun_in_g/vol/protonmass_in_g
            timel=np.append(timel,i)
            denl = np.append(denl,den)
        plt.plot(timel,denl,label=labelneed,lw=2,ls=lsn,color=color)
        plt.yscale('log')
        plt.xlabel(r'$ t\; [{\rm Myr}]$',fontsize=14)
        plt.ylabel(r'$\bar{n}_{\rm ISM}\; [{\rm cm^{-3}}]$',fontsize=14)
        filename='CRplot/denrtime/denrtime_'+runtodo+'.pdf'
        plt.legend(loc='best',fontsize=14)
        plt.tight_layout()
        ptitle = 'tmax = '+str(timel[np.argmax(denl)])+' Myr'
        plt.title(ptitle,fontsize=16)
        plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
        plt.savefig(filename)
        plt.clf()



if wanted=='gasdenz' or wanted == 'gasdenx' or wanted == 'gasdeny':
    withinr=10.
        for runtodo in dirneed:
                for i in [Nsnap]:
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        maindir=info['maindir']
                        halostr=info['halostr']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        cosmo=info['cosmo']
                        color=info['color']
                        withinRv=info['withinRv']
                        usepep=info['usepep']
                        beginno=info['beginno']
                        finalno=info['finalno']
                        firever=info['firever']
                        initsnap=info['initsnap']
                        haveB=info['haveB']
                        M1speed=info['M1speed']
                        Rvirguess=info['Rvir']
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gp = G['p']
                        Gu = G['u']
                        Gm = G['m']
                        cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                        Neb = G['ne']
            from crtestfunction import *
                        Gz = Gp[:,2]
                        Gx = Gp[:,0]
                        Gy = Gp[:,1]
            datasup=1
                        if wanted=='gasdenz':
                                dir='z'
                        if wanted=='gasdeny':
                                dir='y'
                        if wanted=='gasdenx':
                                dir='x'
            cen = findcenz(runtodo,Nsnap,withinr=withinr,dir=dir,datasup=datasup,Gx=Gx,Gy=Gy,Gz=Gz,Gm=Gm)
                        dz = maxlength/nogrid
                        Gnism_in_cm_3l=[]
                        zl =[]
            if wanted=='gasdenz':
                px=Gx;py=Gy;pz=Gz-cen
            if wanted=='gasdeny':
                px=Gz;py=Gx;pz=Gy-cen
            if wanted=='gasdenx':
                px=Gy;py=Gz;pz=Gx-cen
                        for iz in range(nogrid):
                                cutxy = px*px+py*py < withinr*withinr
                                cutz = (pz> dz*iz-maxlength/2.) & (pz<dz*(iz+1)-maxlength/2.)
                                cut = cutxy*cutz
                                Gmcut = Gm[cut]
                                Gm_in_g = Gmcut*1e10*2e33
                                shellvol_in_cm3 = np.pi*dz*np.power(withinr,2)*3.086e21*3.086e21*3.086e21
                                Grho_in_g_cm_3 = Gm_in_g/shellvol_in_cm3
                                Gnism_in_cm_3 = np.sum(1.0/protonmass_in_g*Grho_in_g_cm_3)
                                Gnism_in_cm_3l = np.append(Gnism_in_cm_3l, Gnism_in_cm_3)
                                zl = np.append(zl, dz*(iz+0.5)-maxlength/2.)
                        plt.plot(zl, Gnism_in_cm_3l, label=runtodo)
    if wanted=='gasdenz':
        plt.xlabel('z [kpc]')
        filename='CRplot/gasdenz/gasdensity_'+fmeat+'_z.pdf'
        if wanted=='gasdeny':
                plt.xlabel('y [kpc]')
                filename='CRplot/gasdeny/gasdensity_'+fmeat+'_y.pdf'
        if wanted=='gasdenx':
                plt.xlabel('x [kpc]')
                filename='CRplot/gasdenx/gasdensity_'+fmeat+'_x.pdf'
        plt.ylabel(r'$n_{\rm ISM} [{\rm cm^{-3}}]$')
        plt.legend(loc='best')
        #plt.title('cylinder centered at disk with a radius ' + str(withinr) + ' kpc')
    print filename
        plt.savefig(filename)
        plt.clf()


if wanted=='gmz' or wanted=='crecumz' or wanted=='crecumr' or wanted=='Begyz' or wanted=='ethz' or wanted=='kez':
    rcParams['figure.figsize'] = 5,4
    withinr=100.0
    #withinr=20.
    #maxlength=10.0
    minlength=0.5
    maxlength=100.
    nogrid=500
    normalized =1
    newlabelneed=0
    strlabelneed=0
    legendneed=0
    massloss=0.3
        for runtodo in dirneed:
        info=outdirname(runtodo, Nsnap)
        havecr = info['havecr']
        haveB = info['haveB']
        Nsnapstring = info['Nsnapstring']
        the_prefix = info['the_prefix']
        the_suffix = info['the_suffix']
        the_snapdir = info['the_snapdir']
        if wanted=='crecumz':
            if havecr==0:
                continue
        if wanted=='Begyz':
            if haveB==0:
                continue
                zlist=np.linspace(minlength,maxlength,num=nogrid)
                if wanted=='crecumz' or wanted=='crecumr':
                        crel=zlist*0.
        if wanted=='Begyz':
            Begyl=zlist*0.
        if wanted=='ethz':
            ethl = zlist*0.
        if wanted=='kez':
            kel = zlist*0.
        Gm_in_sunl=zlist*0.
        nooftimes=0
        Esncr=1.
        if (wanted=='crecumz' or wanted=='crecumr') and normalized == 1:
            Nsnapstring0 = '000'
            S0 = readsnapcr(the_snapdir, Nsnapstring0, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
            S0mass = np.sum(S0['m']*1e10)
            print 'Nsnapstring', Nsnapstring
            S = readsnapcr(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
            Smass = np.sum(S['m']*1e10)
            Esncr = 0.1*(Smass-S0mass)*1e51*0.01*(1.+massloss) #total CR energy from SNII in erg
            print 'Esncr', Esncr
                for i in range(startno,Nsnap,snapsep):
                        info=outdirname(runtodo, i)
                        M1speed=info['M1speed']
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        haveB=info['haveB']
            havecr = info['havecr']
            color=info['color']
            newlabel=info['newlabel']
            strlabel=info['strlabel']
            ptitle=title
            if runtitle=='SMC':
                ptitle='Dwarf'
            elif runtitle=='SBC':
                ptitle='Starburst'
            elif runtitle=='MW':
                ptitle=r'$L\star$ Galaxy'
            labelneed=dclabel
            if newlabelneed==1:
                labelneed="\n".join(wrap(newlabel,17)) 
                        if strlabelneed==1:
                                labelneed="\n".join(wrap(strlabel,17)) 
            print 'labelneed,newlabel', labelneed,newlabel
                        G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gp = G['p']
            Gv = G['v']
                        Gu = G['u']
                        Gm = G['m']
            Grho = G['rho']
            if wanted=='crecumz' or wanted=='crecumr':
                cregy = G['cregy']
                print 'CR energy in erg', np.sum(cregy)*1e10*Msun_in_g*km_in_cm*km_in_cm
            if wanted=='Begyz':
                                Bfield = G['B']
                                Bx = Bfield[:,0]
                                By = Bfield[:,1]
                                Bz = Bfield[:,2]
                                B2 = Bx*Bx+By*By+Bz*Bz
                                Begy = B2/8./np.pi*Gm/Grho*kpc_in_cm*kpc_in_cm*kpc_in_cm
            if wanted=='kez':
                Gvx = Gv[:,0]
                Gvy = Gv[:,1]
                Gvz = Gv[:,2]
                ke = 0.5*Gm*(Gvx*Gvx+Gvy*Gvy+Gvz*Gvz)*1e10*Msun_in_g*km_in_cm*km_in_cm
                        Gz = Gp[:,2]
                        Gx = Gp[:,0]
                        Gy = Gp[:,1]
            from crtestfunction import *
                        datasup=1
                        xcen = findcenz(runtodo,Nsnap,withinr=withinr,dir='x',datasup=datasup,Gx=Gx,Gy=Gy,Gz=Gz,Gm=Gm)
                        ycen = findcenz(runtodo,Nsnap,withinr=withinr,dir='y',datasup=datasup,Gx=Gx,Gy=Gy,Gz=Gz,Gm=Gm)
                        zcen = findcenz(runtodo,Nsnap,withinr=withinr,dir='z',datasup=datasup,Gx=Gx,Gy=Gy,Gz=Gz,Gm=Gm)
                        Gz = Gz-zcen;Gy = Gy-ycen; Gx = Gx-xcen
                        dz = maxlength/nogrid
            if normalized==1:
                                cutxy = Gx*Gx+Gy*Gy < withinr
                                cutz = np.absolute(Gz)<dz*nogrid
                cut= cutxy*cutz
                Gmtot = np.sum(Gm[cut]*1e10)
                        for iz in range(nogrid):
                if wanted=='crecumr':
                    cut = np.sqrt(Gx*Gx+Gy*Gy+Gz*Gz)<zlist[iz]
                else:
                    cutxy = Gx*Gx+Gy*Gy < withinr
                    cutz = np.absolute(Gz)<zlist[iz]
                    cut = cutxy*cutz
                                Gmcut = Gm[cut]
                                Gm_in_sun = Gmcut*1e10
                if wanted=='crecumz' or wanted=='crecumr':
                    cregycut = cregy[cut]*1e10*Msun_in_g*km_in_cm*km_in_cm
                    crel[iz] += np.sum(cregycut)
                if wanted=='Begyz':
                    Begycut = Begy[cut]
                    Begyl[iz] += np.sum(Begycut)
                if wanted=='ethz':
                    ethcut = Gu[cut]*km_in_cm*km_in_cm*Gm_in_sun*Msun_in_g
                    ethl[iz] += np.sum(ethcut)
                if wanted=='kez':
                    kecut = ke[cut]
                    kel[iz] += np.sum(kecut)
                                Gm_in_sunl[iz] += np.sum(Gm_in_sun)
            nooftimes+=1
        if wanted=='crecumz' or wanted=='crecumr':
            crel=crel/nooftimes
        if wanted=='Begyz':
            Begyl=Begyl/nooftimes
        if wanted=='ethz':
            ethl=ethl/nooftimes
        if wanted=='kez':
            kel=kel/nooftimes
        Gm_in_sunl=Gm_in_sunl/nooftimes
        if haveB==0:
            lsn='solid'
        else:
            lsn='dashed'
        if wanted == 'gmz':
            if normalized==0:
                plt.plot(zlist, Gm_in_sunl, label=labelneed, ls=lsn,lw=2,color=color)
            else:
                plt.plot(zlist, Gm_in_sunl/Gmtot/(1.+massloss), label=labelneed, ls=lsn,lw=2,color=color)
        if wanted == 'crecumz' or wanted=='crecumr':
            plt.plot(zlist, crel/Esncr, label=labelneed, ls=lsn,lw=2,color=color)
        if wanted == 'Begyz':
            plt.plot(zlist, Begyl, label=labelneed, ls=lsn,lw=2,color=color)
                if wanted == 'ethz':
                        plt.plot(zlist, ethl, label=labelneed, ls=lsn,lw=2,color=color)
        if wanted == 'kez':
            plt.plot(zlist, kel, label=labelneed, ls=lsn,lw=2,color=color)
    plt.yscale('log')
    plt.xscale('log')
    if wanted=='crecumr':
                plt.xlabel(r'$r\; [{\rm kpc}]$',fontsize=18)
    else:
        plt.xlabel(r'$z\; [{\rm kpc}]$',fontsize=18)
    if wanted == 'gmz':
        if normalized==1:
            plt.ylabel(r'$M_{\rm gas}(<z)/M_0$',fontsize=18)
        else:
            plt.ylabel(r'$M_{\rm gas}(<z) [ M_\odot]$',fontsize=18)
        if runtitle=='SMC' and legendneed==1:
            #plt.legend(loc=2,fontsize=8,ncol=3,frameon=False)
            plt.legend(loc=4,fontsize=8,ncol=3)
    if wanted == 'crecumz' or wanted=='crecumr':
        if normalized==1:
            plt.ylabel(r'$E_{\rm cr}(<z)/E_{\rm cr,SN}(tot)$',fontsize=18)
        else:
            plt.ylabel(r'$E_{\rm cr}(<z)\;[{\rm erg}]$',fontsize=18)
        if wanted=='crecumr':
            if normalized==1:
                plt.ylabel(r'$E_{\rm cr}(<r)/E_{\rm cr,SN}(tot)$',fontsize=18)
            else:   
                plt.ylabel(r'$E_{\rm cr}(<r)\;[{\rm erg}]$',fontsize=18)
                if runtitle=='SMC' and legendneed==1:
                        #plt.legend(loc=2,fontsize=8,ncol=3,frameon=False)
                        plt.legend(loc=4,fontsize=8,ncol=3)
    if wanted == 'Begyz':
        plt.ylabel(r'$E_{\rm B}(<z)\;[{\rm erg}]$',fontsize=18)
        if wanted == 'ethz':
                plt.ylabel(r'$E_{\rm th}(<z)\;[{\rm erg}]$',fontsize=18)
        if wanted == 'kez':
                plt.ylabel(r'$E_{\rm KE}(<z)\;[{\rm erg}]$',fontsize=18)
    plt.tight_layout()
    if wanted == 'gmz':
        if normalized==0:
            filename = 'CRplot/gmz/gasmass_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        else:
            filename = 'CRplot/gmz/gasmass_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'_normal.pdf'
    if wanted == 'crecumz':
        filename = 'CRplot/crecumz/crecumz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted == 'crecumr':
                filename = 'CRplot/crecumr/crecumr_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted == 'Begyz':
                filename = 'CRplot/Begyz/Begyz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted == 'ethz':
                filename = 'CRplot/ethz/ethz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted == 'kez':
                filename = 'CRplot/kez/kez_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    print 'filename', filename
    plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
    plt.title(ptitle)
    plt.savefig(filename)
        plt.clf()


if wanted=='crez' or wanted=='dpcrz' or wanted=='dpcrz_rhog':
    rcParams['figure.figsize'] = 5,4
    resoneed=0
    for runtodo in dirneed:
        maxlength = 40
        nogrid=20
        withinr=5
        zl=np.linspace(2,maxlength,num=nogrid)
        dz = maxlength/nogrid
        crden=zl*0.0
        Gmden=zl*0.0
        mtotl=zl*0.0
        numoftimes=0
        info=outdirname(runtodo, Nsnap)
        havecr=info['havecr']
        if havecr==0:
            continue
        for i in range(startno,Nsnap,snapsep):
            info=outdirname(runtodo, i)
            M1speed=info['M1speed']
            rundir=info['rundir']
            runtitle=info['runtitle']
            slabel=info['slabel']
            snlabel=info['snlabel']
            dclabel=info['dclabel']
            resolabel=info['resolabel']
            the_snapdir=info['the_snapdir']
            Nsnapstring=info['Nsnapstring']
            Fcal=info['Fcal']
            iavesfr=info['iavesfr']
            timestep=info['timestep']
            haveB=info['haveB']
            G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
            Gp = G['p']
            #Grho = G['rho']
            Gu = G['u']
            Gm = G['m']
            cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
            cregy *= 1e10*Msun_in_g*km_in_cm*km_in_cm
            #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
            Gz = Gp[:,2]
            Gx = Gp[:,0]
            Gy = Gp[:,1]
            if haveB==0:
                lsn='solid'
            else:
                lsn='dashed'
            #if M1speed>999:
            #   lsn='dotted'
            #if M1speed>1999:
            #   lsn='dashdot'
            for iz in range(len(zl)-1):
                cutxy = (Gx*Gx+Gy*Gy < withinr*withinr)
                cutzu = Gz<zl[iz+1]
                cutzd = Gz>zl[iz]
                cutz = cutzu*cutzd
                cut = cutxy*cutz
                crecut = cregy[cut]
                crden[iz] += np.sum(crecut)/(zl[iz+1]-zl[iz])/kpc_in_cm/(np.pi*withinr*withinr*kpc_in_cm*kpc_in_cm)
                if wanted=='dpcrz_rhog':
                    Gmcut=Gm[cut]*1e10*Msun_in_g
                    Gmden[iz] += np.sum(Gmcut)/(zl[iz+1]-zl[iz])/kpc_in_cm/(np.pi*withinr*withinr*kpc_in_cm*kpc_in_cm)
                    mintot = 0
                    for ipa in [1,2,3,4]:
                        Pa = readsnap(the_snapdir, Nsnapstring, ipa, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Pp = Pa['p']
                        Pm = Pa['m']
                        Pr = np.sqrt(Pp[:,0]*Pp[:,0]+Pp[:,1]*Pp[:,1]+Pp[:,2]*Pp[:,2])
                        cutpr = Pr<zl[iz]
                        Pmtot = np.sum(Pm[cutpr]*1e10*Msun_in_g) #in g
                        mintot += Pmtot
                    mtotl[iz] += mintot
            numoftimes+=1
        crden = crden/numoftimes
        Gmden = Gmden/numoftimes
        mtotl = mtotl/numoftimes
        labelneed=dclabel
                if resoneed==1:
                        if resolabel=='llr':
                                labelneed='low res'
                                lsn = 'solid'
                        if resolabel=='lr':
                                labelneed='std res'
                                lsn = 'dashed'
                        if resolabel=='mr':
                                labelneed='high res'
                                lsn = 'dashdot'
        if wanted=='crez':
            plt.plot(zl[:-1], crden[:-1]*erg_in_eV, label=labelneed,lw=2,ls=lsn)
        if wanted=='dpcrz':
            pcrden = (CRgamma-1.0)*crden
            zlm = (zl[1:]+zl[:-1])/2.
            dpcr = (pcrden[1:]-pcrden[:-1])/(zl[1:]-zl[:-1])/kpc_in_cm
            plt.plot(zlm[:-1], -dpcr[:-1], label=labelneed,lw=2,ls=lsn)
        if wanted=='dpcrz_rhog':
                        pcrden = (CRgamma-1.0)*crden
            print 'Gmden, mtotl', Gmden, mtotl
            rhog = Gmden*mtotl*NewtonG_cgs/(zl*kpc_in_cm)/(zl*kpc_in_cm)
                        zlm = (zl[1:]+zl[:-1])/2.
                        dpcr = (pcrden[1:]-pcrden[:-1])/(zl[1:]-zl[:-1])/kpc_in_cm
                        plt.plot(zlm, -dpcr/rhog[:-1], label=labelneed,lw=2,ls=lsn)
    if wanted=='crez':
        plt.yscale('log')
        plt.xlabel(r'$z [{\mr kpc}]$', fontsize=18)
        plt.ylabel(r'$e_{\rm cr} [{\rm eV/cm^3}]$', fontsize=18)
        plt.legend(loc='best')
        #plt.title('cylinder centered at disk with a radius ' + str(withinr) + ' kpc')
        filename='CRplot/crez/crez_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    if wanted=='dpcrz':
        plt.yscale('log')
        plt.xlabel('z [kpc]')
        plt.ylabel(r'${\rm d}p_{\rm cr}/{\rm d}z [{\rm erg/cm^4}]$')
        plt.legend(loc='best')
        plt.title('cylinder centered at disk with a radius ' + str(withinr) + ' kpc')
        filename='CRplot/dpcrz/dpcrz_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted=='dpcrz_rhog':
                plt.yscale('log')
                plt.xlabel('z [kpc]')
                plt.ylabel(r'${\rm d}p_{\rm cr}/{\rm d}z/(\rho g)$')
                plt.legend(loc='best')
                plt.title('cylinder centered at disk with a radius ' + str(withinr) + ' kpc')
                filename='CRplot/dpcrz_rhog/dpcrz_rhog_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    print 'filename', filename
    plt.savefig(filename)
    plt.clf()




if wanted=='vz' or wanted=='pz' or wanted=='dpz' or wanted=='vzvwed' or wanted=='Tzmwed':
    rcParams['figure.figsize'] = 5,4
    nogrid=20
    maxlength=40.
    minlength=2.
    withinr=20.
    needTcut = 0
    needvcut = 0
    Tlabel=''
    Tno=1
    lw=2
    vcut=0.
    if needTcut==1:
        Tlabel='Tcut'
        Tcut = 1.0e5 #K
        Tno=2
    if needvcut==1:
        Tlabel+='vcut'
        vcut = 30. #km/s
        
        for runtodo in dirneed:
        for tcount in range(Tno):
            zl=np.linspace(minlength,maxlength,num=nogrid)
            pzl=zl*0.0
            Gml=zl*0.0
            Tml=zl*0.0
            numoftimes=0
            for i in range(startno,Nsnap,snapsep):
                info=outdirname(runtodo, i)
                rundir=info['rundir']
                runtitle=info['runtitle']
                slabel=info['slabel']
                snlabel=info['snlabel']
                dclabel=info['dclabel']
                resolabel=info['resolabel']
                the_snapdir=info['the_snapdir']
                Nsnapstring=info['Nsnapstring']
                havecr=info['havecr']
                Fcal=info['Fcal']
                iavesfr=info['iavesfr']
                timestep=info['timestep']
                haveB=info['haveB']
                color=info['color']
                ptitle=title
                if runtitle=='SMC':
                    ptitle='Dwarf'
                elif runtitle=='SBC':
                    ptitle='Starburst'
                elif runtitle=='MW':
                    ptitle=r'$L\star$ Galaxy'
                G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                Gp = G['p']
                Gv = G['v']
                #Grho = G['rho']
                Gu = G['u']
                if needTcut==1 or wanted=='Tzmwed':
                    Tb = Gu
                    Neb = G['ne']
                    rho = G['rho']
                Gm = G['m']*1e10
                #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                Gz = Gp[:,2]
                Gx = Gp[:,0]
                Gy = Gp[:,1]
                Gvz = Gv[:,2]
                if haveB>0:
                    lsn='dashed'
                else:
                    lsn='solid'
                for iz in range(nogrid-1):
                    cutxy = (Gx*Gx+Gy*Gy < withinr*withinr)
                    cutzu = Gz<zl[iz+1]
                    cutzd = Gz>zl[iz]
                    Gvzcut = Gvz*Gz > vcut #outflow gas only
                    cut = cutxy*cutzu*cutzd*Gvzcut
                    if needTcut==1:
                        TrueTemp, converted_rho  = SF.convertTemp(Tb, Neb, rho)
                        if tcount==0:
                            cutT = TrueTemp>Tcut
                            lw = 1
                        elif tcount==1:
                            cutT = TrueTemp<Tcut
                            lw = 2
                        cut = cut*cutT
                    pzcut = Gvz[cut]*Gm[cut]
                    if wanted=='Tzmwed':
                        TrueTemp, converted_rho  = SF.convertTemp(Tb, Neb, rho)
                        Tml[iz] += np.sum(TrueTemp[cut]*Gm[cut])
                    #print 'maxv', np.amax(Gvz[cut])    
                    if wanted=='pz' or wanted=='vz':
                        pzl[iz] += np.sum(np.absolute(pzcut))
                    if wanted=='vz' or wanted=='Tzmwed':
                        Gml[iz] += np.sum(Gm[cut])
                numoftimes+=1
            if wanted=='vz':
                vzl=pzl/Gml
                print 'zl', zl
                print 'vzl', vzl
                print 'pzl', pzl
                print 'Gml', Gml
                plt.plot(zl[:-1], vzl[:-1], label=dclabel,ls=lsn,lw=lw,color=color)
                        if wanted=='Tzmwed':
                                Tzl=Tml/Gml
                                plt.plot(zl[:-1], Tzl[:-1], label=dclabel,ls=lsn,lw=lw,color=color)
                        if wanted=='vzvwed':
                dz = zl[1]-zl[0]
                dV = dz*np.pi*withinr*withinr
                                vzl= pzl*dV
                                plt.plot(zl[:-1], vzl[:-1], label=dclabel,ls=lsn,lw=lw,color=color)
            if wanted=='pz':
                pzl=pzl/numoftimes
                plt.plot(zl[:-1], pzl[:-1], label=dclabel,ls=lsn,lw=lw,color=color)
            if wanted=='dpz':
                pzl=pzl/numoftimes
                dpzl = (pzl[:-2]-pzl[1:-1])*Msun_in_g/(zl[:-2]-zl[1:-1])/kpc_in_cm*km_in_cm
                zlm = (zl[:-2]+zl[1:-1])/2.0
                plt.plot(zlm, dpzl, label=dclabel,ls=lsn,lw=lw,color=color)
        plt.xlabel(r'${\rm z [kpc]}$',fontsize=18)
    if wanted=='vz': 
        #plt.yscale('log')
        plt.ylabel(r'$v_{\rm z} [{\rm km/s}]$',fontsize=18)
        if wanted=='vzvwed':
                #plt.yscale('log')
                plt.ylabel(r'$v_{\rm z} [{\rm km/s}]$',fontsize=18)
    if wanted=='Tzmwed':
        plt.yscale('log')
        plt.ylabel(r'$T [{\rm K}]$',fontsize=18)
    if wanted=='pz':
        plt.yscale('log')
        plt.ylabel(r'$p_{\rm z} [{\rm M_\odot km/s}]$', fontsize=18)
    if wanted=='dpz':
        plt.ylabel(r'${\rm d}p_{\rm z}/{\rm d}z [{\rm g/s}]$', fontsize=18)
        #plt.legend(loc='best')
        #plt.title('cylinder centered at disk with a radius ' + str(15) + ' kpc')
    plt.tight_layout()
    if wanted=='vz':
        totalname='CRplot/vz/vz_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+Tlabel+'.pdf'
        if wanted=='vzvwed':
                totalname='CRplot/vzvwed/vzvwed_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+Tlabel+'.pdf'
        if wanted=='Tzmwed':
                totalname='CRplot/Tzmwed/Tzmwed_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+Tlabel+'.pdf'
        if wanted=='pz':
                totalname='CRplot/pz/pz_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+Tlabel+'.pdf'
        if wanted=='dpz':
                totalname='CRplot/dpz/dpz_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+Tlabel+'.pdf'
    print 'totalname', totalname
    plt.title(ptitle)
    plt.savefig(totalname)
        plt.clf()






if wanted=='nismcumcre':
        for runtodo in dirneed:
                timel=[]
                for i in [Nsnap]:
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        cosmo=info['cosmo']
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
            if cosmo==1:
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,cosmological=1,h0=1)
            else:
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        timeneed=i*0.98*1e9 #in yr
                        Grho = G['rho'] #1e10Msun per kpc^3
                        cregy = G['cregy']*1e10*solar_mass_in_g*km_in_cm*km_in_cm #original cosmic ray energy in 1e10Msun km^2/sec^2
                        Neb = G['ne']
                        timel=np.append(timel,timeneed)
                Gnism_in_cm_3 = (0.78+0.22*Neb*0.76)/proton_mass_in_g*Grho*1e10*solar_mass_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm
                #tpi_in_yr = 2e5/Gnism_in_cm_3*250.0 #pi decay time in yr
                #Lgammagev = cregy_in_erg/tpi_in_yr/sec_in_yr*betapi/nopi_per_gamma #in erg/s
                LogGnism = np.log10(Gnism_in_cm_3)
                LogGnxaxis = np.linspace(-4,4,num=nogrid)
                dx =LogGnxaxis[1]-LogGnxaxis[0]
                Ecrl = []
                for inism in range(nogrid-1):
                        cutg = (LogGnism > LogGnxaxis[inism]) 
                        Ecrcut = cregy[cutg]
                        Ecrl = np.append(Ecrl, np.sum(Ecrcut))
                plt.plot((LogGnxaxis[1:]+LogGnxaxis[:-1])/2.,Ecrl, label=runtodo)
        plt.legend(loc='best')
        plt.yscale('log')
        plt.xlabel(r'$\log (n_{\rm ISM}[{\rm cm^{-3}}])$')
        plt.ylabel(r'$E_{cr} (>n_{\rm ISM})[{\rm erg}]$')
        plt.savefig('CRplot/nismcumcr_'+fmeat+'.pdf')
        plt.clf()




if wanted=='credecaytime':
        for runtodo in dirneed:
                for i in [Nsnap]:
                        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr,timestep=outdirname(runtodo, i)
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Grho = G['rho']
                        cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                        Neb = G['ne']
                        Gnism_in_cm_3 = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21
                        tpi_in_yr = 2e5/Gnism_in_cm_3*250.0 #pi decay time in yr
                        cregy_in_erg = cregy*solar_mass_in_g*1e10*km_in_cm*km_in_cm
                        #Lgammagev = cregy_in_erg/tpi_in_yr/sec_in_yr*betapi/nopi_per_gamma #in erg/s
                        logtpi = np.log10(tpi_in_yr)
                        logtpixaxis = np.linspace(5,14,num=nogrid)
            dx = logtpixaxis[1]-logtpixaxis[0]
                        cregy_in_ergl = []
                        for inism in range(nogrid-1):
                                cutt = (logtpi > logtpixaxis[inism]) & (logtpi < logtpixaxis[inism+1])
                                cregy_in_ergcut = cregy_in_erg[cutt]
                                cregy_in_ergl = np.append(cregy_in_ergl, np.sum(cregy_in_ergcut))
                        plt.plot((logtpixaxis[1:]+logtpixaxis[:-1])/2.,cregy_in_ergl/dx, label=runtodo)
        plt.legend(loc='best')
    plt.yscale('log')
        plt.xlabel(r'$\log (t_{\rm \pi}[{\rm yr}])$')
        #plt.ylabel(r'$\frac{\mathrm{d} E_{\rm cr}}{\mathrm{d} n_{\rm ISM}}\Delta \log n_{\rm ISM}[{\rm erg/cm^3}]$')
        plt.ylabel(r'$\mathrm{d} E_{\rm cr}/\mathrm{d} \log (t_{\rm \pi})[{\rm erg}]$')
        plt.savefig('credecaytime_'+fmeat+'.pdf')
        plt.clf()


if wanted=='crecumdecaytime':
        for runtodo in dirneed:
                for i in [Nsnap]:
                        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr,timestep=outdirname(runtodo, i)
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Grho = G['rho']
                        cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                        Neb = G['ne']
                        Gnism_in_cm_3 = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21
                        tpi_in_yr = 2e5/Gnism_in_cm_3*250.0 #pi decay time in yr
                        cregy_in_erg = cregy*solar_mass_in_g*1e10*km_in_cm*km_in_cm
                        #Lgammagev = cregy_in_erg/tpi_in_yr/sec_in_yr*betapi/nopi_per_gamma #in erg/s
                        logtpi = np.log10(tpi_in_yr)
                        logtpixaxis = np.linspace(5,14,num=nogrid)
                        cregy_in_ergl = []
                        for inism in range(nogrid-1):
                                cutt =  (logtpi > logtpixaxis[inism])
                                cregy_in_ergcut = cregy_in_erg[cutt]
                                cregy_in_ergl = np.append(cregy_in_ergl, np.sum(cregy_in_ergcut))
                        plt.plot((logtpixaxis[1:]+logtpixaxis[:-1])/2.,cregy_in_ergl, label=runtodo)
        plt.legend(loc='best')
        plt.yscale('log')
        plt.xlabel(r'$\log (t_{\rm \pi}[{\rm yr}])$')
        #plt.ylabel(r'$\frac{\mathrm{d} E_{\rm cr}}{\mathrm{d} n_{\rm ISM}}\Delta \log n_{\rm ISM}[{\rm erg/cm^3}]$')
        plt.ylabel(r'$E_{\rm cr}(>t_{\rm \pi})[{\rm erg}]$')
        plt.savefig('crecumdecaytime_'+fmeat+'.pdf')
        plt.clf()


if wanted=='cregamma':
        for runtodo in dirneed:
                for i in [Nsnap]:
                        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr,timestep=outdirname(runtodo, i)
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Grho = G['rho']
                        cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                        Neb = G['ne']
                        Gnism_in_cm_3 = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21
                        tpi_in_yr = 2e5/Gnism_in_cm_3*250.0 #pi decay time in yr
                        cregy_in_erg = cregy*solar_mass_in_g*1e10*km_in_cm*km_in_cm
                        Lgammagev = cregy_in_erg/tpi_in_yr/yr_in_sec*betapi/nopi_per_gamma #in erg/s
            print 'np.sum(Lgammagev)', np.sum(Lgammagev)
                        logLgamma = np.log10(Lgammagev)
                        loggammaxaxis = np.linspace(25,42,num=nogrid)
            dx = loggammaxaxis[1]-loggammaxaxis[0]
                        cregy_in_ergl = []
                        for inism in range(nogrid-1):
                                cutt = ( logLgamma> loggammaxaxis[inism]) & (logLgamma < loggammaxaxis[inism+1])
                                cregy_in_ergcut = cregy_in_erg[cutt]
                                cregy_in_ergl = np.append(cregy_in_ergl, np.sum(cregy_in_ergcut))
                        plt.plot((loggammaxaxis[1:]+loggammaxaxis[:-1])/2.,cregy_in_ergl/dx, label=runtodo)
        plt.legend(loc='best')
    plt.yscale('log')
        plt.xlabel(r'$\log (L_{\rm \gamma}[{\rm erg/s}])$')
        #plt.ylabel(r'$\frac{\mathrm{d} E_{\rm cr}}{\mathrm{d} n_{\rm ISM}}\Delta \log n_{\rm ISM}[{\rm erg/cm^3}]$')
        plt.ylabel(r'$\mathrm{d} E_{\rm cr}\mathrm{d} \log (L_{\gamma})[{\rm erg}]$')
        plt.savefig('cregamma_'+fmeat+'.pdf')
        plt.clf()

if wanted=='nismgamma' or wanted=='nismcre' or wanted=='nismcrad' or wanted=='nismcrg' or wanted=='nismcrl' or wanted == 'rcrad':
    rcParams['figure.figsize'] = 5,4
    nogrid=15
    Rvcut = 0
    atstarburst=0
    trackneed=0
    newlabelneed=1
        for runtodo in dirneed:
                timel=[]
        #enllist=[]
        #precount=0
        info=outdirname(runtodo, Nsnap)
        havecr=info['havecr']
        if havecr==0:
            continue
                if atstarburst==1:
                        if runtodo=='bwsbclr':
                                Nsnap=523
                        if runtodo=='bwsbclrmhd':
                                Nsnap=589
                        if runtodo=='bwsbclrdc0':
                                Nsnap=590
                        if runtodo=='bwsbclrdc27':
                                Nsnap=138
                        if runtodo=='bwsbclrdc28':
                                Nsnap=520
                        if runtodo=='bwsbclrdc29':
                                Nsnap=433
                        if runtodo=='bwsbclrstr':
                                Nsnap=245
                        if runtodo=='bwsbclrdc28mhd':
                                Nsnap=558
                        if runtodo=='bwsbclrdc28str':
                                Nsnap=370
                for i in [Nsnap]:
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        cosmo=info['cosmo']
            color=info['color']
            haveB=info['haveB']
            Rvir = info['Rvir']
            newlabel=info['newlabel']
            if runtitle=='SMC':
                labelneed=dclabel
                if newlabelneed==1:
                    labelneed="\n".join(wrap(newlabel,17))
            else:
                labelneed=''
                        ptitle=title
                        if runtitle=='SMC':
                                ptitle='Dwarf'
                        elif runtitle=='SBC':
                                ptitle='Starburst'
                        elif runtitle=='MW':
                                ptitle=r'$L\star$ Galaxy'
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
            header = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix,\
             extension=the_suffix, havecr=havecr, header_only=1)
            print 'time', header['time']
                        timeneed=i*0.98*1e6 #in yr
                        Grho_t = G['rho'] #1e10Msun per kpc^3
            Gid_t = G['id']
                        cregy_codeunit = G['cregy']
                        cregy = cregy_codeunit*1e10*solar_mass_in_g*km_in_cm*km_in_cm #original cosmic ray energy in 1e10Msun km^2/sec^2
                        Neb = G['ne']
                        cregyl = G['cregyl']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
            if wanted == 'nismcrad' or wanted == 'rcrad':
                cregyad_t = G['cregya']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
            if wanted == 'nismcrg':
                                cregyg_t = G['cregyg']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
            Gp_t = G['p']
            Gx_t = Gp_t[:,0]
            Gy_t = Gp_t[:,1]
            Gz_t = Gp_t[:,2]
            Gr_t = np.sqrt(Gx_t*Gx_t+Gy_t*Gy_t+Gz_t*Gz_t)
            if Rvcut==1:
                Grtcut = Gr_t<1.0*Rvir
                Grho_t=Grho_t[Grtcut]
                cregy_codeunit=cregy_codeunit[Grtcut]
                cregy=cregy[Grtcut]
                Neb = Neb[Grtcut]
                cregyl = cregyl[Grtcut]
                Gid_t = Gid_t[Grtcut]
                Gr_t=Gr_t[Grtcut]
                if wanted == 'nismcrad' or wanted == 'rcrad':
                    cregyad_t = cregyad_t[Grcut]
                if wanted == 'nismcrg':
                                        cregyg_t = cregyg_t[Grcut]
            if wanted == 'nismcrad' or wanted == 'nismcrg' or wanted == 'nismcrl' or wanted == 'rcrad':
                snaptrack=i-snapsep
                info=outdirname(runtodo, snaptrack)
                Nsnapstring = info['Nsnapstring']
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                timetrack=snaptrack*0.98*1e6
                if trackneed==1:
                    Gp= G['p']
                    Gx= Gp[:,0]
                    Gy= Gp[:,1]
                    Gz= Gp[:,2]
                    Gr= np.sqrt(Gx*Gx+Gy*Gy+Gz*Gz)
                    Grho = G['rho']
                    Gid = G['id']
                    Gid, unindex = np.unique(Gid, return_index=True)
                    Gid_t, unindex_t = np.unique(Gid_t, return_index=True)
                    idint0=np.in1d(Gid,Gid_t)
                    # there are some duplicate IDs, from particle splitting? I drop those repetitions.
                    Gid = Gid[idint0]
                    Gidinds = Gid.argsort()
                    Gid_tinds = Gid_t.argsort()
                    Grho = Grho[unindex]
                    Grho = Grho[idint0]
                    Grho = Grho[Gidinds]
                    Grho_t = Grho_t[unindex_t]
                    Grho_t = Grho_t[Gid_tinds]
                                        Gr = Gr[unindex]
                                        Gr = Gr[idint0]
                                        Gr = Gr[Gidinds]
                                        Gr_t = Gr_t[unindex_t]
                                        Gr_t = Gr_t[Gid_tinds]
                if wanted == 'nismcrad' or wanted == 'rcrad':
                    cregyad = G['cregya']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                    if trackneed==1:
                        cregyad = cregyad[unindex]
                        cregyad = cregyad[idint0]
                        cregyad_t = cregyad_t[unindex_t]
                        dcrad = (cregyad_t[Gid_tinds]-cregyad[Gidinds])/(timeneed-timetrack)/yr_in_sec
                if wanted == 'nismcrg':
                    cregyg = G['cregyg']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                                        cregyg = cregyg[idint0]
                                        dcrg = (cregyg_t[Gid_tinds]-cregyg[Gidinds])/(timeneed-timetrack)/yr_in_sec
                                if wanted == 'nismcrl':
                                        cregyl_p = G['cregyl']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                                        cregyl_p = cregyl_p[idint0]
                                        dcrl = (cregyl[Gid_tinds]-cregyl_p[Gidinds])/(timeneed-timetrack)/yr_in_sec
                Gnism_in_cm_3 = 1.0/proton_mass_in_g*Grho_t*1e10*Msun_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm
        if wanted=='nismgamma':
            Lout = outLgamma_nism(Grho,Neb,cregy_codeunit)
            Lgamma = Lout['Lgamma']
                LogGnism = np.log10(Gnism_in_cm_3)
        LogGr = np.log10(Gr_t)
        if trackneed==1:
            LogGnxaxis = np.linspace(-6,2.5,num=nogrid)
        else:
            LogGnxaxis = np.linspace(-4,4.0,num=nogrid)
        if wanted == 'rcrad':
            LogGnxaxis = np.linspace(-1,2.0,num=nogrid)
        dx =LogGnxaxis[1]-LogGnxaxis[0]
                Lgammal = []
        crel=[]
        dcradl=[]
        dcrgl=[]
        dcrll=[]
        for inism in range(nogrid-1):
            if wanted == 'rcrad':
                cutg = (LogGr > LogGnxaxis[inism]) & (LogGr < LogGnxaxis[inism+1])
            else:
                cutg = (LogGnism > LogGnxaxis[inism]) & (LogGnism < LogGnxaxis[inism+1])
            if wanted=='nismgamma':
                Lgammacut = Lgamma[cutg]
                Lgammal = np.append(Lgammal, np.sum(Lgammacut))
            if wanted=='nismcre':
                crecut = cregy[cutg]
                crel = np.append(crel, np.sum(crecut))  
            if wanted=='nismcrad' or wanted == 'rcrad':
                dcradcut = dcrad[cutg]
                dcradl = np.append(dcradl, np.sum(dcradcut)) 
                        if wanted=='nismcrg':  
                                dcrgcut = dcrg[cutg] 
                                dcrgl = np.append(dcrgl, np.sum(dcrgcut))
            if wanted=='nismcrl':
                                dcrlcut = dcrl[cutg]
                                dcrll = np.append(dcrll, np.sum(dcrlcut))
        if haveB>0:
            ls='dashed'
        else:
            ls='solid'
        if wanted == 'rcrad':
            plt.plot((LogGnxaxis[1:]+LogGnxaxis[:-1])/2.,dcradl/dx, label=labelneed, color=color,lw=2,ls=ls)
        if wanted=='nismgamma':
            plt.plot((LogGnxaxis[1:]+LogGnxaxis[:-1])/2.,Lgammal/dx, label=labelneed, color=color,lw=2,ls=ls)
        if wanted=='nismcre':
            plt.plot((LogGnxaxis[1:]+LogGnxaxis[:-1])/2.,crel/dx, label=labelneed, color=color,lw=2,ls=ls)
                if wanted=='nismcrad':
                        plt.plot((LogGnxaxis[1:]+LogGnxaxis[:-1])/2.,dcradl/dx, label=labelneed, color=color,lw=2,ls=ls)
                if wanted=='nismcrg':
                        plt.plot((LogGnxaxis[1:]+LogGnxaxis[:-1])/2.,dcrgl/dx, label=labelneed, color=color,lw=2,ls=ls)
                if wanted=='nismcrl':
                        plt.plot((LogGnxaxis[1:]+LogGnxaxis[:-1])/2.,dcrll/dx, label=labelneed, color=color,lw=2,ls=ls)
    if runtitle=='SMC':
        plt.legend(loc='best',fontsize=9,ncol=1)
    if trackneed==0:
        plt.yscale('log')

        if wanted=='rcrad':
                plt.ylabel(r'$\mathrm{d} \dot{E}_{\rm Ad}/\mathrm{d} \log r\;[{\rm erg/s}]$',fontsize=16)
                filename = 'CRplot/rcrad/rcrad_'+fmeat+'.pdf'
    if wanted=='nismgamma':
        plt.ylabel(r'$\mathrm{d} L_{\gamma}/\mathrm{d} \log (n)[{\rm erg/s}]$',fontsize=16)
        filename = 'CRplot/nismgamma/nismgamma_'+fmeat+'.pdf'
    if wanted=='nismcre':
        plt.ylim(ymin=1e48)
        plt.ylabel(r'$\mathrm{d} E_{\rm cr}/\mathrm{d} \log (n)[{\rm erg}]$',fontsize=16)
        filename = 'CRplot/nismcre/nismcre_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
        if wanted=='nismcrad':
                plt.ylabel(r'$\mathrm{d} \dot{E}_{\rm Ad}/\mathrm{d} \log (n)[{\rm erg/s}]$',fontsize=16)
                filename = 'CRplot/nismcrad/nismcrad_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
        if wanted=='nismcrg':
                plt.ylabel(r'$\mathrm{d} \dot{E}_{\rm SNe}/\mathrm{d} \log (n)[{\rm erg/s}]$',fontsize=16)
                filename = 'CRplot/nismcrg/nismcrg_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
        if wanted=='nismcrl':
                plt.ylabel(r'$\mathrm{d} \dot{E}_{\rm Loss}/\mathrm{d} \log (n)[{\rm erg/s}]$',fontsize=16)
                filename = 'CRplot/nismcrl/nismcrl_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
    print 'filename', filename
    if wanted=='rcrad':
        plt.xlabel(r'$\log (r\;[{\rm kpc}])$',fontsize=16)
    else:
        plt.xlabel(r'$\log (n[{\rm cm^{-3}}])$',fontsize=16)
    plt.title(ptitle,fontsize=21)
        plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
        plt.savefig(filename,bbox_inches='tight')
        plt.clf()

if wanted=='nismcumgamma':
        for runtodo in dirneed:
        timel=[]
                for i in [Nsnap]:
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        cosmo=info['cosmo']
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
            if cosmo==1:
                                        G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
            else:
                                        G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
            timeneed=i*0.98*1e9 #in yr
            Grho = G['rho'] #1e10Msun per kpc^3
            Gm = G['m'] #1e10Msun
            cregy = G['cregy'] # cosmic ray energy in 1e10Msun km^2/sec^2
            cregy_in_erg = cregy*1e10*solar_mass_in_g*km_in_cm*km_in_cm
            Neb = G['ne']
            if havecr>1:
                cregyl = G['cregyl']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
            timel=np.append(timel,timeneed)
        Gnism_in_cm_3 = (0.78+0.22*Neb*0.76)/proton_mass_in_g*Grho*1e10*solar_mass_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm
        #tpi_in_yr = 2e5/Gnism_in_cm_3*250.0 #pi decay time in yr
        #Lgammagev = cregy_in_erg/tpi_in_yr/sec_in_yr*betapi/nopi_per_gamma #in erg/s
        #print 'count nism>1e3', np.count_nonzero(Gnism_in_cm_3[Gnism_in_cm_3>1e3])
        if havecr==1:
            Lout=outLgamma_nism(Grho,Neb,cregy)
            Gnism_in_cm_3=Lout['nism']
            Lgamma=Lout['Lgamma']
        else:
            Eloss=cregyl/7.51e-16/pidecay_fac*betapi/nopi_per_gamma
            Lgamma=np.absolute((Eloss)/float(Nsnap)/0.98e6/yr_in_sec)
        LogGnism = np.log10(Gnism_in_cm_3)
        LogGnxaxis = np.linspace(-4,4,num=nogrid)
        Lgammal = []
        #print 'Eloss', Eloss
        for inism in range(nogrid-1):
            cutg = (LogGnism > LogGnxaxis[inism])
            Lgammacut = Lgamma[cutg]
                        Lgammal = np.append(Lgammal, np.sum(Lgammacut))
        print 'Lgammal', Lgammal
        plt.plot((LogGnxaxis[1:]+LogGnxaxis[:-1])/2.,Lgammal, label=runtodo)
        plt.legend(loc='best')
        plt.yscale('log')
        plt.xlabel(r'$\log (n_{\rm ISM}[{\rm cm^{-3}}])$')
        plt.ylabel(r'$L_{\gamma}(>n_{\rm ISM})[{\rm erg/s}]$')
        plt.savefig('CRplot/nismcumgamma_'+fmeat+'.pdf')
        plt.clf()


if wanted=='gammadecaytime':
        for runtodo in dirneed:
                for i in [Nsnap]:
                        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr,timestep=outdirname(runtodo, i)
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Grho = G['rho']
                        cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                        Neb = G['ne']
                        Gnism_in_cm_3 = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21
                        tpi_in_yr = 2e5/Gnism_in_cm_3*250.0 #pi decay time in yr
                        cregy_in_erg = cregy*solar_mass_in_g*1e10*km_in_cm*km_in_cm
                        Lgammagev = cregy_in_erg/tpi_in_yr/yr_in_sec*betapi/nopi_per_gamma #in erg/s
                        logtpi = np.log10(tpi_in_yr)
                        logtpixaxis = np.linspace(5,14,num=nogrid)
                        dx = logtpixaxis[1]-logtpixaxis[0]
                        Lgammagevl = []
                        for inism in range(nogrid-1):
                                cutt = (logtpi > logtpixaxis[inism]) & (logtpi < logtpixaxis[inism+1])
                                Lgammagevcut = Lgammagev[cutt]
                                Lgammagevl = np.append(Lgammagevl, np.sum(Lgammagevcut))
                        plt.plot((logtpixaxis[1:]+logtpixaxis[:-1])/2.,Lgammagevl/dx, label=runtodo)
        plt.legend(loc='best')
        plt.yscale('log')
        plt.xlabel(r'$\log (t_{ \pi}[{\rm yr}])$')
        #plt.ylabel(r'$\frac{\mathrm{d} E_{\rm cr}}{\mathrm{d} n_{\rm ISM}}\Delta \log n_{\rm ISM}[{\rm erg/cm^3}]$')
        plt.ylabel(r'$\mathrm{d} L_{\gamma} /\mathrm{d} \log (t_{ \pi})[{\rm erg/s}]$')
        plt.savefig('gammadecaytime_'+fmeat+'.pdf')
        plt.clf()


if wanted=='gammacumdecaytime':
        for runtodo in dirneed:
                for i in [Nsnap]:
                        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr,timestep=outdirname(runtodo, i)
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Grho = G['rho']
                        cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                        Neb = G['ne']
                        Gnism_in_cm_3 = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21
                        tpi_in_yr = 2e5/Gnism_in_cm_3*250.0 #pi decay time in yr
                        cregy_in_erg = cregy*solar_mass_in_g*1e10*km_in_cm*km_in_cm
                        Lgammagev = cregy_in_erg/tpi_in_yr/yr_in_sec*betapi/nopi_per_gamma #in erg/s
                        logtpi = np.log10(tpi_in_yr)
                        logtpixaxis = np.linspace(5,14,num=nogrid)
                        Lgammagevl = []
                        for inism in range(nogrid-1):
                                cutt = (logtpi > logtpixaxis[inism])
                                Lgammagevcut = Lgammagev[cutt]
                                Lgammagevl = np.append(Lgammagevl, np.sum(Lgammagevcut))
                        plt.plot((logtpixaxis[1:]+logtpixaxis[:-1])/2.,Lgammagevl, label=runtodo)
        plt.legend(loc='best')
        plt.yscale('log')
        plt.xlabel(r'$\log (t_{ \pi}[{\rm yr}])$')
        #plt.ylabel(r'$\frac{\mathrm{d} E_{\rm cr}}{\mathrm{d} n_{\rm ISM}}\Delta \log n_{\rm ISM}[{\rm erg/cm^3}]$')
        plt.ylabel(r'$L_{\gamma}(>t_{ \pi})[{\rm erg/s}]$')
        plt.savefig('gammacumdecaytime_'+fmeat+'.pdf')
        plt.clf()



if wanted=='gmnism':
        for runtodo in dirneed:
                for i in [Nsnap]:
                        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr,timestep=outdirname(runtodo, i)
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Grho = G['rho']
            Gm = G['m']
                        #cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                        Neb = G['ne']
                        Gnism_in_cm_3 = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21
                        #tpi_in_yr = 2e5/Gnism_in_cm_3*250.0 #pi decay time in yr
                        #betapi = 0.7 ##should be between 0.5-0.9 from Lacki 2011 #fraction of pi has enough energy (>1GeV)
                        #nopi_per_gamma = 3.0
            Gm_in_sun = Gm*1e10 
                        #cm_in_km = 1e5
                        #cregy_in_erg = cregy*solar_mass_in_g*1e10*cm_in_km*cm_in_km
                        #print 'np.sum(Lgammagev)', np.sum(Lgammagev)

                        LogGnism = np.log10(Gnism_in_cm_3)
                        LogGnxaxis = np.linspace(-4,4,num=nogrid)
                        Gml = []
                        dx = LogGnxaxis[1]-LogGnxaxis[0]
                        for inism in range(nogrid-1):
                                cutg = (LogGnism > LogGnxaxis[inism]) & (LogGnism < LogGnxaxis[inism+1])
                                Gmcut = Gm_in_sun[cutg]
                                Gml = np.append(Gml, np.sum(Gmcut))
                        plt.plot((LogGnxaxis[1:]+LogGnxaxis[:-1])/2.,Gml/dx, label=runtodo)
        plt.legend(loc='best')
        plt.yscale('log')
        plt.xlabel(r'$\log (n_{\rm ISM}[{\rm cm^{-3}}])$')
        plt.ylabel(r'$\mathrm{d} M_{\rm gas}/\mathrm{d} \log (n_{\rm ISM})[M_\odot]$')
        plt.savefig('gmnism_'+fmeat+'.pdf')
        plt.clf()


if wanted=='nismcumgm':
    withindisk=1
        for runtodo in dirneed:
                for i in [Nsnap]:
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        cosmo=info['cosmo']
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
            if cosmo==1:
                G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1, cosmological=1)
            else:
                G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
            if withindisk==1:
                if cosmo==1:
                    header=readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1,header_only=1)
                    h0 = header['hubble']
                    atime = header['time']
                    print 'halostr', halostr
                    halosA = SF.read_halo_history_pep(rundir, i, singlesnap=1, firever=firever,halonostr=halostr, hubble=h0, comoving=0, maindir=maindir)
                    redlist = halosA['redshift']
                    haloid = halosA['ID']
                    a_scale = 1.0/(1.0+redlist)
                    xcen = halosA['x']
                    ycen = halosA['y']
                    zcen = halosA['z']
                    Rvir = halosA['R']
                    Lxstar = halosA['Lxstar']
                    Lystar = halosA['Lystar']
                    Lzstar = halosA['Lzstar']
                else:
                    xcen=ycen=zcen=0.
                Gp = G['p']
                Gz = Gp[:,2]-zcen
                Gx = Gp[:,0]-xcen
                Gy = Gp[:,1]-ycen
                if cosmo==1 and rotface==1:
                    Gx, Gy, Gz = SF.rotateL_to_z(Gx,Gy,Gz,Lxstar,Lystar,Lzstar)
                Grxy = np.sqrt(Gx*Gx+Gy*Gy)
                cutxy = Grxy<10.
                cutz = np.absolute(Gz)<0.5
                cutxyz = cutxy*cutz
                        Grho = G['rho']
                        Gm = G['m']
            if withindisk==1:
                Grho=Grho[cutxyz]
                Gm=Gm[cutxyz]
                        Gnism_in_cm_3 = 1.0/protonmass_in_g*Grho*1e10*Msun_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm
                        Gm_in_sun = Gm*1e10
            Neb = np.zeros(len(Grho))
            cregy = np.zeros(len(Grho))
            Lout=outLgamma_nism(Grho,Neb,cregy)
                        Gnism_in_cm_3=Lout['nism']
            print 'count nism>1e3', np.count_nonzero(Gnism_in_cm_3[Gnism_in_cm_3>1e3])
                        LogGnism = np.log10(Gnism_in_cm_3)
                        LogGnxaxis = np.linspace(-4,4,num=nogrid)
                        Gml = []
                        dx = LogGnxaxis[1]-LogGnxaxis[0]
            print 'np.amin(Gm_in_sun)', np.amin(Gm_in_sun)
                        for inism in range(nogrid-1):
                                cutg = (LogGnism > LogGnxaxis[inism]) 
                                Gmcut = Gm_in_sun[cutg]
                                Gml = np.append(Gml, np.sum(Gmcut))
                        plt.plot((LogGnxaxis[1:]+LogGnxaxis[:-1])/2.,Gml, label=runtodo)
        plt.legend(loc='best')
        plt.yscale('log')
        plt.xlabel(r'$\log (n_{\rm ISM}[{\rm cm^{-3}}])$')
        plt.ylabel(r'$M_{\rm gas}(>n_{\rm ISM})[M_\odot]$')
    if withindisk==1:
                figname = 'CRplot/nismcumgm/nismcumgm_'+fmeat+'_sn'+str(Nsnap)+'_withindisk.pdf'
    else:
        figname = 'CRplot/nismcumgm/nismcumgm_'+fmeat+'_sn'+str(Nsnap)+'.pdf'
    print 'figname', figname
        plt.savefig(figname)
        plt.clf()


if wanted=='vzwindm' or wanted=='Twindm':
    rcParams['figure.figsize'] = 5,4
    vcut=0.
    zup=30.
    zdown=20.
    withinr=20.
        for runtodo in dirneed:
                for i in [Nsnap]: 
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep'] 
                        cosmo=info['cosmo']
                        haveB=info['haveB']
            color=info['color']
            ptitle=title
            if runtitle=='SMC':
                ptitle='Dwarf'
            elif runtitle=='SBC':
                ptitle='Starburst'
            elif runtitle=='MW':
                ptitle=r'$L\star$ Galaxy'
            lsn='solid'
            if haveB>0:
                lsn='dashed'
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        if cosmo==1:
                                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1, cosmological=1)
                        else:
                                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gm = G['m']
            Gp = G['p']
            Gv = G['v']
            Tb = G['u']
            rho = G['rho']
            Neb = G['ne']
                        Gm_in_sun = Gm*1e10
            Gx = Gp[:,0]
            Gy = Gp[:,1]
            Gz = Gp[:,2]
            Gvz = Gv[:,2]
            cutv = Gvz*Gz/np.absolute(Gz)>vcut
            cutz = (np.absolute(Gz)>zdown) & (np.absolute(Gz)<zup)
            cutr = np.sqrt(Gx*Gx+Gy*Gy)<withinr
            cut = cutv*cutz*cutr
            if wanted=='Twindm':
                TrueTemp, converted_rho  = SF.convertTemp(Tb, Neb, rho)
                Logvar = np.log10(TrueTemp[cut])
                                xmin=3;xmax=7.0;
                if runtitle=='SMC':
                    xmax=6.0;
            elif wanted=='vzwindm':
                Logvar = np.log10(Gvz[cut])
                xmin=0.5;xmax=3.0;
                if runtitle=='SMC':
                    xmax=2.3;
            Gm_in_sun=Gm_in_sun[cut]
                        Logxaxis = np.linspace(xmin,xmax,num=nogrid)
                        Gml = []
                        for ib in range(nogrid-1):
                                cutg = (Logvar > Logxaxis[ib]) & (Logvar < Logxaxis[ib+1])
                                Gmcut = Gm_in_sun[cutg]
                                Gml = np.append(Gml, np.sum(Gmcut))
                        plt.plot((Logxaxis[1:]+Logxaxis[:-1])/2.,Gml, label=dclabel,color=color,ls=lsn,lw=2)
    if runtitle=='SMC':
        plt.legend(loc='best',fontsize=9,ncol=3)
        plt.yscale('log')
    if wanted=='Twindm':
                plt.xlabel(r'$\log (T[{\rm K}])$',fontsize=18)
                plt.ylabel(r'${\rm d}M_{\rm wind}/{\rm d}\log T\;[M_\odot]$',fontsize=18)
        figname='CRplot/Twindm/Twindm_'+runtodo+'_sn'+str(Nsnap)+'.pdf'
    elif wanted=='vzwindm':
        plt.xlabel(r'$\log (v_z[{\rm km/s}])$',fontsize=18)
        plt.ylabel(r'${\rm d}M_{\rm wind}/{\rm d}\log v_z\;[M_\odot]$',fontsize=18)
        figname='CRplot/vzwindm/vzwindm_'+runtodo+'_sn'+str(Nsnap)+'.pdf'
    plt.title(ptitle,fontsize=18)
    plt.tight_layout()
    print 'figname', figname
        plt.savefig(figname)
        plt.clf()


if wanted=='Bcumgm':
        for runtodo in dirneed:
                for i in [Nsnap]:
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        cosmo=info['cosmo']
            haveB=info['haveB']
            if haveB<1:
                'no B field'
                exit()
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        if cosmo==1:
                                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1, cosmological=1)
                        else:
                                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        GB = G['B']
                        Gm = G['m']
                        Gm_in_sun = Gm*1e10
                        LogGB = np.log10(np.sqrt(GB[:,0]*GB[:,0]+GB[:,1]*GB[:,1]+GB[:,2]*GB[:,2]))
                        LogGnxaxis = np.linspace(-10,0,num=nogrid)
                        Gml = []
                        for ib in range(nogrid-1):
                                cutg = (LogGB > LogGnxaxis[ib])
                                Gmcut = Gm_in_sun[cutg]
                                Gml = np.append(Gml, np.sum(Gmcut))
                        plt.plot((LogGnxaxis[1:]+LogGnxaxis[:-1])/2.,Gml, label=dclabel)
        plt.legend(loc='best')
        plt.yscale('log')
        plt.xlabel(r'$\log (B[Gauss])$')
        plt.ylabel(r'$M_{\rm gas}(>B)[M_\odot]$')
        plt.savefig('CRplot/Bcumgm/Bcumgm_'+fmeat+'_sn'+str(Nsnap)+'.pdf')
        plt.clf()


if wanted=='Begytime' or wanted=='Begy_mtime' or wanted=='Brmstime':
        for runtodo in dirneed:
        snapl=[]
        Begyl=[]
        Gml=[]
        Gvoll=[]
        Begydencutl=[]
        Gmdencutl=[]
        Gvoldencutl=[]
        cregylist=[]
        info=outdirname(runtodo, Nsnap)
        haveB=info['haveB']
        withinr=10 #kpc
        maxlength = 1#kpc
        if haveB<1:
            'no B field'
            continue
                for i in range(startno,Nsnap,snapsep):
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
            maindir=info['maindir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        cosmo=info['cosmo']
            color=info['color']
                        haveB=info['haveB']
                        if haveB<1:
                                'no B field'
                exit()
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
            try:
                if cosmo==1:
                    G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1, cosmological=1)
                    h0=1
                else:
                    G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                    h0=0
                Gpos = G['p']
            except KeyError:
                print 'KeyError'
                break
            header = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo, header_only=1)
            ascale = header['time']
            #print 'this time', ascale
            thisred = 1./ascale-1.
            hubble = header['hubble']
            if cosmo==1:
                halosingle = read_halo_history(rundir, halonostr='00',hubble=hubble, comoving=0, maindir=maindir, singlesnap=1, atime=ascale)
                xcen = halosingle['x']
                ycen = halosingle['y']
                zcen = halosingle['z']
                xvcen = halosingle['xv']
                yvcen = halosingle['yv']
                zvcen = halosingle['zv']
                Rvirnow = halosingle['R']
                MgAHF = halosingle['Mg']
            #       print 'xcenl[0]', xcenl[0]
            #       print 'thisred', thisred
                #print 'cen', xcen, ycen, zcen
                #print 'MgAHF', MgAHF
                #print 'Rvir', Rvirnow
            else:
                xcen=0
                ycen=0
                zcen=0
                xvcen=0
                yvcen=0
                zvcen=0
            Gx = Gpos[:,0]-xcen
            Gy = Gpos[:,1]-ycen
            Gz = Gpos[:,2]-zcen
            Grxy = np.sqrt(Gx*Gx+Gy*Gy)
                        GB = G['B'] # in Gauss
                        Gm = G['m'] #in 1e10 Msun
            Grho = G['rho'] #in 1e10 Msun/kpc^3
            Neb = G['ne']
            Bx = GB[:,0]
            By = GB[:,1]
            Bz = GB[:,2]
            B2 = Bx*Bx+By*By+Bz*Bz
            if havecr>0:
                cregy = G['cregy']*1e10*Msun_in_g*km_in_cm*km_in_cm
            print 'np.amax(np.sqrt(B2))', np.amax(np.sqrt(B2))
                        print 'np.amin(np.sqrt(B2))', np.amin(np.sqrt(B2))
            cutr = Grxy<withinr
            cutz = np.absolute(Gz)<maxlength
            cut = cutr*cutz
            Grho = Grho[cut]
            Gm = Gm[cut]
            B2 = B2[cut]
            Neb = Neb[cut]
            if havecr>0:
                cregyt=np.sum(cregy[cut])
                        Gnism = (0.78+0.22*Neb*0.76)/protonmass_in_g*Grho*1e10*Msun_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm
            Gvol = Gm/Grho*kpc_in_cm*kpc_in_cm*kpc_in_cm
            Begy = B2/8./np.pi*Gvol
            Gtotm_in_g = np.sum(Gm)*1e10*Msun_in_g
            dencut = Gnism>1.0 #cm^-3
            Gmdencut = np.sum(Gm[dencut])*1e10*Msun_in_g
            Begydencut = np.sum(Begy[dencut])
            sumBegy = np.sum(Begy)  
            sumGvol = np.sum(Gvol)
            sumGvoldencut = np.sum(Gvol[dencut])
            Begyl = np.append(sumBegy,Begyl)
            Gml = np.append(Gtotm_in_g,Gml)
            Begydencutl = np.append(Begydencut,Begydencutl)
            Gmdencutl = np.append(Gmdencut,Gmdencutl)
            print 'Gmdencutl',Gmdencutl
            Gvoldencutl = np.append(sumGvoldencut,Gvoldencutl)
            Gvoll = np.append(sumGvol,Gvoll)
            snapl = np.append(i,snapl)
            if havecr>0:
                cregylist = np.append(cregyt,cregylist)
        Begyl=np.array(Begyl)
        Gvoll=np.array(Gvoll)
        Gml  =np.array(Gml)
        Begydencutl = np.array(Begydencutl)
        Gvoldencutl = np.array(Gvoldencutl)
        if wanted=='Begytime':
            totV= 2.*maxlength*np.pi*withinr*withinr*kpc_in_cm*kpc_in_cm*kpc_in_cm #in cm3
            plt.plot(snapl,Begyl*erg_in_eV/totV,color=color,lw=2,ls='dashed', label=dclabel)
                        if havecr>0:
                                plt.plot(snapl,cregylist*erg_in_eV/totV,lw=2,ls='dashdot',color=color)
        if wanted=='Begy_mtime':
            B_ml = Begyl/Gml
            plt.plot(snapl,B_ml,color=color,lw=2,ls='dashed', label=dclabel)    
            plt.plot(snapl,Begydencutl/Gmdencutl,color=color,lw=1,ls='solid')
            if havecr>0:
                plt.plot(snapl,cregylist/Gml,lw=2,ls='dashdot',color=color)
        if wanted=='Brmstime':
            Brms = np.sqrt(Begyl/Gvoll)
            Brmsdencut = np.sqrt(Begydencutl/Gvoldencutl)
            plt.plot(snapl,Brms*1.e6,color=color,lw=2,ls='dashed', label=dclabel)
            plt.plot(snapl,1.e6*Brmsdencut,color=color,lw=1,ls='solid')
    if wanted=='Begytime':
        plt.legend(loc='best')
        plt.yscale('log')
        plt.ylabel('Magnetic energy density (eV/cm^3)')
        plt.xlabel('Myr')
        plt.savefig('CRplot/Begytime/Begytime_'+fmeat+'.pdf')
        plt.clf()
        if wanted=='Begy_mtime':
                plt.legend(loc='best')
                plt.yscale('log')
                plt.ylabel('Magnetic energy per mass (erg/g)')
                plt.xlabel('Myr')
                plt.savefig('CRplot/Begy_mtime/Begy_mtime_'+fmeat+'.pdf')
                plt.clf()
    if wanted=='Brmstime':
                plt.legend(loc='best')
                plt.yscale('log')
                plt.ylabel(r'RMS Magnetic field ($\mu$G)')
                plt.xlabel('Myr')
                plt.savefig('CRplot/Brmstime/Brmstime_'+fmeat+'.pdf')
                plt.clf()





if wanted=='pvolcumN':
        for runtodo in dirneed:
                for i in [Nsnap]:
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        cosmo=info['cosmo']
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        if cosmo==1:
                                G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1, cosmological=1)
                        else:
                                G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Grho = G['rho']
                        Gm = G['m']
            Gvol_in_kpc3 = Gm/Grho
                        LogGvol = np.log10(1/Gvol_in_kpc3)
                        LogGvolxaxis = np.linspace(-4,4,num=nogrid)
                        GNl = []
                        dx = LogGvolxaxis[1]-LogGvolxaxis[0]
                        for inism in range(nogrid-1):
                                cutg = (LogGvol > LogGvolxaxis[inism])
                                GNcut = Gm[cutg]
                                GNl = np.append(GNl, len(GNcut))
                        plt.plot((LogGvolxaxis[1:]+LogGvolxaxis[:-1])/2.,GNl, label=runtodo)
        plt.legend(loc='best')
        plt.yscale('log')
        plt.xlabel(r'$\log (1/V_{\rm gas}[{\rm kpc^{-3}}])$')
        plt.ylabel(r'$N(>1/V_{\rm gas})$')
        plt.savefig('CRplot/pvolcumN_'+fmeat+'.pdf')
        plt.clf()



if wanted=='pdxcumN':
        for runtodo in dirneed:
                for i in [Nsnap]:
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        cosmo=info['cosmo']
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        if cosmo==1:
                                G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1, cosmological=1)
                        else:
                                G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gh = G['h']
                        Gdx_in_pc = Gh*1e3*1.61199/32. #(4/3*pi)^(1/2) and 32 is NumNgb
                        LogGdx = np.log10(Gdx_in_pc)
                        LogGxaxis = np.linspace(-5,5,num=nogrid)
                        GNl = []
                        dx = LogGxaxis[1]-LogGxaxis[0]
                        for inism in range(nogrid-1):
                                cutg = (LogGdx < LogGxaxis[inism])
                                GNcut = Gh[cutg]
                                GNl = np.append(GNl, len(GNcut))
                        plt.plot((LogGxaxis[1:]+LogGxaxis[:-1])/2.,GNl, label=runtodo)
        plt.legend(loc='best')
        plt.yscale('log')
    plt.xlim([5,-5])
        plt.xlabel(r'$\log ({\rm d} x_{\rm gas}[{\rm pc}])$')
        plt.ylabel(r'$N(>{\rm d} x_{\rm gas})$')
        plt.savefig('CRplot/pdxcumN_'+fmeat+'.pdf')
        plt.clf()


if wanted=='pdxcumcre':
        for runtodo in dirneed:
                for i in [Nsnap]:
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        cosmo=info['cosmo']
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        if cosmo==1:
                                G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1, cosmological=1)
                        else:
                                G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gh = G['h']
            cregy_in_erg = G['cregy']*1e10*Msun_in_g*km_in_cm*km_in_cm
                        Gdx_in_pc = Gh*1e3*1.61199/32. #(4/3*pi)^(1/2) and 32 is NumNgb
                        LogGdx = np.log10(Gdx_in_pc)
                        LogGxaxis = np.linspace(-5,5,num=nogrid)
                        Gcrel = []
                        dx = LogGxaxis[1]-LogGxaxis[0]
                        for inism in range(nogrid-1):
                                cutg = (LogGdx < LogGxaxis[inism])
                                Gcrecut = cregy_in_erg[cutg]
                                Gcrel = np.append(Gcrel, np.sum(Gcrecut))
                        plt.plot((LogGxaxis[1:]+LogGxaxis[:-1])/2.,Gcrel, label=runtodo)
        plt.legend(loc='best')
        plt.yscale('log')
        plt.xlim([5,-5])
        plt.xlabel(r'$\log ({\rm d} x_{\rm gas}[{\rm pc}])$')
        plt.ylabel(r'$E_{\rm CR}(>{\rm d} x_{\rm gas})(erg)$')
        plt.savefig('CRplot/pdxcumcre_'+fmeat+'.pdf')
        plt.clf()



if wanted=='pvolcumcre':
        for runtodo in dirneed:
                for i in [Nsnap]:
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        cosmo=info['cosmo']
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        if cosmo==1:
                                G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1, cosmological=1)
                        else:
                                G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Grho = G['rho']
                        Gm = G['m']
            cregy = G['cregy']*1e10*Msun_in_g*km_in_cm*km_in_cm
                        Gvol_in_kpc3 = Gm/Grho
                        LogGvol = np.log10(1/Gvol_in_kpc3)
                        LogGvolxaxis = np.linspace(-4,4,num=nogrid)
                        Gcrel = []
                        dx = LogGvolxaxis[1]-LogGvolxaxis[0]
                        for inism in range(nogrid-1):
                                cutg = (LogGvol > LogGvolxaxis[inism])
                                Gcrecut = cregy[cutg]
                                Gcrel = np.append(Gcrel, np.sum(Gcrecut))
                        plt.plot((LogGvolxaxis[1:]+LogGvolxaxis[:-1])/2.,Gcrel, label=runtodo)
        plt.legend(loc='best')
        plt.yscale('log')
        plt.xlabel(r'$\log (1/V_{\rm gas}[{\rm kpc^{-3}}])$')
        plt.ylabel(r'$E_{\rm CR}(>1/V_{\rm gas})({\rm erg})$')
        plt.savefig('CRplot/pvolcumcre_'+fmeat+'.pdf')
        plt.clf()




if wanted=='gasdensph':
    withinr=100.0
        for runtodo in dirneed:
        Gnism_in_cm_3l=[]
                for i in [Nsnap]:
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        cosmo=info['cosmo']
                        usepep=info['usepep']
                        finalno=info['finalno']
                        beginno=info['beginno']
                        firever=info['firever']
                        halostr=info['halostr']
                        maindir=info['maindir']

                        print 'Nsnapstring', Nsnapstring
                        if cosmo==1:
                                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                        else:
                                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gp = G['p']
            Gm = G['m'] 
            print 'Gm', Gm
                        #Grho = G['rho']
                        Gz = Gp[:,2]
                        Gx = Gp[:,0]
                        Gy = Gp[:,1]
                        dr = withinr/nogrid
                        radl =[]
                        for irad in range(nogrid):
                                if cosmo ==1:
                                        Gp = G['p']
                                        Gx = Gp[:,0]
                                        Gy = Gp[:,1]
                                        Gz = Gp[:,2]
                                        header=readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1,header_only=1)
                                        h0 = header['hubble']
                                        atime = header['time']
                                        if usepep==1:
                                                halosA = read_halo_history_pep(rundir, finalno, beginno=beginno,\
                         singlesnap=0, firever=firever,halonostr=halostr, hubble=h0, comoving=1, maindir=maindir)
                                                afactor=atime
                                        else:
                                                halosA = read_halo_history(rundir, maindir=maindir, halonostr=halostr, hubble=h0, comoving=0)
                                                afactor=1.0
                                        redlist = halosA['redshift']
                                        haloid = halosA['ID']
                                        a_scale = 1.0/(1.0+redlist)
                                        xcenl = halosA['x']*afactor
                                        ycenl = halosA['y']*afactor
                                        zcenl = halosA['z']*afactor
                                        Rvirl = halosA['R']*afactor
                                        xcen = np.interp(atime,a_scale,xcenl)
                                        ycen = np.interp(atime,a_scale,ycenl)
                                        zcen = np.interp(atime,a_scale,zcenl)
                                        Rvir = np.interp(atime,a_scale,Rvirl)
                    print 'center', xcen, ycen, zcen
                    print 'Rvir', Rvir
                                        Gxrel = Gx-xcen
                                        Gyrel = Gy-ycen
                                        Gzrel = Gz-zcen
                                        Gr = np.sqrt(Gxrel*Gxrel+Gyrel*Gyrel+Gzrel*Gzrel)
                                else:
                                        Gr = np.sqrt(Gx*Gx+Gy*Gy+Gz*Gz)
                                cut = (Gr > dr*irad) & (Gr < dr*(irad+1))
                                Gmcut = Gm[cut]
                                Gm_in_g = Gmcut*1e10*Msun_in_g
                                shellvol_in_cm3 = 4.0/3.0*np.pi*(-np.power(dr*irad,3)+np.power(dr*(irad+1),3))*kpc_in_cm*kpc_in_cm*kpc_in_cm
                                Grho_in_g_cm_3 = Gm_in_g/shellvol_in_cm3
                                Gnism_in_cm_3 = np.sum(Grho_in_g_cm_3/proton_mass_in_g)
                                Gnism_in_cm_3l = np.append(Gnism_in_cm_3l, Gnism_in_cm_3)
                                radl = np.append(radl, dr*(irad+0.5))
            print 'radl', radl
            print 'Gnism_in_cm_3l', Gnism_in_cm_3l
                        plt.plot(radl, Gnism_in_cm_3l, label=dclabel)
        plt.xlabel('r [kpc]')
        plt.ylabel(r'$n_{\rm ISM} [{\rm cm^{-3}}]$')
        plt.legend(loc='best')
        plt.yscale('log')
        plt.title('spherical shell centered at the center of the disk')
        plt.savefig('CRplot/gasdensph/gasdensity_'+fmeat+'_sph.pdf')
        plt.clf()


if wanted=='credensph':
        for runtodo in dirneed:
                for i in [Nsnap]:
                        info=outdirname(runtodo, Nsnap=Nsnap)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        cosmo=info['cosmo']
            usepep=info['usepep']
            finalno=info['finalno']
            beginno=info['beginno']
            firever=info['firever']
            halostr=info['halostr']
            maindir=info['maindir'] 

            print 'Nsnapstring', Nsnapstring
                        if cosmo==1:
                                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                        else:
                                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gp = G['p']
            Gm = G['m']
                        #Grho = G['rho']
                        cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                        #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                        Gz = Gp[:,2]
                        Gx = Gp[:,0]
                        Gy = Gp[:,1]
                        dr = withinr/nogrid
                        cregydenl=[]
                        radl =[]
            print 'cregy', np.sum(cregy)
            print 'Gm', Gm
                        for irad in range(nogrid):
                if cosmo ==1:
                    Gp = G['p']
                    Gx = Gp[:,0]
                    Gy = Gp[:,1]
                    Gz = Gp[:,2]
                    header=readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1,header_only=1)
                    h0 = header['hubble']
                    atime = header['time']
                    if usepep==1:
                        halosA = read_halo_history_pep(rundir, finalno, beginno=beginno,\
                         singlesnap=1, firever=firever,halonostr=halostr, hubble=h0, comoving=0, maindir=maindir)
                                                afactor=atime
                    else:
                        halosA = read_halo_history(rundir, maindir=maindir, halonostr=halostr, hubble=h0, comoving=0)
                                                afactor=1.0
                    redlist = halosA['redshift']
                    haloid = halosA['ID']
                    a_scale = 1.0/(1.0+redlist)
                    xcenl = halosA['x']*afactor
                    ycenl = halosA['y']*afactor
                    zcenl = halosA['z']*afactor
                    Rvirl = halosA['R']
                    xcen = np.interp(atime,a_scale,xcenl)
                    ycen = np.interp(atime,a_scale,ycenl)
                    zcen = np.interp(atime,a_scale,zcenl)
                    Rvir = np.interp(atime,a_scale,Rvirl)
                    Gxrel = Gx-xcen
                    Gyrel = Gy-ycen
                    Gzrel = Gz-zcen
                    Gr = np.sqrt(Gxrel*Gxrel+Gyrel*Gyrel+Gzrel*Gzrel)
                else:
                    Gr = np.sqrt(Gx*Gx+Gy*Gy+Gz*Gz)
                                cut = (Gr > dr*irad) & (Gr < dr*(irad+1))
                                cregycut = cregy[cut]
                print 'cregycut', cregycut
                                shellvol_in_kpc3 = 4.0/3.0*np.pi*(-np.power(dr*irad,3)+np.power(dr*(irad+1),3))
                                cregyden = cregycut/shellvol_in_kpc3
                cregyden *= 1e10*Msun_in_g*km_in_cm*km_in_cm/kpc_in_cm/kpc_in_cm/kpc_in_cm*erg_in_eV
                                cregydenl = np.append(cregydenl, np.sum(cregyden))
                #print 'unit conversion constant', 1e10*Msun_in_g*km_in_cm*km_in_cm/kpc_in_cm/kpc_in_cm/kpc_in_cm*erg_in_eV
                                radl = np.append(radl, dr*(irad+0.5))
                        plt.plot(radl, cregydenl, label=runtodo)
        plt.xlabel('r [kpc]')
#   plt.ylim([0.1,100])
#        plt.ylabel(r'$e_{\rm CR} [{\rm 10^{10}M_{\odot}km^2/s^2/kpc^3}]$')
        plt.ylabel(r'$e_{\rm CR} [{\rm eV/cm^3}]$')
        plt.legend(loc='best')
        plt.yscale('log')
        plt.title('spherical shell centered at the center of the disk')
        plt.savefig('CRplot/credensity_'+fmeat+'_sph.pdf')
        plt.clf()

if wanted=='gammadensph':
        for runtodo in dirneed:
                for i in [Nsnap]:
                        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr,timestep=outdirname(runtodo, i)
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gp = G['p']
                        Grho = G['rho']
                        Gu = G['u']
                        Gm = G['m']
                        cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                        Neb = G['ne']
                        #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                        Gz = Gp[:,2]
                        Gx = Gp[:,0]
                        Gy = Gp[:,1]
                        Gnism_in_cm_3 = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21
                        tpi_in_yr = 2e5/Gnism_in_cm_3*250.0 #pi decay time in yr
                        cregy_in_erg = cregy*solar_mass_in_g*1e10*km_in_cm*km_in_cm
                        Lgammagev = cregy_in_erg/tpi_in_yr/yr_in_sec*betapi/nopi_per_gamma #in erg/s
                        dr = withinr/nogrid
                        Lgdenl=[]
                        radl =[]
                        for irad in range(nogrid):
                                cut = (Gx*Gx+Gy*Gy+Gz*Gz > dr*irad*dr*irad) & (Gx*Gx+Gy*Gy+Gz*Gz < dr*(irad+1)*dr*(irad+1))
                                Nebcut = Neb[cut]
                                Lgammacut = Lgammagev[cut]
                                shellvol_in_kpc3 = 4.0/3.0*np.pi*(-np.power(dr*irad,3)+np.power(dr*(irad+1),3))
                                Lgden = Lgammacut/shellvol_in_kpc3
                Lgdenl = np.append(Lgdenl, np.sum(Lgden))
                                radl = np.append(radl, dr*(irad+0.5))
                        plt.plot(radl, Lgdenl, label=runtodo)
        plt.xlabel('r [kpc]')
        plt.ylabel(r'$l_{\gamma} [{\rm erg/s/kpc^{3}}]$')
    plt.yscale('log')
        plt.legend(loc='best')
        plt.title('spherical shell centered at the center of the disk')
        plt.savefig('gammadensity_'+fmeat+'_sph.pdf')
        plt.clf()

if wanted=='Gmencsph':
        for runtodo in dirneed:
                for i in [Nsnap]:
                        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr, timestep=outdirname(runtodo, i)
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gp = G['p']
                        #Grho = G['rho']
                        Gu = G['u']
                        Gm = G['m']
                        cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                        Neb = G['ne']
                        #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                        Gz = Gp[:,2]
                        Gx = Gp[:,0]
                        Gy = Gp[:,1]
                        dr = withinr/nogrid
                        Gmenc_in_sunl=[]
                        radl =[]
                        for irad in range(nogrid):
                                cut = (Gx*Gx+Gy*Gy+Gz*Gz < dr*(irad+1)*dr*(irad+1))
                                Nebcut = Neb[cut]
                                Gmcut = Gm[cut]
                                Gmenc_in_sun = Gmcut*1e10
                                Gmenc_in_sunl = np.append(Gmenc_in_sunl, np.sum(Gmenc_in_sun))
                                radl = np.append(radl, dr*(irad+0.5))
                        plt.plot(radl, Gmenc_in_sunl, label='kappa='+dclabel+'; '+'time step= '+timestep)
        plt.xlabel('r [kpc]')
        plt.ylabel(r'$M_{\rm gas} [ M_{\odot}]$')
        plt.legend(loc='best')
        plt.yscale('log')
        plt.title('Acummulative gas mass within a sphere centered at the center of the disk')
        plt.savefig('CRplot/Gmenc_'+fmeat+'_sph.pdf')
        plt.clf()


if wanted=='cresph':
        for runtodo in dirneed:
                for i in [Nsnap]:
                        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr,timestep=outdirname(runtodo, i)
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gp = G['p']
                        Grho = G['rho']
                        Gu = G['u']
                        Gm = G['m']
                        cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                        Neb = G['ne']
                        #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                        Gz = Gp[:,2]
                        Gx = Gp[:,0]
                        Gy = Gp[:,1]
                        dr = withinr/nogrid
                        crel=[]
                        radl =[]
                        for irad in range(nogrid):
                                cut = (Gx*Gx+Gy*Gy+Gz*Gz < dr*(irad)*dr*(irad))
                                crecut = cregy[cut]
                                crel = np.append(crel, np.sum(crecut))
                                radl = np.append(radl, dr*(irad))
                        plt.plot(radl, crel, label=runtodo)
        plt.xlabel('r [kpc]')
        plt.ylabel(r'$E_{rm cr} [{\rm 10^{10}M_{\odot}km^2/s^2}]$')
        plt.yscale('log')
        plt.legend(loc='best')
        plt.title('Cosmic ray energy within a sphere centered at the center of the disk')
        plt.savefig('cre_'+fmeat+'_sph.pdf')
        plt.clf()



if wanted=='gammasph':
    rcParams['figure.figsize'] = 5,4
    withinr = 15. #kpc
    nogrid=50
    usecen=1
    newlabelneed=1
    dr = withinr/nogrid
    from crtestfunction import *
        for runtodo in dirneed:
        info=outdirname(runtodo, Nsnap)
        havecr=info['havecr']
        if havecr==0:
            continue
        print 'runtodo', runtodo
        radl = np.linspace(dr,withinr,num=nogrid)
        Lgammalist=radl*0.
        nooftimes=0
                for i in range(startno,Nsnap,snapsep):
            info=outdirname(runtodo, i)
            rundir=info['rundir']
            runtitle=info['runtitle']
            slabel=info['slabel']
            snlabel=info['snlabel']
            dclabel=info['dclabel']
            resolabel=info['resolabel']
            the_snapdir=info['the_snapdir']
            Nsnapstring=info['Nsnapstring']
            havecr=info['havecr']
            Fcal=info['Fcal']
            iavesfr=info['iavesfr']
            timestep=info['timestep']
            maindir=info['maindir']
            haveB=info['haveB']
            color=info['color']
            newlabel=info['newlabel']
            ptitle=title
            if runtitle=='SMC':
                labelneed=dclabel
                if newlabelneed==1:
                    labelneed="\n".join(wrap(newlabel,17))
            else:
                labelneed=''
            if runtitle=='SMC':
                ptitle='Dwarf'
            elif runtitle=='SBC':
                ptitle='Starburst'
            elif runtitle=='MW':
                ptitle=r'$L\star$ Galaxy'
            if haveB==0:
                lsn='solid'
            else:
                lsn='dashed'
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gp = G['p']
                        Grho = G['rho']
                        Gu = G['u']
                        Gm = G['m']
                        cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                        Neb = G['ne']
                        #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                        Gz = Gp[:,2]
                        Gx = Gp[:,0]
                        Gy = Gp[:,1]
            if usecen==1:
                datasup=1
                xcen = findcenz(runtodo,Nsnap,withinr=withinr,dir='x',datasup=datasup,Gx=Gx,Gy=Gy,Gz=Gz,Gm=Gm)
                ycen = findcenz(runtodo,Nsnap,withinr=withinr,dir='y',datasup=datasup,Gx=Gx,Gy=Gy,Gz=Gz,Gm=Gm)
                zcen = findcenz(runtodo,Nsnap,withinr=withinr,dir='z',datasup=datasup,Gx=Gx,Gy=Gy,Gz=Gz,Gm=Gm)
                Gx=Gx-xcen;Gy=Gy-ycen;Gz=Gz-zcen
            outdata=outLgamma_nism(Grho,Neb,cregy)
            Lgamma=outdata['Lgamma']
            print 'Lgamma', Lgamma
                        for irad in range(nogrid):
                                cut = (Gx*Gx+Gy*Gy+Gz*Gz < radl[irad])
                                Lgammacut = Lgamma[cut]
                Lgammalist[irad]+=np.sum(Lgammacut)
            nooftimes+=1
        Lgammalist=Lgammalist/nooftimes
        plt.plot(radl, Lgammalist, label=labelneed,color=color,ls=lsn,lw=2)
        plt.xlabel(r'${\rm r [kpc]}$',fontsize=18)
        plt.ylabel(r'$L_{\gamma} [{\rm erg/s}]$',fontsize=18)
        plt.yscale('log')
    plt.title(ptitle)
    if runtitle=='SMC':
        plt.legend(loc='best',fontsize=10,ncol=2)
        #plt.title('Gamma ray luminosity within a sphere centered at the center of the disk')
    #plt.tight_layout()
        plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
        plt.savefig('CRplot/gammasph/gammasph_'+fmeat+'.pdf',bbox_inches='tight')
        plt.clf()

if wanted=='decaytimesph':
    icolor=0
        for runtodo in dirneed:
                for i in [Nsnap]:
                        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr,timestep=outdirname(runtodo, i)
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gp = G['p']
                        #Grho = G['rho']
                        Gu = G['u']
                        Gm = G['m']
                        cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                        Neb = G['ne']
                        #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                        Gz = Gp[:,2]
                        Gx = Gp[:,0]
                        Gy = Gp[:,1]
                        dr = withinr/nogrid
                        tpi_in_yrl=[]
                        radl =[]
                        for irad in range(nogrid):
                                cut = (Gx*Gx+Gy*Gy+Gz*Gz < dr*(irad+1)*dr*(irad+1))
                                Nebcut = Neb[cut]
                                Gmcut = Gm[cut]
                                Gm_in_g = Gmcut*1e10*2e33
                                shellvol_in_cm3 = 4.0/3.0*np.pi*(np.power(dr*irad,3))*3.086e21*3.086e21*3.086e21
                                Grho_in_g_cm_3 = Gm_in_g/shellvol_in_cm3
                                Gnism_in_cm_3 = np.sum((0.78+0.22*Nebcut*0.76)/protonmass_in_g*Grho_in_g_cm_3)
                tpi_in_yr = 2e5/Gnism_in_cm_3*250.0 #pi decay time in yr
                                tpi_in_yrl = np.append(tpi_in_yrl, np.sum(tpi_in_yr))
                                radl = np.append(radl, dr*irad)
                        plt.plot(radl, tpi_in_yrl, label=runtodo, color=colortable[icolor])
            print 'icolor', icolor
        icolor +=1
    diffusioncoefficient = 3e28 #in cm^2/s
    esctime = radl*radl*3.086e21*3.086e21/diffusioncoefficient/3.2e7
    plt.plot(radl, esctime*10., label=r'Escape time for $\kappa_{di} = 3\times 10^{27}$', ls='dashed', color=colortable[1])
    plt.plot(radl, esctime, label=r'Escape time for $\kappa_{di} = 3\times 10^{28}$', ls='dashed', color=colortable[2])
        plt.xlabel('r [kpc]')
        plt.ylabel(r'$t_\pi [{\rm yr}]$')
        plt.legend(loc='best')
        plt.yscale('log')
        plt.title('Within sphere centered at the center of the disk')
        plt.savefig('decaytime_'+fmeat+'_sph.pdf')
        plt.clf()


if wanted=='decayratiosph':
        icolor=0
    diffusiontable=[0, 3e27, 3e28]
        for runtodo in dirneed:
                for i in [Nsnap]:
                        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr=outdirname(runtodo, i)
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gp = G['p']
                        #Grho = G['rho']
                        Gu = G['u']
                        Gm = G['m']
                        cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                        Neb = G['ne']
                        #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                        Gz = Gp[:,2]
                        Gx = Gp[:,0]
                        Gy = Gp[:,1]
                        dr = withinr/nogrid
                        tpi_in_yrl=[]
                        radl =[]
                        for irad in range(nogrid):
                                cut = (Gx*Gx+Gy*Gy+Gz*Gz < dr*(irad+1)*dr*(irad+1))
                                Nebcut = Neb[cut]
                                Gmcut = Gm[cut]
                                Gm_in_g = Gmcut*1e10*2e33
                                shellvol_in_cm3 = 4.0/3.0*np.pi*(np.power(dr*irad,3))*3.086e21*3.086e21*3.086e21
                                Grho_in_g_cm_3 = Gm_in_g/shellvol_in_cm3
                                Gnism_in_cm_3 = np.sum((0.78+0.22*Nebcut*0.76)/protonmass_in_g*Grho_in_g_cm_3)
                                tpi_in_yr = 2e5/Gnism_in_cm_3*250.0 #pi decay time in yr
                                tpi_in_yrl = np.append(tpi_in_yrl, np.sum(tpi_in_yr))
                                radl = np.append(radl, dr*irad)
            diffusioncoefficient = diffusiontable[icolor] #in cm^2/s
            esctime = radl*radl*3.086e21*3.086e21/diffusioncoefficient/3.2e7
                    if icolor>0:
                   plt.plot(radl, esctime/(tpi_in_yrl+esctime), label=runtodo, color=colortable[icolor])
                        print 'icolor', icolor
                icolor +=1
        plt.xlabel('r [kpc]')
        plt.ylabel(r'fraction of decayed CR $\sim 1/(1+t_{\pi}/t_{\rm esc})$')
        plt.legend(loc='best')
        plt.yscale('log')
        plt.title('Within sphere centered at the center of the disk')
        plt.savefig('decayratio_'+fmeat+'_sph.pdf')
    plt.clf()


if wanted=='gammar':
        for runtodo in dirneed:
                for i in [Nsnap]:
                        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr=outdirname(runtodo, i)
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gp = G['p']
                        Grho = G['rho']
                        Gu = G['u']
                        Gm = G['m']
                        cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                        Neb = G['ne']
                        #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                        Gz = Gp[:,2]
                        Gx = Gp[:,0]
                        Gy = Gp[:,1]
                        Gnism_in_cm_3 = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21
                        tpi_in_yr = 2e5/Gnism_in_cm_3*250.0 #pi decay time in yr
                        cregy_in_erg = cregy*solar_mass_in_g*1e10*km_in_cm*km_in_cm
                        Lgammagev = cregy_in_erg/tpi_in_yr/yr_in_sec*betapi/nopi_per_gamma #in erg/s
                        dr = withinr/nogrid
                        Lgammal=[]
                        radl =[]
                        for irad in range(nogrid):
                                cutxy = (Gx*Gx+Gy*Gy < dr*irad*dr*irad)
                                cutz = np.absolute(Gz-med)< maxlength/2.
                                cut = cutxy*cutz
                                Nebcut = Neb[cut]
                                Lgammacut = Lgammagev[cut]
                                Lgammal = np.append(Lgammal, np.sum(Lgammacut))
                                radl = np.append(radl, dr*(irad))
                        plt.plot(radl, Lgammal, label=runtodo)
        plt.xlabel('r [kpc]')
        plt.ylabel(r'enclosed $L_{\gamma} [{\rm erg/s}]$')
        plt.yscale('log')
        plt.legend(loc='best')
        plt.title('cylinder centered at disk '+str(med)+' kpc above z plane with a height ' + str(maxlength) + ' kpc')
        plt.savefig('gamma_'+fmeat+'_r_zmed'+str(med)+'.pdf')
        plt.clf()

if wanted=='gammaz':
        for runtodo in dirneed:
                for i in [Nsnap]:
                        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr=outdirname(runtodo, i)
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gp = G['p']
                        Grho = G['rho']
                        Gu = G['u']
                        Gm = G['m']
                        cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                        Neb = G['ne']
                        #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                        Gz = Gp[:,2]
                        Gx = Gp[:,0]
                        Gy = Gp[:,1]
                        Gnism_in_cm_3 = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21
                        tpi_in_yr = 2e5/Gnism_in_cm_3*250.0 #pi decay time in yr
                        cregy_in_erg = cregy*solar_mass_in_g*1e10*km_in_cm*km_in_cm
                        Lgammagev = cregy_in_erg/tpi_in_yr/yr_in_sec*betapi/nopi_per_gamma #in erg/s
                        dz = maxlength/nogrid
                        Lgammal=[]
                        zl =[]
                        for iz in range(nogrid):
                                cutxy = Gx*Gx+Gy*Gy < withinr
                                cutz = np.absolute(Gz)<dz*(iz+1)
                                cut = cutxy*cutz
                                Nebcut = Neb[cut]
                                Lgammacut = Lgammagev[cut]
                                Lgammal = np.append(Lgammal, np.sum(Lgammacut))
                                zl = np.append(zl, dz*(iz))
                        plt.plot(zl, Lgammal, label=runtodo)
        plt.xlabel('z [kpc]')
        plt.ylabel(r'enclosed $L_{\gamma} [{\rm erg/s}]$')
        plt.yscale('log')
        plt.legend(loc='best')
        plt.title('cylinder centered at disk with a radius ' + str(withinr) + ' kpc')
        plt.savefig('gamma_'+fmeat+'_z.pdf')
        plt.clf()

if wanted=='crdensol': #cosmic ray energy density at solar circle (observation ~ 1 eV/cm^3)
    rotface=1
        for runtodo in dirneed:
        snaplist=[]
        credenlist=[]
                for i in range(startno, Nsnap, snapsep):
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
            snumadd = info['snumadd']
            halostr = info['halostr']
            cosmo = info['cosmo']
                        if cosmo==1:
                                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                        else:
                                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        header = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=cosmo,cosmological=cosmo, header_only=1)
                        ascale = header['time']
                        #print 'this time', ascale
                        thisred = 1./ascale-1.
                        hubble = header['hubble']
                        if cosmo==1:
                                if usepep==1:
                                        halosingle = SF.read_halo_history_pep(rundir, Nsnap, singlesnap=1, firever=firever,halonostr=halostr, hubble=hubble, comoving=0, maindir=maindir)
                                else:
                                        halosingle = read_halo_history(rundir, halonostr=halostr,hubble=hubble, comoving=0, maindir=maindir, singlesnap=1, atime=ascale,snumadd=snumadd)
                                xcen = halosingle['x']
                                ycen = halosingle['y']
                                zcen = halosingle['z']
                                xvcen = halosingle['xv']
                                yvcen = halosingle['yv']
                                zvcen = halosingle['zv']
                                Rvirnow = halosingle['R']
                                MgAHF = halosingle['Mg']
                        else:
                                xcen=0
                                ycen=0
                                zcen=0
                                xvcen=0
                                yvcen=0
                                zvcen=0
                        Gp = G['p']
            Gv = G['v']
                        Grho = G['rho']
                        Gu = G['u']
                        Gm = G['m']
                        cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                        Neb = G['ne']
                        #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                        Gz = Gp[:,2]-zcen
                        Gx = Gp[:,0]-xcen
                        Gy = Gp[:,1]-ycen
            Gvz = Gv[:,2]-zvcen
            Gvx = Gv[:,0]-xvcen
            Gvy = Gv[:,1]-yvcen
                        if rotface==1:
                                cutr = Gr < 5. #kpc
                                Gmcutr = Gm[cutr]; #to avoid overflow
                                Gxcutr = Gx[cutr]; Gycutr = Gy[cutr]; Gzcutr = Gz[cutr];
                                Gvxcutr = Gvx[cutr]; Gvycutr = Gvy[cutr]; Gvzcutr = Gvz[cutr];
                                Lang = [0.,0.,0.]
                                for i in range(len(Gxcutr)):
                                        Lang += Gmcutr[i]*np.cross([Gxcutr[i],Gycutr[i],Gzcutr[i]],[Gvxcutr[i],Gvycutr[i],Gvzcutr[i]])
                                #test
                                Gx, Gy, Gz = SF.rotateL_to_z(Gx,Gy,Gz,Lang[0],Lang[1],Lang[2])
                                Gvx, Gvy, Gvz = SF.rotateL_to_z(Gvx,Gvy,Gvz,Lang[0],Lang[1],Lang[2])
                        cregy_in_erg = cregy*solar_mass_in_g*1e10*km_in_cm*km_in_cm
            cregy_in_eV = cregy_in_erg*erg_in_eV
            solar_radius = 8.
            wdisk=2.
            insol = solar_radius+wdisk/2.
            outsol = solar_radius-wdisk/2.
            hdisk = 2.
            cutxy = (Gx*Gx+Gy*Gy < insol*insol) & (Gx*Gx+Gy*Gy > outsol*outsol) 
            cutz = np.absolute(Gz) < hdisk/2.
            cut = cutxy*cutz
            cregy_in_eV_cut = np.sum(cregy_in_eV[cut])
            creden_in_eV_per_cm3 = cregy_in_eV_cut/(np.pi*(insol*insol-outsol*outsol))/hdisk/kpc_in_cm/kpc_in_cm/kpc_in_cm
            print 'cosmic ray energy density (eV/cm^3) around solar circle', creden_in_eV_per_cm3
            snaplist=np.append(snaplist,i)
            credenlist=np.append(credenlist,creden_in_eV_per_cm3)
                #plt.plot(snaplist,credenlist,label='kappa='+dclabel+'; '+'time step= '+timestep)
        #plt.plot(snaplist,credenlist,label='kappa='+dclabel)
        labelneed=runtodo
        if dclabelneed==1:
            labelneed=r'$\kappa=$'+dclabel
        plt.plot(snaplist,credenlist,label=labelneed)
        plt.axhline(y=1.0,ls='--',color='k')
        plt.xlabel('Myr')
        plt.xlim(xmax=Nsnap)
        plt.title(runtitle+' CR around solar circle')
        plt.ylabel(r'$E_{\rm CR}({\rm eV/cm^3})$')
        plt.yscale('log')
        plt.legend(loc='best')
        plt.savefig('CRplot/CR_solarcircle_'+fmeat+'.pdf')
        plt.clf()

if wanted=='crdenmidplane' or wanted=='gasdenmidplane': #cosmic ray energy density at solar circle (observation ~ 1 eV/cm^3)
    resoneed=0
    rotface=1
    newlabelneed=1
        rcParams['figure.figsize'] = 5,4
        for runtodo in dirneed:
        print 'runtodo', runtodo
        if wanted=='crdenmidplane':
            info=outdirname(runtodo, Nsnap)
            havecr=info['havecr']
            dclabel=info['dclabel']
            haveB=info['haveB']
            if havecr==0:
                continue
        withinr=15.
        nogrid = 15
        maxlength=0.1 #thickness
        dr = withinr/nogrid
        radlist=np.linspace(0.001,withinr,num=nogrid)
        if wanted=='crdenmidplane':
            credenlist=radlist*0.
        if wanted=='gasdenmidplane':
            gasdenlist=radlist*0.
        numoftimes=0

        for i in range(startno,Nsnap,snapsep):
            snaplist=[]
            info=outdirname(runtodo, i)
            rundir=info['rundir']
            runtitle=info['runtitle']
            slabel=info['slabel']
            snlabel=info['snlabel']
            dclabel=info['dclabel']
            resolabel=info['resolabel']
            the_snapdir=info['the_snapdir']
            Nsnapstring=info['Nsnapstring']
            havecr=info['havecr']
            Fcal=info['Fcal']
            iavesfr=info['iavesfr']
            timestep=info['timestep']
            cosmo=info['cosmo']
            maindir=info['maindir']
            color=info['color']
            haveB=info['haveB']
            M1speed=info['M1speed']
            newlabel=info['newlabel']
            snumadd=info['snumadd']
            usepep=info['usepep']
            halostr=info['halostr']
            ptitle=title
            if runtitle=='SMC':
                ptitle='Dwarf'
            elif runtitle=='SBC':
                ptitle='Starburst'
            elif runtitle=='MW':
                ptitle=r'$L\star$ Galaxy'
            labelneed=dclabel
            if newlabelneed==1:
                labelneed="\n".join(wrap(newlabel,17))
            if havecr==0 and wanted=='crdenmidplane':
                continue
            if cosmo==1:
                h0=1
            else:
                h0=0
            if cosmo==1:
                datasup=0;
            else:
                datasup=1;
            Gextra = readsnapwcen(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix,\
             havecr=havecr,h0=h0,cosmo=cosmo, usepep=usepep, maindir=maindir,snumadd=snumadd,rotface=rotface,\
             datasup=datasup,runtodo=runtodo,rundir=rundir,halostr=halostr)
            Gx = Gextra['x']; Gy = Gextra['y']; Gz = Gextra['z'];
            Gvx = Gextra['vx']; Gvy = Gextra['vy']; Gvz = Gextra['vz'];
                        Grho = Gextra['rho']; Gu = Gextra['u']; Gm = Gextra['m']
                        if wanted=='crdenmidplane':
                                cregy = Gextra['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
            if wanted=='crdenmidplane':
                cregy_in_erg = cregy*solar_mass_in_g*1e10*km_in_cm*km_in_cm
                cregy_in_eV = cregy_in_erg*erg_in_eV
            for irad in range(len(radlist)-1):
                cutxy = ((Gx)*(Gx)+(Gy)*(Gy) > radlist[irad]*radlist[irad]) & ((Gx)*(Gx)+(Gy)*(Gy) < radlist[irad+1]*radlist[irad+1])
                cutz = (Gz)*(Gz) < maxlength*maxlength/4.
                cut = cutxy*cutz
                shellvol_in_cm3 = np.pi*(-np.power(radlist[irad],2)+np.power(radlist[irad+1],2))*kpc_in_cm*kpc_in_cm*kpc_in_cm*maxlength
                if wanted=='crdenmidplane':
                    cregy_in_eV_cut = np.sum(cregy_in_eV[cut])
                    creden_in_eV_per_cm3 = cregy_in_eV_cut/shellvol_in_cm3
                    credenlist[irad]+=creden_in_eV_per_cm3
                if wanted=='gasdenmidplane':
                    Gm_in_g=Gm[cut]*1e10*Msun_in_g
                    Grho_in_g_cm_3 = np.sum(Gm_in_g)/shellvol_in_cm3
                    try:
                        Nebave = np.average(Neb[cut],weights=Gm[cut])
                    except ZeroDivisionError:
                        Nebave = 0
                    Gnism_in_cm_3 = 1.0/proton_mass_in_g*Grho_in_g_cm_3
                    #print 'radlist, Gnism_in_cm_3', radlist[irad], Gnism_in_cm_3
                    gasdenlist[irad]+=Gnism_in_cm_3
            numoftimes+=1
        if haveB>0:
            lsn='dashed'
        else:
            lsn='solid'
                if resoneed==1:
                        if resolabel=='llr':
                                labelneed='low res'
                                lsn = 'solid'
                        if resolabel=='lr':
                                labelneed='std res'
                                lsn = 'dashed'
                        if resolabel=='mr':
                                labelneed='high res'
                                lsn = 'dashdot'
        #if M1speed>999:
        #   lsn='dotted'
        #if M1speed>1999:
    #       lsn='dashdot'
        if wanted=='crdenmidplane': 
            plt.plot(radlist[:-1],credenlist[:-1]/numoftimes,label=labelneed,lw=2,ls=lsn,color=color)
        if wanted=='gasdenmidplane':
            plt.plot(radlist[:-1],gasdenlist[:-1]/numoftimes,label=labelneed,lw=2,ls=lsn,color=color)
    if wanted=='crdenmidplane':
        #plt.plot(8,1.8,ls='none', marker='*', markersize=3, color=0.5)
        plt.xlabel(r'${\rm r [kpc]}$', fontsize=18)
        #plt.axhline(y=1.0,ls='--',color='k')
        plt.xlim(xmax=np.amax(radlist))
        plt.ylabel(r'$<e_{\rm CR}>[{\rm eV/cm^3}]$', fontsize=18)
        plt.yscale('log')
        #plt.legend(loc='best',frameon=False)
        filename='CRplot/crdenmidplane/CR_midplane_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
        if wanted=='gasdenmidplane':
                plt.xlabel(r'${\rm r [kpc]}$', fontsize=18)
                plt.xlim(xmax=np.amax(radlist))
                #plt.title(runtitle+' time averaged gas density at midplane (width='+str(2*maxlength)+'kpc)')
                plt.ylabel(r'$<n_{\rm ISM}>[{\rm cm^{-3}}]$',fontsize=18)
        if runtitle=='SMC':
            plt.legend(loc='best',fontsize=10,ncol=2)
                filename='CRplot/gasdenmidplane/gas_midplane_'+fmeat+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    print 'filename', filename
        plt.title(ptitle,fontsize=18)
        plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
    plt.savefig(filename,bbox_inches='tight')
    plt.clf()



if wanted=='crdenv': #cosmic ray energy density at solar circle (observation ~ 1 eV/cm^3)
        for runtodo in dirneed:
                snaplist=[]
                credenlist=[]
                info=outdirname(runtodo, Nsnap)
                rundir=info['rundir']
                runtitle=info['runtitle']
                slabel=info['slabel']
                snlabel=info['snlabel']
                dclabel=info['dclabel']
                resolabel=info['resolabel']
                the_snapdir=info['the_snapdir']
                Nsnapstring=info['Nsnapstring']
                havecr=info['havecr']
                Fcal=info['Fcal']
                iavesfr=info['iavesfr']
                timestep=info['timestep']
        cosmo=info['cosmo']
        if cosmo==1:
            h0=1
        else:
            h0=0
        header = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo, header_only=1)
        ascale = header['time']
        print 'this time', ascale
        thisred = 1./ascale-1.
        hubble = header['hubble']
        print 'hubble', hubble
        print 'the_snapdir', the_snapdir
        G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo)
                Gp = G['p']
                Grho = G['rho']
                Gu = G['u']
                Gm = G['m']
                cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                Neb = G['ne']
                #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                Gz = Gp[:,2]
                Gx = Gp[:,0]
                Gy = Gp[:,1]
        if cosmo==1:
            halosingle = read_halo_history(rundir, halonostr='00',hubble=hubble, comoving=0, maindir=maindir, singlesnap=1, atime=ascale)
            xcen = halosingle['x']
            ycen = halosingle['y']
            zcen = halosingle['z']
            xvcen = halosingle['xv']
            yvcen = halosingle['yv']
            zvcen = halosingle['zv']
            Rvirnow = halosingle['R']
            MgAHF = halosingle['Mg']
        #       print 'xcenl[0]', xcenl[0]
        #       print 'thisred', thisred
            print 'cen', xcen, ycen, zcen
            print 'MgAHF', MgAHF
            #print 'Rvir', Rvirnow
        else:
            xcen=0
            ycen=0
            zcen=0
            xvcen=0
            yvcen=0
            zvcen=0
                cregy_in_erg = cregy*solar_mass_in_g*1e10*km_in_cm*km_in_cm
                cregy_in_eV = cregy_in_erg*erg_in_eV
                hdisk=5.
                rdisk = 3.
        dz = hdisk/nogrid
        zl=[]
                for iv in range(nogrid):
                        cutxy = ((Gx-xcen)*(Gx-xcen)+(Gy-ycen)*(Gy-ycen) < rdisk*rdisk)
                        cutz = (np.absolute(Gz-zcen) < dz*(iv+1.0)) & (np.absolute(Gz-zcen) > dz*(iv))
                        cut = cutxy*cutz
                        cregy_in_eV_cut = np.sum(cregy_in_eV[cut])
                        creden_in_eV_per_cm3 = cregy_in_eV_cut/(np.pi*rdisk*rdisk*dz*2.0)/kpc_in_cm/kpc_in_cm/kpc_in_cm
                        credenlist=np.append(credenlist,creden_in_eV_per_cm3)
            zl = np.append(zl,dz*(iv+0.5))
                plt.plot(zl,credenlist,label=runtodo)
        plt.axhline(y=1.0,ls='--',color='k')
        plt.xlabel('z (kpc)')
        plt.xlim(xmax=np.amax(zl))
        plt.title(' CR energy density within radius ='+str(rdisk)+' kpc')
        plt.ylabel(r'$E_{\rm CR}({\rm eV/cm^3})$')
        #plt.yscale('log')
        plt.legend(loc='best')
        plt.savefig('CRplot/CR_vertical_'+fmeat+'.pdf')
        plt.clf()



if wanted=='crdenplanes': #cosmic ray energy density at solar circle (observation ~ 1 eV/cm^3)
        for runtodo in dirneed:
                snaplist=[]
                creden0list=[]
                creden1list=[]
                creden2list=[]
                info=outdirname(runtodo, Nsnap)
                rundir=info['rundir']
                runtitle=info['runtitle']
                slabel=info['slabel']
                snlabel=info['snlabel']
                dclabel=info['dclabel']
                resolabel=info['resolabel']
                the_snapdir=info['the_snapdir']
                Nsnapstring=info['Nsnapstring']
                havecr=info['havecr']
                Fcal=info['Fcal']
                iavesfr=info['iavesfr']
                timestep=info['timestep']
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                Gp = G['p']
                Grho = G['rho']
                Gu = G['u']
                Gm = G['m']
                cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                Neb = G['ne']
                #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                Gz = Gp[:,2]
                Gx = Gp[:,0]
                Gy = Gp[:,1]
                cregy_in_erg = cregy*solar_mass_in_g*1e10*km_in_cm*km_in_cm
                cregy_in_eV = cregy_in_erg*erg_in_eV
                radlist = np.linspace(1.1,20,num=10)
                wdisk=2.
                hdisk = 0.25
                for radius in radlist:
                        insol = radius+wdisk/2.
                        outsol = radius-wdisk/2.
                        cutxy = (Gx*Gx+Gy*Gy < insol*insol) & (Gx*Gx+Gy*Gy > outsol*outsol)
                        cutz0 = np.absolute(Gz) < hdisk/2.
            cutz1 = np.absolute(Gz-1.0) < hdisk/2.
            cutz2 = np.absolute(Gz-2.0) < hdisk/2.
                        cut0 = cutxy*cutz0
            cut1 = cutxy*cutz1
            cut2 = cutxy*cutz2
                        cregy_in_eV_cut0 = np.sum(cregy_in_eV[cut0])
                        cregy_in_eV_cut1 = np.sum(cregy_in_eV[cut1])
                        cregy_in_eV_cut2 = np.sum(cregy_in_eV[cut2])
                        creden0_in_eV_per_cm3 = cregy_in_eV_cut0/(np.pi*(insol*insol-outsol*outsol))/hdisk/kpc_in_cm/kpc_in_cm/kpc_in_cm
                        creden1_in_eV_per_cm3 = cregy_in_eV_cut1/(np.pi*(insol*insol-outsol*outsol))/hdisk/kpc_in_cm/kpc_in_cm/kpc_in_cm
                        creden2_in_eV_per_cm3 = cregy_in_eV_cut2/(np.pi*(insol*insol-outsol*outsol))/hdisk/kpc_in_cm/kpc_in_cm/kpc_in_cm
                        creden0list=np.append(creden0list,creden0_in_eV_per_cm3)
            creden1list=np.append(creden1list,creden1_in_eV_per_cm3)
            creden2list=np.append(creden2list,creden2_in_eV_per_cm3)
                #plt.plot(radlist,credenlist,label='kappa='+dclabel+'; '+'time step= '+timestep)
                #plt.plot(radlist,credenlist,label='kappa='+dclabel)
                plt.plot(radlist,creden0list,label='midplane')
                plt.plot(radlist,creden1list,label='1kpc')
                plt.plot(radlist,creden2list,label='2kpc')
        plt.xlabel('r (kpc)')
        plt.xlim(xmax=np.amax(radlist))
        plt.title(runtitle+' CR energy density')
        plt.ylabel(r'$E_{\rm CR}({\rm eV/cm^3})$')
        plt.yscale('log')
        plt.legend(loc='best')
        plt.savefig('CRplot/CR_planes_'+runtodo+'.pdf')
        plt.clf()




if wanted=='cramap' or wanted=='cramapv':
    withinr=30
    maxlength=30
    rotface=1
    rcParams['figure.figsize'] = 5, 5
        for runtodo in dirneed:
                snaplist=[]
                enclist=[]
                englist=[]
                enllist=[]
                endlist=[]
                enplist=[]
                enalist=[]
                avesfrl=[]
                prel=0
                preg=0
                prec=0
                pred=0
                prep=0
                presm = 0
                presnap = 0
                for i in [Nsnap]:
            info=outdirname(runtodo, i)
            rundir=info['rundir'] 
            runtitle=info['runtitle']
            slabel=info['slabel']
            snlabel=info['snlabel']
            dclabel=info['dclabel']
            resolabel=info['resolabel']
            the_snapdir=info['the_snapdir']
            Nsnapstring=info['Nsnapstring']
            havecr=info['havecr']
            Fcal=info['Fcal']
            iavesfr=info['iavesfr']
            timestep=info['timestep']
            maindir=info['maindir']
            haveB=info['haveB']
            usepep=info['usepep']
            snapsep=info['snapsep']
            halostr=info['halostr']
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
            cosmo=info['cosmo']
            firever = info['firever']
            if cosmo==1:
                h0=1
            else:
                h0=0
            header = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo, header_only=1)
            ascale = header['time']
            print 'this time', ascale
            thisred = 1./ascale-1.
            hubble = header['hubble']
            print 'hubble', hubble
            G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo)
            if cosmo==1:
                if usepep==1:
                    halosingle = SF.read_halo_history_pep(rundir, Nsnap, singlesnap=1,\
                     firever=firever,halonostr=halostr, hubble=hubble, comoving=1, maindir=maindir)
                    afactor=atime
                else:
                    halosingle = read_halo_history(rundir, halonostr='00',hubble=hubble, comoving=0, maindir=maindir, singlesnap=1, atime=ascale)
                    afactor=1.0
                xcen = halosingle['x']*afactor
                ycen = halosingle['y']*afactor
                zcen = halosingle['z']*afactor
                xvcen = halosingle['xv']
                yvcen = halosingle['yv']
                zvcen = halosingle['zv']
                Rvirnow = halosingle['R']*afactor
                MgAHF = halosingle['Mg']
            #       print 'xcenl[0]', xcenl[0]
            #       print 'thisred', thisred
                print 'cen', xcen, ycen, zcen
                print 'MgAHF', MgAHF
                #print 'Rvir', Rvirnow
            else:
                xcen=0
                ycen=0
                zcen=0
                xvcen=0
                yvcen=0
                zvcen=0
            Gpos = G['p']
            Gvel = G['v']
            Gid_t = G['id']
            Gx = Gpos[:,0]-xcen
            Gy = Gpos[:,1]-ycen
            Gz = Gpos[:,2]-zcen
            Gvx = Gvel[:,0]-xvcen   #km/s
            Gvy = Gvel[:,1]-yvcen
            Gvz = Gvel[:,2]-zvcen
            if cosmo==1 and rotface==1:
                Lxstar = halosingle['Lxstar']; Lystar = halosingle['Lystar']; Lzstar = halosingle['Lzstar']
                Gx, Gy, Gz = SF.rotateL_to_z(Gx,Gy,Gz,Lxstar,Lystar,Lzstar)
            Grho = G['rho']*1e10 #Msun/kpc^3
            Gm = G['m']*1e10 #Msun
            GEint = G['u']*km_in_cm*km_in_cm*Gm*solar_mass_in_g
            if havecr>0:
                cregyl = G['cregyl']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                cregyg = G['cregyg']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                cregy  = G['cregy']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                cregyd = G['cregyd']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                        if havecr>4:
                                cregyp = G['cregyp']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                        if havecr>5:
                                cregya = G['cregya']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
            if wanted == 'cramap':
                cutz = np.absolute(Gz)<maxlength
            elif wanted == 'cramapv':
                cutz = np.absolute(Gy)<maxlength
            if haveB>0:
                Bfield = G['B']
                Bx = Bfield[:,0]
                By = Bfield[:,1]
                Bz = Bfield[:,2]
            Gx=Gx[cutz]
            Gy=Gy[cutz]
            Gz=Gz[cutz]
            Gvx=Gvx[cutz]
            Gvy=Gvy[cutz]
            Gvz=Gvz[cutz]
            Gm=Gm[cutz]
            Grho=Grho[cutz]
            GEint=GEint[cutz]
            Gid_t=Gid_t[cutz]
            if havecr>0:
                cregyl=cregyl[cutz]
                cregyg=cregyg[cutz]
                cregy=cregy[cutz]
                cregyd=cregyd[cutz]
            
                        if havecr>4:
                                cregyp =cregyp[cutz]
                        if havecr>5:
                                cregya = cregya[cutz]
                #print 'cregyp', cregyp
                #print 'np.sum(cregyp)', np.sum(cregyp)
            if haveB>0:
                Bx = Bx[cutz]
                By = By[cutz]
                Bz = Bz[cutz] 
                                B2 = Bx*Bx+By*By+Bz*Bz
                Begy = B2/8./np.pi*Gm/Grho*kpc_in_cm*kpc_in_cm*kpc_in_cm
            cellvol = 2.0*withinr*2.0*withinr*2.0*maxlength/100./100.*kpc_in_cm*kpc_in_cm*kpc_in_cm
            snaptrack=i-snapsep
            info=outdirname(runtodo, snaptrack)
            Nsnapstring = info['Nsnapstring']
            G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo)
            Gid = G['id']
            Gid, unindex = np.unique(Gid_t, return_index=True)
            Gid_t, unindex_t = np.unique(Gid_t, return_index=True)
            idint0=np.in1d(Gid,Gid_t)
            # there are some duplicate IDs, from particle splitting? I drop those repetitions.
            Gid = Gid[idint0]
            Gidinds = Gid.argsort()
            Gid_tinds = Gid_t.argsort()
            timeneed = i*0.98*1e6
            timetrack=snaptrack*0.98*1e6
                        if havecr>0:
                                cregyl_p = G['cregyl']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                                cregyl_p = cregyl_p[unindex_t]
                                cregyl = cregyl[unindex]
                                dcrl = (cregyl[Gid_tinds]-cregyl_p[Gidinds])/(timeneed-timetrack)/yr_in_sec
                                cregyg_p = G['cregyg']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                                cregyg_p = cregyg_p[unindex_t]
                                cregyg = cregyg[unindex]
                                dcrg = (cregyg[Gid_tinds]-cregyg_p[Gidinds])/(timeneed-timetrack)/yr_in_sec
                                cregyd_p = G['cregyd']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                                cregyd_p = cregyd_p[unindex_t]
                                cregyd = cregyd[unindex]
                                dcrd = (cregyd[Gid_tinds]-cregyd_p[Gidinds])/(timeneed-timetrack)/yr_in_sec
                        if havecr>4:
                                cregyp_p = G['cregyp']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                                cregyp_p = cregyp_p[unindex_t]
                                cregyp = cregyp[unindex]
                                dcrp = (cregyp[Gid_tinds]-cregyp_p[Gidinds])/(timeneed-timetrack)/yr_in_sec
                        if havecr>5:
                                cregya_p = G['cregya']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                cregya_p = cregya_p[unindex_t]
                cregya = cregya[unindex]
                dcra = (cregya[Gid_tinds]-cregya_p[Gidinds])/(timeneed-timetrack)/yr_in_sec
            Gxd = Gx[unindex_t]
            Gxd = Gxd[Gid_tinds]
            Gyd = Gy[unindex_t]
            Gyd = Gyd[Gid_tinds]
            Gzd = Gz[unindex_t]
            Gzd = Gzd[Gid_tinds]
            if havecr>0:
                if wanted == 'cramap':
                    Hm, xedges, yedges = np.histogram2d(Gyd, Gxd, bins=100,\
                    range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(dcrl/cellvol))
                    plt.xlabel('x (kpc)')
                    plt.ylabel('y (kpc)')
                elif wanted == 'cramapv':
                    Hm, xedges, yedges = np.histogram2d(Gzd, Gxd, bins=100,\
                    range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(dcrl/cellvol))
                    plt.xlabel('x (kpc)')
                    plt.ylabel('z (kpc)')
                im = plt.imshow(np.log10(Hm), interpolation='nearest', origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
                plt.colorbar(im,fraction=0.046, pad=0.04)
                plt.tight_layout()
                if wanted == 'cramap':
                    plt.savefig('CRplot/cramap/'+runtodo+'_cregyl_sn'+str(Nsnap)+'.pdf')
                elif wanted == 'cramapv':
                    plt.savefig('CRplot/cramapv/'+runtodo+'_cregyl_sn'+str(Nsnap)+'_v.pdf')
                plt.clf()
            if havecr>0:
                if wanted == 'cramap':
                    Hm, xedges, yedges = np.histogram2d(Gyd, Gxd, bins=100,\
                    range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(dcrg/cellvol))
                    plt.xlabel('x (kpc)')
                    plt.ylabel('y (kpc)')
                elif wanted == 'cramapv':
                    Hm, xedges, yedges = np.histogram2d(Gzd, Gxd, bins=100,\
                    range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(dcrg/cellvol))
                    plt.xlabel('x (kpc)')
                    plt.ylabel('z (kpc)')
                im = plt.imshow(np.log10(Hm),interpolation='nearest', origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
                plt.colorbar(im,fraction=0.046, pad=0.04)
                plt.tight_layout()
                if wanted == 'cramap':
                    plt.savefig('CRplot/cramap/'+runtodo+'_cregyg_sn'+str(Nsnap)+'.pdf')
                elif wanted == 'cramapv':
                    plt.savefig('CRplot/cramapv/'+runtodo+'_cregyg_sn'+str(Nsnap)+'_v.pdf')
                plt.clf()
            if havecr>4:
                if wanted=='cramap':
                    Hm, xedges, yedges = np.histogram2d(Gyd, Gxd, bins=100,\
                    range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(dcrp/cellvol))
                    plt.xlabel('x (kpc)')
                    plt.ylabel('y (kpc)')
                elif wanted == 'cramapv':
                    Hm, xedges, yedges = np.histogram2d(Gzd, Gxd, bins=100,\
                    range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(dcrp/cellvol))
                    plt.xlabel('x (kpc)')
                    plt.ylabel('z (kpc)')
                im = plt.imshow(np.log10(Hm),interpolation='nearest', origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
                plt.colorbar(im,fraction=0.046, pad=0.04)
                plt.tight_layout()
                if wanted == 'cramap':
                    plt.savefig('CRplot/cramap/'+runtodo+'_cregyp_sn'+str(Nsnap)+'.pdf')
                elif wanted == 'cramapv':
                    plt.savefig('CRplot/cramapv/'+runtodo+'_cregyp_sn'+str(Nsnap)+'_v.pdf')
                plt.clf()
            if havecr>0:
                if wanted=='cramap':
                    Hmc, xedges, yedges = np.histogram2d(Gy, Gx, bins=100,\
                    range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(cregy/cellvol))
                    plt.xlabel('x (kpc)')
                    plt.ylabel('y (kpc)')
                elif wanted == 'cramapv':
                    Hmc, xedges, yedges = np.histogram2d(Gz, Gx, bins=100,\
                    range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(cregy/cellvol))
                    plt.xlabel('x (kpc)')
                    plt.ylabel('z (kpc)')
                im = plt.imshow(np.log10(Hmc), interpolation='nearest', origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
                plt.colorbar(im,fraction=0.046, pad=0.04)
                plt.tight_layout()
                if wanted == 'cramap':
                    plt.savefig('CRplot/cramap/'+runtodo+'_cregy_sn'+str(Nsnap)+'.pdf')
                elif wanted == 'cramapv':
                    plt.savefig('CRplot/cramapv/'+runtodo+'_cregy_sn'+str(Nsnap)+'_v.pdf')
                plt.clf()
                        if wanted=='cramap':
                                Hm, xedges, yedges = np.histogram2d(Gy, Gx, bins=100,range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(Gm/cellvol))
                                plt.xlabel('x (kpc)')
                                plt.ylabel('y (kpc)')
                        elif wanted == 'cramapv':
                                Hm, xedges, yedges = np.histogram2d(Gz, Gx, bins=100,range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(Gm/cellvol))
                                plt.xlabel('x (kpc)')
                                plt.ylabel('z (kpc)')
                        im = plt.imshow(np.log10(Hm), interpolation='nearest', origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
                        plt.colorbar(im,fraction=0.046, pad=0.04)
                        plt.tight_layout()
                        if wanted == 'cramap':
                                plt.savefig('CRplot/cramap/'+runtodo+'_gm_sn'+str(Nsnap)+'.pdf')
                        elif wanted == 'cramapv':
                                plt.savefig('CRplot/cramapv/'+runtodo+'_gm_sn'+str(Nsnap)+'_v.pdf')
                        plt.clf()
                        if wanted == 'cramapv':
                Hm, xedges, yedges = np.histogram2d(Gz, Gx, bins=100,range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(Gm))
                                Hmz, xedges, yedges = np.histogram2d(Gz, Gx, bins=100,range=[[-withinr,withinr],[-withinr,withinr]], weights=Gm*Gvz)
                                plt.xlabel('x (kpc)')
                                plt.ylabel('z (kpc)')
                                im = plt.imshow(np.log10(np.absolute(Hmz))-np.log10(Hm), interpolation='nearest', origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
                                plt.colorbar(im,fraction=0.046, pad=0.04)
                                plt.tight_layout()
                plt.title('mass averaged velocity [log(km/s)]')
                plt.savefig('CRplot/cramapv/'+runtodo+'_vz_sn'+str(Nsnap)+'.pdf')
                plt.clf()
            if havecr>4:
                                cutpos=dcra<0.0
                                Gyn = Gyd[cutpos]
                                Gxn = Gxd[cutpos]
                                Gzn = Gzd[cutpos]
                                dcraneg = dcra[cutpos]
                if wanted=='cramap':
                    Hm, xedges, yedges = np.histogram2d(Gyn, Gxn, bins=100,\
                    range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(dcraneg/cellvol))
                    plt.xlabel('x (kpc)')
                    plt.ylabel('y (kpc)')
                elif wanted == 'cramapv':
                    Hm, xedges, yedges = np.histogram2d(Gzn, Gxn, bins=100,\
                    range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(dcraneg/cellvol))
                    plt.xlabel('x (kpc)')
                    plt.ylabel('z (kpc)')
                im = plt.imshow(np.log10(Hm), interpolation='nearest', origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
                plt.colorbar(im,fraction=0.046, pad=0.04)
                plt.tight_layout()
                if wanted == 'cramap':
                    plt.savefig('CRplot/cramap/'+runtodo+'_cread_neg_sn'+str(Nsnap)+'.pdf')
                elif wanted == 'cramapv':
                    plt.savefig('CRplot/cramapv/'+runtodo+'_cread_neg_sn'+str(Nsnap)+'_v.pdf')
                plt.clf()
                cutneg=dcra>0.0
                Gyp = Gyd[cutneg]
                Gxp = Gxd[cutneg]
                Gzp = Gzd[cutneg]
                dcrapos = dcra[cutneg]
                if wanted=='cramap':
                    Hm, xedges, yedges = np.histogram2d(Gyp, Gxp, bins=100,\
                    range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(dcrapos/cellvol))
                    plt.xlabel('x (kpc)')
                    plt.ylabel('y (kpc)')
                elif wanted == 'cramapv':
                    Hm, xedges, yedges = np.histogram2d(Gzp, Gxp, bins=100,\
                    range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(dcrapos/cellvol))
                    plt.xlabel('x (kpc)')
                    plt.ylabel('z (kpc)')
                im = plt.imshow(np.log10(Hm), interpolation='nearest', origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
                plt.colorbar(im,fraction=0.046, pad=0.04)
                plt.tight_layout()
                if wanted == 'cramap':
                    plt.savefig('CRplot/cramap/'+runtodo+'_cread_pos_sn'+str(Nsnap)+'.pdf')
                elif wanted == 'cramapv':
                    plt.savefig('CRplot/cramapv/'+runtodo+'_cread_pos_sn'+str(Nsnap)+'_v.pdf')
                plt.clf()
            if wanted=='cramap':
                                Hmi, xedges, yedges = np.histogram2d(Gy, Gx, bins=100,range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(GEint/cellvol))
                                plt.xlabel('x (kpc)')
                                plt.ylabel('y (kpc)')
                        elif wanted == 'cramapv':
                                Hmi, xedges, yedges = np.histogram2d(Gz, Gx, bins=100,range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(GEint/cellvol))
                                plt.xlabel('x (kpc)')
                                plt.ylabel('z (kpc)')
                        im = plt.imshow(np.log10(Hmi), interpolation='nearest', origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
                        plt.colorbar(im,fraction=0.046, pad=0.04)
                        plt.tight_layout()
                        if wanted == 'cramap':
                                plt.savefig('CRplot/cramap/'+runtodo+'_eint_sn'+str(Nsnap)+'.pdf')
                        elif wanted == 'cramapv':
                                plt.savefig('CRplot/cramapv/'+runtodo+'_eint_sn'+str(Nsnap)+'_v.pdf')
                        plt.clf()
            if havecr>0:
                if wanted=='cramap':
                    plt.xlabel('x (kpc)')
                    plt.ylabel('y (kpc)')
                elif wanted == 'cramapv':
                    plt.xlabel('x (kpc)')
                    plt.ylabel('z (kpc)')
                im = plt.imshow(np.log10(Hmc)-np.log10(Hmi), interpolation='nearest', origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
                plt.colorbar(im,fraction=0.046, pad=0.04)
                plt.tight_layout()
                if wanted == 'cramap':
                    plt.savefig('CRplot/cramap/'+runtodo+'_cre_eint_sn'+str(Nsnap)+'.pdf')
                elif wanted == 'cramapv':
                    plt.savefig('CRplot/cramapv/'+runtodo+'_cre_eint_sn'+str(Nsnap)+'_v.pdf')
                plt.clf()
                        if haveB>0:
                                if wanted=='cramap':
                    HmB, xedges, yedges = np.histogram2d(Gy, Gx, bins=100,range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(Begy/cellvol))
                                        plt.xlabel('x (kpc)')
                                        plt.ylabel('y (kpc)')
                                elif wanted == 'cramapv':
                    HmB, xedges, yedges = np.histogram2d(Gz, Gx, bins=100,range=[[-withinr,withinr],[-withinr,withinr]], weights=np.absolute(Begy/cellvol))
                                        plt.xlabel('x (kpc)')
                                        plt.ylabel('z (kpc)')
                                im = plt.imshow(np.log10(HmB)/2., interpolation='nearest', origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
                                plt.colorbar(im,fraction=0.046, pad=0.04)
                                plt.tight_layout()
                plt.title('Magnetic energy per pixel')
                                if wanted == 'cramap':
                                        plt.savefig('CRplot/cramap/'+runtodo+'_B_sn'+str(Nsnap)+'.pdf')
                                elif wanted == 'cramapv':
                                        plt.savefig('CRplot/cramapv/'+runtodo+'_B_sn'+str(Nsnap)+'v.pdf')
                                plt.clf()


if wanted=='crarad':
        for runtodo in dirneed:
                snaplist=[]
                enclist=[]
                englist=[]
                enllist=[]
                endlist=[]
                enplist=[]
                enalist=[]
                avesfrl=[]
                prel=0
                preg=0
                prec=0
                pred=0
                prep=0
                presm = 0
                presnap = []
        crecuml = []
        crecumg = []
        crecuma = []
        crecumd = []
        crecum  = []
        crecump = []
                for i in [Nsnap]:
                        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr, timestep=outdirname(runtodo, i)
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gpos = G['p']
                        Gx = Gpos[:,0]
                        Gy = Gpos[:,1]
                        Gz = Gpos[:,2]
                        Grho = G['rho']
                        Gm = G['m']*1e10
            Gr = np.sqrt(np.square(Gx)+np.square(Gy)+np.square(Gz))
                        cregyl = G['cregyl']/1e6/3.2e7*2e53
                        cregyg = G['cregyg']/1e6/3.2e7*2e53
                        cregy  = G['cregy']/1e6/3.2e7*2e53
                        cregyd = G['cregyd']/1e6/3.2e7*2e53
                        if havecr>4:
                                cregyp = G['cregyp']/1e6/3.2e7*2e53
                        if havecr>5:
                                cregya = G['cregya']/1e6/3.2e7*2e53
            rad=np.linspace(0.1,withinr, num=nogrid)        
            for i in range(len(rad)):
                crecuml = np.append(crecuml,np.sum(cregyl[Gr<rad[i]]))
                crecumg = np.append(crecumg, np.sum(cregyg[Gr<rad[i]]))
                crecuma = np.append(crecuma, np.sum(cregya[Gr<rad[i]]))
                crecum = np.append(crecum, np.sum(cregy[Gr<rad[i]]))
                crecumd = np.append(crecumd, np.sum(cregyd[Gr<rad[i]]))
                crecump = np.append(crecump, np.sum(cregyp[Gr<rad[i]]))
            crecumad = crecumd + crecumg
            plt.plot(rad, crecum,label='CR energy')
            plt.plot(rad, crecumg, label='SNe')
            plt.plot(rad, crecuml, label='Loss')
            plt.plot(rad, crecumd-crecuml-crecump, label='Ad')
            plt.plot(rad, crecuma, label='test Ad')
            plt.plot(rad, crecump, label='Flux')
            plt.plot(rad, crecumad, label='SNe+Loss+Ad+Flux')
            plt.legend()
            plt.ylabel(r'$\left \langle \rm{dE/dt}  \right \rangle$ (erg/s)')
            plt.xlabel('r (kpc)')
            plt.savefig('CRplot/'+runtodo+'_rad.pdf')
            plt.clf()

if wanted=='sfrrad' or wanted=='nsrad' or wanted == 'cumnsmrad':
        rcParams['figure.figsize'] = 5, 4
    withinr = 15.
    minr = 0.5
    maxlength=0.25
    usecen=0
    normalized=1
    usenscen=1
    usesphrad=1
    uselog=1
    Myrneed=500
    ageneed=200 #Myr
        for runtodo in dirneed:
        info=outdirname(runtodo, Myrneed)
        rundir=info['rundir']
        runtitle=info['runtitle']
        slabel=info['slabel']
        snlabel=info['snlabel']
        dclabel=info['dclabel']
        resolabel=info['resolabel']
        the_snapdir=info['the_snapdir']
        Nsnapstring=info['Nsnapstring']
        havecr=info['havecr']
        haveB=info['haveB']
        Fcal=info['Fcal']
        iavesfr=info['iavesfr']
        timestep=info['timestep']
        color=info['color']
                ptitle=title
                if runtitle=='SMC':
                        ptitle='Dwarf'
                elif runtitle=='SBC':
                        ptitle='Starburst'
                elif runtitle=='MW':
                        ptitle=r'$L\star$ Galaxy'
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix)
                Gpos = G['p']
                Gx = Gpos[:,0]
                Gy = Gpos[:,1]
                Gz = Gpos[:,2]
        Gm = G['m']
                S = readsnapcr(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix)
                Sage = S['age']
        Spos = S['p']
        Sx = Spos[:,0]
        Sy = Spos[:,1]
        Sz = Spos[:,2]
                Sm = S['m']
        if usecen==1 and cosmo==0:
                        from crtestfunction import findcenz
                        datasup=1
                        xcen = findcenz(runtodo,Nsnap,withinr=withinr,dir='x',datasup=datasup,Gx=Gx,Gy=Gy,Gz=Gz,Gm=Gm)
                        ycen = findcenz(runtodo,Nsnap,withinr=withinr,dir='y',datasup=datasup,Gx=Gx,Gy=Gy,Gz=Gz,Gm=Gm)
                        zcen = findcenz(runtodo,Nsnap,withinr=withinr,dir='z',datasup=datasup,Gx=Gx,Gy=Gy,Gz=Gz,Gm=Gm)
            Sx=Sx-xcen;Sy=Sy-ycen;Sz=Sz-zcen
            if cosmo==1:
                Sp = S['p']
                Sx = Sp[:,0]
                Sy = Sp[:,1]
                Sz = Sp[:,2]
                header=readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1,header_only=1)
                h0 = header['hubble']
                atime = header['time']
                if usepep==1:
                    halosA = SF.read_halo_history_pep(rundir, finalno, beginno=beginno, singlesnap=0, firever=firever,halonostr=halostr, hubble=h0, comoving=0, maindir=maindir)
                    afactor=atime
                else:
                    halosA = SF.read_halo_history(rundir, maindir=maindir, halonostr=halostr, hubble=h0, comoving=0)
                    afactor=1.0
                redlist = halosA['redshift']
                haloid = halosA['ID']
                a_scale = 1.0/(1.0+redlist)
                xcenl = halosA['x']*afactor
                ycenl = halosA['y']*afactor
                zcenl = halosA['z']*afactor
                Rvirl = halosA['R']*afactor
                Mstarl = halosA['Ms']
                xcen = np.interp(atime,a_scale,xcenl)
                ycen = np.interp(atime,a_scale,ycenl)
                zcen = np.interp(atime,a_scale,zcenl)
                Rvir = np.interp(atime,a_scale,Rvirl)
                Mstar = np.interp(atime,a_scale,Mstarl)
                Sxrel = Sx-xcen
                Syrel = Sy-ycen
                Szrel = Sz-zcen
                Sr = np.sqrt(Sxrel*Sxrel+Syrel*Syrel+Szrel*Szrel)
                cutrvs = Sr<Rvir*Rfrac
                Smi = Smi[cutrvs]
                Sage = Sage[cutrvs]
        print 'np.amax(Sage)', np.amax(Sage)
        print 'np.amin(Sage)', np.amin(Sage)
        if uselog==1:
            radlist = np.logspace(np.log10(minr),np.log10(withinr),num=30)
        else:
            radlist = np.linspace(minr,withinr,num=30)
                wdisk=radlist[1]-radlist[0]
                hdisk = maxlength
        agecutlow = np.amax(Sage)-0.001*ageneed
        agecuthigh = np.amax(Sage)
        Smdenlist = []
        Smlist = []
        cuta = (Sage> agecutlow) & (Sage <agecuthigh)
        Sx=Sx[cuta]; Sy=Sy[cuta]; Sz=Sz[cuta]; Sm=Sm[cuta];
        if usenscen==1 and cosmo==0:
            from crtestfunction import findcenz
            datasup=1
            xcen = findcenz(runtodo,Nsnap,withinr=withinr,dir='x',datasup=datasup,Gx=Sx,Gy=Sy,Gz=Sz,Gm=Sm)
            ycen = findcenz(runtodo,Nsnap,withinr=withinr,dir='y',datasup=datasup,Gx=Sx,Gy=Sy,Gz=Sz,Gm=Sm)
            zcen = findcenz(runtodo,Nsnap,withinr=withinr,dir='z',datasup=datasup,Gx=Sx,Gy=Sy,Gz=Sz,Gm=Sm)
            Sx=Sx-xcen;Sy=Sy-ycen;Sz=Sz-zcen;
                for radius in radlist:
            if usesphrad==1:
                cutxyz = (Sx*Sx+Sy*Sy+Sz*Sz < radius*radius)
            else:
                cutxy = (Sx*Sx+Sy*Sy < radius*radius)
                cutz = np.absolute(Sz) < hdisk/2.
                cutxyz =cutxy*cutz  
            Smlist = np.append(Smlist,np.sum(Sm[cutxyz])*1e10)
        Smlist = np.array(Smlist)
        if usesphrad==1:
            vollist = 4./3.*np.pi*(radlist[1:]*radlist[1:]*radlist[1:]-radlist[:-1]*radlist[:-1]*radlist[:-1])
        else:
            vollist = np.pi*(radlist[1:]*radlist[1:]-radlist[:-1]*radlist[:-1])*hdisk
        Smdenlist = (Smlist[1:]-Smlist[:-1])/vollist
        radhlist = (radlist[1:]+radlist[:-1])/2.
        sfrden = Smdenlist/(agecuthigh-agecutlow)/1.0e9
                #plt.plot(radlist,sfrden,label='kappa='+dclabel+'; '+'time step= '+timestep)
        if haveB==0:
            lsn='solid'
        else:
            lsn='dashed'
        if wanted=='sfrrad':
            plt.plot(radhlist,sfrden,label=dclabel,color=color,ls=lsn,lw=2)
        if wanted=='nsrad':
            plt.plot(radhlist,Smdenlist,label=dclabel,color=color,ls=lsn,lw=2)
        if wanted=='cumnsmrad':
            if normalized==1:
                plt.plot(radlist,Smlist/Smlist[-1],label=dclabel,color=color,ls=lsn,lw=2)
            else:
                                plt.plot(radlist,Smlist,label=dclabel,color=color,ls=lsn,lw=2)
        print 'radlist', radlist
        print 'Smlist', Smlist
        print 'halfrad', np.power(10., np.interp(Smlist[-1]/2.,Smlist,np.log10(radlist)))
        plt.xlabel(r'${\rm r\;[kpc]}$',fontsize=18)
        #plt.xlim(xmax=np.amax(radlist))
    if uselog==1:
        plt.xscale('log')
        #plt.title(runtitle+' SFR density at midplane')
    if wanted == 'sfrrad':
        plt.ylabel(r'$\rho_{\rm SFR}({\rm M_{\odot}/yr/kpc^3})$',fontsize=18)
    if wanted == 'nsrad':
        plt.ylabel(r'$\rho_{\rm *,new}({\rm M_{\odot}/kpc^3})$',fontsize=18)
        if runtitle=='SMC':
            plt.legend(loc=1,fontsize=10,ncol=3)
        if wanted == 'cumnsmrad':
        if normalized==1:
            plt.ylabel(r'$M_{\rm *,new}(<r)/M_{\rm *,new}(r_{\rm last})$',fontsize=18)
        else:
            plt.ylabel(r'$M_{\rm *,new}(<r)[{\rm M_{\odot}}]$',fontsize=18)
                if runtitle=='SMC':
                        plt.legend(loc=1,fontsize=10,ncol=3)
        plt.yscale('log')
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.title(ptitle,fontsize=18)
    plt.tight_layout()
    if wanted == 'sfrrad':
        if usesphrad==1:    
                        figname='CRplot/sfrrad/SFR_sph_'+fmeat+'.pdf'
        else:
            figname='CRplot/sfrrad/SFR_midplane_'+fmeat+'.pdf'
    if wanted == 'nsrad':
                if usesphrad==1:
                        figname='CRplot/nsrad/newstar_sph_'+fmeat+'.pdf'
        else:
            figname='CRplot/nsrad/newstar_midplane_'+fmeat+'.pdf'
        if wanted == 'cumnsmrad':
        if normalized==1:
            if usesphrad==1:
                figname='CRplot/cumnsmrad/newstar_sph_'+fmeat+'_normalized.pdf'
            else:
                figname='CRplot/cumnsmrad/newstar_midplane_'+fmeat+'_normalized.pdf'
        else:
                        if usesphrad==1:
                                figname='CRplot/cumnsmrad/newstar_sph_'+fmeat+'.pdf'
                        else:
                                figname='CRplot/cumnsmrad/newstar_midplane_'+fmeat+'.pdf'
    print 'figname', figname
    plt.savefig(figname)
    plt.clf()


if wanted=='sfrv':
        for runtodo in dirneed:
                info=outdirname(runtodo, Nsnap)
                rundir=info['rundir']
                runtitle=info['runtitle']
                slabel=info['slabel']
                snlabel=info['snlabel']
                dclabel=info['dclabel']
                resolabel=info['resolabel']
                the_snapdir=info['the_snapdir']
                Nsnapstring=info['Nsnapstring']
                havecr=info['havecr']
                Fcal=info['Fcal']
                iavesfr=info['iavesfr']
                timestep=info['timestep']
                S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix)
                Sage = S['age']
                Spos = S['p']
                Sx = Spos[:,0]
                Sy = Spos[:,1]
                Sz = Spos[:,2]
                Sm = S['m']
                #print 'Sage', Sage, np.amax(Sage), np.amin(Sage)
                rdisk=3.
                height=2.
        dz=height/nogrid
                agecut = 0.2 #in Gyr
                Smdenlist = []
        zl =[]
        for iv in range(nogrid):
            cutxy = (Sx*Sx+Sy*Sy < rdisk*rdisk)
            cutz = (Sz*Sz<(iv+1.0)*dz*(iv+1.0)*dz)&(Sz*Sz > iv*dz*iv*dz)
                        cuta = Sage< agecut
                        cut = cutxy*cutz*cuta
            vol_in_kpc3 = np.pi*(rdisk*rdisk*dz*2)
                        Smdencut = np.sum(Sm[cut])*1e10/vol_in_kpc3 #in Msun/kpc^3
                        Smdenlist=np.append(Smdenlist,Smdencut)
            zl = np.append(zl,(iv+0.5)*dz)
                Smdenlist = np.array(Smdenlist)
                sfrden = Smdenlist/agecut/1.0e9
                plt.plot(zl,sfrden,label=runtodo)
                print 'sfrden', sfrden
        plt.xlabel('z (kpc)')
        plt.xlim(xmax=np.amax(zl))
        plt.title(runtitle+' SFR density at within = '+str(rdisk)+'kpc')
        plt.ylabel(r'$\rho_{\rm SFR}({\rm M_{\odot}/yr/kpc^3})$')
        plt.yscale('log')
        plt.legend(loc='best')
        plt.savefig('CRplot/SFR_vertical_'+fmeat+'.pdf')


if wanted=='sfrarad':
        for runtodo in dirneed:
                rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr, timestep=outdirname(runtodo, Nsnap)
                S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix)
                Sage = S['age']
                Spos = S['p']
                Sx = Spos[:,0]
                Sy = Spos[:,1]
                Sz = Spos[:,2]
                Sm = S['m']
                #print 'Sage', Sage, np.amax(Sage), np.amin(Sage)
                radlist = np.linspace(1.1,20,num=10)
                wdisk=2.
                hdisk = 1.
                agecut = 0.2 #in Gyr
                Smlist = []
                for radius in radlist:
                        insol = radius+wdisk/2.
                        outsol = radius-wdisk/2.
                        cutxy = (Sx*Sx+Sy*Sy < insol*insol)
                        cutz = np.absolute(Sz) < hdisk/2.
                        cuta = Sage< agecut
                        cut = cutxy*cutz*cuta
                        Smcut = np.sum(Sm[cut])*1e10 #in Msun
                        Smlist=np.append(Smlist,Smcut)
                Smlist = np.array(Smlist)
                sfrcum = Smlist/agecut/1.0e9 #in Msun/yr
                #plt.plot(radlist,sfrden,label='kappa='+dclabel+'; '+'time step= '+timestep)
                plt.plot(radlist,sfrcum,label=runtodo)
                print 'sfrden', sfrcum
        plt.xlabel('r (kpc)')
        plt.xlim(xmax=np.amax(radlist))
        plt.title(runtitle+' cumulative SFR at midplane')
        plt.ylabel(r'$M_{\rm *,new}({\rm M_{\odot}/yr})$')
        #plt.yscale('log')
        plt.legend(loc='best')
        plt.savefig('CRplot/SFR_cum_midplane_'+fmeat+'.pdf')
        plt.clf()


if wanted=='crcpsnap':
        for runtodo in dirneed:
                snaplist=[]
                enclist=[]
                englist=[]
                enllist=[]
                endlist=[]
                enplist=[]
                avesfrl=[]
                prel=0
                preg=0
                prec=0
                pred=0
                prep=0
                presm = 0
                presnap = 0
                for i in range(startno,Nsnap, snapsep):
                        rundir, runtitle, slabel, snlabel, dclabel, resolabel, the_snapdir, Nsnapstring, havecr, Fcal, iavesfr, timestep=outdirname(runtodo, i)
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        cregyl = G['cregyl']
                        cregyg = G['cregyg']
                        cregy  = G['cregy']
                        cregyd = G['cregyd']
                        if havecr>4:
                                cregyp = G['cregyp']
                        try:
                                Sm = np.sum(S['m'])
                        except KeyError:
                                Sm = 0.
                        cregygt=np.sum(cregyg)
                        cregyt =np.sum(cregy)
                        cregylt=np.sum(cregyl)
                        cregydt=np.sum(cregyd)
                        if havecr>4:
                                cregydtp=np.sum(cregyp)
                        eng = (cregygt-preg)/1e6/3.2e7*2e53/float(snapsep)
                        enc = (cregyt-prec)/1e6/3.2e7*2e53/float(snapsep)
                        enl = (cregylt-prel)*2e53/1e6/3.2e7/float(snapsep)
                        end = (cregydt-pred)*2e53/1e6/3.2e7/float(snapsep)
                        if havecr>4:
                                enp = (cregydtp-prep)*2e53/1e6/3.2e7/float(snapsep)
                                print 'CR energy dtp', enp
                        print 'CR energy loss rate (erg/s)', enl
                        print 'CR energy gain rate (erg/s)', eng
                        print 'CR energy change rate (erg/s)', enc
                        print 'CR energy dt rate (erg/s)', end  #including adiabatic heating and streaming
                        snaplist.append(i)
                        enclist.append(enc)
                        englist.append(eng)
                        enllist.append(enl)
                        endlist.append(end)
                        if havecr>4:
                                enplist.append(enp)
                        preg = cregygt
                        prec = cregyt
                        prel = cregylt
                        pred = cregydt
                        if havecr>4:
                                prep = cregydtp
                        if i>presnap:
                                avesfr=(Sm-presm)*1e4/0.98/(i-presnap)
                        else:
                                avesfr=0
                        presnap = i
                        presm = Sm
                        print 'avesfr', avesfr
                        avesfrl.append(avesfr)
                avesfrl=np.array(avesfrl)
                enclist=np.array(enclist)
                englist=np.array(englist)
                endlist=np.array(endlist)
                enllist=np.array(enllist)
                enplist=np.array(enplist)
                plt.plot(snaplist[1:], enclist[1:], label='Change')
                plt.plot(snaplist[1:], englist[1:], label='SNe')
                if havecr > 4:
                        plt.plot(snaplist[1:], endlist[1:], label='Ad')
                #plt.plot(snaplist, endlist, label='CR energy dt')
                plt.plot(snaplist[1:], enllist[1:], label='Cool')
                #plt.plot(snaplist, englist+endlist+enllist, label='CR energy estimate')
                if havecr>4:
                        plt.plot(snaplist[1:], englist[1:]+endlist[1:]+enllist[1:], label='SNe+Ad+Cool')
                #plt.yscale('log')
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                plt.xlabel('Myr', fontsize=25)
                plt.subplots_adjust(left=0.2,bottom=0.2, right=0.75)
                plt.ylabel('dE/dt (erg/s)', fontsize=25)
                plt.savefig('CRplot/cresnap_'+fmeat+'.pdf')
                plt.clf()
                #Lsfr = 3.8e-4 *avesfr*2e33/3.2e7*9e20
                #above is the coefficient for Salpeter only; for Kroupa, the coefficient is 50% larger:
                Lsfr = 3.8e-4 *1.5*avesfrl*2e33/3.2e7*9e20
                Lsfrideal = 3.8e-4 * 1.5 * 1.0 * 2e33 / 3.2e7 *9e20
                if havecr >4:
                        Lgammae = (enllist+endlist)/7.51e-16/2e5/3.2e7/250.0*0.7/3.0
                Lgamma = (enllist)/7.51e-16/2e5/3.2e7/250.0*0.7/3.0
                #Lgammae_sfr = Lgammae/Lsfr
                Lgamma_sfr = Lgamma/Lsfr
                plt.plot(snaplist[1:], np.absolute(Lgamma_sfr[1:]), label=r'Cool')
                #plt.plot(snaplist[1:], np.absolute(Lgammae_sfr[1:]), label='Cool+ad')
                #plt.plot(snaplist[1:], np.absolute(Lgamma[1:]/Lsfrideal), label='SFR=1')
                plt.yscale('log')
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                plt.xlabel('Myr', fontsize=25)
                plt.ylabel(r'$\frac{L_{\gamma}}{L_{\rm SF}}$', fontsize=30)
                plt.savefig('CRplot/gammasfrsnap_'+fmeat+'.pdf')
                plt.clf()

if wanted=='crdrad':
        for runtodo in dirneed:
                snaplist=[]
                enclist=[]
                englist=[]
                enllist=[]
                endlist=[]
                enplist=[]
                enalist=[]
                avesfrl=[]
                prel=0
                preg=0
                prec=0
                pred=0
                prep=0
                presm = 0
                presnap = []
                crecuml = []
                crecumg = []
                crecuma = []
                crecumd = []
                crecum  = []
                crecump = []
                for i in [Nsnap]:
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gpos = G['p']
                        Gx = Gpos[:,0]
                        Gy = Gpos[:,1]
                        Gz = Gpos[:,2]
                        Grho = G['rho']
                        Gm = G['m']*1e10
                        Gr = np.sqrt(np.square(Gx)+np.square(Gy)+np.square(Gz))
                        cregyl = G['cregyl']/1e6/3.2e7*2e53
                        cregyg = G['cregyg']/1e6/3.2e7*2e53
                        cregy  = G['cregy']/1e6/3.2e7*2e53
                        cregyd = G['cregyd']/1e6/3.2e7*2e53
                        if havecr>4:
                                cregyp = G['cregyp']/1e6/3.2e7*2e53
                        if havecr>5:
                                cregya = G['cregya']/1e6/3.2e7*2e53
                        rad=np.linspace(0.1,withinr, num=nogrid)
                        for i in range(len(rad)-1):
                cutu=Gr<rad[i+1]
                cutd=Gr>rad[i]
                cutz=np.absolute(Gz)<maxlength
                cut=cutu*cutd*cutz
                vol = 2.0*maxlength*np.pi*(np.power(rad[i+1],2)-np.power(rad[i],2))
                                crecuml = np.append(crecuml,np.sum(cregyl[cut])/vol)
                                crecumg = np.append(crecumg, np.sum(cregyg[cut])/vol)
                                crecuma = np.append(crecuma, np.sum(cregya[cut])/vol)
                                crecum = np.append(crecum, np.sum(cregy[cut])/vol)
                                crecumd = np.append(crecumd, np.sum(cregyd[cut])/vol)
                                crecump = np.append(crecump, np.sum(cregyp[cut])/vol)
                        crecumad = crecumd + crecumg
                        plt.plot(rad[1:], crecum,label='CR energy')
                        plt.plot(rad[1:], crecumg, label='SNe')
                        plt.plot(rad[1:], crecuml, label='Loss')
                        plt.plot(rad[1:], crecumd-crecuml-crecump, label='Ad')
                        plt.plot(rad[1:], crecuma, label='test Ad')
                        plt.plot(rad[1:], crecump, label='Flux')
                        plt.plot(rad[1:], crecumad, label='SNe+Loss+Ad+Flux')
                        plt.legend()
            plt.title('CR density within the disk')
                        plt.ylabel(r'$\left \langle \rm{de/dt}  \right \rangle {\rm (erg/s/kpc^3)}$')
                        plt.xlabel('r (kpc)')
                        plt.savefig('CRplot/'+runtodo+'_crdrad.pdf')
                        plt.clf()


if wanted=='cratime':
        for runtodo in dirneed:
                snaplist=[]
                enclist=[]
                englist=[]
                enllist=[]
                endlist=[]
                enplist=[]
                enalist=[]
                avesfrl=[]
                prel=0
                preg=0
                prec=0
                pred=0
                prep=0
        prea=0
                presm = 0
                presnap = []
                crecuml = []
                crecumg = []
                crecuma = []
                crecumd = []
                crecum  = []
                crecump = []
        timel = []
                for i in range(startno,Nsnap,snapsep):
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
            cosmo=info['cosmo']
            exceptcool=info['exceptcool']
            havemetal=info['havemetal']
            if havecr==0:
                continue
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
            print 'exceptcool', exceptcool
            try:
                if cosmo==1:
                    G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr, cosmological=1,h0=1,exceptcool=exceptcool,havemetal=havemetal)
                    header=G['header']
                    timeneed=header[2]
                else:
                    G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,exceptcool=exceptcool, havemetal=havemetal)
                Gpos = G['p']
                Gx = Gpos[:,0]
                Gy = Gpos[:,1]
                Gz = Gpos[:,2]
                Gv = G['v']
                Gvx = Gv[:,0]
                Gvy = Gv[:,1]
                Gvz = Gv[:,2]
                Grho = G['rho']
                Gm = G['m']*1e10
                Gu = G['u']
                Gr = np.sqrt(np.square(Gx)+np.square(Gy)+np.square(Gz))
                cutr = Gr<1
                print 'Gucutr', np.sum(Gu[cutr])
                print 'Gvxcutr', np.sum(Gvx[cutr])
                print 'Gmcutr', np.sum(Gm[cutr])
                cregyl = G['cregyl']*2e53
                cregyg = G['cregyg']*2e53
                cregy  = G['cregy']*2e53
                cregyd = G['cregyd']*2e53
                print 'cregycutr', np.sum(cregy[cutr])
            except KeyError:
                continue
                        if havecr>4:
                                cregyp = G['cregyp']*2e53
                        if havecr>5:
                                cregya = G['cregya']*2e53
                        crecum = np.append(crecum, np.sum(cregy))
            #print 'ad min max', np.amin(cregya), np.amax(cregya)
            crecuml = np.append(crecuml,np.sum(cregyl))
            crecumg = np.append(crecumg, np.sum(cregyg))
            crecuma = np.append(crecuma, np.sum(cregya))
            crecumd = np.append(crecumd, np.sum(cregyd))
            crecump = np.append(crecump, np.sum(cregyp))
            print 'cregy', np.sum(cregy)
            print 'cregyp', np.sum(cregyp)
            if cosmo==1:
                readtimelist=readtime(firever=2)
                snap2list=readtimelist['snaplist']
                time2list=readtimelist['timelist']
                a2list=readtimelist['alist']
                timenow = np.interp(timeneed,a2list,time2list)*1e3
            else:
                timenow=i*0.98
            timel = np.append(timel, timenow)
            prec=np.sum(cregy)
            prel=np.sum(cregyl)
            preg=np.sum(cregyg)
            prea=np.sum(cregya)
            pred=np.sum(cregyd)
            prep=np.sum(cregyp)
        if useM1==1:
            crecumad = crecumd + crecumg + crecuml
        else:
            crecumad = crecumd + crecumg
        if needlog == 1:
            plt.plot(timel, np.log10(np.absolute(crecum)), marker='s',label='CR energy')
            plt.plot(timel, np.log10(np.absolute(crecumg)), marker='s',label='SNe')
            plt.plot(timel, np.log10(np.absolute(crecuml)), marker='s',label='Loss')
            plt.plot(timel, np.log10(np.absolute(crecuma)),marker='s', label='Ad')
            plt.plot(timel, np.log10(np.absolute(crecump)), marker='s',label='Flux')
            if nosum==0:
                plt.plot(timel, np.log10(np.absolute(crecumad)), marker='s',label='Sum')
        else:
            plt.plot(timel, crecum, marker='s',label='CR energy')
            plt.plot(timel, crecumg, marker='s',label='SNe')
            plt.plot(timel, crecuml, marker='s',label='Loss')
            plt.plot(timel, crecuma,marker='s', label='Ad')
            plt.plot(timel, crecump, marker='s',label='Flux')
            if nosum==0:
                plt.plot(timel, crecumad, marker='s',label='Sum')
        plt.legend(loc='best')
        if needlog == 1:
                        plt.ylabel(r'${\rm log}\left \langle \rm{E}  \right \rangle$ (erg)')
        else:
            plt.ylabel(r'$\left \langle \rm{E}  \right \rangle$ (erg)')
        plt.xlabel('t (Myr)')
        #plt.ylim([-3.0e43,4.0e43])
        plt.title('diffusion coefficient = ' +dclabel) 
        plt.savefig('CRplot/cratime/'+runtodo+'_time.pdf')
        plt.clf()

if wanted=='crdtime':
        for runtodo in dirneed:
                snaplist=[]
                enclist=[]
                englist=[]
                enllist=[]
                endlist=[]
                enplist=[]
                enalist=[]
                avesfrl=[]
        prea=0
                prel=0
                preg=0
                prec=0
                pred=0
                prep=0
                presm = 0
                presnap = []
                crecuml = []
                crecumg = []
                crecuma = []
                crecumd = []
                crecum  = []
                crecump = []
                timel = []
                for i in range(0,Nsnap,snapsep):
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gpos = G['p']
                        Gx = Gpos[:,0]
                        Gy = Gpos[:,1]
                        Gz = Gpos[:,2]
                        Grho = G['rho']
                        Gm = G['m']*1e10
                        Gr = np.sqrt(np.square(Gx)+np.square(Gy)+np.square(Gz))
                        cregyl = G['cregyl']/1e6/3.2e7*2e53
                        cregyg = G['cregyg']/1e6/3.2e7*2e53
                        cregy  = G['cregy']/1e6/3.2e7*2e53
                        cregyd = G['cregyd']/1e6/3.2e7*2e53
                        if havecr>4:
                                cregyp = G['cregyp']/1e6/3.2e7*2e53
                        if havecr>5:
                                cregya = G['cregya']/1e6/3.2e7*2e53

                        crecuml = np.append(crecuml,(np.sum(cregyl)-prel)/(i*0.98*snapsep))
                        crecumg = np.append(crecumg, (np.sum(cregyg)-preg)/(i*0.98*snapsep))
                        crecuma = np.append(crecuma, (np.sum(cregya)-prea)/(i*0.98*snapsep))
                        crecum = np.append(crecum, (np.sum(cregy)-prec)/(i*0.98*snapsep))
                        crecumd = np.append(crecumd, (np.sum(cregyd)-pred)/(i*0.98*snapsep))
                        crecump = np.append(crecump, (np.sum(cregyp)-prep)/(i*0.98*snapsep))
                        timel = np.append(timel, i*0.98)
            prel=np.sum(cregyl)
            preg=np.sum(cregyg)
            prea=np.sum(cregya)
            prec=np.sum(cregy)
            pred=np.sum(cregyd)
            prep=np.sum(cregyp)
                if useM1==1:
                        crecumad = crecumd + crecumg + crecuml
                else:
                        crecumad = crecumd + crecumg
                plt.plot(timel, crecum,label='CR energy')
                plt.plot(timel, crecumg, label='SNe')
                plt.plot(timel, crecuml, label='Loss')
                plt.plot(timel, crecumd-crecuml-crecump, label='Ad')
                #plt.plot(rad, crecuma, label='test Ad')
                plt.plot(timel, crecump, label='Flux')
                plt.plot(timel, crecumad, label='SNe+Loss+Ad+Flux')
                plt.legend()
                plt.ylabel(r'$\left \langle \rm{dE/dt}  \right \rangle$ (erg/s)')
                plt.xlabel('t (Myr)')
                plt.savefig('CRplot/'+runtodo+'_dtime.pdf')




if wanted=='gasdenv':
        for runtodo in dirneed:
                for i in [Nsnap]:
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gp = G['p']
                        #Grho = G['rho']
                        Gu = G['u']
                        Gm = G['m']
                        cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                        Neb = G['ne']
                        #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                        Gz = Gp[:,2]
                        Gx = Gp[:,0]
                        Gy = Gp[:,1]
            withinz=5.0
            rdisk=3.0
                        dz = withinz/nogrid
                        Gnism_in_cm_3l=[]
                        zl =[]
                        for iv in range(nogrid):
                                cutxy = (Gx*Gx+Gy*Gy < rdisk*rdisk)
                                cutz = (Gz*Gz<(iv+1.0)*dz*(iv+1.0)*dz)&(Gz*Gz > iv*dz*iv*dz)
                                cut=cutxy*cutz
                                Nebcut = Neb[cut]
                                Gmcut = Gm[cut]
                                Gm_in_g = Gmcut*1e10*2e33
                                vol_in_cm3 = np.pi*(rdisk*rdisk*3.086e21*3.086e21*3.086e21*dz*2)
                                Grho_in_g_cm_3 = Gm_in_g/vol_in_cm3
                                Gnism_in_cm_3 = np.sum((0.78+0.22*Nebcut*0.76)/protonmass_in_g*Grho_in_g_cm_3)
                                Gnism_in_cm_3l = np.append(Gnism_in_cm_3l, Gnism_in_cm_3)
                                zl = np.append(zl, dz*(iv+0.5))
                        plt.plot(zl, Gnism_in_cm_3l, label=runtodo)
        plt.xlabel('z [kpc]')
        plt.ylabel(r'$n_{\rm ISM} [{\rm cm^{-3}}]$')
        plt.legend(loc='best')
        plt.yscale('log')
        plt.title('vertical gas density within radius = '+str(rdisk)+'kpc')
        plt.savefig('CRplot/gasdensity_'+fmeat+'_v.pdf')
        plt.clf()



if wanted=='lgratio':
        for runtodo in dirneed:
                snaplist=[]
                enclist=[]
                englist=[]
                enllist=[]
                endlist=[]
                enplist=[]
                enalist=[]
                avesfrl=[]
                prel=0
                preg=0
                prec=0
                pred=0
                prep=0
                presm = 0
                presnap = []
                crecuml = []
                crecumg = []
                crecuma = []
                crecumd = []
                crecum  = []
                crecump = []
                timel = []
                for i in range(0,Nsnap,snapsep):
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gpos = G['p']
                        Gx = Gpos[:,0]
                        Gy = Gpos[:,1]
                        Gz = Gpos[:,2]
                        Grho = G['rho']
                        Gm = G['m']*1e10
                        Gr = np.sqrt(np.square(Gx)+np.square(Gy)+np.square(Gz))
                        cregyl = G['cregyl']/1e6/3.2e7*2e53
                        cregyg = G['cregyg']/1e6/3.2e7*2e53
                        cregy  = G['cregy']/1e6/3.2e7*2e53
                        cregyd = G['cregyd']/1e6/3.2e7*2e53
                        if havecr>4:
                                cregyp = G['cregyp']/1e6/3.2e7*2e53
                        if havecr>5:
                                cregya = G['cregya']/1e6/3.2e7*2e53
                        crecuml = np.append(crecuml,np.sum(cregyl))
                        crecumg = np.append(crecumg, np.sum(cregyg))
                        crecuma = np.append(crecuma, np.sum(cregya))
                        crecum = np.append(crecum, np.sum(cregy))
                        crecumd = np.append(crecumd, np.sum(cregyd))
                        crecump = np.append(crecump, np.sum(cregyp))
                        timel = np.append(timel, i*0.98)
        crecumg=np.array(crecumg)
        crecuml=np.array(crecuml)
                plt.plot(timel, crecuml/crecumg, label=runtodo)
                plt.legend(bbox_to_anchor=(1.1, 1.05))
                plt.ylabel(r'$\left \langle \rm{E_{\bf Loss}/E_{\bf SNe}}  \right \rangle$')
                plt.xlabel('t (Myr)')
                plt.title('diffusion coefficient = ' +dclabel)
    plt.savefig('CRplot/lgratio_'+fmeat+'.pdf')


if wanted=='dlgratio':
        for runtodo in dirneed:
                snaplist=[]
                enclist=[]
                englist=[]
                enllist=[]
                endlist=[]
                enplist=[]
                enalist=[]
                avesfrl=[]
                prel=0
                preg=0
                prec=0
                pred=0
                prep=0
                presm = 0
                presnap = []
                crecuml = []
                crecumg = []
                crecuma = []
                crecumd = []
                crecum  = []
                crecump = []
                timel = []
                for i in range(0,Nsnap,snapsep):
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gpos = G['p']
                        Gx = Gpos[:,0]
                        Gy = Gpos[:,1]
                        Gz = Gpos[:,2]
                        Grho = G['rho']
                        Gm = G['m']*1e10
                        Gr = np.sqrt(np.square(Gx)+np.square(Gy)+np.square(Gz))
                        cregyl = G['cregyl']/1e6/3.2e7*2e53
                        cregyg = G['cregyg']/1e6/3.2e7*2e53
                        cregy  = G['cregy']/1e6/3.2e7*2e53
                        cregyd = G['cregyd']/1e6/3.2e7*2e53
                        if havecr>4:
                                cregyp = G['cregyp']/1e6/3.2e7*2e53
                        if havecr>5:
                                cregya = G['cregya']/1e6/3.2e7*2e53
                        crecuml = np.append(crecuml,np.sum(cregyl))
                        crecumg = np.append(crecumg, np.sum(cregyg))
                        crecuma = np.append(crecuma, np.sum(cregya))
                        crecum = np.append(crecum, np.sum(cregy))
                        crecumd = np.append(crecumd, np.sum(cregyd))
                        crecump = np.append(crecump, np.sum(cregyp))
                        timel = np.append(timel, i*0.98)
                crecumg=np.array(crecumg)
                crecuml=np.array(crecuml)
        plt.ylim([0.01,1.5])
                plt.plot(timel[1:], (-crecuml[1:]+crecuml[:-1])/(crecumg[1:]-crecumg[:-1]), label=runtodo)
                plt.legend(loc='best',fontsize=16)
        plt.yscale('log')
                plt.ylabel(r'$\left \langle \rm{\Delta E_{Loss}/\Delta E_{SNe}}  \right \rangle$')
                plt.xlabel('t (Myr)')
                plt.title('diffusion coefficient = ' +dclabel)
        plt.savefig('CRplot/dlgratio_'+fmeat+'.pdf')

if wanted=='gsmratio':
        for runtodo in dirneed:
                snaplist=[]
                enclist=[]
                englist=[]
                enllist=[]
                endlist=[]
                enplist=[]
                enalist=[]
                avesfrl=[]
                prel=0
                preg=0
                prec=0
                pred=0
                prep=0
                presm = 0
                presnap = []
                crecuml = []
                crecumg = []
                crecuma = []
                crecumd = []
                crecum  = []
                crecump = []
                timel = []
        smen = []
                for i in range(startno,Nsnap,snapsep):
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        G = readsnap(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gpos = G['p']
                        Gx = Gpos[:,0]
                        Gy = Gpos[:,1]
                        Gz = Gpos[:,2]
                        Grho = G['rho']
                        Gm = G['m']*1e10
            SMtot = np.sum(S['m'])*1e10
            Energysm = SMtot*1.989e30*3.0e8*3.0e8/1e6/3.2e7*1e7*6.2e-4
                        Gr = np.sqrt(np.square(Gx)+np.square(Gy)+np.square(Gz))
                        cregyl = G['cregyl']/1e6/3.2e7*2e53
                        cregyg = G['cregyg']/1e6/3.2e7*2e53
                        cregy  = G['cregy']/1e6/3.2e7*2e53
                        cregyd = G['cregyd']/1e6/3.2e7*2e53
                        if havecr>4:
                                cregyp = G['cregyp']/1e6/3.2e7*2e53
                        if havecr>5:
                                cregya = G['cregya']/1e6/3.2e7*2e53
                        crecuml = np.append(crecuml,np.sum(cregyl))
                        crecumg = np.append(crecumg, np.sum(cregyg))
                        crecuma = np.append(crecuma, np.sum(cregya))
                        crecum = np.append(crecum, np.sum(cregy))
                        crecumd = np.append(crecumd, np.sum(cregyd))
                        crecump = np.append(crecump, np.sum(cregyp))
            smen = np.append(smen, Energysm)
                        timel = np.append(timel, i*0.98)
                crecumg=np.array(crecumg)
                smen=np.array(smen)
                plt.plot(timel, crecumg/smen, label=runtodo)
                plt.legend(bbox_to_anchor=(1.1, 1.05))
                plt.ylabel(r'$\left \langle \rm{E_{SNe}/E_{SFR_to_SNe}}  \right \rangle$')
                plt.xlabel('t (Myr)')
                plt.title('diffusion coefficient = ' +dclabel)
        plt.savefig('CRplot/gsmratio_'+fmeat+'.pdf')



if wanted=='dirage':
        for runtodo in dirneed:
                snaplist=[]
                avesfrl=[]
                presm = 0
                presnap = 0
                for i in [startno]:
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        cosmo=info['cosmo']
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        if cosmo==1:
                                S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                        else:
                                S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
            print 'Sage', S['age']
            print 'mass ave age', np.sum(S['age']*S['m'])/np.sum(S['m'])



if wanted=='dg_sfr':
        for runtodo in dirneed:
                snaplist=[]
                enclist=[]
                englist=[]
                enllist=[]
                endlist=[]
                enplist=[]
                enalist=[]
                sml=[]
        nsml=[]
                prel=0
                preg=0
                prec=0
                pred=0
                prep=0
                presm = 0
                presnap = 0
        pretime = 0
                crecuml = []
                crecumg = []
                crecuma = []
                crecumd = []
                crecum  = []
                crecump = []
                timel = []
                for i in range(0,Nsnap,snapsep):
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
            cosmo=info['cosmo']
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        if cosmo==1:
                                S = readsnapcr(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                        else:
                                S = readsnapcr(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        try:
                                Sm = np.sum(S['m'])
                                Smi = S['m']
                                Sage = S['age']
                                Sm = np.sum(Smi)
                                header=S['header']
                                timeneed=header[2]
                                tcut=Sage>timeneed-0.001
                                nsm = np.sum(Smi[tcut])
                        except KeyError:
                                Sm = 0.
                timeneed=0.
                nsm=0.
                        snaplist.append(i)
            if cosmo==1:
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
            else:
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gpos = G['p']
                        Gx = Gpos[:,0]
                        Gy = Gpos[:,1]
                        Gz = Gpos[:,2]
                        Grho = G['rho']
                        Gm = G['m']*1e10
                        Gr = np.sqrt(np.square(Gx)+np.square(Gy)+np.square(Gz))
                        cregyl = G['cregyl']*2e53
                        cregyg = G['cregyg']*2e53
                        cregy  = G['cregy']*2e53
                        cregyd = G['cregyd']*2e53
                        if havecr>4:
                                cregyp = G['cregyp']*2e53
                        if havecr>5:
                                cregya = G['cregya']*2e53
                        crecuml = np.append(crecuml,np.sum(cregyl))
                        crecumg = np.append(crecumg, np.sum(cregyg))
                        crecuma = np.append(crecuma, np.sum(cregya))
                        crecum = np.append(crecum, np.sum(cregy))
                        crecumd = np.append(crecumd, np.sum(cregyd))
                        crecump = np.append(crecump, np.sum(cregyp))
                        timel = np.append(timel, i*0.98)
            sml=np.append(sml,Sm)
            nsml=np.append(nsml,nsm)
            pretime=timeneed
                crecumg=np.array(crecumg)
                crecuml=np.array(crecuml)
        prefac = 1e10/0.4*0.0037*1.0e51
        snesfr = sml*1e10/0.4*0.0037*1.0e51
        dg_sfr=(crecumg[1:]-crecumg[:-1])/(nsml[1:]*prefac)/snapsep
        print 'dg_sfr', dg_sfr
                plt.plot(timel[1:], dg_sfr, label=runtodo)
        plt.yscale('log')
                plt.legend(loc='best')
                plt.ylabel(r'$\left \langle \rm{\Delta E_{SNe,CR}/\Delta E_{SNe,SFR}}  \right \rangle$')
                plt.xlabel('t (Myr)')
                plt.title('diffusion coefficient = ' +dclabel)
        plt.savefig('CRplot/dg_sfr_'+fmeat+'.pdf')


if wanted=='gamma_aveCR':
        for runtodo in dirneed:
        Lgammal=[]
        Lgammaavel=[]
        Lgammapl=[]
        enllist=[]
        timel=[]
        Gmcutwl=[]
        cregywl=[]
        Nebcutal=[]
                for i in range(startno,Nsnap,snapsep):
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        cosmo=info['cosmo']
                        G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gp = G['p']
                        Grho = G['rho']
                        Gu = G['u']
                        Gm = G['m']
                        cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                        Neb = G['ne']
            cregyl = G['cregyl']*1e10*solar_mass_in_g*km_in_cm*km_in_cm #in erg
                        #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                        Gz = Gp[:,2]
                        Gx = Gp[:,0]
                        Gy = Gp[:,1]
                        Gnism_in_cm_3 = (0.78+0.22*Neb*0.76)/proton_mass_in_g*Grho*1e10*solar_mass_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm
                        #tpi_in_s = pidecay_fac/Gnism_in_cm_3 #pi decay time in s
            tpi_in_s = 1.0/hadronicdecayrate/Gnism_in_cm_3 #pi decay time in s from Guo 2007
                        cregy_in_erg = cregy*solar_mass_in_g*1e10*km_in_cm*km_in_cm
            #tcut = tpi_in_s/betapi*nopi_per_gamma > 3e10
            #print 'total number of gas particles', np.count_nonzero(tpi_in_s)
            #print 'no of gas particles that have too much CR',np.count_nonzero(tpi_in_s[tcut])
            # here we account for the exponential in the code at each timestep. Here the time step is 3e10s
            #timestepcode=3e10
                        #Lgammagev = cregy_in_erg*(1.-np.exp(-timestepcode/tpi_in_s*betapi/nopi_per_gamma))/timestepcode #in erg
            Lgammagev = cregy_in_erg/tpi_in_s*betapi/nopi_per_gamma
            cutxy = Gx*Gx+Gy*Gy < withinr*withinr
            cutz = np.absolute(Gz)<maxlength
            cut = cutxy*cutz
            Lgammacut = Lgammagev[cut]
            gasvol = Gm[cut]/Grho[cut] #in kpc^3
            cylindervol = withinr*withinr*2.0*maxlength*np.pi #in kpc^3
            crden = np.sum(cregy_in_erg[cut])/cylindervol #in erg/kpc^3
            Lgammaave = crden*gasvol/tpi_in_s[cut]*betapi/nopi_per_gamma
            #print 'gas volume', gasvol
            #print 'cylinder volume', withinr*withinr*2.0*maxlength*np.pi
            Gmcutw = np.sum(Gm[cut])
            cregyw = np.sum(cregy[cut])
            try:
                Nebcuta = np.average(Neb[cut],weights=Gm[cut])
            except ZeroDivisionError:
                Nebcuta = 0
                        dr = withinr/nogrid
            Lgammap=0.0 
                        for irad in range(nogrid):
                                cutxy = (Gx*Gx+Gy*Gy > dr*irad*dr*irad) & (Gx*Gx+Gy*Gy < dr*(irad+1)*dr*(irad+1))
                                cutz = Gz*Gz < maxlength*maxlength
                                cutr=cutxy*cutz
                                Nebcut = Neb[cutr]
                                Gmcut = Gm[cutr]
                                try:
                                        Nebave = np.average(Nebcut,weights=Gmcut)
                                except ZeroDivisionError:
                                        Nebave = 0
                                Gm_in_g = Gmcut*1e10*2e33
                                shellvol_in_cm3 = np.pi*(-np.power(dr*irad,2)+np.power(dr*(irad+1),2))*kpc_in_cm*kpc_in_cm*kpc_in_cm*2.0*maxlength
                                Grho_in_g_cm_3 = np.sum(Gm_in_g)/shellvol_in_cm3
                                Gnism_in_cm_3p = (0.78+0.22*Nebave*0.76)/proton_mass_in_g*Grho_in_g_cm_3
                tpi_in_sp = pidecay_fac/Gnism_in_cm_3p #pi decay time in s
                cregycut = np.sum(cregy_in_erg[cutr])
                Lgammap += cregycut/tpi_in_sp*betapi/nopi_per_gamma #in erg/s
                #print 'Gnism_in_cm_3p', Gnism_in_cm_3p
            Lgammapl = np.append(Lgammapl,Lgammap)
            Lgammal=np.append(Lgammal,np.sum(Lgammacut))
            Lgammaavel=np.append(Lgammaavel,np.sum(Lgammaave))
            enllist=np.append(enllist,np.sum(cregyl[cut]))
            Gmcutwl=np.append(Gmcutwl, Gmcutw)
            cregywl=np.append(cregywl,cregyw)
            Nebcutal=np.append(Nebcutal,Nebcuta)
            timel=np.append(timel,float(i)*0.98*1e6)
        Lgamma_l = (enllist[1:]-enllist[:-1])/((timel[1:]-timel[:-1])*yr_in_sec)/(hadronicdecayrate+coulombdecayrate)*hadronicdecayrate*betapi/nopi_per_gamma
        #Consider CR density ~ 1eV/cm^3 and nsim ~ 1cm^-3 
        tpi_in_s_ideal=pidecay_fac
        vol_in_cm3 = np.pi*np.power(withinr,2)*kpc_in_cm*kpc_in_cm*kpc_in_cm*2.0*maxlength
        cregy_in_erg_ideal = 1.0*vol_in_cm3/erg_in_eV
        Lgamma_ideal = cregy_in_erg_ideal/tpi_in_s_ideal*betapi/nopi_per_gamma
        cregyaden_in_eV_cm_3 = cregywl*1e10*solar_mass_in_g*km_in_cm*km_in_cm*erg_in_eV/vol_in_cm3
        nism_in_cm_3 = Gmcutwl*1e10*2e33/vol_in_cm3*(0.78+0.22*Nebcutal*0.76)/proton_mass_in_g
        tpi_in_s_ave = pidecay_fac/nism_in_cm_3
        Lgammaave_nism_cre = cregywl*1e10*solar_mass_in_g*km_in_cm*km_in_cm/tpi_in_s_ave*betapi/nopi_per_gamma
                Lsfr = Kroupa_Lsf*1.0*solar_mass_in_g/yr_in_sec*cspeed_in_cm_s*cspeed_in_cm_s
        print 'Lgammal', Lgammal
        print 'Lgamma_l', Lgamma_l
        print 'Lgammapl', Lgammapl
        print 'Lgamma_ideal', Lgamma_ideal
        print 'nism_in_cm_3', nism_in_cm_3
        print 'cregyaden_in_eV_cm_3', cregyaden_in_eV_cm_3
        print 'Lgammaave_nism_cre', Lgammaave_nism_cre
        print 'Lsfr*2e-4', Lsfr*2e-4
        print 'hadronic fraction', hadronicdecayrate/(hadronicdecayrate+coulombdecayrate)
        plt.plot(timel[1:]/1e6, Lgammal[1:], label='from particles')
        plt.plot(timel[1:]/1e6, Lgammaavel[1:], label='from particles (ave CR)')
        plt.plot(timel[1:]/1e6, Lgammapl[1:], label='ave cylinders')
        plt.plot(timel[1:]/1e6, Lgammaave_nism_cre[1:], label='ave whole')
        plt.plot(timel[1:]/1e6, np.absolute(Lgamma_l), label='from code')
    #   plt.axhline(y=Lgamma_ideal, label='ideal',ls='dashed')  
        plt.title(runtodo)
                plt.legend(loc='best')
                plt.ylabel(r'$L_{\gamma} {\rm (erg/s)}$')
                plt.xlabel('t (Myr)')
        plt.savefig('CRplot/Lgamma_aveCR_'+runtodo+'.pdf')
        plt.clf()






if wanted=='gamma_partit':
        for runtodo in dirneed:
                Lgammal=[]
                Lgammapl=[]
                enllist=[]
                timel=[]
                Gmcutwl=[]
                cregywl=[]
                Nebcutal=[]
                for i in range(startno,Nsnap,snapsep):
            try:
                info=outdirname(runtodo, i)
                rundir=info['rundir']
                maindir=info['maindir']
                runtitle=info['runtitle']
                slabel=info['slabel']
                snlabel=info['snlabel']
                dclabel=info['dclabel']
                resolabel=info['resolabel']
                the_snapdir=info['the_snapdir']
                Nsnapstring=info['Nsnapstring']
                havecr=info['havecr']
                Fcal=info['Fcal']
                iavesfr=info['iavesfr']
                timestep=info['timestep']
                cosmo=info['cosmo']
                if cosmo==1:
                    h0=1
                else:
                    h0=0
                header = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo, header_only=1)
                ascale = header['time']
                print 'this time', ascale
                thisred = 1./ascale-1.
                hubble = header['hubble']
                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr, h0=h0,cosmological=cosmo)
                Gp = G['p']
                Grho = G['rho']
                Gu = G['u']
                Gm = G['m']
                cregy = G['cregy'] #cosmic ray energy in 1e10Msun km^2/sec^2
                #print 'np.sum(cregy)', np.sum(cregy)
                Neb = G['ne']
                cregyl = G['cregyl']*1e10*solar_mass_in_g*km_in_cm*km_in_cm #in erg
                #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                Gz = Gp[:,2]
                Gx = Gp[:,0]
                Gy = Gp[:,1]
                #print 'Gcm', np.sum(Gx*Gm)/np.sum(Gm)
                Gnism_in_cm_3 = (0.78+0.22*Neb*0.76)/proton_mass_in_g*Grho*1e10*solar_mass_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm
                #tpi_in_s = pidecay_fac/Gnism_in_cm_3 #pi decay time in s
                tpi_in_s = 1.0/hadronicdecayrate/Gnism_in_cm_3 #pi decay time in s from Guo 2007
                cregy_in_erg = cregy*solar_mass_in_g*1e10*km_in_cm*km_in_cm
                #tcut = tpi_in_s/betapi*nopi_per_gamma > 3e10
                #print 'total number of gas particles', np.count_nonzero(tpi_in_s)
                #print 'no of gas particles that have too much CR',np.count_nonzero(tpi_in_s[tcut])
                # here we account for the exponential in the code at each timestep. Here the time step is 3e10s
                timestepcode=3e10
                Lgammagev = cregy_in_erg*(1.-np.exp(-timestepcode/tpi_in_s*betapi/nopi_per_gamma))/timestepcode #in erg 
                if cosmo==1:
                    halosingle = read_halo_history(rundir, halonostr='00',hubble=hubble, comoving=0, maindir=maindir,singlesnap=1,atime=ascale)
                    xcen = halosingle['x'] #center of the halo physical kpc
                    ycen = halosingle['y']
                    zcen = halosingle['z']
                    Rvirnow = halosingle['R'] #physical kpc
                    Gr = np.sqrt((Gx-xcen)*(Gx-xcen)+(Gy-ycen)*(Gy-ycen)+(Gz-zcen)*(Gz-zcen))
                    #print 'hubble', hubble
                    print 'cen', xcen, ycen, zcen
                else:
                    xcen=0
                    ycen=0
                    zcen=0
                    xvcen=0
                    yvcen=0
                    zvcen=0 
                Gxr = Gx-xcen
                Gyr = Gy-ycen
                Gzr = Gz-zcen   
                cutxy = Gxr*Gxr+Gyr*Gyr < withinr*withinr
                cutz = np.absolute(Gzr)<maxlength
                cut = cutxy*cutz
                Lgammacut = Lgammagev[cut]
                
                Gmcutw = np.sum(Gm[cut])
                cregyw = np.sum(cregy[cut])
                try:
                    Nebcuta = np.average(Neb[cut],weights=Gm[cut])
                except ZeroDivisionError:
                    Nebcuta = 0
                dr = withinr/nogrid
                Lgammap=0.0
                for irad in range(nogrid):
                    cutxy = (Gxr*Gxr+Gyr*Gyr > dr*irad*dr*irad) & (Gxr*Gxr+Gyr*Gyr < dr*(irad+1)*dr*(irad+1))
                    cutz = Gzr*Gzr < maxlength*maxlength
                    cutr=cutxy*cutz
                    Nebcut = Neb[cutr]
                    Gmcut = Gm[cutr]
                    try:
                        Nebave = np.average(Nebcut,weights=Gmcut)
                    except ZeroDivisionError:
                        Nebave = 0
                    Gm_in_g = Gmcut*1e10*2e33
                    shellvol_in_cm3 = np.pi*(-np.power(dr*irad,2)+np.power(dr*(irad+1),2))*kpc_in_cm*kpc_in_cm*kpc_in_cm*2.0*maxlength
                    Grho_in_g_cm_3 = np.sum(Gm_in_g)/shellvol_in_cm3
                    Gnism_in_cm_3p = (0.78+0.22*Nebave*0.76)/proton_mass_in_g*Grho_in_g_cm_3
                    tpi_in_sp = pidecay_fac/Gnism_in_cm_3p #pi decay time in s
                    cregycut = np.sum(cregy_in_erg[cutr])
                    Lgammap += cregycut/tpi_in_sp*betapi/nopi_per_gamma #in erg/s
                    #print 'Gnism_in_cm_3p', Gnism_in_cm_3p
                Lgammapl = np.append(Lgammapl,Lgammap)
                Lgammal=np.append(Lgammal,np.sum(Lgammacut))
                enllist=np.append(enllist,np.sum(cregyl[cut]))
                Gmcutwl=np.append(Gmcutwl, Gmcutw)
                cregywl=np.append(cregywl,cregyw)
                Nebcutal=np.append(Nebcutal,Nebcuta)
                timel=np.append(timel,float(i)*0.98*1e6)
            except (KeyError,IOError):
                                print 'no snapshot'
                Lgamma_l = (enllist[1:]-enllist[:-1])/((timel[1:]-timel[:-1])*yr_in_sec)/(hadronicdecayrate+coulombdecayrate)*hadronicdecayrate*betapi/nopi_per_gamma
                #Consider CR density ~ 1eV/cm^3 and nsim ~ 1cm^-3 
                tpi_in_s_ideal=pidecay_fac
                vol_in_cm3 = np.pi*np.power(withinr,2)*kpc_in_cm*kpc_in_cm*kpc_in_cm*2.0*maxlength
                cregy_in_erg_ideal = 1.0*vol_in_cm3/erg_in_eV
                Lgamma_ideal = cregy_in_erg_ideal/tpi_in_s_ideal*betapi/nopi_per_gamma
                cregyaden_in_eV_cm_3 = cregywl*1e10*solar_mass_in_g*km_in_cm*km_in_cm*erg_in_eV/vol_in_cm3
                nism_in_cm_3 = Gmcutwl*1e10*2e33/vol_in_cm3*(0.78+0.22*Nebcutal*0.76)/proton_mass_in_g
                tpi_in_s_ave = pidecay_fac/nism_in_cm_3
                Lgammaave_nism_cre = cregywl*1e10*solar_mass_in_g*km_in_cm*km_in_cm/tpi_in_s_ave*betapi/nopi_per_gamma
                Lsfr = Kroupa_Lsf*1.0*solar_mass_in_g/yr_in_sec*cspeed_in_cm_s*cspeed_in_cm_s
                print 'Lgammal', Lgammal
                print 'Lgamma_l', Lgamma_l
                print 'Lgammapl', Lgammapl
                print 'Lgamma_ideal', Lgamma_ideal
                print 'nism_in_cm_3', nism_in_cm_3
                print 'cregyaden_in_eV_cm_3', cregyaden_in_eV_cm_3
                print 'Lgammaave_nism_cre', Lgammaave_nism_cre
                print 'Lsfr*2e-4', Lsfr*2e-4
                print 'hadronic fraction', hadronicdecayrate/(hadronicdecayrate+coulombdecayrate)
                plt.plot(timel[1:]/1e6, Lgammal[1:], label='from particles')
                plt.plot(timel[1:]/1e6, Lgammapl[1:], label='ave cylinders')
                plt.plot(timel[1:]/1e6, Lgammaave_nism_cre[1:], label='ave whole')
                plt.plot(timel[1:]/1e6, np.absolute(Lgamma_l), label='from code')
        #       plt.axhline(y=Lgamma_ideal, label='ideal',ls='dashed')  
        plt.yscale('log')
                plt.title(runtodo)
                plt.legend(loc='best')
                plt.ylabel(r'$L_{\gamma} {\rm (erg/s)}$')
                plt.xlabel('t (Myr)')
                plt.savefig('CRplot/Lgamma_ave_'+runtodo+'.pdf')
                plt.clf()



if wanted=='gaske': 
        for runtodo in dirneed:
                snaplist=[]
                gkelist=[]
                for i in range(startno, Nsnap, snapsep):
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gp = G['p']
                        Grho = G['rho']
                        Gu = G['u']
                        Gm = G['m']
            Gv = G['v']
                        #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                        Gz = Gp[:,2]
                        Gx = Gp[:,0]
                        Gy = Gp[:,1]
            Gvx = Gv[:,0]
            Gvy = Gv[:,1]
            Gvz = Gv[:,2]
            gke_in_Msun_km2_s2 = np.sum(0.5*Gm*1e10*(Gvx*Gvx+Gvy*Gvy+Gvz*Gvz))
            gke_in_erg = gke_in_Msun_km2_s2*2e33*1e5*1e5
                        snaplist=np.append(snaplist,i)
                        gkelist=np.append(gkelist,gke_in_erg)
        print 'len', len(snaplist), len(gkelist)
                plt.plot(snaplist,gkelist,label=runtodo)
        plt.xlabel('Myr')
        plt.xlim(xmax=Nsnap)
        plt.ylabel(r'$E_{\rm KE}({\rm erg})$')
        plt.yscale('log')
        plt.legend(loc='best')
        plt.savefig('CRplot/gasKE_'+fmeat+'.pdf')
        plt.clf()


if wanted=='gaskechange':
        for runtodo in dirneed:
                snaplist=[]
                gkelist=[]
        avesfrl=[]
        avesfrnewl=[]
        presnap=startno
                for i in range(startno, Nsnap, snapsep):
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
            cosmo=info['cosmo']
                        G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gp = G['p']
                        Grho = G['rho']
                        Gu = G['u']
                        Gm = G['m']
                        Gv = G['v']
                        #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                        Gz = Gp[:,2]
                        Gx = Gp[:,0]
                        Gy = Gp[:,1]
                        Gvx = Gv[:,0]
                        Gvy = Gv[:,1]
                        Gvz = Gv[:,2]
                        gke_in_Msun_km2_s2 = np.sum(0.5*Gm*1e10*(Gvx*Gvx+Gvy*Gvy+Gvz*Gvz))
                        gke_in_erg = gke_in_Msun_km2_s2*2e33*1e5*1e5
                        snaplist=np.append(snaplist,i)
                        gkelist=np.append(gkelist,gke_in_erg)
                        if cosmo==1:
                                S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                        else:
                                S = readsnap(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        try:
                                Smi = S['m']
                                Sage = S['age']
                                Sm = np.sum(Smi)
                                header=S['header']
                                timeneed=header[2]
                                #print 'time', header[2]
                                #print 'Sage', Sage
                        except KeyError:
                                Sm = 0.
                                Smi= 0.
                                Sage = 0.
                                timeneed=0.
                        if i>presnap:
                                if cosmo==1:
                                        readtimelist=readtime(firever=2)
                                        snap2list=readtimelist['snaplist']
                                        time2list=readtimelist['timelist']
                                        a2list=readtimelist['alist']
                                        tnow = np.interp(timeneed,a2list,time2list)*1e9
                                        pret = np.interp(pretime,a2list,time2list)*1e9
                                        avesfr=(Sm-presm)*1e10/(tnow-pret)
                                else:
                                        avesfr=(Sm-presm)*1e4/0.98/(i-presnap)
                        else:
                                avesfr=0
                        if Sm>1e-9:
                                tcut=Sage>pretime
                                Smnew = np.sum(Smi[tcut])
                                avesfrnew = Smnew*10./(timeneed-pretime)
                        else:
                                avesfrnew=0.
                        presnap = i
                        presm = Sm
                        pretime=timeneed
                        print 'avesfr', avesfr
                        print 'avesfrnew', avesfrnew
                        print 'Sm', Sm
                        avesfrl.append(avesfr)
                        avesfrnewl.append(avesfrnew)
                print 'len', len(snaplist), len(gkelist)
        avesfrnewl=np.array(avesfrnewl)
        avesfrl=np.array(avesfrl)
        SNerate = avesfrnewl*1e51/0.4*0.0037/3.2e7
        plt.plot(snaplist[1:],(gkelist[1:]-gkelist[:-1])/(snaplist[1:]-snaplist[:-1])/1e6/3.2e7,label=runtodo)
        plt.plot(snaplist, SNerate, ls='--')
        plt.xlabel('Myr')
        plt.xlim(xmax=Nsnap)
        plt.ylabel(r'${\rm d} E_{\rm KE}/{\rm d} t({\rm erg/s})$')
        plt.yscale('log')
        plt.legend(loc='best')
        plt.savefig('CRplot/gasKEchange_'+fmeat+'.pdf')
        plt.clf()

if wanted=='dirheader':
        for runtodo in dirneed:
                snaplist=[]
                avesfrl=[]
                avesfrnewl=[]
                presm = 0
                presnap = 0
                pretime=0
                for i in range(startno,Nsnap, snapsep):
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        cosmo=info['cosmo']
                        color=info['color']
                        print 'runtodo, color', runtodo, color
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        if cosmo==1:
                                S = readsnapcr(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1,header_only=1)
                        else:
                                S = readsnapcr(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,header_only=1)
            print 'time', S
            print 'boxsize', S['boxsize']


if wanted=='outflowall':
        for runtodo in dirneed:
                snaplist=[]
                ofalist=[]
                for i in range(startno, Nsnap, snapsep):
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        Gp = G['p']
                        Grho = G['rho']
                        Gu = G['u']
                        Gm = G['m']
                        Gv = G['v']
                        #Gnism = (0.78+0.22*Neb*0.76)/1.67e-24*Grho*1e10*1.99e33/3.086e21/3.086e21/3.086e21 #gas number density in ISM 
                        Gz = Gp[:,2]
                        Gx = Gp[:,0]
                        Gy = Gp[:,1]
                        Gvx = Gv[:,0]
                        Gvy = Gv[:,1]
                        Gvz = Gv[:,2]
                        outflowrate_in_Msun_yr = np.sum(Gm*1e10*(Gvx*Gx+Gvy*Gy+Gvz*Gz)/(Gx*Gx+Gy*Gy+Gz*Gz)/kpc_in_cm*km_in_cm*yr_in_sec)
                        snaplist=np.append(snaplist,i)
                        ofalist=np.append(ofalist,outflowrate_in_Msun_yr)
        ofalist=np.array(ofalist)
                plt.plot(snaplist,ofalist,label=runtodo)
        plt.plot(snaplist,np.absolute(ofalist),ls='dashed')
        print 'ofalist', ofalist
        plt.xlabel('Myr')
        plt.xlim(xmax=Nsnap)
        plt.ylabel(r'$\dot{M}_{\rm out} ({\rm M_\odot/yr})$')
        plt.yscale('log')
        plt.legend(loc='best')
        plt.savefig('CRplot/outflowall_'+fmeat+'.pdf')
        plt.clf()


if wanted=='outflowSasha' or wanted=='massloadingSasha' or wanted=='outflowBooth' or wanted=='massloadingBooth':
    # in Sasha definition, the outflow rate is defined with a shell 0.2-0.3Rvir
    vcut=0 #to avoid counting inflow
        rcParams['figure.figsize'] = 5,4
    M1labelneed=0
    resoneed=0
    newlabelneed=1
    diffusionsolverneed=0
    legendneed=0
        for runtodo in dirneed:
        snaplist=[]
                ofalist=[]
        nsml=[]
        timel=[]
        pretime=0.
        info=outdirname(runtodo, Nsnap)
        M1speed=info['M1speed']
                for i in range(startno, Nsnap, snapsep):
            try:
                info=outdirname(runtodo, i)
                maindir=info['maindir']
                rundir=info['rundir']
                runtitle=info['runtitle']
                slabel=info['slabel']
                snlabel=info['snlabel']
                dclabel=info['dclabel']
                resolabel=info['resolabel']
                the_snapdir=info['the_snapdir']
                Nsnapstring=info['Nsnapstring']
                havecr=info['havecr']
                Fcal=info['Fcal']
                iavesfr=info['iavesfr']
                timestep=info['timestep']
                color=info['color']
                cosmo=info['cosmo']
                withinRv=info['withinRv']
                Rvirguess =info['Rvir']
                haveB=info['haveB']
                M1speed=info['M1speed']
                ptitle=title
                if runtitle=='SMC':
                    ptitle='Dwarf'
                elif runtitle=='SBC':
                    ptitle='Starburst'
                elif runtitle=='MW':
                    ptitle=r'$L\star$ Galaxy'
                if cosmo==1:
                    h0=1
                else:
                    h0=0
                header = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo, header_only=1)
                ascale = header['time']
                print 'this time', ascale
                thisred = 1./ascale-1.
                hubble = header['hubble']
                print 'hubble', hubble
                if cosmo==1:
                    halosingle = read_halo_history(rundir, halonostr='00',hubble=hubble, comoving=0, maindir=maindir, singlesnap=1, atime=ascale)
                                        xcen = halosingle['x']
                                        ycen = halosingle['y']
                                        zcen = halosingle['z']
                    xvcen = halosingle['xv']
                    yvcen = halosingle['yv']
                    zvcen = halosingle['zv']
                    Rvirnow = halosingle['R'] 
                    MgAHF = halosingle['Mg']
                else:
                    xcen=0
                    ycen=0
                    zcen=0
                    xvcen=0
                    yvcen=0
                    zvcen=0
                                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo)
                                Gp = G['p']
                                Grho = G['rho']
                                Gu = G['u']
                                Gm = G['m']
                                Gv = G['v']
                                Gz = Gp[:,2]-zcen
                                Gx = Gp[:,0]-xcen
                                Gy = Gp[:,1]-ycen
                                Gvx = Gv[:,0]-xvcen
                                Gvy = Gv[:,1]-yvcen
                                Gvz = Gv[:,2]-zvcen
                Gr = np.sqrt((Gx)*(Gx)+(Gy)*(Gy)+(Gz)*(Gz))
                if cosmo==1:
                    Rvirguess=Rvirnow
                if wanted=='massloadingSasha' or wanted=='outflowSasha':
                    upbound=0.3*Rvirguess #kpc
                    lowbound=0.2*Rvirguess
                    dL = upbound-lowbound
                    cutu = Gr<upbound
                    cutl = Gr>lowbound
                    cut = cutu*cutl
                    Gmc = Gm[cut]
                    Gvxc = Gvx[cut]
                    Gvyc = Gvy[cut]
                    Gvzc = Gvz[cut]
                    Gxc = Gx[cut]
                    Gyc = Gy[cut]
                    Gzc = Gz[cut]
                    Grc = Gr[cut]
                    Gvdotr_r=(Gvxc*Gxc+Gvyc*Gyc+Gvzc*Gzc)/Grc
                    cutv=Gvdotr_r>vcut
                    Gvdotr_r_out=Gvdotr_r[cutv]
                    outflowrate_in_Msun_yr = np.sum(Gmc[cutv]*1e10*Gvdotr_r_out/dL/kpc_in_cm*km_in_cm*yr_in_sec)
                if wanted=='massloadingBooth' or wanted=='outflowBooth':
                    #upbound=55 #kpc
                    #lowbound=45
                    #upbound=32 #kpc
                    #lowbound=28 
                                        upbound=23 #kpc
                                        lowbound=17
                    rbound=15
                    #upbound=12 #kpc
                    #lowbound=8 #kpc
                                        dL = upbound-lowbound
                                        cutu = np.absolute(Gz)<upbound
                                        cutl = np.absolute(Gz)>lowbound
                    cutr = np.sqrt(Gx*Gx+Gy*Gy)<rbound
                                        cut = cutu*cutl*cutr
                                        Gmc = Gm[cut]
                                        Gvxc = Gvx[cut]
                                        Gvyc = Gvy[cut]
                                        Gvzc = Gvz[cut]
                                        Gxc = Gx[cut]
                                        Gyc = Gy[cut]
                                        Gzc = Gz[cut]
                                        Grc = Gr[cut]
                                        cutv=Gvzc*Gzc/np.absolute(Gzc)>vcut
                                        Gvzout=np.absolute(Gvzc[cutv])
                                        outflowrate_in_Msun_yr = np.sum(Gmc[cutv]*1e10*Gvzout/dL/kpc_in_cm*km_in_cm*yr_in_sec)
                snaplist=np.append(snaplist,i)
                ofalist=np.append(ofalist,outflowrate_in_Msun_yr)
                if wanted=='massloadingSasha' or wanted=='massloadingBooth':
                    S = readsnapcr(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo)
                    header=S['header']
                    timeneed=header[2]
                    Smi=S['m']
                    Sage=S['age']
                    Sp = S['p']
                    Sx = Sp[:,0]
                    Sy = Sp[:,1]
                    Sz = Sp[:,2]
                    header=readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo,header_only=1)
                    h0 = header['hubble']
                    atime = header['time']
                    Sxrel = Sx-xcen
                    Syrel = Sy-ycen
                    Szrel = Sz-zcen
                    Sr = np.sqrt(Sxrel*Sxrel+Syrel*Syrel+Szrel*Szrel)
                    cutrvs = Sr<Rvirguess
                    Smi = Smi[cutrvs]
                    Sage = Sage[cutrvs]
                    Sm = np.sum(Smi)*1e10 #in solar mass
                    tcut=Sage>pretime
                    Nsm = np.sum(Smi[tcut])*1e10
                    pretime=timeneed
                    if cosmo==1:
                        readtimelist=readtime(firever=2)
                        snap2list=readtimelist['snaplist']
                        time2list=readtimelist['timelist']
                        a2list=readtimelist['alist']
                        tnow = np.interp(timeneed,a2list,time2list)*1e9
                        pret = np.interp(pretime,a2list,time2list)*1e9
                    if cosmo==1:
                        timel.append(tnow)
                    else:
                        timel.append(float(i)*0.98*1e6)
                    nsml=np.append(nsml,Nsm)

            except (KeyError,IOError): 
                print 'no snapshot'
                continue
                #if wanted=='massloadingSasha' or wanted=='massloadingBooth':
                #   timel.append(pretime)
                #   nsml=np.append(nsml,0.)
        ofalist=np.array(ofalist)
        if wanted=='massloadingSasha' or wanted=='massloadingBooth': 
            nsml=np.array(nsml)
            timel=np.array(timel)
            avesfrl=(nsml[1:])/(timel[1:]-timel[:-1]) #in Msun/yr
        labelneed=dclabel
        if haveB>0:
            lsn = 'dashed'
        else:
            lsn = 'solid'
        if M1labelneed==1:
            if M1speed>499:
                lsn = 'solid'
                labelneed='M1=500'
            if M1speed>999:
                lsn = 'dashed'
                labelneed='M1=1000'
            if M1speed>1999:
                lsn = 'dashdot'
                labelneed='M1=2000'
            if M1speed>3999:
                lsn = 'dotted'
                labelneed='M1=4000'
        if resoneed==1:
            if resolabel=='llr':
                lsn= 'solid'
                labelneed='Lowest res'
            if resolabel=='lr':
                lsn= 'dashed'
                labelneed='Lower res'
            if resolabel=='mr':
                lsn='dashdot'
                labelneed='Standard res'
                if diffusionsolverneed==1:
                        if runtodo=='bwmwlrdc27ds':
                                labelneed='Zeroth moment'
                                lsn = 'dashed'
                        else:
                                labelneed='Two moment'
                                lsn = 'solid'
                
        if wanted=='outflowSasha' or wanted=='outflowBooth':
            plt.plot(snaplist,ofalist,label=labelneed,color=color,lw=2,ls=lsn)
        if wanted=='massloadingSasha' or wanted=='massloadingBooth':
            if M1labelneed==1 and color=='r':
                plt.plot(snaplist[1:],ofalist[1:]/avesfrl,color=color,lw=2,ls=lsn)
            else:
                plt.plot(snaplist[1:],ofalist[1:]/avesfrl,label=labelneed,color=color,lw=2,ls=lsn)
            print 'ofalist[1:]/avesfrl', ofalist[1:]/avesfrl
            print 'ofalist[1:]', ofalist[1:]
            print 'avesfrl', avesfrl
            print 'nsml', nsml
            print 'timel', timel
    if wanted=='outflowSasha':
        plt.xlabel(r'${\rm Myr}$',fontsize=18)
        plt.xlim(xmax=Nsnap)
        plt.ylabel(r'$\dot{M}_{\rm out,r} ({\rm M_\odot/yr})$',fontsize=18)
        plt.yscale('log')
        plt.legend(loc='best')
        filename='CRplot/outflowSasha/outflowSasha_'+fmeat+'.pdf'
        if wanted=='outflowBooth':
                plt.xlabel(r'${\rm Myr}$', fontsize=18)
                plt.xlim(xmax=Nsnap)
                plt.ylabel(r'$\dot{M}_{\rm out,z} ({\rm M_\odot/yr})$', fontsize=18)
                plt.yscale('log')
                plt.legend(loc='best')
                filename='CRplot/outflowBooth/outflowBooth_'+fmeat+'.pdf'
    if wanted=='massloadingSasha':
        plt.xlabel(r'${\rm Myr}$', fontsize=18)
        plt.xlim(xmax=Nsnap)
        plt.ylabel(r'$\dot{M}_{\rm out,r}/{\rm SFR}$', fontsize=18)
        plt.yscale('log')
        plt.legend(loc='best')
        filename='CRplot/massloadingSasha/massloadingSasha_'+fmeat+'.pdf'
        if wanted=='massloadingBooth':
                plt.xlabel(r'${\rm Myr}$', fontsize=18)
                plt.xlim(xmax=Nsnap)
                plt.ylabel(r'$\dot{M}_{\rm out,z}/{\rm SFR}$', fontsize=18)
                plt.yscale('log')
        plt.title(ptitle,fontsize=18)
        if runtitle=='SMC' and legendneed==1:
            plt.legend(loc='best',fontsize=8,ncol=3)
        filename='CRplot/massloadingBooth/massloadingBooth_'+fmeat+'.pdf'
    print 'filename', filename
        plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=16)
        plt.savefig(filename,bbox_inches='tight')
    plt.clf()




if wanted=='testcr':
        rcParams['figure.figsize'] = 5, 5
        for runtodo in dirneed:
                snaplist=[]
                enclist=[]
                englist=[]
                enllist=[]
                endlist=[]
                enplist=[]
                enalist=[]
                avesfrl=[]
                prel=0
                preg=0
                prec=0
                pred=0
                prep=0
                presm = 0
                presnap = 0
        precregyp = 0
        precregyl = 0
                for i in [Nsnap-1,Nsnap]:
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        maindir=info['maindir']
                        haveB=info['haveB']
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        cosmo=info['cosmo']
                        if cosmo==1:
                                h0=1
                        else:
                                h0=0
                        header = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo, header_only=1)
                        ascale = header['time']
                        print 'this time', ascale
                        thisred = 1./ascale-1.
                        hubble = header['hubble']
                        print 'hubble', hubble
                        G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=h0,cosmological=cosmo)
                        if cosmo==1:
                                halosingle = read_halo_history(rundir, halonostr='00',hubble=hubble, comoving=0, maindir=maindir, singlesnap=1, atime=ascale)
                                xcen = halosingle['x']
                                ycen = halosingle['y']
                                zcen = halosingle['z']
                                xvcen = halosingle['xv']
                                yvcen = halosingle['yv']
                                zvcen = halosingle['zv']
                                Rvirnow = halosingle['R']
                                MgAHF = halosingle['Mg']
                        #       print 'xcenl[0]', xcenl[0]
                        #       print 'thisred', thisred
                                print 'cen', xcen, ycen, zcen
                                print 'MgAHF', MgAHF
                                #print 'Rvir', Rvirnow
                        else:
                                xcen=0
                                ycen=0
                                zcen=0
                                xvcen=0
                                yvcen=0
                                zvcen=0
                        Gpos = G['p']
                        Gvel = G['v']
                        Gx = Gpos[:,0]-xcen
                        Gy = Gpos[:,1]-ycen
                        Gz = Gpos[:,2]-zcen
                        Gvx = Gvel[:,0]-xvcen   #km/s
                        Gvy = Gvel[:,1]-yvcen
                        Gvz = Gvel[:,2]-zvcen
                        Grho = G['rho']*1e10 #Msun/kpc^3
                        Gm = G['m']*1e10 #Msun
                        GEint = G['u']*km_in_cm*km_in_cm*Gm*solar_mass_in_g
                        if havecr>0:
                                cregyl = G['cregyl']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                                cregyg = G['cregyg']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                                cregy  = G['cregy']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                                cregyd = G['cregyd']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                        if havecr>4:
                                cregyp = G['cregyp']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
                        if havecr>5:
                                cregya = G['cregya']*1e10*solar_mass_in_g*km_in_cm*km_in_cm
            diffp=np.sum(cregyp-precregyp)
            diffl=np.sum(cregyl-precregyl)
            precregyp = cregyp
            precregyl = cregyl
            print 'diffp, diffl', diffp, diffl

    


if wanted=='nismB':
    nogrid=15
        for runtodo in dirneed:
        LogGnxaxis = np.linspace(-7,3,num=nogrid)
        aveB2l = 0.0*LogGnxaxis
        errl = 0.0*LogGnxaxis+1.0
        numoftimes=0
        info=outdirname(runtodo, Nsnap)
        haveB=info['haveB']
        if haveB<1:
            'no B field'
            continue
                for i in range(startno,Nsnap,snapsep):
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        cosmo=info['cosmo']
                        haveB=info['haveB']
            color=info['color']
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        if cosmo==1:
                                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1, cosmological=1)
                        else:
                                G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        GB = G['B']
                        Gm = G['m']
            Gpos = G['p']
                        Grho = G['rho']
                        Neb = G['ne']
            Gx = Gpos[:,0]
            Gy = Gpos[:,1]
            Gz = Gpos[:,2]
            Gr = np.sqrt(Gx*Gx+Gy*Gy)
            cutr = Gr<10.0 #kpc
            cutz = np.absolute(Gz)<1.0 #kpc
            cut = cutr*cutz
            Grho = Grho[cut]
            Neb = Neb[cut]
            Gm = Gm[cut]
                        Gnism_in_cm_3 = 1.0/proton_mass_in_g*Grho*1e10*Msun_in_g/kpc_in_cm/kpc_in_cm/kpc_in_cm
                        GB2 = GB[:,0]*GB[:,0]+GB[:,1]*GB[:,1]+GB[:,2]*GB[:,2]
            GB2 = GB2[cut]
                        LogGnism = np.log10(Gnism_in_cm_3)
            print 'LogGnism', LogGnism
                        for ib in range(nogrid-1):
                try:
                    cutgu = (LogGnism > LogGnxaxis[ib])
                    cutgd = (LogGnism < LogGnxaxis[ib+1])
                    cutg = cutgu*cutgd
                    GB2cut = GB2[cutg]
                    #print 'GB2cut', GB2cut
                    aveB2 = np.average(GB2cut, weights=Gm[cutg]/Grho[cutg])
                    print 'aveB', np.sqrt(aveB2)*1.0e6
                    aveB2l[ib] += aveB2
                except ZeroDivisionError:
                    print 'ZeroDivisionError'
                    errl[ib]=-1.0
            numoftimes+=1
        cuterr = errl>0.
        xdata = LogGnxaxis[cuterr]
        ydata = np.log10(aveB2l[cuterr]/numoftimes)/2.+6.
        xdata = xdata[:-1]
        ydata = ydata[:-1]
        print 'xdata', xdata
        print 'ydata', ydata
        x0=np.array([0.0,0.0])
        def func(x, a, b):
            return a+b*x
        try:
            outfit=optimize.curve_fit(func, xdata, ydata, x0)
            alpha=outfit[0][1]
            offset = outfit[0][0]
            print 'alpha', alpha
        except (RuntimeError, TypeError, NameError):
            print "alpha not exist"
            alpha=-3
            offset=-3
            continue
        plt.plot(xdata,offset+alpha*xdata,color=color, label=dclabel+'; slope = '+str(alpha))
        plt.plot(xdata,ydata, ls='none',marker='s', mfc=color)
        plt.legend(loc='best')
        plt.xlabel(r'$\log (n_{\rm ISM}[{\rm cm^{-3}}])$')
        plt.ylabel(r'$\log (B_{\rm rms} [\mu G])$')
    filename = 'CRplot/nismB/nismB_'+fmeat+'_sn'+str(startno)+'_sn'+str(Nsnap)+'.pdf'
    print filename
        plt.savefig(filename)
        plt.clf()

if wanted=='dirage' or wanted=='diragestd' or wanted=='diragemax' or wanted=='dirage10_200Myr' or\
 wanted=='diragestd10_200Myr' or wanted=='dirage200Myr' or wanted=='diragehis' or wanted=='dirssfr10_200Myr':
        rcParams['figure.figsize'] = 8,4
        Rfrac=0.25
        normalizedsm=1
    startt = 6
    endt = 13.7
    withinRv=1
    agestd10Myrl = []
        agemax10Myrl = []
        agestd200Myrl = []
        agemax200Myrl = []
    Msl = []
    crlabl = []
        for runtodo in dirneed:
        time10Myrl = np.linspace(startt,endt,int((endt-startt)/0.01))
        time200Myrl = np.linspace(startt,endt,int((endt-startt)/0.2))
                for i in [Nsnap]:
            print 'runtodo', runtodo
                        info=outdirname(runtodo, i)
                        rundir=info['rundir']
                        maindir=info['maindir']
            print 'rundir, maindir', rundir, maindir
                        halostr=info['halostr']
                        runtitle=info['runtitle']
                        slabel=info['slabel']
                        snlabel=info['snlabel']
                        dclabel=info['dclabel']
            print 'dclabel', dclabel
                        resolabel=info['resolabel']
                        the_snapdir=info['the_snapdir']
                        Nsnapstring=info['Nsnapstring']
                        havecr=info['havecr']
                        Fcal=info['Fcal']
                        iavesfr=info['iavesfr']
                        timestep=info['timestep']
                        cosmo=info['cosmo']
                        color=info['color']
                        withinRv=info['withinRv']
                        usepep=info['usepep']
                        beginno=info['beginno']
                        finalno=info['finalno']
                        firever=info['firever']
                        initsnap=info['initsnap']
                        haveB=info['haveB']
                        M1speed=info['M1speed']
                        Rvirguess=info['Rvir']
            crlabel=info['crlabel']
            snumadd=info['snumadd']
                        if cosmo==0:
                                inittime=initsnap
                        else:
                                inittime=0
                                print 'set initial time'
                        print 'the_snapdir', the_snapdir
                        print 'Nsnapstring', Nsnapstring
                        print 'havecr', havecr
                        print 'withinRv', withinRv
                        print 'cosmo', cosmo
                        try:
                                if cosmo==1:
                                        S = readsnapcr(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1)
                                else:
                                        S = readsnapcr(the_snapdir, Nsnapstring, 4, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                        except KeyError:
                                print 'Keyerror'
                                break
                        try:
                                header=S['header']
                                timeneed=header[2]
                                print 'timeneed', timeneed
                                Smi=S['m']
                                Sage=S['age']
                print 'Smi', Smi
                print 'Sage', Sage
                                if withinRv ==1 and cosmo==1:
                                        Sp = S['p']
                                        Sx = Sp[:,0]
                                        Sy = Sp[:,1]
                                        Sz = Sp[:,2]
                                        header=readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1,cosmological=1,header_only=1)
                                        h0 = header['hubble']
                                        atime = header['time']
                                        if usepep==1:
                                                halosA = SF.read_halo_history_pep(rundir, finalno, beginno=beginno,\
                        singlesnap=0, firever=firever,halonostr=halostr, hubble=h0, comoving=1, maindir=maindir)
                        afactor=atime
                                        else:
                                                halosA = SF.read_halo_history(rundir, maindir=maindir, halonostr=halostr,\
                         hubble=h0, comoving=0, snumadd=snumadd)
                        afactor=1.0
                                        redlist = halosA['redshift']
                                        haloid = halosA['ID']
                                        a_scale = 1.0/(1.0+redlist)
                                        xcenl = halosA['x']*atime
                                        ycenl = halosA['y']*atime
                                        zcenl = halosA['z']*atime
                                        Rvirl = halosA['R']*atime
                    Mstarl = halosA['Ms']
                                        xcen = np.interp(atime,a_scale,xcenl)
                                        ycen = np.interp(atime,a_scale,ycenl)
                                        zcen = np.interp(atime,a_scale,zcenl)
                                        Rvir = np.interp(atime,a_scale,Rvirl)
                    Mstar = np.interp(atime,a_scale,Mstarl)
                                        Sxrel = Sx-xcen
                                        Syrel = Sy-ycen
                                        Szrel = Sz-zcen
                                        Sr = np.sqrt(Sxrel*Sxrel+Syrel*Syrel+Szrel*Szrel)
                                        cutrvs = Sr<Rvir*Rfrac
                                        Smi = Smi[cutrvs]
                                        Sage = Sage[cutrvs]
                    print 'Sage', Sage
                if cosmo==1:
                    readtimelist=readtime(firever=2)
                    snap2list=readtimelist['snaplist']
                    time2list=readtimelist['timelist']
                    a2list=readtimelist['alist']
                    Sage = np.interp(Sage,a2list,time2list)
                histage, bin_edges = np.histogram(Sage, weights=1e10*Smi,bins=1000)
                histage10Myr, bin10Myr = np.histogram(Sage[Sage>6.01], weights=1e10*Smi[Sage>6.01], bins=time10Myrl)
                histage200Myr, bin200Myr = np.histogram(Sage[Sage>6.01], weights=1e10*Smi[Sage>6.01], bins=time200Myrl)
                ha10Myr = histage10Myr/10./1e6+1e-10
                ha200Myr = histage200Myr/200./1e6
                print 'median 10 Myr', np.median(np.log10(ha10Myr))
                print 'median 200 Myr', np.median(np.log10(ha200Myr))
                #print 'ha10Myr[~np.isfinite(ha10Myr)]', ha10Myr[~np.isfinite(ha10Myr)]
                print 'max sfr 10Myr', np.amax(np.log10(ha10Myr))
                print 'std 10Myr', np.std(np.log10(ha10Myr))
                print 'std 10Myr (>1e-5)', np.std(np.log10(ha10Myr[ha10Myr>1e-5]))
                print 'std 200Myr', np.std(np.log10(ha200Myr))
                agestd10Myrl = np.append(agestd10Myrl,np.std(np.log10(ha10Myr)))
                agemax10Myrl = np.append(agemax10Myrl,np.max(np.log10(ha10Myr)))
                                agestd200Myrl = np.append(agestd200Myrl,np.std(np.log10(ha200Myr)))
                                agemax200Myrl = np.append(agemax200Myrl,np.max(np.log10(ha200Myr)))
                Msl = np.append(Msl,Mstar)      
                crlabl = np.append(crlabl,crlabel)
                        except KeyError:
                                print 'key error'
                exit()
            #print 'histage', histage
        #plt.plot(np.absolute(bin_edges[1:]+bin_edges[:-1])/2.,histage/1e9/np.absolute(bin_edges[1:]-bin_edges[:-1]),label=dclabel)
            if wanted=='diragehis':
                sfr10Myr = np.log10(ha10Myr)
                histsfr, binsfr = np.histogram(sfr10Myr, bins=15)
                plt.plot(binsfr[:-1],histsfr,label=dclabel,lw=1)
            if wanted=='dirage':
                plt.plot(np.absolute(bin10Myr[1:]+bin10Myr[:-1])/2.,\
                ha10Myr/np.absolute(bin10Myr[1:]-bin10Myr[:-1]),label=dclabel,lw=1)
                plt.ylim(ymin=1e-3)
                        if wanted=='dirage200Myr':
                                plt.plot(np.absolute(bin200Myr[1:]+bin200Myr[:-1])/2.,\
                                ha200Myr/np.absolute(bin200Myr[1:]-bin200Myr[:-1]),label=dclabel,lw=1)
                                plt.ylim(ymin=1e-3)
                        if wanted=='dirage10_200Myr':
                                plt.plot(ha10Myr[::21],ha200Myr,label=dclabel,ls='none',marker='o')
                        if wanted=='dirssfr10_200Myr':
                                plt.plot(ha10Myr[::21]/Msl[-1],ha200Myr/Msl[-1],label=dclabel,ls='none',marker='o') 

    if wanted=='diragestd':
        plt.plot(Msl[crlabl==''],agestd10Myrl[crlabl==''],ls='none',marker='o',label='Hydro')
        plt.plot(Msl[crlabl=='mhdcv'],agestd10Myrl[crlabl=='mhdcv'],ls='none',marker='s',label='MHD')
        plt.plot(Msl[crlabl=='cr70'],agestd10Myrl[crlabl=='cr70'],ls='none',marker='s',label=r'$\kappa$=3e28')
        plt.plot(Msl[crlabl=='cr700'],agestd10Myrl[crlabl=='cr700'],ls='none',marker='^',label=r'$\kappa$=3e29')        
        plt.xscale('log')
        plt.xlabel(r'${\rm M}_{\rm *}[{\rm M}_\odot]$')
        plt.ylabel('Scatter in SFR (10 Myr) [dex]')
        figname='CRplot/diragestd/diragestd_'+fmeat+'.pdf'
        if wanted=='diragestd10_200Myr':
        rat = agestd10Myrl-agestd200Myrl
                plt.plot(Msl[crlabl==''],rat[crlabl==''],ls='none',marker='o',label='Hydro')
                plt.plot(Msl[crlabl=='mhdcv'],rat[crlabl=='mhdcv'],ls='none',marker='s',label='MHD')
                plt.plot(Msl[crlabl=='cr70'],rat[crlabl=='cr70'],ls='none',marker='s',label=r'$\kappa$=3e28')
                plt.plot(Msl[crlabl=='cr700'],rat[crlabl=='cr700'],ls='none',marker='^',label=r'$\kappa$=3e29')
                plt.xscale('log')
                plt.xlabel(r'${\rm M}_{\rm *}[{\rm M}_\odot]$')
                plt.ylabel(r' SFR (10 Myr)/ SFR (200Myr) Scatter [dex]')
                figname='CRplot/diragestd10_200Myr/diragestd10_200Myr_'+fmeat+'.pdf'
        if wanted=='diragemax':
        ssfrmax = agemax10Myrl-np.log10(Msl)
                plt.plot(Msl[crlabl==''],ssfrmax[crlabl==''],ls='none',marker='o',label='Hydro')
                plt.plot(Msl[crlabl=='mhdcv'],ssfrmax[crlabl=='mhdcv'],ls='none',marker='s',label='MHD')
                plt.plot(Msl[crlabl=='cr70'],ssfrmax[crlabl=='cr70'],ls='none',marker='s',label=r'$\kappa$=3e28')
                plt.plot(Msl[crlabl=='cr700'],ssfrmax[crlabl=='cr700'],ls='none',marker='^',label=r'$\kappa$=3e29')
                plt.xscale('log')
                plt.xlabel(r'${\rm M}_{\rm *}[{\rm M}_\odot]$')
                plt.ylabel(r'Log Maximum SSFR (10 Myr) [${\rm yr^{-1}}$]')
                figname='CRplot/diragemax/diragemax_'+fmeat+'.pdf'
    if wanted=='dirage':
        plt.xlabel('stellar age [Gyr]',fontsize=16)
        plt.ylabel(r'${\rm d\;\hat{M}/d\;age [M_\odot/yr]}$',fontsize=18)
        plt.yscale('log')
        plt.legend(loc='best')
        figname='CRplot/dirage/dirage_'+fmeat+'.pdf'
        if wanted=='dirage200Myr':
                plt.xlabel('stellar age [Gyr]',fontsize=16)
                plt.ylabel(r'${\rm d\;\hat{M}/d\;age [M_\odot/yr]}$',fontsize=18)
                plt.yscale('log')
                plt.legend(loc='best')
                figname='CRplot/dirage/dirage200Myr_'+fmeat+'.pdf'
        if wanted=='dirage10_200Myr':
        sfrline = np.power(10.,np.linspace(-3,2,num=20))
        print 'sfrline', sfrline
        plt.plot(sfrline,sfrline,ls='dashed')
                plt.xlabel(r'${\rm {\rm SFR}\;(10Myr)[M_\odot/yr]}$',fontsize=16)
                plt.ylabel(r'${\rm {\rm SFR}\;(200Myr)[M_\odot/yr]}$',fontsize=16)
        plt.xscale('log')
                plt.yscale('log')
        plt.xlim(xmin=1e-4)
                plt.legend(loc='best')
                figname='CRplot/dirage/dirage10_200Myr_'+fmeat+'.pdf'
        if wanted=='dirssfr10_200Myr':
                sfrline = np.power(10.,np.linspace(-13,-9,num=20))
                plt.plot(sfrline,sfrline,ls='dashed')
                plt.xlabel(r'${\rm {\rm sSFR}\;(10Myr)[yr^{-1}]}$',fontsize=16)
                plt.ylabel(r'${\rm {\rm sSFR}\;(200Myr)[yr^{-1}]}$',fontsize=16)
        plt.xlim(xmin=1e-13)
                plt.xscale('log')
                plt.yscale('log')
                plt.legend(loc='best')
                figname='CRplot/dirage/dirssfr10_200Myr_'+fmeat+'.pdf'
        if wanted=='diragehis':
                plt.xlabel('Log (SFR[10Myr])',fontsize=16)
                plt.ylabel('Frequency',fontsize=18)
                plt.legend(loc='best')
                figname='CRplot/dirage/diragehis_'+fmeat+'.pdf'
    plt.legend()
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.savefig(figname,bbox_inches='tight')
    plt.clf()



if wanted=='Tz' or wanted=='rhoz' or wanted=='vzz' or wanted=='Tvz' or wanted=='pcrpth':
    #startno=450
    #Nsnap=500
    from crtestfunction import *    
        rcParams['figure.figsize'] = 5, 4
    Tcut=-1.0 #K
    highTcut = 1e10
    vcut=0.0 #outflow velocity cut
    vhcut = 1e10 #upper outflow velocity cut
    withinr=20.0
    zup=100.0
    zdown=1.0
    trackgas=0
    tlabel=''
    userad=0
    if trackgas==1:
        snaptrack = Nsnap
        Tcut_t=1.0 #K
        highTcut_t = 1e10
        vcut_t=0.0 #outflow velocity cut
        vhcut_t = 1e10 #upper outflow velocity cut
        withinr_t=20.0
        zup_t=30.0
        zdown_t=20.0
        tlabel='track'
    if userad==1:
        tlabel+='rad'
    if wanted=='Tz':
        extent = [zdown,zup,1,9]
    if wanted=='rhoz':
        extent = [zdown,zup,-7,4] 
    if wanted=='vzz':
        extent = [zdown,zup,0,3.2] 
    if wanted=='Tvz':
        extent = [0,3.2,1,9] 
    if wanted=='pcrpth':
        extent = [-20,-11,-20,-11]
        nobin=51
        needcontour=1
        for runtodo in dirneed:
        info=outdirname(runtodo, Nsnap)
        cosmo=info['cosmo']
        runtitle=info['runtitle']
        havecr=info['havecr']
        if wanted=='pcrpth' and havecr==0:
            continue
                ptitle=title
                if runtitle=='SMC':
                        ptitle='Dwarf'
                elif runtitle=='SBC':
                        ptitle='Starburst'
                elif runtitle=='MW':
                        ptitle=r'$L\star$ Galaxy'

        #if cosmo==1:
    #       userad=1
        Hadd = np.zeros((nobin-1,nobin-1))
        if trackgas==1:
            data = gaswindphase(runtodo,snaptrack,Tcut=Tcut_t,highTcut=highTcut_t,\
            vcut=vcut_t,vhcut=vhcut_t,withinr=withinr_t,zup=zup_t,zdown=zdown_t,userad=userad)
            Gid_t = data['Gid']
            nooftimes=0
            for i in range(startno,Nsnap,snapsep):
            if trackgas==1:
                vcut = -1e10 #if we track gas particles, we should consider all velocities
                if nooftimes==1:
                    continue
            data = gaswindphase(runtodo,i,Tcut=Tcut,highTcut=highTcut,\
            vcut=vcut,vhcut=vhcut,withinr=withinr,zup=zup,zdown=zdown,userad=userad)
            partZ = data['partZ']
            partR = data['partR']
            vz = data['vz']
            vr = data['vr']
            TrueTemp = data['TrueTemp']
            Gu = data['Gu']
            rho = data['rho']
            cregy = data['cregy']
            converted_rho = data['convertedrho']
            Gmass = data['Gmass']   
            vmax = data['vmax']
            #vmax=8
            vmin = data['vmin']
            #vmin=0
            Gid = data['Gid']
            if trackgas==1:
                idint0=np.in1d(Gid,Gid_t)
                if userad==1:
                    partR=partR[idint0]
                    vr=vr[idint0]
                else:
                    partZ=partZ[idint0]
                    vz=vz[idint0]
                TrueTemp=TrueTemp[idint0]
                converted_rho=converted_rho[idint0]
                Gmass = Gmass[idint0]
            if userad==1:
                x = partR
            else:
                x = partZ
            if wanted=='Tvz':
                if userad==1:
                    x = np.log10(vr)
                else:
                    x = np.log10(vz)
            if wanted=='pcrpth':
                pth = (GAMMA-1.0)*Gu*u_codetocgs*(rho*rho_codetocgs) #cgs
                x = np.log10(pth)
                print 'x', x
            if wanted=='Tz' or wanted=='Tvz':
                y = np.log10(TrueTemp)
            if wanted=='rhoz':
                y = np.log10(converted_rho)
            if wanted=='vzz':
                if userad==1:
                    y = np.log10(vr)
                else:
                    y = np.log10(vz)
            if wanted=='pcrpth':
                pcr = (CRgamma-1.0)*cregy*cregy_codetocgs*(rho*rho_codetocgs)/(Gmass*m_codetocgs) #cgs
                y = np.log10(pcr)   
                print 'y', y
                gridxx = np.linspace(extent[0],extent[1],nobin)
                gridyy = np.linspace(extent[2],extent[3],nobin)
                H, xedges, yedges = np.histogram2d(x, y, bins=[gridxx, gridyy],weights=Gmass*1e10)
                H=H.T
                Hadd += H
                nooftimes += 1
                Hadd=Hadd/nooftimes
                if needcontour==1:
                        levels = np.linspace(vmin,vmax,num=9)
                        plt.contourf(xedges[:-1], yedges[:-1], np.log10(Hadd), levels=levels, extent=extent,origin='lower',aspect='auto')
                else:
                        plt.imshow(np.log10(Hadd), extent=extent, interpolation='nearest',origin='lower',aspect='auto',vmax=vmax,vmin=vmin)
                cbar = plt.colorbar(extend='both', norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax))
                cbar.set_label(r'$\mathrm{Log}_{10} (M {\rm [M_\odot]})$', rotation=270, labelpad=20)
        if userad==1:
            plt.xlabel(r"$r {\rm [kpc]}$")
        else:
            plt.xlabel(r"$z {\rm [kpc]}$")
        if wanted=='Tvz':
            if userad==1:
                plt.xlabel(r"$\mathrm{Log}_{10} (v_r {\rm [km/s]})$")
            else:
                plt.xlabel(r"$\mathrm{Log}_{10} (v_z {\rm [km/s]})$")
        if wanted=='pcrpth':
            plt.xlabel(r"$\mathrm{Log}_{10} (P_{\rm th} {\rm [erg/cm^3]})$")
            plt.ylabel(r"$\mathrm{Log}_{10} (P_{\rm cr} {\rm [erg/cm^3]})$")
        if wanted=='rhoz':
            plt.ylabel(r"$\mathrm{Log}_{10} (n {\rm [cm^{-3}]})$")
        if wanted=='Tz' or wanted=='Tvz':
            plt.ylabel(r"$\mathrm{Log}_{10} (T {\rm [K]})$")
        if wanted=='vzz':
            if userad==1:
                plt.ylabel(r"$\mathrm{Log}_{10} (v_r {\rm [km/s]})$")
            else:
                plt.ylabel(r"$\mathrm{Log}_{10} (v_z {\rm [km/s]})$")
                if wanted == 'Tz':
                        if needcontour==1:
                                totalname = 'CRplot/Tz/Tz_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+'_contour'+tlabel+'.pdf'
                        else:
                                totalname = 'CRplot/Tz/Tz_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+tlabel+'.pdf'
                if wanted == 'rhoz':
                        if needcontour==1:
                                totalname = 'CRplot/rhoz/rhoz_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+'_contour'+tlabel+'.pdf'
                        else:
                                totalname = 'CRplot/rhoz/rhoz_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+tlabel+'.pdf'
                if wanted == 'vzz':
                        totalname = 'CRplot/vzz/vzz_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+tlabel+'.pdf'
                if wanted == 'Tvz':
                        totalname = 'CRplot/Tvz/Tvz_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+tlabel+'.pdf'
        if wanted == 'pcrpth':
            totalname = 'CRplot/pcrpth/pcrpth_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+tlabel+'.pdf'
            plt.title(ptitle,fontsize=18)
            plt.tight_layout()
            print totalname
            plt.savefig(totalname)
            plt.clf()






if wanted=='phaseTz' or wanted=='phaseTwindz':
    from crtestfunction import gaswindphase
    rcParams['figure.figsize'] = 5, 4
    diffon=1
    userad=1
    vcut=30.0 #outflow velocity cut
    if wanted=='phaseTz':
        vcut=-1e10
        withinr=20.0
        zup=30.0
        zdown=10.0
        nooftimes=0
        nobin=51
        zl = np.linspace(zdown,zup,num=nobin)
        if diffon==1:
                filesuffix = 'diff'
        else:
                filesuffix = 'cum'
        for runtodo in dirneed:
        hotmass = zl*0.0
        warmmass = zl*0.0
        lowmass = zl*0.0
        coldmass = zl*0.0
        for i in range(startno,Nsnap,snapsep):
            data=gaswindphase(runtodo,i,vcut=vcut,withinr=withinr,zup=zup,zdown=zdown,userad=userad)
            partX=data['partX']
            partY=data['partY']
            partZ=data['partZ']
            partR=data['partR']
            TrueTemp=data['TrueTemp']
            Gmass=data['Gmass'] 
            for iz in range(nobin-1):
                if diffon==1:
                    if userad==1:
                        cutr = (np.absolute(partR)>zl[iz]) & (np.absolute(partR)<zl[iz+1])
                        cut=cutr
                    else: 
                        cutz = (np.absolute(partZ)>zl[iz]) & (np.absolute(partZ)<zl[iz+1])
                        cut=cutz
                else:
                    if userad==1:
                        cutr = (np.absolute(partR)>zdown) & (np.absolute(partR)<zl[iz+1])
                        cut=cutr
                    else:   
                        cutz = (np.absolute(partZ)>zdown) & (np.absolute(partZ)<zl[iz+1])
                        cut=cutz
                                TrueTempcut = TrueTemp[cut]
                                Gmasscut = Gmass[cut]
                                warmcut = np.log10(TrueTempcut) < 5.3
                                lowcut = np.log10(TrueTempcut) < 4.7
                                coldcut = np.log10(TrueTempcut) < 4.0
                                hotmass[iz+1] += np.sum(Gmasscut)
                                warmmass[iz+1] += np.sum(Gmasscut[warmcut])
                                lowmass[iz+1] += np.sum(Gmasscut[lowcut])
                                coldmass[iz+1] += np.sum(Gmasscut[coldcut])
                        nooftimes+=1
                hotmass = hotmass/nooftimes
                warmmass = warmmass/nooftimes
                lowmass = lowmass/nooftimes
                coldmass = coldmass/nooftimes
                if diffon==1:
                        dz = (zup-zdown)/nobin
                        hotmass=hotmass/dz
                        warmmass=warmmass/dz
                        lowmass=lowmass/dz
                        coldmass=coldmass/dz
                plt.plot(zl,hotmass,color='r',label=r'$T\;[{\rm K}]>10^{5.3}$')
                plt.plot(zl,warmmass,color='y',label=r'$10^{4.7}<T\;[{\rm K}]<10^{5.3}$')
                plt.plot(zl,lowmass,color='b', label=r'$10^{4.0}<T\;[{\rm K}]<10^{4.7}$')
                plt.plot(zl,coldmass,color='m', label=r'$T\;[{\rm K}]<10^{4.0}$')
                plt.fill_between(zl, hotmass, warmmass, where=hotmass >= warmmass, facecolor='r', interpolate=True)
                plt.fill_between(zl, warmmass, lowmass, where=warmmass >= lowmass, facecolor='y', interpolate=True)
                plt.fill_between(zl, lowmass, coldmass, where=lowmass >= coldmass, facecolor='b', interpolate=True)
                plt.fill_between(zl, coldmass, 0.0*coldmass, where=coldmass >= 0.0*coldmass, facecolor='m', interpolate=True)
                if wanted == 'phaseTz':
                        totalname = 'CRplot/phaseTz/phaseTz_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+filesuffix+'.pdf'
                elif wanted == 'phaseTwindz':
                        totalname = 'CRplot/phaseTwindz/phaseTwindz_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+filesuffix+'.pdf'
        if userad==1:
            plt.xlabel(r"$r\;[{\rm kpc}]$",fontsize=18)
        else:
            plt.xlabel(r"$z\;[{\rm kpc}]$",fontsize=18)
        if userad==1:
            if diffon==1:
                                if wanted == 'phaseTwindz':
                                        plt.ylabel(r"${\rm d}M_{\rm w}/{\rm d}r\;[{\rm M}_\odot/{\rm kpc}]$",fontsize=18)
                                else:
                                        plt.ylabel(r"${\rm d}M/{\rm d}r\;[{\rm M}_\odot/{\rm kpc}]$",fontsize=18)
                        else:
                                if wanted == 'phaseTwindz':
                                        plt.ylabel(r"$M_{\rm w}(<r)\;[{\rm M}_\odot]$",fontsize=18)
                                else:
                                        plt.ylabel(r"$M(<r)\;[{\rm M}_\odot]$",fontsize=18)
            if diffon==1:
                if wanted == 'phaseTwindz':
                    plt.ylabel(r"${\rm d}M_{\rm w}/{\rm d}z\;[{\rm M}_\odot/{\rm kpc}]$",fontsize=18)
                else:
                    plt.ylabel(r"${\rm d}M/{\rm d}z\;[{\rm M}_\odot/{\rm kpc}]$",fontsize=18)
            else:
                if wanted == 'phaseTwindz':
                    plt.ylabel(r"$M_{\rm w}(<z)\;[{\rm M}_\odot]$",fontsize=18)
                else:
                    plt.ylabel(r"$M(<z)\;[{\rm M}_\odot]$",fontsize=18)
                if wanted == 'phaseTwindz':
                        plt.legend(loc='best',frameon=False,fontsize=12)
                else:
                        plt.legend(loc=1,frameon=False,fontsize=12)
                plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                plt.tight_layout()
                print totalname
                plt.title(title,fontsize=18)
                plt.savefig(totalname)
                plt.clf()







if wanted=='Tbar' or wanted=='Twindbar':
        from crtestfunction import gaswindphase
        rcParams['figure.figsize'] = 5, 5
        diffon=0
    userad=1
    newlabelneed=0
    strlabelneed=0
    useRv=0
    if wanted=='Twindbar':
        vcut=30.0 #km/s
        filesuffix='wind'
    else:
        vcut=-1e10 #outflow velocity cut
        filesuffix=''
        zup=200.0
        zdown=20.0
        nooftimes=0
        nobin=2
        zl = np.linspace(zdown,zup,num=nobin)
        hl = []
    wl = []
    ll = []
    cl = []
    lablist=[]
        for runtodo in dirneed:
        hotmass = zl*0.0
        warmmass = zl*0.0
        lowmass = zl*0.0
        coldmass = zl*0.0
        info=outdirname(runtodo, Nsnap)
        dclabel=info['dclabel']
        cosmo=info['cosmo']
        firever=info['firever']
        halostr=info['halostr']
        h0=info['h0']
        maindir=info['maindir']
        Rvir=info['Rvir']
        rundir=info['rundir']
        newlabel=info['newlabel']
        runtitle=info['runtitle']
        if cosmo==1:
            userad=1
            useRv=1
        if userad==1:
            withinr=1e5
        else:
            withinr=20.0
        ptitle=title
        labelneed=dclabel
        if newlabelneed==1:
            labelneed="\n ".join(wrap(newlabel,29))
        if strlabelneed==1:
            labelneed="\n ".join(wrap(strlabel,40))
                lablist=np.append(lablist,labelneed)
        if runtitle=='SMC':
            ptitle='Dwarf'
        elif runtitle=='SBC':
            ptitle='Starburst'
        elif runtitle=='MW':
            ptitle=r'$L\star$ Galaxy'
                for i in range(startno,Nsnap,snapsep):
            print 'in the loop'
            if cosmo==1 and userad==1 and useRv==1:
                halosA = SF.read_halo_history_pep(rundir, i, singlesnap=1, firever=firever,halonostr=halostr, hubble=h0, comoving=0, maindir=maindir)
                Rvir = halosA['R']
                zup = Rvir
                zdown = Rvir*0.1
            if cosmo==0 and useRv==1:
                                zup = Rvir
                                zdown = Rvir*0.1
                        data=gaswindphase(runtodo,i,rotface=0,vcut=vcut,withinr=withinr,zup=zup,zdown=zdown,userad=userad)
            partX=data['partX']
            partY=data['partY']
            partZ=data['partZ']
            partp=partZ
            if userad==1:
                partp = np.sqrt(partX*partX+partY*partY+partZ*partZ)
                        TrueTemp=data['TrueTemp']
                        Gmass=data['Gmass']*1e10
                        for iz in range(nobin-1):
                cutz = (np.absolute(partp)>zl[iz]) & (np.absolute(partp)<zl[iz+1])
                cut=cutz
                                TrueTempcut = TrueTemp[cut]
                                Gmasscut = Gmass[cut]
                                warmcut = np.log10(TrueTempcut) < 5.3
                                lowcut = np.log10(TrueTempcut) < 4.7
                                coldcut = np.log10(TrueTempcut) < 4.0
                                hotmass[iz+1] += np.sum(Gmasscut)
                                warmmass[iz+1] += np.sum(Gmasscut[warmcut])
                                lowmass[iz+1] += np.sum(Gmasscut[lowcut])
                                coldmass[iz+1] += np.sum(Gmasscut[coldcut])
                        nooftimes+=1
        hotm = hotmass[-1]/nooftimes
        warmm = warmmass[-1]/nooftimes
        lowm = lowmass[-1]/nooftimes
        coldm = coldmass[-1]/nooftimes
        dz = (zup-zdown)/nobin
        hotm=hotm/dz
        warmm=warmm/dz
        lowm=lowm/dz
        coldm=coldm/dz
        perh = 1.
        perw = warmm/hotm
        perl = lowm/hotm
        perc = coldm/hotm
        hl = np.append(hl,perh)
        wl = np.append(wl,perw)
        ll = np.append(ll,perl)
        cl = np.append(cl,perc)
    ind = [x for x, _ in enumerate(dirneed)]
    plt.bar(ind, hl, width=0.8, color='r',label=r'$T\;[{\rm K}]>10^{5.3}$')
    plt.bar(ind, wl, width=0.8, color='y',label=r'$10^{4.7}<T\;[{\rm K}]<10^{5.3}$')
    plt.bar(ind, ll, width=0.8, color='b', label=r'$10^{4.0}<T\;[{\rm K}]<10^{4.7}$')
    plt.bar(ind, cl, width=0.8, color='m', label=r'$T\;[{\rm K}]<10^{4.0}$')
    plt.ylim(ymax=1.2)
    plt.legend(loc=2, fontsize=10,ncol=2,frameon=False)
    print 'lablist', lablist
    plt.xticks(ind, lablist, rotation='vertical',ha="left",fontsize=10)
    plt.tick_params(axis='both', which='both',direction='in',bottom=True,top=True,left=True,right=True,labelsize=10)
    if wanted=='Tbar':
        totalname = 'CRplot/Tbar/Tbar_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+filesuffix+'.pdf'
    if wanted=='Twindbar':
        totalname = 'CRplot/Twindbar/Twindbar_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+filesuffix+'.pdf'
    plt.ylabel("Mass Fraction",fontsize=14)
    plt.tight_layout()
    print totalname
    plt.title(ptitle,fontsize=18)
    plt.savefig(totalname)
    plt.clf()




if wanted=='Tztrack' or wanted=='rhoztrack' or wanted=='vzztrack' or wanted=='Tvztrack':
    rcParams['figure.figsize'] = 5, 4
    from crtestfunction import *
    from matplotlib.patches import FancyArrowPatch 
    def arrow(x,y,ax,n,color):
        d = len(x)//(n+1)    
        ind = np.arange(d,len(x),d)
        for i in ind:
        ar = FancyArrowPatch ((x[i-1],y[i-1]),(x[i],y[i]), 
                      arrowstyle='->', mutation_scale=20,color=color)
        ax.add_patch(ar)
    rcParams['figure.figsize'] = 5, 4
    userad=0
    snaptrack = Nsnap
    Tcut_t=1.0 #K
    highTcut_t = 1e3
    vcut_t=30.0 #outflow velocity cut
    vhcut_t = 1e10 #upper outflow velocity cut
    withinr_t=20.0
    zup_t=30.0
    zdown_t=20.0
    zup=zup_t
    zdown=0.
    if wanted=='Tztrack':
            extent = [zdown,zup,3.5,6.5]
    if wanted=='rhoztrack':
            extent = [zdown,zup,-7,4]
    if wanted=='vzztrack':
            extent = [zdown,zup,1.0,2.8]
    if wanted=='Tvztrack':
            extent = [1,3.2,3.5,7]
    fig, ax = plt.subplots()
        for runtodo in dirneed:
        try:
            info = outdirname(runtodo, Nsnap)
            runtitle=info['runtitle']
            ptitle=title
            ymin=extent[2]
            ymax = extent[3]
            if runtitle=='SMC':
                ptitle='Dwarf'
                if wanted=='Tztrack' or wanted=='Tvztrack':
                    ymax = 5.5; ymin = 3.5
                if wanted=='vzztrack':
                    ymax = 2.2; ymin = 1.0
            elif runtitle=='SBC':
                if wanted=='Tztrack' or wanted=='Tvztrack':
                    ymax = 6.0; ymin = 3.5
                if wanted=='vzztrack':
                    ymax = 2.5; ymin = 1.0
                ptitle='Starburst'
            elif runtitle=='MW':
                ptitle=r'$L\star$ Galaxy'
            xpoints=[]; ypoints=[];
            #ax.set_aspect("equal")
            for i in range(startno,Nsnap,snapsep):
                info = outdirname(runtodo, i)
                color=info['color']
                dclabel=info['dclabel']
                haveB=info['haveB']
                if haveB==1:
                    lsn='dashed'
                else:
                    lsn='solid'
                data = gastrack(runtodo,Nsnap,i,Tcut_t=Tcut_t,highTcut_t=highTcut_t,\
                    vcut_t=vcut_t,vhcut_t=vhcut_t,withinr_t=withinr_t,zup_t=zup_t,zdown_t=zdown_t,userad=userad)
                Gmass = data['Gmass']*1e10  
                if wanted=='Tvztrack':
                    vz = np.absolute(data['vz'])
                    TrueTemp = data['TrueTemp']
                    print 'vz', vz
                    print 'Gmass', Gmass
                    
                    xmed = np.log10(weighted_quantile(vz,50.,sample_weight=Gmass))
                    ymed = np.log10(weighted_quantile(TrueTemp,50.,sample_weight=Gmass))
                if wanted=='Tztrack':
                    TrueTemp = data['TrueTemp']
                    partZ = np.absolute(data['partZ'])
                    xmed = weighted_quantile(partZ,50.,sample_weight=Gmass)
                    ymed = np.log10(weighted_quantile(TrueTemp,50.,sample_weight=Gmass))
                if wanted=='vzztrack':
                    vz = np.absolute(data['vz'])
                    partZ = np.absolute(data['partZ'])
                    xmed = weighted_quantile(partZ,50.,sample_weight=Gmass)
                    ymed = np.log10(weighted_quantile(vz,50.,sample_weight=Gmass))
                if wanted=='rhoztrack':
                        convertedrho = data['convertedrho']
                        partZ = np.absolute(data['partZ'])
                        xmed = weighted_quantile(partZ,50.,sample_weight=Gmass)
                        ymed = np.log10(weighted_quantile(convertedrho,50.,sample_weight=Gmass))
                xpoints = np.append(xpoints,xmed)
                ypoints = np.append(ypoints,ymed)
            ax.plot(xpoints, ypoints,color=color,ls=lsn,label=dclabel)
            ax.scatter(xpoints[:-1],ypoints[:-1],color=color)
            ax.scatter(xpoints[-1],ypoints[-1],color=color,marker='>')
        except ValueError:
            continue
        #if wanted=='Tvztrack':
        #   arrow(xpoints,ypoints,ax,3,color)
    plt.xlim(xmin=extent[0],xmax=extent[1])
    plt.ylim(ymin=ymin,ymax=ymax)
    if wanted=='Tztrack' or wanted=='Tvztrack':
        plt.ylabel(r"$\mathrm{Log}_{10} (T {\rm [K]})$",fontsize=18)
    if userad==1:
        plt.xlabel(r"$r {\rm [kpc]}$",fontsize=18)
    else:
        plt.xlabel(r"$z {\rm [kpc]}$",fontsize=18)
    if wanted=='Tvztrack':
        if userad==1:
            plt.xlabel(r"$\mathrm{Log}_{10} (v_r {\rm [km/s]})$",fontsize=18)
        else:
            plt.xlabel(r"$\mathrm{Log}_{10} (v_z {\rm [km/s]})$",fontsize=18)   
    if wanted=='rhoztrack':
        plt.ylabel(r"$\mathrm{Log}_{10} (n {\rm [cm^{-3}]})$",fontsize=18)
    if wanted=='Tztrack' or wanted=='Tvztrack':
        plt.ylabel(r"$\mathrm{Log}_{10} (T {\rm [K]})$",fontsize=18)
        if runtitle=='MW':
            plt.legend(loc='best', fontsize=9, ncol=3,scatterpoints=1)
    if wanted=='vzztrack':
        if userad==1:
            plt.ylabel(r"$\mathrm{Log}_{10} (v_r {\rm [km/s]})$",fontsize=18)
        else:
            plt.ylabel(r"$\mathrm{Log}_{10} (v_z {\rm [km/s]})$",fontsize=18)
    if wanted == 'Tvztrack':
        totalname = 'CRplot/Tvztrack/Tvztrack_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    if wanted == 'Tztrack':
        totalname = 'CRplot/Tztrack/Tztrack_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    if wanted == 'rhoztrack':
        totalname = 'CRplot/rhoztrack/rhoztrack_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    if wanted == 'vzztrack':
        totalname = 'CRplot/vzztrack/vzztrack_'+runtodo+'_sn'+str(startno)+'_'+str(Nsnap)+'.pdf'
    plt.title(ptitle,fontsize=18)
    plt.tight_layout()
    print totalname
    plt.savefig(totalname)
    plt.clf()

