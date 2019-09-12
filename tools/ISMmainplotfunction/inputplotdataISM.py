import collections
from crdenmidplanetestinput import crdenmidplanetestinput
from nHmidplanetestinput import nHmidplanetestinput
from eturtimetestinput import eturtimetestinput
from samson_functions import checkkey
from rhoTgridtestinput import rhoTgridtestinput
from CRerzdistestinput import CRerzdistestinput
from CRerzdis_volave_testinput import CRerzdis_volave_testinput
from CRerzdis_single_testinput import CRerzdis_single_testinput
from rhogzdistestinput import rhogzdistestinput
from star_midplane_outtestinput import star_midplane_outtestinput
from kslawtestinput import kslawtestinput
from synchrotrontestinput import synchrotrontestinput
from Faradayrotationtestinput import Faradayrotationtestinput
from vzcylouttestinput import vzcylouttestinput
from hgashalftestinput import hgashalftestinput
from sfrturtestinput import sfrturtestinput
from sfrtestinput import sfrtestinput
from rhogzdispartestinput import rhogzdispartestinput
from gasdenTtestinput import gasdenTtestinput
from gasdenThistestinput import gasdenThistestinput
from outputsnipshot import outputsnipshot
from generatepressure import generatepressure
from generategfield import generategfield
from crdenfromdata_testinput import crdenfromdata_testinput
from sfrmockturtestinput import sfrmockturtestinput
from Blitzmidplanepressure_testinput import Blitzmidplanepressure_testinput
from turdifscale_testinput import turdifscale_testinput
from sfrhotgasmass_testinput import sfrhotgasmass_testinput
from Tztrack_outtestinput import Tztrack_outtestinput



nested_dict = lambda: collections.defaultdict(nested_dict)
inputplotdict = nested_dict()



inputplotdict['figTztrackstart_m12i_580']={
    '_plotfunction':Tztrack_outtestinput,
    'dirneed':[['m12imhdcvhr','m12icr_700hr']],
    'wanted':'Tztrack',
    'startno':580,
    'Nsnap':600,
    'snapsep':2,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12i_trackstart_590',
    'Tcut_t':1e5,
    'highTcut_t':1e10,
    'vcut_t':2e2,
    'vhcut_t':1e10,
    'withinr_t':10.0,
    'zup_t':2.0,
    'zdown_t':1.0,
    'zup':0.01,
    'zdown':1e2,
    'trackstart':1    
    }


inputplotdict['figTztrackstart_m12i']={
    '_plotfunction':Tztrack_outtestinput,
    'dirneed':[['m12imhdcvhr','m12icr_700hr']],
    'wanted':'Tztrack',
    'startno':510,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12i_trackstart',
    'Tcut_t':1e5,
    'highTcut_t':1e10,
    'vcut_t':2e2,
    'vhcut_t':1e10,
    'withinr_t':10.0,
    'zup_t':2.0,
    'zdown_t':1.0,
    'zup':0.01,
    'zdown':1e2,
    'trackstart':1    
    }


inputplotdict['figTztrack_m12i']={
    '_plotfunction':Tztrack_outtestinput,
    'dirneed':[['m12imhdcvhr','m12icr_700hr']],
    'wanted':'Tztrack',
    'startno':510,
    'Nsnap':600,
    'snapsep':90,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12i_mult',
    'Tcut_t':1e5,
    'highTcut_t':1e10,
    'vcut_t':-1e10,
    'vhcut_t':1e10,
    'withinr_t':10.0,
    'zup_t':10.0,
    'zdown_t':5.0,
    'zup':1.0,
    'zdown':1e2,
    'hotwarmmode':1    
    }

inputplotdict['figTztrack_m12i_580']={
    '_plotfunction':Tztrack_outtestinput,
    'dirneed':[['m12imhdcvhr','m12icr_700hr']],
    'wanted':'Tztrack',
    'startno':580,
    'Nsnap':600,
    'snapsep':5,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12i_mult',
    'Tcut_t':1e5,
    'highTcut_t':1e10,
    'vcut_t':-1e10,
    'vhcut_t':1e10,
    'withinr_t':10.0,
    'zup_t':10.0,
    'zdown_t':5.0,
    'zup':1.0,
    'zdown':1e2,
    'hotwarmmode':1    
    }





inputplotdict['figer_m12ihr0_25kpc']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':[['m12imhdcvhr'],['m12icr_700hr']],
    'wanted':'er',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12grid0_25kpc_r10kpc',
    'withinr':10.0,
    'maxlength':1.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':0,
    'outHI':0,
    'griddir':'grid0_25kpc'
    }

inputplotdict['figrealpr_m12ihr0_25kpc']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':[['m12imhdcvhr'],['m12icr_700hr']],
    'wanted':'realpr',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12grid0_25kpc_r10kpc',
    'withinr':10.0,
    'maxlength':1.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':0,
    'outHI':0,
    'griddir':'grid0_25kpc'
    }


inputplotdict['gf11bgrid0_5kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m11bmhdcv','m11bcr_700'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40
    }

inputplotdict['gf11fgrid0_5kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m11fmhdcv','m11fcr_700'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40
    }

inputplotdict['gf11hgrid0_5kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m11hmhdcv','m11hcr_700'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40
    }

inputplotdict['gf11dgrid0_5kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m11dmhdcv','m11dcr_700'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40
    }

inputplotdict['gf11ggrid0_5kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m11gmhdcv','m11gcr_700'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40
    }


inputplotdict['gpm11bgrid0_5kpcHI']={
    '_plotfunction':generatepressure,
    'dirneed':['m11bmhdcv','m11bcr_700'],
    'wanted':'',
    'startno':580,
    'Nsnap':601,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40,
    'outHI':1
    }

inputplotdict['gpm11fgrid0_5kpcHI']={
    '_plotfunction':generatepressure,
    'dirneed':['m11fmhdcv','m11fcr_700'],
    'wanted':'',
    'startno':580,
    'Nsnap':601,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40,
    'outHI':1
    }

inputplotdict['gpm11hgrid0_5kpcHI']={
    '_plotfunction':generatepressure,
    'dirneed':['m11hmhdcv','m11hcr_700'],
    'wanted':'',
    'startno':580,
    'Nsnap':601,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40,
    'outHI':1
    }

inputplotdict['gpm11dgrid0_5kpcHI']={
    '_plotfunction':generatepressure,
    'dirneed':['m11dmhdcv','m11dcr_700'],
    'wanted':'',
    'startno':580,
    'Nsnap':601,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40,
    'outHI':1
    }


inputplotdict['gpm11ggrid0_5kpcHI']={
    '_plotfunction':generatepressure,
    'dirneed':['m11gmhdcv','m11gcr_700'],
    'wanted':'',
    'startno':580,
    'Nsnap':601,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40,
    'outHI':1
    }


inputplotdict['gf12imhdcvgrid0_125kpc600']={
    '_plotfunction':generategfield,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':600,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gf12icr_700grid0_125kpc600']={
    '_plotfunction':generategfield,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':600,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gpm12imhdcvgrid0_125kpc600']={
    '_plotfunction':generatepressure,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':600,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gpm12imhdcvgrid0_125kpc600HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':600,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160,
    'outHI':1
    }


inputplotdict['gpm12icr_700grid0_125kpc600']={
    '_plotfunction':generatepressure,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':600,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }

inputplotdict['gpm12icr_700grid0_125kpc600HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':600,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160,
    'outHI':1
    }


inputplotdict['gf12fmhdcvgrid0_125kpc600']={
    '_plotfunction':generategfield,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':600,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gf12fcr_700grid0_125kpc600']={
    '_plotfunction':generategfield,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':600,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gpm12fmhdcvgrid0_125kpc600']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':600,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }

inputplotdict['gpm12fmhdcvgrid0_125kpc600HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':600,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160,
    'outHI':1
    }



inputplotdict['gpm12fcr_700grid0_125kpc600']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':600,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12f',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }

inputplotdict['gpm12fcr_700grid0_125kpc600HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':600,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12f',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160,
    'outHI':1
    }


inputplotdict['gf12mmhdcvgrid0_125kpc600']={
    '_plotfunction':generategfield,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':600,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gf12mcr_700grid0_125kpc600']={
    '_plotfunction':generategfield,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':600,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gpm12mmhdcvgrid0_125kpc600']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':600,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gpm12mmhdcvgrid0_125kpc600HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':600,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160,
    'outHI':1
    }


inputplotdict['gpm12mcr_700grid0_125kpc600']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':600,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12m',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }

inputplotdict['gpm12mcr_700grid0_125kpc600HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':600,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12m',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160,
    'outHI':1
    }


inputplotdict['gf12imhdcvgrid0_125kpc590']={
    '_plotfunction':generategfield,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':590,
    'Nsnap':590,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gf12icr_700grid0_125kpc590']={
    '_plotfunction':generategfield,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':590,
    'Nsnap':590,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gpm12imhdcvgrid0_125kpc590']={
    '_plotfunction':generatepressure,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':590,
    'Nsnap':590,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }

inputplotdict['gpm12imhdcvgrid0_125kpc590HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':590,
    'Nsnap':590,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160,
    'outHI':1
    }


inputplotdict['gpm12icr_700grid0_125kpc590']={
    '_plotfunction':generatepressure,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':590,
    'Nsnap':590,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }

inputplotdict['gpm12icr_700grid0_125kpc590HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':590,
    'Nsnap':590,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160,
    'outHI':1
    }


inputplotdict['gf12fmhdcvgrid0_125kpc590']={
    '_plotfunction':generategfield,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':590,
    'Nsnap':590,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gf12fcr_700grid0_125kpc590']={
    '_plotfunction':generategfield,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':590,
    'Nsnap':590,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gpm12fmhdcvgrid0_125kpc590']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':590,
    'Nsnap':590,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }

inputplotdict['gpm12fmhdcvgrid0_125kpc590HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':590,
    'Nsnap':590,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160,
    'outHI':1
    }


inputplotdict['gpm12fcr_700grid0_125kpc590']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':590,
    'Nsnap':590,
    'snapsep':10,
    'fmeat':'m12f',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }

inputplotdict['gpm12fcr_700grid0_125kpc590HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':590,
    'Nsnap':590,
    'snapsep':10,
    'fmeat':'m12f',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160,
    'outHI1':1
    }


inputplotdict['gf12mmhdcvgrid0_125kpc590']={
    '_plotfunction':generategfield,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':590,
    'Nsnap':590,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gf12mcr_700grid0_125kpc590']={
    '_plotfunction':generategfield,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':590,
    'Nsnap':590,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gpm12mmhdcvgrid0_125kpc590']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':590,
    'Nsnap':590,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }

inputplotdict['gpm12mmhdcvgrid0_125kpc590HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':590,
    'Nsnap':590,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160,
    'outHI':1
    }


inputplotdict['gpm12mcr_700grid0_125kpc590']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':590,
    'Nsnap':590,
    'snapsep':10,
    'fmeat':'m12m',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }

inputplotdict['gpm12mcr_700grid0_125kpc590HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':590,
    'Nsnap':590,
    'snapsep':10,
    'fmeat':'m12m',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160,
    'outHI':1
    }



inputplotdict['gf12imhdcvgrid0_125kpc580']={
    '_plotfunction':generategfield,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gf12icr_700grid0_125kpc580']={
    '_plotfunction':generategfield,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gpm12imhdcvgrid0_125kpc580']={
    '_plotfunction':generatepressure,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }

inputplotdict['gpm12imhdcvgrid0_125kpc580HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160,
    'outHI':1
    }


inputplotdict['gpm12icr_700grid0_125kpc580']={
    '_plotfunction':generatepressure,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }

inputplotdict['gpm12icr_700grid0_125kpc580HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160,
    'outHI':1
    }


inputplotdict['gf12fmhdcvgrid0_125kpc580']={
    '_plotfunction':generategfield,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gf12fcr_700grid0_125kpc580']={
    '_plotfunction':generategfield,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gpm12fmhdcvgrid0_125kpc580']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }

inputplotdict['gpm12fmhdcvgrid0_125kpc580HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160,
    'outHI':1
    }


inputplotdict['gpm12fcr_700grid0_125kpc580']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'fmeat':'m12f',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }

inputplotdict['gpm12fcr_700grid0_125kpc580HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'fmeat':'m12f',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160,
    'outHI':1
    }


inputplotdict['gf12mmhdcvgrid0_125kpc580']={
    '_plotfunction':generategfield,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gf12mcr_700grid0_125kpc580']={
    '_plotfunction':generategfield,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gpm12mmhdcvgrid0_125kpc580']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gpm12mmhdcvgrid0_125kpc580HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160,
    'outHI':1
    }


inputplotdict['gpm12mcr_700grid0_125kpc580']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'fmeat':'m12m',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gpm12mcr_700grid0_125kpc580HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'fmeat':'m12m',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160,
    'outHI':1
    }


inputplotdict['gf12imhdcvgrid0_125kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gf12icr_700grid0_125kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gpm12imhdcvgrid0_125kpc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gpm12icr_700grid0_125kpc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gf12fmhdcvgrid0_125kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gf12fcr_700grid0_125kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gpm12fmhdcvgrid0_125kpc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gpm12fcr_700grid0_125kpc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'fmeat':'m12f',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gf12mmhdcvgrid0_125kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }




inputplotdict['gf12mcr_700grid0_125kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gpm12mmhdcvgrid0_125kpc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_125kpc',
    'parallel':1,
    'nogrid':160
    }


inputplotdict['gpm12mcr_700grid0_1kpc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'fmeat':'m12m',
    'maxlength':20.0,
    'griddir':'grid0_1kpc',
    'parallel':1,
    'nogrid':200
    }



inputplotdict['gpm12fmhdcvgrid1kpccc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12f',
    'maxlength':20.0,
    'griddir':'grid1kpc/cutcold',
    'nogrid':21,
    'parallel':1,
    'cutcold':1
    }


inputplotdict['gpm12fcr_700grid1kpccc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12f',
    'maxlength':20.0,
    'griddir':'grid1kpc/cutcold',
    'nogrid':21,
    'parallel':1,
    'cutcold':1
    }


inputplotdict['gpm12mmhdcvgrid1kpccc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12m',
    'maxlength':20.0,
    'griddir':'grid1kpc/cutcold',
    'nogrid':21,
    'parallel':1,
    'cutcold':1
    }


inputplotdict['gpm12mcr_700grid1kpccc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12m',
    'maxlength':20.0,
    'griddir':'grid1kpc/cutcold',
    'nogrid':21,
    'parallel':1,
    'cutcold':1
    }



inputplotdict['gpm12imhdcvgrid1kpccc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid1kpc/cutcold',
    'nogrid':21,
    'parallel':1,
    'cutcold':1
    }


inputplotdict['gpm12fmhdcvgrid1kpccc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12f',
    'maxlength':20.0,
    'griddir':'grid1kpc/cutcold',
    'nogrid':21,
    'parallel':1,
    'cutcold':1
    }


inputplotdict['gpm12fcr_700grid1kpccc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12f',
    'maxlength':20.0,
    'griddir':'grid1kpc/cutcold',
    'nogrid':21,
    'parallel':1,
    'cutcold':1
    }


inputplotdict['gpm12mmhdcvgrid1kpccc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12m',
    'maxlength':20.0,
    'griddir':'grid1kpc/cutcold',
    'nogrid':21,
    'parallel':1,
    'cutcold':1
    }


inputplotdict['gpm12mcr_700grid1kpccc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12m',
    'maxlength':20.0,
    'griddir':'grid1kpc/cutcold',
    'nogrid':21,
    'parallel':1,
    'cutcold':1
    }



inputplotdict['gpm12imhdcvgrid2kpccc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid2kpc/cutcold',
    'nogrid':10,
    'cutcold':1,
    'parallel':1
    }


inputplotdict['gpm12icr_700grid2kpccc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid2kpc/cutcold',
    'nogrid':10,
    'parallel':1,
    'cutcold':1
    }


inputplotdict['gf12imhdcvgrid2kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid2kpc',
    'parallel':1,
    'nogrid':10
    }


inputplotdict['gf12icr_700grid2kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid2kpc',
    'parallel':1,
    'nogrid':10
    }


inputplotdict['gpm12imhdcvgrid2kpc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid2kpc',
    'parallel':1,
    'nogrid':10
    }


inputplotdict['gpm12icr_700grid2kpc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid2kpc',
    'parallel':1,
    'nogrid':10
    }


inputplotdict['gf12fmhdcvgrid2kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid2kpc',
    'parallel':1,
    'nogrid':10
    }


inputplotdict['gf12fcr_700grid2kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid2kpc',
    'parallel':1,
    'nogrid':10
    }


inputplotdict['gpm12fmhdcvgrid2kpc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid2kpc',
    'parallel':1,
    'nogrid':10
    }


inputplotdict['gpm12fcr_700grid2kpc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12f',
    'maxlength':20.0,
    'griddir':'grid2kpc',
    'parallel':1,
    'nogrid':10
    }


inputplotdict['gf12mmhdcvgrid2kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid2kpc',
    'parallel':1,
    'nogrid':10
    }


inputplotdict['gf12mcr_700grid2kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid2kpc',
    'parallel':1,
    'nogrid':10
    }


inputplotdict['gpm12mmhdcvgrid2kpc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid2kpc',
    'parallel':1,
    'nogrid':10
    }

inputplotdict['gpm12mmhdcvgrid2kpcHI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid2kpc',
    'parallel':1,
    'nogrid':10,
    'outHI':1
    }

inputplotdict['gpm12mcr_700grid2kpcHI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid2kpc',
    'parallel':1,
    'nogrid':10,
    'outHI':1
    }



inputplotdict['gpm12mcr_700grid2kpc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12m',
    'maxlength':20.0,
    'griddir':'grid2kpc',
    'parallel':1,
    'nogrid':10
    }

inputplotdict['gf12imhdcvgrid0_25kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80
    }

inputplotdict['gf12imhdcvgrid0_25kpc598']={
    '_plotfunction':generategfield,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':598,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80
    }

inputplotdict['gf12imhdcvgrid0_25kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80
    }


inputplotdict['gf12icr_700grid0_25kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80
    }


inputplotdict['gpm12imhdcvgrid0_25kpc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80
    }

inputplotdict['gpm12imhdcvgrid0_25kpcHI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80,
    'outHI':1
    }


inputplotdict['gpm12imhdcvgrid0_25kpc590HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':590,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80,
    'outHI':1
    }


inputplotdict['gpm12icr_700grid0_25kpc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80
    }


inputplotdict['gpm12icr_700grid0_25kpcHI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80,
    'outHI':1
    }


inputplotdict['gpm12icr_700grid0_25kpc590HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':590,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80,
    'outHI':1
    }



inputplotdict['gf12fmhdcvgrid0_25kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80
    }


inputplotdict['gf12fmhdcvgrid0_25kpc590']={
    '_plotfunction':generategfield,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':590,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80
    }

inputplotdict['gf12fmhdcvgrid0_25kpc595']={
    '_plotfunction':generategfield,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':595,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80
    }


inputplotdict['gf12fcr_700grid0_25kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80
    }


inputplotdict['gf12fcr_700grid0_25kpc590']={
    '_plotfunction':generategfield,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':590,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80
    }

inputplotdict['gf12fcr_700grid0_25kpc595']={
    '_plotfunction':generategfield,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':595,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80
    }

inputplotdict['gpm12fmhdcvgrid0_25kpc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80
    }


inputplotdict['gpm12fmhdcvgrid0_25kpc580HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':586,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80,
    'outHI':1
    }

inputplotdict['gpm12fmhdcvgrid0_25kpc585HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':585,
    'Nsnap':591,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80,
    'outHI':1
    }


inputplotdict['gpm12fmhdcvgrid0_25kpc590HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':590,
    'Nsnap':596,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80,
    'outHI':1
    }

inputplotdict['gpm12fmhdcvgrid0_25kpc595HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':595,
    'Nsnap':601,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80,
    'outHI':1
    }




inputplotdict['gpm12fcr_700grid0_25kpc580HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':586,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80,
    'outHI':1
    }



inputplotdict['gpm12fcr_700grid0_25kpc585HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':585,
    'Nsnap':591,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80,
    'outHI':1
    }

inputplotdict['gpm12fcr_700grid0_25kpc590HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':590,
    'Nsnap':596,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80,
    'outHI':1
    }


inputplotdict['gpm12fcr_700grid0_25kpc595HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':595,
    'Nsnap':601,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80,
    'outHI':1
    }

inputplotdict['gf12mmhdcvgrid0_25kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80
    }


inputplotdict['gf12mmhdcvgrid0_25kpc590']={
    '_plotfunction':generategfield,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':590,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80
    }

inputplotdict['gf12mcr_700grid0_25kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80
    }


inputplotdict['gf12mcr_700grid0_25kpc590']={
    '_plotfunction':generategfield,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':590,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80
    }

inputplotdict['gf12mcr_700grid0_25kpc595']={
    '_plotfunction':generategfield,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':595,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80
    }

inputplotdict['gf12mmhdcvgrid0_25kpc593']={
    '_plotfunction':generategfield,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':593,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80
    }

inputplotdict['gf12mmhdcvgrid0_25kpc596']={
    '_plotfunction':generategfield,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':596,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80
    }


inputplotdict['gpm12mmhdcvgrid0_25kpc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80
    }


inputplotdict['gpm12mmhdcvgrid0_25kpc580HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':586,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80,
    'outHI':1
    }


inputplotdict['gpm12mmhdcvgrid0_25kpc585HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':585,
    'Nsnap':591,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80,
    'outHI':1
    }


inputplotdict['gpm12mmhdcvgrid0_25kpc590HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':590,
    'Nsnap':596,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80,
    'outHI':1
    }


inputplotdict['gpm12mmhdcvgrid0_25kpc595HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':595,
    'Nsnap':601,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80,
    'outHI':1
    }


inputplotdict['gpm12mcr_700grid0_25kpc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80
    }


inputplotdict['gpm12mcr_700grid0_25kpc580HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':586,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80,
    'outHI':1
    }


inputplotdict['gpm12mcr_700grid0_25kpc585HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':585,
    'Nsnap':591,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80,
    'outHI':1
    }


inputplotdict['gpm12mcr_700grid0_25kpc590HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':590,
    'Nsnap':596,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80,
    'outHI':1
    }


inputplotdict['gpm12mcr_700grid0_25kpc595HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':595,
    'Nsnap':601,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_25kpc',
    'parallel':1,
    'nogrid':80,
    'outHI':1
    }


inputplotdict['gf12imhdcvgrid0_5kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40
    }


inputplotdict['gf12imhdcvgrid0_5kpc590']={
    '_plotfunction':generategfield,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':590,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40
    }


inputplotdict['gf12icr_700grid0_5kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40
    }


inputplotdict['gf12icr_700grid0_5kpc590']={
    '_plotfunction':generategfield,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':590,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40
    }


inputplotdict['gpm12imhdcvgrid0_5kpc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40
    }

inputplotdict['gpm12imhdcvgrid0_5kpcHI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40,
    'outHI':1
    }


inputplotdict['gpm12icr_700grid0_5kpc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40
    }


inputplotdict['gpm12icr_700grid0_5kpcHI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40,
    'outHI':1
    }


inputplotdict['gf12fmhdcvgrid0_5kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40
    }


inputplotdict['gf12fcr_700grid0_5kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40
    }


inputplotdict['gpm12fmhdcvgrid0_5kpc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40
    }

inputplotdict['gpm12fmhdcvgrid0_5kpcHI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40,
    'outHI':1
    }


inputplotdict['gpm12fcr_700grid0_5kpc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40
    }

inputplotdict['gpm12fcr_700grid0_5kpcHI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40,
    'outHI':1
    }


inputplotdict['gf12mmhdcvgrid0_5kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40
    }


inputplotdict['gf12mcr_700grid0_5kpc']={
    '_plotfunction':generategfield,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40
    }


inputplotdict['gpm12mmhdcvgrid0_5kpc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40
    }


inputplotdict['gpm12mmhdcvgrid0_5kpcHI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40,
    'outHI':1
    }

inputplotdict['gpm12mcr_700grid0_5kpc']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40
    }


inputplotdict['gpm12mcr_700grid0_5kpcHI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'griddir':'grid0_5kpc',
    'parallel':1,
    'nogrid':40,
    'outHI':1
    }

inputplotdict['gf12imhdcv']={
    '_plotfunction':generategfield,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'parallel':1,
    'nogrid':21
    }


inputplotdict['gf12icr_700']={
    '_plotfunction':generategfield,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'parallel':1,
    'nogrid':21
    }


inputplotdict['gpm12imhdcv']={
    '_plotfunction':generatepressure,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'parallel':1,
    'nogrid':21
    }


inputplotdict['gpm12imhdcvHI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'parallel':1,
    'nogrid':21,
    'outHI':1
    }



inputplotdict['gpm12icr_700']={
    '_plotfunction':generatepressure,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'parallel':1,
    'nogrid':21
    }

inputplotdict['gpm12icr_700HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12icr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'parallel':1,
    'nogrid':21,
    'outHI':1
    }


inputplotdict['gf12fmhdcv']={
    '_plotfunction':generategfield,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'parallel':1,
    'nogrid':21
    }


inputplotdict['gf12fcr_700']={
    '_plotfunction':generategfield,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'parallel':1,
    'nogrid':21
    }


inputplotdict['gpm12fmhdcv']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'parallel':1,
    'nogrid':21
    }

inputplotdict['gpm12fmhdcvHI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'parallel':1,
    'nogrid':21,
    'outHI':1
    }


inputplotdict['gpm12fcr_700']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'parallel':1,
    'nogrid':21
    }

inputplotdict['gpm12fcr_700HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'parallel':1,
    'nogrid':21,
    'outHI':1
    }


inputplotdict['gf12mmhdcv']={
    '_plotfunction':generategfield,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'parallel':1,
    'nogrid':21
    }


inputplotdict['gf12mcr_700']={
    '_plotfunction':generategfield,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'parallel':1,
    'nogrid':21
    }


inputplotdict['gpm12mmhdcv']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'parallel':1,
    'nogrid':21
    }

inputplotdict['gpm12mmhdcvHI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'parallel':1,
    'nogrid':21,
    'outHI':1
    }

inputplotdict['gpm12mmhdcvHI10']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mmhdcvhr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'parallel':1,
    'nogrid':21,
    'outHI':1
    }


inputplotdict['gpm12mcr_700']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'parallel':1,
    'nogrid':21
    }



inputplotdict['gpm12mcr_700HI']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'fmeat':'m12i',
    'maxlength':20.0,
    'parallel':1,
    'nogrid':21,
    'outHI':1
    }

inputplotdict['gpm12mcr_700HI10']={
    '_plotfunction':generatepressure,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'fmeat':'m12i',
    'maxlength':20.0,
    'parallel':1,
    'nogrid':21,
    'outHI':1
    }

inputplotdict['outssm12i580']={
    '_plotfunction':outputsnipshot,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':591,
    'snapsep':1,
    'fmeat':'m12i'
    }

inputplotdict['outssm12i590']={
    '_plotfunction':outputsnipshot,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'',
    'startno':590,
    'Nsnap':601,
    'snapsep':1,
    'fmeat':'m12i'
    }


inputplotdict['outssm12f580']={
    '_plotfunction':outputsnipshot,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':591,
    'snapsep':1,
    'fmeat':'m12f'
    }


inputplotdict['outssm12fcr_700590']={
    '_plotfunction':outputsnipshot,
    'dirneed':['m12fcr_700hr'],
    'wanted':'',
    'startno':590,
    'Nsnap':601,
    'snapsep':1,
    'fmeat':'m12f'
    }

inputplotdict['outssm12f590']={
    '_plotfunction':outputsnipshot,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'',
    'startno':590,
    'Nsnap':601,
    'snapsep':1,
    'fmeat':'m12f'
    }


inputplotdict['outssm12mcr_700580']={
    '_plotfunction':outputsnipshot,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':591,
    'snapsep':1,
    'fmeat':'m12m'
    }

inputplotdict['outssm12mcr_700590']={
    '_plotfunction':outputsnipshot,
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':590,
    'Nsnap':601,
    'snapsep':1,
    'fmeat':'m12m'
    }

inputplotdict['outssm12m580']={
    '_plotfunction':outputsnipshot,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':591,
    'snapsep':1,
    'fmeat':'m12m'
    }


inputplotdict['outssm12m590']={
    '_plotfunction':outputsnipshot,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'',
    'startno':590,
    'Nsnap':601,
    'snapsep':1,
    'fmeat':'m12m'
    }

inputplotdict['outssm11g580']={
    '_plotfunction':outputsnipshot,
    'dirneed':['m11gmhdcv','m11gcr_700'],
    'wanted':'',
    'startno':580,
    'Nsnap':591,
    'snapsep':10,
    'fmeat':'m11g'
    }


inputplotdict['outssm11g590']={
    '_plotfunction':outputsnipshot,
    'dirneed':['m11gmhdcv','m11gcr_700'],
    'wanted':'',
    'startno':590,
    'Nsnap':601,
    'snapsep':10,
    'fmeat':'m11g'
    }


inputplotdict['outssm11d580']={
    '_plotfunction':outputsnipshot,
    'dirneed':['m11dmhdcv','m11dcr_700'],
    'wanted':'',
    'startno':580,
    'Nsnap':591,
    'snapsep':10,
    'fmeat':'m11d'
    }


inputplotdict['outssm11d590']={
    '_plotfunction':outputsnipshot,
    'dirneed':['m11dmhdcv','m11dcr_700'],
    'wanted':'',
    'startno':590,
    'Nsnap':601,
    'snapsep':10,
    'fmeat':'m11g'
    }

inputplotdict['outssm11h580']={
    '_plotfunction':outputsnipshot,
    'dirneed':['m11hmhdcv','m11hcr_700'],
    'wanted':'',
    'startno':580,
    'Nsnap':591,
    'snapsep':10,
    'fmeat':'m11h'
    }


inputplotdict['outssm11h590']={
    '_plotfunction':outputsnipshot,
    'dirneed':['m11hmhdcv','m11hcr_700'],
    'wanted':'',
    'startno':590,
    'Nsnap':601,
    'snapsep':10,
    'fmeat':'m11h'
    }

inputplotdict['outssm11f580']={
    '_plotfunction':outputsnipshot,
    'dirneed':['m11fmhdcv','m11fcr_700'],
    'wanted':'',
    'startno':580,
    'Nsnap':591,
    'snapsep':10,
    'fmeat':'m11f'
    }


inputplotdict['outssm11f590']={
    '_plotfunction':outputsnipshot,
    'dirneed':['m11fmhdcv','m11fcr_700'],
    'wanted':'',
    'startno':590,
    'Nsnap':601,
    'snapsep':10,
    'fmeat':'m11f'
    }

inputplotdict['outssm11b580']={
    '_plotfunction':outputsnipshot,
    'dirneed':['m11bmhdcv','m11bcr_700'],
    'wanted':'',
    'startno':580,
    'Nsnap':591,
    'snapsep':10,
    'fmeat':'m11b'
    }


inputplotdict['outssm11b590']={
    '_plotfunction':outputsnipshot,
    'dirneed':['m11bmhdcv','m11bcr_700'],
    'wanted':'',
    'startno':590,
    'Nsnap':601,
    'snapsep':10,
    'fmeat':'m11b'
    }


inputplotdict['outssm12']={
    '_plotfunction':outputsnipshot,
    #'dirneed':['m12imhdcvhr','m12icr_700hr','m12fmhdcvhr','m12fcr_700hr','m12mmhdcvhr','m12mcr_700hr'],
    'dirneed':['m12mcr_700hr'],
    'wanted':'',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'fmeat':'m12'
    }

inputplotdict['outssm11m12noCR']={
    '_plotfunction':outputsnipshot,
    'dirneed':['fm11b','f46','fm11h','fm11d','f61','fm12i','fm12f','fm12m'],
    'wanted':'',
    'startno':580,
    'Nsnap':601,
    'snapsep':10,
    'fmeat':'m12'
    }

inputplotdict['outssm12test']={
    '_plotfunction':outputsnipshot,
    'dirneed':['m12imhdcvhr'],
    'wanted':'',
    'startno':590,
    'Nsnap':591,
    'snapsep':10,
    'fmeat':'m12'
    }

inputplotdict['figrhogfrompardist_m11ghr0kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m11gmhdcv','m11gcr_700'],
    'wanted':'dpdz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fcr700hr0kpc',
    'withoutr':0.0,
    'withinr':0.0,
    'maxlength':10.0,
    'nogrid':10,
    'usehalfz':1
    }


inputplotdict['figrhogfrompardist_m12fhr0kpchighz']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'dpdz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fcr700hr0kpchighz',
    'withoutr':0.0,
    'withinr':0.0,
    'maxlength':100.0,
    'nogrid':10,
    'usehalfz':1
    }

inputplotdict['figrhogfrompardist_m12fhr0kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'dpdz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fcr700hr0kpc',
    'withoutr':0.0,
    'withinr':0.0,
    'maxlength':10.0,
    'nogrid':10,
    'usehalfz':1
    }


inputplotdict['figintrhogfrompardist_m12fhr10kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fcr700hr10kpc',
    'withoutr':9.0,
    'withinr':10.0,
    'maxlength':10.0,
    'nogrid':10,
    'usehalfz':1
    }


inputplotdict['figrhogfrompardist_m12fhr10kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'dpdz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fcr700hr10kpc',
    'withoutr':9.0,
    'withinr':11.0,
    'maxlength':10.0,
    'nogrid':10,
    'usehalfz':1
    }


inputplotdict['figrhogfrompardist_m12mhr0kpchighz']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'dpdz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mcr700hr0kpchighz',
    'withoutr':0.0,
    'withinr':0.0,
    'maxlength':100.0,
    'nogrid':10,
    'usehalfz':1
    }


inputplotdict['figrhogfrompardist_m12mhr0kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'dpdz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mcr700hr0kpc',
    'withoutr':0.0,
    'withinr':0.0,
    'maxlength':10.0,
    'nogrid':10,
    'usehalfz':1
    }


inputplotdict['figrhogfrompardist_m12mhr']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'dpdz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mcr700hr7kpc',
    'withoutr':6.0,
    'withinr':8.0,
    'maxlength':10.0,
    'nogrid':10,
    'usehalfz':1
    }


inputplotdict['figintrhogfrompardist_m12mhr10kpcz10kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'pz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mcr700hr10kpcz10kpc',
    'withoutr':9.0,
    'withinr':10.0,
    'maxlength':20.0,
    'nogrid':10,
    'usehalfz':1
    }


inputplotdict['figintrhogfrompardist_m12mhr10kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'pz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mcr700hr10kpc',
    'withoutr':9.0,
    'withinr':10.0,
    'maxlength':10.0,
    'nogrid':10,
    'usehalfz':1
    }



inputplotdict['figrhogfrompardist_m12mhr10kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'dpdz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mcr700hr10kpc',
    'withoutr':9.0,
    'withinr':10.0,
    'maxlength':10.0,
    'nogrid':10,
    'usehalfz':1
    }



inputplotdict['figrhogfrompardist_m12ihr0kpchighz']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'dpdz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hr0kpchighz',
    'withoutr':0.0,
    'withinr':0.0,
    'maxlength':100.0,
    'nogrid':10,
    'usehalfz':1
    }


inputplotdict['figrhogfrompardist_m12ihr0kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'dpdz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hr0kpc',
    'withoutr':0.0,
    'withinr':0.0,
    'maxlength':10.0,
    'nogrid':10,
    'usehalfz':1
    }


inputplotdict['figrhogfrompardist_m12fhr1kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'dpdz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fcr700hr1kpc',
    'withoutr':0.0,
    'withinr':1.0,
    'maxlength':10.0,
    'nogrid':5,
    'usehalfz':1
    }



inputplotdict['figrhogfrompardist_m12ihr1kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'dpdz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hr1kpc',
    'withoutr':0.0,
    'withinr':1.0,
    'maxlength':10.0,
    'nogrid':5,
    'usehalfz':1
    }


inputplotdict['figrhogfrompardist_m12ihr']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'dpdz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hrfromz10kpc',
    'withoutr':6.0,
    'withinr':8.0,
    'maxlength':10.0,
    'nogrid':5,
    'usehalfz':1
    }


inputplotdict['figintrhogfrompardist_m12ihr']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'pz',
    'startno':550,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hr8kpc',
    'withoutr':6.0,
    'withinr':10.0,
    'maxlength':10.0,
    'nogrid':5,
    'usehalfz':1
    }


inputplotdict['figrhogfrompardist_m12ihrz10kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'dpdz',
    'startno':550,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hr10kpc',
    'withoutr':6.0,
    'withinr':10.0,
    'maxlength':20.0,
    'nogrid':10,
    'usehalfz':1
    }


inputplotdict['figrhogfrompardist_m12ihrtest']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'dpdz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hr10kpctest',
    'withoutr':6.0,
    'withinr':10.0,
    'maxlength':10.0,
    'nogrid':11,
    'usehalfz':0
    }



inputplotdict['figintrhogfrompardist_m12mhrtest']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mcr700hr10kpctest',
    'withoutr':6.0,
    'withinr':10.0,
    'maxlength':20.0,
    'nogrid':20,
    'usehalfz':0
    }



inputplotdict['figintrhogfrompardist_m12fhrtest']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fcr700hr10kpctest',
    'withoutr':6.0,
    'withinr':10.0,
    'maxlength':20.0,
    'nogrid':20,
    'usehalfz':0
    }


inputplotdict['figintrhogfrompardist_m12imhdcvhrtest']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12imhdcvtest',
    'withoutr':6.0,
    'withinr':10.0,
    'maxlength':20.0,
    'nogrid':20,
    'usehalfz':0
    }




inputplotdict['figintrhogfrompardist_m12ihrtest']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12iusepredata',
    'withoutr':6.0,
    'withinr':10.0,
    'maxlength':20.0,
    'nogrid':21,
    'usehalfz':0,
    'usepredata':1
    }


inputplotdict['figintrhogfrompardist_m12fhrtest']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fusepredata',
    'withoutr':6.0,
    'withinr':10.0,
    'maxlength':20.0,
    'nogrid':21,
    'usehalfz':0,
    'usepredata':1
    }

inputplotdict['figintrhogfrompardist_m12mhrtest']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12musepredata',
    'withoutr':6.0,
    'withinr':10.0,
    'maxlength':20.0,
    'nogrid':21,
    'usehalfz':0,
    'usepredata':1
    }


inputplotdict['figintrhogfrompardist_m12fhrceng2kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fgrid2kpc',
    'withoutr':6.0,
    'withinr':10.0,
    'maxlength':20.0,
    'nogrid':10,
    'usehalfz':0,
    'usepredata':1,
    'griddir':'grid2kpc'
    }

inputplotdict['figintrhogfrompardist_m12mhrceng2kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid2kpctktest',
    'withoutr':6.0,
    'withinr':10.0,
    'maxlength':20.0,
    'nogrid':10,
    'usehalfz':0,
    'usepredata':1,
    'griddir':'grid2kpc'
    }



inputplotdict['figintrhogfrompardist_m12fhrceng1kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fgrid1kpc',
    'withoutr':6.0,
    'withinr':10.0,
    'maxlength':20.0,
    'nogrid':21,
    'usehalfz':0,
    'usepredata':1,
    'griddir':'grid1kpc'
    }

inputplotdict['figintrhogfrompardist_m12mhrceng1kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid1kpc',
    'withoutr':6.0,
    'withinr':10.0,
    'maxlength':20.0,
    'nogrid':21,
    'usehalfz':0,
    'usepredata':1,
    'griddir':'grid1kpc'
    }


inputplotdict['figintrhogfrompardist_m12fhrceng0_5kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fgrid0_5kpc',
    'withoutr':6.0,
    'withinr':10.0,
    'maxlength':20.0,
    'nogrid':40,
    'usehalfz':0,
    'usepredata':1,
    'griddir':'grid0_5kpc'
    }



inputplotdict['figintrhogfrompardist_m12mhrceng0_5kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid0_5kpc',
    'withoutr':6.0,
    'withinr':10.0,
    'maxlength':20.0,
    'nogrid':40,
    'usehalfz':0,
    'usepredata':1,
    'griddir':'grid0_5kpc'
    }

inputplotdict['figintrhogfrompardist_m12ihrceng0_5kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12igrid0_5kpc',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':40,
    'usehalfz':0,
    'usepredata':1,
    'griddir':'grid0_5kpc'
    }


inputplotdict['figintrhogfrompardist_m12ihrceng2kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12igrid2kpc',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':10,
    'usehalfz':0,
    'usepredata':1,
    'griddir':'grid2kpc'
    }

inputplotdict['figintrhogfrompardist_m12ihr2kpccc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':580,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12igrid2kpccc',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':10,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':1,
    'griddir':'grid2kpc'
    }


inputplotdict['figintrhogfrompardist_m12mhr2kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid2kpc_r3kpc',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':10,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'griddir':'grid2kpc'
    }


inputplotdict['figintrhogfrompardist_m12fhr2kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fgrid2kpc_r3kpc',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':10,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'griddir':'grid2kpc'
    }


inputplotdict['figintrhogfrompardist_m12ihr2kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12igrid2kpc_r3kpc',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':10,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'griddir':'grid2kpc'
    }

inputplotdict['figintrhogfrompardist_m12mhr2kpc_usekez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid2kpc_r3kpc_usekez',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':10,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid2kpc'
    }


inputplotdict['figintrhogfrompardist_m12fhr2kpc_usekez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fgrid2kpc_r3kpc_usekez',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':10,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid2kpc'
    }


inputplotdict['figintrhogfrompardist_m12ihr2kpc_ueskez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12igrid2kpc_r3kpc_usekez',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':10,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid2kpc'
    }



inputplotdict['figintrhogfrompardist_m12mhr1kpc_usekez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid1kpc_r3kpc_usekez',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':21,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid1kpc'
    }


inputplotdict['figintrhogfrompardist_m12fhr1kpc_usekez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fgrid1kpc_r3kpc_usekez',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':21,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid1kpc'
    }


inputplotdict['figintrhogfrompardist_m12ihr1kpc_usekez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12igrid1kpc_r3kpc_usekez',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':21,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid1kpc'
    }



inputplotdict['figintrhogfrompardist_m12mhr0_5kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid0_5kpc_r3kpc',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':40,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'griddir':'grid0_5kpc'
    }


inputplotdict['figintrhogfrompardist_m12fhr0_5kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fgrid0_5kpc_r3kpc',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':40,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'griddir':'grid0_5kpc'
    }


inputplotdict['figintrhogfrompardist_m12ihr0_5kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12igrid0_5kpc_r3kpc',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':40,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'griddir':'grid0_5kpc'
    }



inputplotdict['figintrhogfrompardist_m12mhr0_5kpc_usekez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid0_5kpc_r3kpc_usekez',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':40,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_5kpc'
    }


inputplotdict['figintrhogfrompardist_m12fhr0_5kpc_usekez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fgrid0_5kpc_r3kpc_usekez',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':40,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_5kpc'
    }


inputplotdict['figintrhogfrompardist_m12ihr0_5kpc_usekez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12igrid0_5kpc_r3kpc_usekez',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':40,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_5kpc'
    }


inputplotdict['figintrhogfrompardist_m12mhr0_5kpcr9']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid0_5kpc_r9kpc',
    'withoutr':6.0,
    'withinr':9.0,
    'maxlength':20.0,
    'nogrid':40,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'griddir':'grid0_5kpc'
    }


inputplotdict['figintrhogfrompardist_m12fhr0_5kpcr9']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fgrid0_5kpc_r9kpc',
    'withoutr':6.0,
    'withinr':9.0,
    'maxlength':20.0,
    'nogrid':40,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'griddir':'grid0_5kpc'
    }

inputplotdict['figintrhogfrompardist_m12fhr0_5kpcr9kez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fgrid0_5kpc_r9kpc',
    'withoutr':6.0,
    'withinr':9.0,
    'maxlength':20.0,
    'nogrid':40,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'griddir':'grid0_5kpc',
    'usekez':1
    }


inputplotdict['figintrhogfrompardist_m12ihr0_5kpcr9']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12igrid0_5kpc_r9kpc',
    'withoutr':6.0,
    'withinr':9.0,
    'maxlength':20.0,
    'nogrid':40,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'griddir':'grid0_5kpc'
    }


inputplotdict['figintrhogfrompardist_m12mhr0_25kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid0_25kpc_r3kpc',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'griddir':'grid0_25kpc'
    }


inputplotdict['figintrhogfrompardist_m12fhr0_25kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fgrid0_25kpc_r3kpc',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'griddir':'grid0_25kpc'
    }


inputplotdict['figintrhogfrompardist_m12ihr0_25kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12igrid0_25kpc_r3kpc',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'griddir':'grid0_25kpc'
    }


inputplotdict['figintrhogfrompardist_m12mhr0_25kpcr9']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid0_25kpc_r9kpc',
    'withoutr':6.0,
    'withinr':9.0,
    'maxlength':20.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'griddir':'grid0_25kpc'
    }


inputplotdict['figintrhogfrompardist_m12fhr0_25kpcr9']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fgrid0_25kpc_r9kpc',
    'withoutr':6.0,
    'withinr':9.0,
    'maxlength':20.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'griddir':'grid0_25kpc'
    }


inputplotdict['figintrhogfrompardist_m12ihr0_25kpcr9']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12igrid0_25kpc_r9kpc',
    'withoutr':6.0,
    'withinr':9.0,
    'maxlength':20.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'griddir':'grid0_25kpc'
    }



inputplotdict['figintrhogfrompardist_m12mhr0_25kpc_usekez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid0_25kpc_r3kpc_usekez',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_25kpc'
    }


inputplotdict['figintrhogfrompardist_m12fhr0_25kpc_usekez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fgrid0_25kpc_r3kpc_usekez',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_25kpc'
    }


inputplotdict['figintrhogfrompardist_m12ihr0_25kpc_usekez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12igrid0_25kpc_r3kpc_usekez',
    'withoutr':6.0,
    'withinr':3.0,
    'maxlength':20.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_25kpc'
    }



inputplotdict['figintrhogfrompardist_m12mhr0_25kpcr9_usekez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid0_25kpc_r9kpc_usekez',
    'withoutr':6.0,
    'withinr':9.0,
    'maxlength':20.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_25kpc'
    }


inputplotdict['figintrhogfrompardist_m12fhr0_25kpcr9_usekez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fgrid0_25kpc_r9kpc_usekez',
    'withoutr':6.0,
    'withinr':9.0,
    'maxlength':20.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_25kpc'
    }


inputplotdict['figintrhogfrompardist_m12ihr0_25kpcr9_usekez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12igrid0_25kpc_r9kpc_usekez',
    'withoutr':6.0,
    'withinr':9.0,
    'maxlength':20.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_25kpc'
    }


inputplotdict['figintrhogfrompardist_m12ihr0_25kpcr0_9_usekez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12igrid0_25kpc_r0_9kpc_usekez',
    'withoutrnew':0.1,
    'withinr':9.0,
    'maxlength':20.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_25kpc'
    }


inputplotdict['figintrhogfrompardist_m12fhr0_25kpcr0_9_usekez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fgrid0_25kpc_r0_9kpc_usekez',
    'withoutrnew':0.1,
    'withinr':9.0,
    'maxlength':20.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_25kpc'
    }



inputplotdict['figintrhogfrompardist_m12mhr0_25kpcr0_9_usekez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid0_25kpc_r0_9kpc_usekez',
    'withoutrnew':0.1,
    'withinr':9.0,
    'maxlength':20.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_25kpc'
    }

inputplotdict['figintrhogfrompardist_m12mhr0_5kpcr0_9_usekez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid0_5kpc_r0_9kpc_usekez',
    'withoutrnew':0.1,
    'withinr':9.0,
    'withoutr':0.1,
    'maxlength':2.0,
    'nogrid':40,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_5kpc'
    }

inputplotdict['figintrhogfrompardist_m12hr0_5kpcr0_9_usekez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':[['m12imhdcvhr','m12fmhdcvhr','m12mmhdcvhr'],['m12icr_700hr','m12fcr_700hr','m12mcr_700hr']],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12grid0_5kpc_r0_9kpc_usekez',
    'withoutrnew':0.1,
    'withinr':9.0,
    'withoutr':0.1,
    'maxlength':2.0,
    'nogrid':40,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_5kpc'
    }


inputplotdict['figintrhogfrompardist_m12hr0_25kpcr0_9_usekez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':[['m12imhdcvhr','m12fmhdcvhr','m12mmhdcvhr'],['m12icr_700hr','m12fcr_700hr','m12mcr_700hr']],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12grid0_25kpc_r0_9kpc_usekez',
    'withoutrnew':0.1,
    'withinr':9.0,
    'withoutr':0.1,
    'maxlength':2.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_25kpc'
    }


inputplotdict['figintrhogfrompardist_m11g0_5kpcr0_9_usekez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m11gmhdcv','m11gcr_700'],
    #'dirneed':['m11gmhdcv'],
    'wanted':'pz',
    'startno':600,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m11ggrid0_5kpc_r0_9kpc_usekez',
    'withoutrnew':0.1,
    'withinr':9.0,
    'withoutr':0.1,
    'maxlength':2.0,
    'nogrid':40,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':0,
    'griddir':'grid0_5kpc'
    }

inputplotdict['figintrhogfrompardist_m11b0_5kpcr0_9_usekez']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m11bmhdcv','m11bcr_700'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m11bgrid0_5kpc_r0_9kpc_usekez',
    'withoutrnew':0.1,
    'withinr':9.0,
    'withoutr':0.1,
    'maxlength':2.0,
    'nogrid':40,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_5kpc'
    }

inputplotdict['figpz_m12ihr0_25kpcr9_usekez']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12igrid0_25kpc_r9kpc_usekez',
    'withoutr':6.0,
    'withinr':9.0,
    'maxlength':20.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_25kpc'
    }

inputplotdict['figpz_m12ihr0_25kpc_usekez']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12igrid0_25kpc_r9kpc_usekez',
    'withoutr':0.1,
    'withinr':9.0,
    'maxlength':20.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_25kpc'
    }


inputplotdict['figpr_m12ihr0_5kpc_usekez']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'pr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12igrid0_5kpc_r9kpc_usekez',
    'withoutr':0.1,
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_5kpc'
    }


inputplotdict['figtdspr_m12hr0_5kpc_usekez']={
    '_plotfunction':turdifscale_testinput,
    'dirneed':[['m12imhdcvhr','m12fmhdcvhr','m12mmhdcvhr'],['m12icr_700hr','m12fcr_700hr','m12mcr_700hr']],
    'wanted':'pr',
    'startno':600,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12grid0_5kpc_r9kpc_usekez',
    'withoutr':0.1,
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1
    }

inputplotdict['figtdsvz_m12hr0_5kpc_usekez']={
    '_plotfunction':turdifscale_testinput,
    'dirneed':[['m12imhdcvhr','m12fmhdcvhr','m12mmhdcvhr'],['m12icr_700hr','m12fcr_700hr','m12mcr_700hr']],
    'wanted':'vz',
    'startno':600,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12grid0_5kpc_r9kpc_usekez',
    'withoutr':0.1,
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1
    }

inputplotdict['figtdspz_m12hr0_5kpc_usekez']={
    '_plotfunction':turdifscale_testinput,
    'dirneed':[['m12imhdcvhr','m12fmhdcvhr','m12mmhdcvhr'],['m12icr_700hr','m12fcr_700hr','m12mcr_700hr']],
    'wanted':'pz',
    'startno':600,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12grid0_5kpc_r9kpc_usekez',
    'withoutr':0.1,
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1
    }

inputplotdict['figpr_m12hr0_5kpc_usekez']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':[['m12imhdcvhr','m12fmhdcvhr','m12mmhdcvhr'],['m12icr_700hr','m12fcr_700hr','m12mcr_700hr']],
    'wanted':'pr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12grid0_5kpc_r9kpc_usekez',
    'withoutr':0.1,
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_5kpc'
    }

inputplotdict['figvr_m12hr0_5kpc_usekez']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':[['m12imhdcvhr','m12fmhdcvhr','m12mmhdcvhr'],['m12icr_700hr','m12fcr_700hr','m12mcr_700hr']],
    'wanted':'vr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12grid0_5kpc_r9kpc_usekez',
    'withoutr':0.1,
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_5kpc'
    }


inputplotdict['figpr_m12ihr0_25kpc_usekez']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'pr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12igrid0_25kpc_r9kpc_usekez',
    'withoutr':0.1,
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_25kpc'
    }


inputplotdict['figvr_m12ihr0_5kpc_usekez']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'vr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12igrid0_5kpc_r9kpc_usekez',
    'withoutr':0.1,
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_5kpc'
    }


inputplotdict['figvr_m12ihr0_25kpc_usekez']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'vr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12igrid0_25kpc_r9kpc_usekez',
    'withoutr':0.1,
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_25kpc'
    }

inputplotdict['figpz_m12ihr0_125kpcr9_usekez']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'pz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12igrid0_125kpc_r9kpc_usekez',
    'withoutr':6.0,
    'withinr':9.0,
    'maxlength':20.0,
    'nogrid':160,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_125kpc'
    }

inputplotdict['figvr_m12fhr0_5kpc_usekez']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'vr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fgrid0_5kpc_r9kpc_usekez',
    'withoutr':0.1,
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_5kpc'
    }


inputplotdict['figvr_m12fhr0_25kpc_usekez']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'vr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fgrid0_25kpc_r9kpc_usekez',
    'withoutr':0.1,
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_25kpc'
    }


inputplotdict['figvr_m12mhr0_5kpc_usekez']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'vr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid0_5kpc_r9kpc_usekez',
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_5kpc'
    }


inputplotdict['figvr_m12mhr0_25kpc_usekez']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'vr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid0_25kpc_r9kpc_usekez',
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_25kpc'
    }





inputplotdict['figvr_m12mhr0_125kpc_usekez']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'vr',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid0_125kpc_r9kpc_usekez',
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':160,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_125kpc'
    }

inputplotdict['figvz_m12mhr0_5kpc_usekez']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'vz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid0_5kpc_r9kpc_usekez',
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':40,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_5kpc'
    }


inputplotdict['figvz_m12mhr0_25kpc_usekez']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'vz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid0_25kpc_r9kpc_usekez',
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_25kpc'
    }

inputplotdict['figvz_m12mhr2kpc_usekezHI']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'vz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid2kpc_r9kpc_usekez',
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':10,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'outHI':1,
    'griddir':'grid2kpc'
    }

inputplotdict['figvturr_m12mhr2kpc_usekezHI']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'vturr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid2kpc_r9kpc_usekez',
    'withinr':9.0,
    'maxlength':3.0,
    'nogrid':10,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'outHI':1,
    'griddir':'grid2kpc'
    }


inputplotdict['figvz_m12mhr1kpc_usekezHI']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'vz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid2kpc_r9kpc_usekez',
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':21,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'outHI':1,
    'griddir':'grid1kpc'
    }


inputplotdict['figvturr_m12mhr1kpc_usekezHI']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'vturr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid1kpc_r9kpc_usekez',
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':21,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'outHI':1,
    'griddir':'grid1kpc'
    }

inputplotdict['figvz_m12fhr1kpc_usekezHI']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'vz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fgrid2kpc_r9kpc_usekez',
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':21,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'outHI':1,
    'griddir':'grid1kpc'
    }

inputplotdict['figvz_m12hr0_25kpc_usekezHI']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':[['m12imhdcvhr','m12fmhdcvhr','m12mmhdcvhr'],['m12icr_700hr','m12fcr_700hr','m12mcr_700hr']],
    'wanted':'vz',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12grid0_25kpc_r9kpc_usekez',
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'outHI':1,
    'griddir':'grid0_25kpc'
    }


inputplotdict['figvturr_m12ihr1kpc_usekezHI']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'vturr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12igrid1kpc_r9kpc_usekez',
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':21,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'outHI':1,
    'griddir':'grid1kpc'
    }




inputplotdict['figvturr_m12mhr1kpc_usekezHI']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'vturr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid1kpc_r9kpc_usekez',
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':21,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'outHI':1,
    'griddir':'grid1kpc'
    }


inputplotdict['figvturr_m12mhr0_5kpc_usekezHI']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'vturr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid0_5kpc_r9kpc_usekez',
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':40,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'outHI':1,
    'griddir':'grid0_5kpc'
    }


inputplotdict['figvturr_m12ihr0_5kpc_usekezHI']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'vturr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12igrid0_5kpc_r9kpc_usekez',
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':40,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'outHI':1,
    'griddir':'grid0_5kpc'
    }


inputplotdict['figvturr_m12hr0_5kpc_usekezHI']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':[['m12imhdcvhr','m12fmhdcvhr','m12mmhdcvhr'],['m12icr_700hr','m12fcr_700hr','m12mcr_700hr']],
    'wanted':'vturr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12grid0_5kpc_r9kpc_usekez',
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':40,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'outHI':1,
    'griddir':'grid0_5kpc'
    }

inputplotdict['figvturr_m12hr0_25kpc_usekezHI']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':[['m12imhdcvhr','m12fmhdcvhr','m12mmhdcvhr'],['m12icr_700hr','m12fcr_700hr','m12mcr_700hr']],
    'wanted':'vturr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12grid0_25kpc_r9kpc_usekez',
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'outHI':1,
    'griddir':'grid0_25kpc'
    }

inputplotdict['figvturr_m12mhr0_25kpc_usekezHI']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'vturr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mgrid0_25kpc_r9kpc_usekez',
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'outHI':1,
    'griddir':'grid0_25kpc'
    }


inputplotdict['figvz_m12fhr0_5kpc_usekez']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'vz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fgrid0_5kpc_r9kpc_usekez',
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':40,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_5kpc'
    }


inputplotdict['figvz_m12fhr0_25kpc_usekez']={
    '_plotfunction':crdenfromdata_testinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'vz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fgrid0_25kpc_r9kpc_usekez',
    'withinr':9.0,
    'maxlength':1.0,
    'nogrid':80,
    'usehalfz':0,
    'usepredata':1,
    'cutcold':0,
    'usekez':1,
    'griddir':'grid0_25kpc'
    }





inputplotdict['figrhogfrompardist_m12ihr10kpc']={
    '_plotfunction':rhogzdispartestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'dpdz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hr10kpc',
    'withoutr':-1.0,
    'withinr':10.0,
    'maxlength':10.0,
    'nogrid':10,
    'usehalfz':1
    }


inputplotdict['figsfrtur_m_m11']={
    '_plotfunction':sfrturtestinput,
    #'dirneed':['m12icr_700hr'],
    'dirneed':['m11bmhdcv','f476'],
    #,'m12fmhdcvhr','m12fcr_700hr','m12mmhdcvhr','m12mcr_700hr','m12imhdcvhr','m12icr_700hr'
    'wanted':'sfrtur_m',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':0.2,
    'withinr': 20.0,
    'withoutr': 0.1
    }



inputplotdict['figsfrtur_m11']={
    '_plotfunction':sfrturtestinput,
    #'dirneed':['m12icr_700hr'],
    'dirneed':['m11bmhdcv','f476'],
    #,'m12fmhdcvhr','m12fcr_700hr','m12mmhdcvhr','m12mcr_700hr','m12imhdcvhr','m12icr_700hr'
    'wanted':'sfrtur',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':0.2,
    'withinr': 20.0,
    'withoutr': 0.1
    }

inputplotdict['figsfrtur_m11_m12']={
    '_plotfunction':sfrturtestinput,
    #'dirneed':['m12icr_700hr'],
    'dirneed':['m11gcr_700','m11gmhdcv','m11dmhdcv','m11hmhdcv','m12bcr_700hr','m12wcr_700hr','m12rcr_700hr','m11bmhdcv',\
               'm11bcr_700','m11fmhdcv','m11fcr_700','fm11q',\
               'f476','f61','f46','fm12f','fm12i','fm12m','m12fmhdcvhr',\
               'm12fcr_700hr','m12mmhdcvhr','m12mcr_700hr',\
               'm12imhdcvhr','m12icr_700hr','f573',],
    #,'m12fmhdcvhr','m12fcr_700hr','m12mmhdcvhr','m12mcr_700hr','m12imhdcvhr','m12icr_700hr'
    'wanted':'sfrtur',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':0.2,
    'withinr': 20.0,
    'withoutr': 0.1
    }


inputplotdict['fig_sigmaz_m11_m12']={
    '_plotfunction':sfrturtestinput,
    'dirneed':['m12icr_700hr'],
    #'dirneed':['m11gcr_700','m11gmhdcv','m11dmhdcv','m11hmhdcv','m12bcr_700hr','m12wcr_700hr','m12rcr_700hr','m11bmhdcv',\
    #           'm11bcr_700','m11fmhdcv','m11fcr_700','fm11q',\
    #           'f476','f61','f46','fm12f','fm12i','fm12m','m12fmhdcvhr',\
    #           'm12fcr_700hr','m12mmhdcvhr','m12mcr_700hr',\
    #           'm12imhdcvhr','m12icr_700hr','f573',],
    #,'m12fmhdcvhr','m12fcr_700hr','m12mmhdcvhr','m12mcr_700hr','m12imhdcvhr','m12icr_700hr'
    'wanted':'sigmatur_mweighted',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':1.0,
    'withinr': 10.0,
    'withoutr': 0.1
    }


inputplotdict['fig_sigmaz_m12']={
    '_plotfunction':sfrturtestinput,
    #'dirneed':['m12icr_700hr'],
    #'dirneed':['m12bcr_700hr','m12wcr_700hr','m12rcr_700hr','fm12f','fm12i','fm12m','m12fmhdcvhr',\
    #           'm12fcr_700hr','m12mmhdcvhr','m12mcr_700hr',\
    #           'm12imhdcvhr','m12icr_700hr'],
    'dirneed':['fm12f','fm12i','fm12m',\
              'm12mmhdcvhr','m12imhdcvhr','m12fmhdcvhr',\
              'm12mcr_700hr','m12fcr_700hr','m12icr_700hr'],
    #'dirneed':['m12fcr_700hr','m12mcr_700hr','m12icr_700hr'],
    'wanted':'sigmatur_mweighted',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':1.0,
    'withinr': 10.0,
    'withoutr': 0.1
    }


inputplotdict['figsfrtur_m_m11_m12']={
    '_plotfunction':sfrturtestinput,
    #'dirneed':['m12icr_700hr'],
    'dirneed':['m11gcr_700','m11gmhdcv','m11dmhdcv','m11hmhdcv','m12bcr_700hr','m12wcr_700hr','m12rcr_700hr','m11bmhdcv',\
               'm11bcr_700','m11fmhdcv','m11fcr_700','fm11q',\
               'f476','f61','f46','fm12f','fm12i','fm12m','m12fmhdcvhr',\
               'm12fcr_700hr','m12mmhdcvhr','m12mcr_700hr',\
               'm12imhdcvhr','m12icr_700hr','f573',],
    #,'m12fmhdcvhr','m12fcr_700hr','m12mmhdcvhr','m12mcr_700hr','m12imhdcvhr','m12icr_700hr'
    'wanted':'sfrtur_m',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':0.2,
    'withinr': 20.0,
    'withoutr': 0.1
    }


inputplotdict['figsfrmocktur_m12i']={
    '_plotfunction':sfrmockturtestinput,
    'dirneed':['m12icr_700hr'],
    #'dirneed':['m12fmhdcvhr','m12fcr_700hr','m12mmhdcvhr','m12mcr_700hr','m12imhdcvhr','m12icr_700hr'],
    'wanted':'',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':10.0,
    'withinr': 20.0,
    'withoutr': 0.1,
    'usesnipshot':1
    }


inputplotdict['figsfrmocktur_m11gcr_700']={
    '_plotfunction':sfrmockturtestinput,
    'dirneed':['m11gcr_700'],
    #'dirneed':['m12fmhdcvhr','m12fcr_700hr','m12mmhdcvhr','m12mcr_700hr','m12imhdcvhr','m12icr_700hr'],
    'wanted':'',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':10.0,
    'withinr': 20.0,
    'withoutr': 0.1,
    'usesnipshot':0
    }


inputplotdict['figsfrmocktur_m11bcr_700']={
    '_plotfunction':sfrmockturtestinput,
    'dirneed':['m11bcr_700'],
    #'dirneed':['m12fmhdcvhr','m12fcr_700hr','m12mmhdcvhr','m12mcr_700hr','m12imhdcvhr','m12icr_700hr'],
    'wanted':'',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':10.0,
    'withinr': 20.0,
    'withoutr': 0.1,
    'usesnipshot':0
    }


inputplotdict['figsfrmocktur_m12imhdcv']={
    '_plotfunction':sfrmockturtestinput,
    'dirneed':['m12imhdcvhr'],
    #'dirneed':['m12fmhdcvhr','m12fcr_700hr','m12mmhdcvhr','m12mcr_700hr','m12imhdcvhr','m12icr_700hr'],
    'wanted':'SFRturrad',
    'startno':590,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':10.0,
    'withinr': 20.0,
    'withoutr': 0.1,
    'usesnipshot':1
    }


inputplotdict['figPextBlitz_m12']={
    '_plotfunction':Blitzmidplanepressure_testinput,
    #'dirneed':['m12icr_700hr','m12imhdcvhr'],
    #'dirneed':['m12fmhdcvhr'],
    'dirneed':['m12fmhdcvhr','m12fcr_700hr','m12mmhdcvhr','m12mcr_700hr','m12imhdcvhr','m12icr_700hr'],
    'wanted':'',
    'startno':580,
    'Nsnap':600,
    'snapsep':1,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fmhdcv',
    'maxlength':2.0,
    'withinr': 10.0,
    'withoutr': 0.1,
    'usesnipshot':1
    }


inputplotdict['figsfrmockturtot_m12']={
    '_plotfunction':sfrmockturtestinput,
    #'dirneed':['m12icr_700hr','m12imhdcvhr'],
    'dirneed':['m12fmhdcvhr','m12fcr_700hr','m12mmhdcvhr','m12mcr_700hr','m12imhdcvhr','m12icr_700hr'],
    #'dirneed':['m12fmhdcvhr','m12fcr_700hr','m12mmhdcvhr','m12mcr_700hr','m12imhdcvhr','m12icr_700hr'],
    'wanted':'SFRturtot',
    'startno':580,
    'Nsnap':601,
    'snapsep':1,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12',
    'maxlength':2.0,
    'withinr': 10.0,
    'withoutr': 0.1,
    'usesnipshot':1
    }


inputplotdict['figsfrmockturtot_m11_m12']={
    '_plotfunction':sfrmockturtestinput,
    #'dirneed':['m12icr_700hr','m12imhdcvhr'],
    'dirneed':['m11gmhdcv','m11gcr_700','m11dmhdcv','m11dcr_700','m11hmhdcv','m11hcr_700',\
               'm11bmhdcv','m11bcr_700','m11fmhdcv','m11fcr_700',\
               'm12fmhdcvhr','m12fcr_700hr','m12mmhdcvhr','m12mcr_700hr','m12imhdcvhr','m12icr_700hr'],
    #'dirneed':['m12fmhdcvhr','m12fcr_700hr','m12mmhdcvhr','m12mcr_700hr','m12imhdcvhr','m12icr_700hr'],
    'wanted':'SFRturtot',
    'startno':580,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m11_m12',
    'maxlength':2.0,
    'withinr': 10.0,
    'withoutr': 0.1,
    'usesnipshot':1
    }



inputplotdict['figsfrtur_m12i']={
    '_plotfunction':sfrturtestinput,
    #'dirneed':['m12icr_700hr'],
    'dirneed':['m12fmhdcvhr','m12fcr_700hr','m12mmhdcvhr','m12mcr_700hr','m12imhdcvhr','m12icr_700hr'],
    'wanted':'sfrtur',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':0.2,
    'withinr': 20.0,
    'withoutr': 0.1
    }

inputplotdict['figsfrvth_m12i']={
    '_plotfunction':sfrturtestinput,
    #'dirneed':['m12icr_700hr'],
    'dirneed':['fm12f','fm12i','fm12m',
               'm12fmhdcvhr','m12fcr_700hr','m12mmhdcvhr',
               'm12mcr_700hr','m12imhdcvhr','m12icr_700hr'],
    'wanted':'vth_mweighted',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':1.0,
    'withinr': 10.0,
    'withoutr': 0.1
    }



inputplotdict['figsfrv_m12']={
    '_plotfunction':sfrturtestinput,
    #'dirneed':['m12icr_700hr'],
    'dirneed':['fm12f','fm12i','fm12m',
               'm12fmhdcvhr','m12fcr_700hr','m12mmhdcvhr',
               'm12mcr_700hr','m12imhdcvhr','m12icr_700hr'],
    'wanted':'v_mweighted',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':1.0,
    'withinr': 10.0,
    'withoutr': 0.1
    }

inputplotdict['figsfr_m12']={
    '_plotfunction':sfrtestinput,
    #'dirneed':['m12icr_700hr'],
    'dirneed':['m12fmhdcvhr','m12fcr_700hr','m12mmhdcvhr',
               'm12mcr_700hr','m12imhdcvhr','m12icr_700hr'],
    'wanted':'sfr',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':1.0,
    'withinr': 10.0,
    'withoutr': 0.1
    }

inputplotdict['figsfrhotgasmass_test']={
    '_plotfunction':sfrhotgasmass_testinput,
    'dirneed':['m12fmhdcvhr'],
    'wanted':'sfrhotgasmass',
    'startno':590,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'test',
    'maxlength':1.0,
    'withinr': 10.0,
    'withoutr': 0.1,
    'minsfr':1e-2,
    'maxsfr':20
    }

inputplotdict['figsfrhotgasmass_m12']={
    '_plotfunction':sfrhotgasmass_testinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr','m12mmhdcvhr','m12mcr_700hr',\
               'm12imhdcvhr','m12icr_700hr'],
    'wanted':'sfrhotgasmass',
    'startno':590,
    'Nsnap':601,
    'snapsep':1,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12',
    'maxlength':1.0,
    'withinr': 10.0,
    'withoutr': 0.1,
    'minsfr':1e-2,
    'maxsfr':20
    }

inputplotdict['figsfrhotgasmass_m12_20kpc']={
    '_plotfunction':sfrhotgasmass_testinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr','m12mmhdcvhr','m12mcr_700hr',\
               'm12imhdcvhr','m12icr_700hr'],
    'wanted':'sfrhotgasmass',
    'startno':590,
    'Nsnap':601,
    'snapsep':1,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12_20kpc',
    'maxlength':1.0,
    'withinr': 20.0,
    'withoutr': 0.1,
    'minsfr':1e-2,
    'maxsfr':20
    }


inputplotdict['figsfrv_m11_m12']={
    '_plotfunction':sfrturtestinput,
    'dirneed':['m11gcr_700','m11gmhdcv','m11dmhdcv','m11hmhdcv','m12bcr_700hr','m12wcr_700hr','m12rcr_700hr','m11bmhdcv',\
               'm11bcr_700','m11fmhdcv','m11fcr_700','fm11q',\
               'f476','f61','f46','fm12f','fm12i','fm12m','m12fmhdcvhr',\
               'm12fcr_700hr','m12mmhdcvhr','m12mcr_700hr',\
               'm12imhdcvhr','m12icr_700hr','f573',],
    'wanted':'v_mweighted',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m11_m12',
    'maxlength':1.0,
    'withinr': 10.0,
    'withoutr': 0.1,
    'minsfr':1e-2,
    'maxsfr':20
    }


inputplotdict['figeturm12fhr']={
    '_plotfunction':eturtimetestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'eturtime',
    'startno':580,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fhr'
    }

inputplotdict['figeturm12ihr']={
    '_plotfunction':eturtimetestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'eturtime',
    'startno':580,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12ihr'
    }


inputplotdict['figeturm12mhr']={
    '_plotfunction':eturtimetestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'eturtime',
    'startno':580,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mhr'
    }


inputplotdict['figetur_m12ihr']={
    '_plotfunction':eturtimetestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'etur_mtime',
    'startno':580,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12ihr'
    }

inputplotdict['figetur_m12mhr']={
    '_plotfunction':eturtimetestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'etur_mtime',
    'startno':580,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mhr'
    }



inputplotdict['figetur_smc']={
    '_plotfunction':eturtimetestinput,
    'dirneed':['bwsmclr','bwsmclrmhd','bwsmclrdc28mhd'],
    'wanted':'etur_mtime',
    'startno':400,
    'Nsnap':500,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':''
    }




inputplotdict['figCRz_m12m']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['m12mmhdcv','m12mcr_b_70','m12mcr_700'],
    'wanted':'crdenz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':''
    }


inputplotdict['figabsvz_m12i']={
    '_plotfunction':vzcylouttestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'absvzout',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':100.0,
    'withinr': 15.0,
    'withoutr': 8.0,
    'title':'m12ires7000'
    }


inputplotdict['figabsvz_m12m']={
    '_plotfunction':vzcylouttestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'absvzout',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':100.0,
    'withinr': 15.0,
    'withoutr': 8.0,
    'title':'m12mres7000'
    }



inputplotdict['figabsvz_m12f']={
    '_plotfunction':vzcylouttestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'absvzout',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':100.0,
    'withinr': 15.0,
    'withoutr': 8.0,
    'title':'m12fres7000'
    }

inputplotdict['figgmz_m12i']={
    '_plotfunction':vzcylouttestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'gmz',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':10.0,
    'withinr': 15.0,
    'withoutr': 8.0,
    'title':'m12ires7000'
    }



inputplotdict['fighg_m12m']={
    '_plotfunction':hgashalftestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'hgas',
    'startno':580,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':20.0,
    'withinr': 15.0,
    'withoutr': 1.0,
    'title':'m12mres7000'
    }


inputplotdict['fighg_m12f']={
    '_plotfunction':hgashalftestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'hgas',
    'startno':580,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':20.0,
    'withinr': 15.0,
    'withoutr': 1.0,
    'title':'m12fres7000'
    }

inputplotdict['fighg_m12i']={
    '_plotfunction':hgashalftestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'hgas',
    'startno':580,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':20.0,
    'withinr': 15.0,
    'withoutr': 1.0,
    'title':'m12ires7000'
    }

inputplotdict['fighcr_m12i']={
    '_plotfunction':hgashalftestinput,
    'dirneed':['m12icr_700hr'],
    'wanted':'hcr',
    'startno':580,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':20.0,
    'withinr': 15.0,
    'withoutr': 1.0,
    'title':'m12ires7000'
    }


inputplotdict['fighcr_m12m']={
    '_plotfunction':hgashalftestinput,
    'dirneed':['m12mcr_700hr'],
    'wanted':'hcr',
    'startno':580,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':20.0,
    'withinr': 15.0,
    'withoutr': 1.0,
    'title':'m12mres7000'
    }


inputplotdict['fighth_m12i']={
    '_plotfunction':hgashalftestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'hth',
    'startno':580,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':20.0,
    'withinr': 15.0,
    'withoutr': 1.0,
    'title':'m12ires7000'
    }

inputplotdict['fighth_m12m']={
    '_plotfunction':hgashalftestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'hth',
    'startno':580,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':20.0,
    'withinr': 15.0,
    'withoutr': 1.0,
    'title':'m12mres7000'
    }

inputplotdict['fighth_m12f']={
    '_plotfunction':hgashalftestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'hth',
    'startno':580,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':20.0,
    'withinr': 15.0,
    'withoutr': 1.0,
    'title':'m12fres7000'
    }

inputplotdict['figmgasz_m12i']={
    '_plotfunction':hgashalftestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'mgasz',
    'startno':580,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':5.0,
    'withinr': 15.0,
    'withoutr': 1.0,
    'title':'m12ires7000'
    }


inputplotdict['figmgasz_m12m']={
    '_plotfunction':hgashalftestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'mgasz',
    'startno':580,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':5.0,
    'withinr': 15.0,
    'withoutr': 1.0,
    'title':'m12mres7000'
    }

inputplotdict['figmgasz_m12f']={
    '_plotfunction':hgashalftestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'mgasz',
    'startno':580,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':5.0,
    'withinr': 15.0,
    'withoutr': 1.0,
    'title':'m12fres7000'
    }


inputplotdict['figgasdenTz_m12m']={
    '_plotfunction':gasdenTtestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'gasdenTz',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':10.0,
    'withinr': 10.0,
    'withoutr': 0.1,
    'title':'m12mres7000'
    }


inputplotdict['figgasdenTz_m12i']={
    '_plotfunction':gasdenTtestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'gasdenTz',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':20.0,
    'withinr': 10.0,
    'withoutr': 0.1,
    'title':'m12ires7000'
    }

inputplotdict['figgasdenTz_m12']={
    '_plotfunction':gasdenTtestinput,
    #'dirneed':[['m12imhdcvhr','m12icr_700hr']],
    'dirneed':[['m12imhdcvhr','m12icr_700hr'],['m12fmhdcvhr','m12fcr_700hr'],['m12mmhdcvhr','m12mcr_700hr']],
    'wanted':'gasdenTz',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':20.0,
    'withinr': 10.0,
    'withoutr': 0.1,
    'title':['m12ires7000','m12fres7000','m12mres7000']
    }

inputplotdict['figgasdenTz_m12_atr8kpc']={
    '_plotfunction':gasdenTtestinput,
    #'dirneed':[['m12imhdcvhr','m12icr_700hr']],
    #'dirneed':[['m12mmhdcvhr','m12mcr_700hr']],
    'dirneed':[['m12imhdcvhr','m12icr_700hr'],['m12fmhdcvhr','m12fcr_700hr'],['m12mmhdcvhr','m12mcr_700hr']],
    'wanted':'gasdenTz',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'at8kpc',
    'maxlength':4.0,
    'withinr': 9.0,
    'withoutr': 7.0,
    'phicut':1,
    'usehz':0,
    'title':['m12ires7000','m12fres7000','m12mres7000']
    }

inputplotdict['figgasdenTz_m12f_atr8kpc']={
    '_plotfunction':gasdenTtestinput,
    #'dirneed':[['m12imhdcvhr','m12icr_700hr']],
    'dirneed':[['m12fmhdcvhr','m12fcr_700hr']],
    #'dirneed':[['m12imhdcvhr','m12icr_700hr'],['m12fmhdcvhr','m12fcr_700hr'],['m12mmhdcvhr','m12mcr_700hr']],
    'wanted':'gasdenTz',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':5.0,
    'withinr': 7.5,
    'withoutr': 5.5,
    'phicut':1,
    'usehz':0,
    'title':['m12ires7000','m12fres7000','m12mres7000']
    }


inputplotdict['figrhoTvol_m12']={
    '_plotfunction':gasdenThistestinput,
    'dirneed':[['m12imhdcvhr','m12icr_700hr'],['m12fmhdcvhr','m12fcr_700hr'],['m12mmhdcvhr','m12mcr_700hr']],
    'wanted':'rhoTvol',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':1.0,
    'withinr': 10.0,
    'withoutr': 0.1,
    'title':['m12ires7000','m12fres7000','m12mres7000']
    }

inputplotdict['figrhoTvoloutz_m12']={
    '_plotfunction':gasdenThistestinput,
    'dirneed':[['m12imhdcvhr','m12icr_700hr'],['m12fmhdcvhr','m12fcr_700hr'],['m12mmhdcvhr','m12mcr_700hr']],
    'wanted':'rhoTvol',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':1.0,
    'withinr': 10.0,
    'withoutr': 0.1,
    'outz':1,
    'zup':10.0,
    'zdown':1.0,
    'title':['m12ires7000','m12fres7000','m12mres7000']
    }


inputplotdict['figrhoTvol_m12out']={
    '_plotfunction':gasdenThistestinput,
    'dirneed':[['m12imhdcvhr','m12icr_700hr'],['m12fmhdcvhr','m12fcr_700hr'],['m12mmhdcvhr','m12mcr_700hr']],
    'wanted':'rhoTvol',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'out',
    'maxlength':400.0,
    'withinr': 200.0,
    'withoutr': 30.0,
    'title':['m12ires7000','m12fres7000','m12mres7000']
    }

inputplotdict['figrhoTmass_m12']={
    '_plotfunction':gasdenThistestinput,
    'dirneed':[['m12imhdcvhr','m12icr_700hr'],['m12fmhdcvhr','m12fcr_700hr'],['m12mmhdcvhr','m12mcr_700hr']],
    'wanted':'rhoTmass',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':1.0,
    'withinr': 10.0,
    'withoutr': 0.1,
    'title':['m12ires7000','m12fres7000','m12mres7000']
    }


inputplotdict['figgasdenTr_m12i']={
    '_plotfunction':gasdenTtestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'gasdenTr',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':1.0,
    'withinr': 10.0,
    'withoutr': 0.1,
    'title':'m12ires7000'
    }

inputplotdict['figgasdenTr_m12']={
    '_plotfunction':gasdenTtestinput,
    'dirneed':[['m12imhdcvhr','m12icr_700hr'],['m12fmhdcvhr','m12fcr_700hr'],['m12mmhdcvhr','m12mcr_700hr']],
    'wanted':'gasdenTr',
    'startno':580,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':1.0,
    'withinr': 10.0,
    'withoutr': 0.1,
    'title':['m12ires7000','m12fres7000','m12mres7000']
    }


inputplotdict['figgmz_m12m']={
    '_plotfunction':vzcylouttestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'gmz',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':10.0,
    'withinr': 15.0,
    'withoutr': 8.0,
    'title':'m12mres7000'
    }


inputplotdict['figanHmid_m12i']={
    '_plotfunction':nHmidplanetestinput,
    'dirneed':['m12imhdcvhr','m12icr_b_70hr','m12icr_700hr'],
    'wanted':'nHmidplane',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':0.2,
    'withinr': 20.0,
    'withoutr': 0.1
    }


inputplotdict['figanHmid_m12m']={
    '_plotfunction':nHmidplanetestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'nHmidplane',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':0.2,
    'withinr': 20.0,
    'withoutr': 0.1
    }

inputplotdict['figanHz_m12i']={
    '_plotfunction':nHmidplanetestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'nHz',
    'startno':580,
    'Nsnap':601,
    'snapsep':1,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':4.0,
    'withinr': 5.0,
    'withoutr': 4.0,
    'phicut':0,
    'useobs':1
    }

inputplotdict['figanHz_m12f']={
    '_plotfunction':nHmidplanetestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'nHz',
    'startno':580,
    'Nsnap':601,
    'snapsep':1,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':4.0,
    'withinr': 5.0,
    'withoutr': 4.0,
    'phicut':0,
    'useobs':1
    }

inputplotdict['figanHz_m12m']={
    '_plotfunction':nHmidplanetestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'nHz',
    'startno':580,
    'Nsnap':601,
    'snapsep':1,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':4.0,
    'withinr': 5.0,
    'withoutr': 4.0,
    'phicut':0,
    'useobs':1
    }



inputplotdict['figCRgasmid_m12mgas']={
    '_plotfunction':crdenmidplanetestinput,
    'dirneed':['m12mcr_b_70','m12mcr_700'],
    'wanted':'crdenmidplane',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'normalizedsm':0,
    'M1labelneed':0,
    'M1runlabelneed':0,
    'resoneed':0,
    'diffusionsolverneed':0,
    'refereelabelneed':0,
    'newlabelneed':1,
    'strlabelneed':0,
    'showstarburst':0,
    'legendneed':1,
    'correctIa':0,
    'runlabelneed':0,
    'withincr':0
    }


inputplotdict['figCRgasmid_m12igas']={
    '_plotfunction':crdenmidplanetestinput,
    'dirneed':['m12imhdcvhr','m12icr_b_70hr','m12icr_700hr'],
    'wanted':'crdenmidplane',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'normalizedsm':0,
    'M1labelneed':0,
    'M1runlabelneed':0,
    'resoneed':0,
    'diffusionsolverneed':0,
    'refereelabelneed':0,
    'newlabelneed':1,
    'strlabelneed':0,
    'showstarburst':0,
    'legendneed':1,
    'correctIa':0,
    'runlabelneed':0,
    'withincr':0
    }

inputplotdict['figrhoTgrid_m12i']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'rhoT',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12i',
    'mtitle':'m12ires7000',
    'withinr':10.0,
    'withoutr':0.01,
    'zup':0.5,
    'zdown':-0.5
    }


inputplotdict['figrhoTgrid_m12f']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'rhoT',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12f',
    'mtitle':'m12fres7000',
    'withinr':10.0,
    'withoutr':0.01,
    'zup':0.5,
    'zdown':-0.5
    }

inputplotdict['figrhoTgrid_m12m']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'rhoT',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12m',
    'mtitle':'m12mres7000',
    'withinr':10.0,
    'withoutr':0.01,
    'zup':0.5,
    'zdown':-0.5
    }

inputplotdict['figrhoTgrid_m12']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12imhdcvhr','m12fmhdcvhr','m12mmhdcvhr','m12icr_700hr','m12fcr_700hr','m12mcr_700hr'],
    'wanted':'rhoT',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12',
    'mtitle':'',
    'withinr':10.0,
    'withoutr':0.01,
    'zup':0.5,
    'zdown':-0.5,
    'addtitleinbox':1
    }

inputplotdict['figvzzgrid_m12']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12imhdcvhr','m12fmhdcvhr','m12mmhdcvhr','m12icr_700hr','m12fcr_700hr','m12mcr_700hr'],
    'wanted':'vzz',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12',
    'mtitle':'',
    'withinr':10.0,
    'withoutr':0.01,
    'zup':10,
    'zdown':-10,
    'addtitleinbox':1,
    'normalized':0
    }

inputplotdict['figTvzgrid_m12']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12imhdcvhr','m12fmhdcvhr','m12mmhdcvhr','m12icr_700hr','m12fcr_700hr','m12mcr_700hr'],
    'wanted':'Tvz',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12',
    'mtitle':'',
    'withinr':10.0,
    'withoutr':0.01,
    'zup':10,
    'zdown':1.0,
    'addtitleinbox':1,
    'normalized':1
    }

inputplotdict['figTvznologgrid_m12']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12imhdcvhr','m12fmhdcvhr','m12mmhdcvhr','m12icr_700hr','m12fcr_700hr','m12mcr_700hr'],
    'wanted':'Tvznolog',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12',
    'mtitle':'',
    'withinr':10.0,
    'withoutr':0.01,
    'zup':10,
    'zdown':1.0,
    'addtitleinbox':1,
    'normalized':1
    }


inputplotdict['figTvznologgrid_m12i']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'Tvznolog',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12iz5_10kpc',
    'mtitle':'',
    'withinr':10.0,
    'withoutr':0.01,
    'zup':10,
    'zdown':5.0,
    'addtitleinbox':1,
    'normalized':0,
    'useverplot':1
    }

inputplotdict['figTvznologgrid_m12i_50_100']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'Tvznolog',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12iz50_100kpc',
    'mtitle':'',
    'withinr':100.0,
    'withoutr':0.01,
    'zup':100,
    'zdown':50.0,
    'addtitleinbox':1,
    'normalized':1
    }


inputplotdict['figrhozgrid_m12']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12imhdcvhr','m12fmhdcvhr','m12mmhdcvhr','m12icr_700hr','m12fcr_700hr','m12mcr_700hr'],
    'wanted':'rhoz',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12',
    'mtitle':'',
    'withinr':10.0,
    'withoutr':0.01,
    'zup':10,
    'zdown':1.0,
    'addtitleinbox':1,
    'normalized':1
    }


inputplotdict['figpcrpthgrid_m12']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12icr_700hr','m12fcr_700hr','m12mcr_700hr'],
    'wanted':'pcrpth',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12',
    'mtitle':'',
    'withinr':10.0,
    'withoutr':0.01,
    'zup':0.5,
    'zdown':-0.5,
    'addtitleinbox':1
    }

inputplotdict['figrhopcrgrid_m12']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12icr_700hr','m12fcr_700hr','m12mcr_700hr'],
    'wanted':'rhoPcr',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12',
    'mtitle':'',
    'withinr':10.0,
    'withoutr':0.01,
    'zup':0.5,
    'zdown':-0.5,
    'addtitleinbox':1
    }

inputplotdict['figrhopthgrid_m12']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12imhdcvhr','m12fmhdcvhr','m12mmhdcvhr','m12icr_700hr','m12fcr_700hr','m12mcr_700hr'],
    'wanted':'rhoPth',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12',
    'mtitle':'',
    'withinr':10.0,
    'withoutr':0.01,
    'zup':0.5,
    'zdown':-0.5,
    'addtitleinbox':1
    }


inputplotdict['figrhozgrid_m12hotgas']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12imhdcvhr','m12fmhdcvhr','m12mmhdcvhr','m12icr_700hr','m12fcr_700hr','m12mcr_700hr'],
    'wanted':'rhoz',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12hotgas',
    'mtitle':'',
    'Tcut':1e5,
    'withinr':10.0,
    'withoutr':0.01,
    'zup':10,
    'zdown':1.0,
    'addtitleinbox':1,
    'normalized':1
    }

inputplotdict['figrhoTgrid_m12out']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12imhdcvhr','m12fmhdcvhr','m12mmhdcvhr','m12icr_700hr','m12fcr_700hr','m12mcr_700hr'],
    'wanted':'rhoT',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12',
    'mtitle':'',
    'withinr':200.0,
    'withoutr':40.0,
    'zup':200.0,
    'zdown':-200.0,
    'addtitleinbox':1,
    'extendedrange':1
    }





inputplotdict['figrhoTgrid_m12out']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12imhdcvhr','m12fmhdcvhr','m12mmhdcvhr','m12icr_700hr','m12fcr_700hr','m12mcr_700hr'],
    'wanted':'rhoT',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12whole',
    'mtitle':'',
    'withinr':200.0,
    'withoutr':0.1,
    'zup':200.0,
    'zdown':-200.0,
    'addtitleinbox':1,
    'extendedrange':1
    }


inputplotdict['figrhoTgrid_m11']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m11bmhdcv','m11bcr_700'],
    'wanted':'rhoT',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m11b',
    'mtitle':'',
    'withinr':10.0,
    'withoutr':0.01,
    'zup':0.5,
    'zdown':-0.5,
    'addtitleinbox':1
    }



inputplotdict['figrhoTgrid_m12iin']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12imhdcvhr','m12icr_b_70hr','m12icr_700hr'],
    'wanted':'rhoT',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12iin',
    'mtitle':r'm12i ($r\;[{\rm kpc}] <$ 4)',
    'withinr':4.0,
    'withoutr':0.01,
    'zup':1.0,
    'zdown':-1.0
    }


inputplotdict['figTcrTgrid_m12']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12icr_700hr','m12fcr_700hr','m12mcr_700hr'],
    'wanted':'TcrT',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12',
    'mtitle':'m12ires7000',
    'withinr':10.0,
    'withoutr':0.01,
    'zup':0.5,
    'zdown':-0.5,
    'addtitleinbox':1
    }


inputplotdict['figrhoTcrgrid_m12']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12icr_700hr','m12fcr_700hr','m12mcr_700hr'],
    'wanted':'rhoTcr',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12',
    'mtitle':'m12ires7000',
    'withinr':10.0,
    'withoutr':0.01,
    'zup':0.5,
    'zdown':-0.5,
    'addtitleinbox':1
    }

inputplotdict['figrhoBgrid_sbcin']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['bwsbclrmhd','bwsbclrdc28mhd'],
    'wanted':'rhoB',
    'startno':500,
    'Nsnap':501,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'sbcin',
    'mtitle':r'sbc (0 $< r\;[{\rm kpc}] <$ 10)',
    'withinr':10.0,
    'withoutr':0.0,
    'zup':1.0,
    'zdown':-1.0
    }


inputplotdict['figrhoBgrid_smcin']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['bwsmclrmhd','bwsmclrdc28mhd'],
    'wanted':'rhoB',
    'startno':500,
    'Nsnap':501,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'smcin',
    'mtitle':r'smc (0 $< r\;[{\rm kpc}] <$ 4)',
    'withinr':10.0,
    'withoutr':0.0,
    'zup':1.0,
    'zdown':-1.0
    }


inputplotdict['figrhoTgrid_smclr']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['bwsmclrmhd','bwsmclrdc28mhd'],
    'wanted':'rhoT',
    'startno':500,
    'Nsnap':501,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'smclr',
    'mtitle':r'smc (0 $< r\;[{\rm kpc}] <$ 10)',
    'withinr':10.0,
    'withoutr':0.0,
    'zup':1.0,
    'zdown':-1.0
    }



inputplotdict['figrhoBgrid_mwin']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['bwmwmrmhd','bwmwmrdc28mhd'],
    'wanted':'rhoB',
    'startno':500,
    'Nsnap':501,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'mwin',
    'mtitle':r'm12i (0 $< r\;[{\rm kpc}] <$ 4)',
    'withinr':4.0,
    'withoutr':0.0,
    'zup':1.0,
    'zdown':-1.0
    }


inputplotdict['figrhoTgrid_mwlr']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['bwmwlrmhd','bwmwlrdc28mhd'],
    'wanted':'rhoT',
    'startno':500,
    'Nsnap':501,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'mwlr',
    'mtitle':r'm12i (0 $< r\;[{\rm kpc}] <$ 10)',
    'withinr':200.0,
    'withoutr':0.0,
    'zup':200.0,
    'zdown':-200.0,
    'extendedrange':1
    }



inputplotdict['figrhoTgrid_mwlrAd']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['bwmwlrmhd','bwmwlrdc28mhd','bwmwlrdc28mhdAd'],
    'wanted':'rhoT',
    'startno':250,
    'Nsnap':251,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'mwlrAd',
    'mtitle':'',
    'withinr':200.0,
    'withoutr':0.0,
    'zup':200.0,
    'zdown':-200.0,
    'extendedrange':1
    }


inputplotdict['figrhoBgrid_m12iin']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12imhdcvhr','m12icr_b_70hr','m12icr_700hr'],
    'wanted':'rhoB',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12iin',
    'mtitle':r'm12i (0 $< r\;[{\rm kpc}] <$ 4)',
    'withinr':4.0,
    'withoutr':0.0,
    'zup':1.0,
    'zdown':-1.0
    }



inputplotdict['figrhoBgrid_m12iout']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12imhdcvhr','m12icr_b_70hr','m12icr_700hr'],
    'wanted':'rhoB',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12iout',
    'mtitle':r'm12i (4 $< r\;[{\rm kpc}] <$ 8)',
    'withinr':8.0,
    'withoutr':4.0,
    'zup':1.0,
    'zdown':-1.0
    }


inputplotdict['figrhoBgrid_m12iRv']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12imhdcvhr','m12icr_b_70hr','m12icr_700hr'],
    'wanted':'rhoB',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12iRv',
    'mtitle':r'm12i ($r\;[{\rm kpc}] < R_{\rm vir}$)',
    'withinr':200.0,
    'withoutr':0.01,
    'zup':200.0,
    'zdown':-200.0
    }


inputplotdict['figrhoBgrid_m12iCGM']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'rhoB',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12CGM',
    'mtitle':r'm12i ($50 {\rm kpc} <r\;[{\rm kpc}] < 200 {\rm kpc}$)',
    'withinr':200.0,
    'withoutr':50.0,
    'zup':200.0,
    'zdown':-200.0
    }


inputplotdict['figrhoTgrid_m12iCGM']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'rhoT',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12iCGM',
    'mtitle':r'm12i ($50 {\rm kpc} <r\;[{\rm kpc}] < 200 {\rm kpc}$)',
    'withinr':200.0,
    'withoutr':50.0,
    'zup':200.0,
    'zdown':-200.0,
    'extendedrange':1
    }



inputplotdict['figrhoBgrid_m12mout']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'rhoB',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mout',
    'mtitle':r'm12i (4 $< r\;[{\rm kpc}] <$ 8)',
    'withinr':8.0,
    'withoutr':4.0,
    'zup':1.0,
    'zdown':-1.0
    }


inputplotdict['figrhoBgrid_m12m']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'rhoB',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mout',
    'mtitle':r'm12m (0 $< r\;[{\rm kpc}] <$ 8)',
    'withinr':8.0,
    'withoutr':0.0,
    'zup':1.0,
    'zdown':-1.0
    }

inputplotdict['figrhoBgrid_m12f']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_b_70hr','m12fcr_700hr'],
    'wanted':'rhoB',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fout',
    'mtitle':r'm12f (0 $< r\;[{\rm kpc}] <$ 8)',
    'withinr':8.0,
    'withoutr':0.0,
    'zup':1.0,
    'zdown':-1.0
    }



inputplotdict['figrhoTgrid_m12iout']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12imhdcvhr','m12icr_b_70hr','m12icr_700hr'],
    'wanted':'rhoT',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12iout',
    'mtitle':r'm12i (4 $< r\;[{\rm kpc}] <$ 8)',
    'withinr':8.0,
    'withoutr':4.0,
    'zup':1.0,
    'zdown':-1.0,
    'extendedrange':1
    }

inputplotdict['figvzzgrid_m12i']={
    '_plotfunction':rhoTgridtestinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'vzz',
    'startno':600,
    'Nsnap':601,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12ihr',
    'mtitle':'m12ires7000',
    'withinr':9.0,
    'withoutr':0.1,
    'zup':10.0,
    'zdown':-10.0,
    'normalized':0
    }

inputplotdict['figCRerzdist_m12m']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['fm12m','m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'crdenr',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mcr700hrrz',
    'withoutr':0.01,
    'withinr':8.0,
    'maxlength':0.2,
    'nogrid':10
    }

inputplotdict['figCRezrdist_m12m']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['fm12m','m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'crdenz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mcr700hrzr',
    'withoutr':0.01,
    'withinr':8.0,
    'maxlength':0.2,
    'nogrid':10
    }


inputplotdict['figCRerdist_m12m']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['fm12m','m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'crdenr',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mcr700hr',
    'withoutr':0.01,
    'withinr':8.0,
    'maxlength':0.2,
    'nogrid':10
    }


inputplotdict['figCRerdist_m12m_1kpc']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['fm12m','m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'crdenr',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mcr700hr1kpc',
    'withoutr':0.01,
    'withinr':15.0,
    'maxlength':1.0,
    'nogrid':10
    }


inputplotdict['figCRerdist_m12f']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['fm12f','m12fmhdcvhr','m12fcr_b_70hr','m12fcr_700hr'],
    'wanted':'crdenr',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fcr700hr',
    'withoutr':0.01,
    'withinr':8.0,
    'maxlength':0.2,
    'nogrid':10
    }


inputplotdict['figCRerdist_m12b']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['m12bcr_700hr'],
    'wanted':'crdenr',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12bcr700hr',
    'withoutr':0.01,
    'withinr':8.0,
    'maxlength':0.2,
    'nogrid':10
    }



inputplotdict['figCRerdist_m12r']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['m12rcr_700hr'],
    'wanted':'crdenr',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12rcr700hr',
    'withoutr':0.01,
    'withinr':8.0,
    'maxlength':0.2,
    'nogrid':10
    }


inputplotdict['figCRerdist_m12w']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['m12wcr_700hr'],
    'wanted':'crdenr',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12wcr700hr',
    'withoutr':0.01,
    'withinr':8.0,
    'maxlength':0.2,
    'nogrid':10
    }


inputplotdict['figCRerdist_m12i']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['m12imhdcvhr','m12icr_b_70hr','m12icr_700hr'],
    'wanted':'crdenr',
    'startno':600,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hr',
    'withoutr':0.01,
    'withinr':8.0,
    'maxlength':0.2,
    'nogrid':10
    }


inputplotdict['figCRerdist_m12i_1kpc']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['m12imhdcvhr','m12icr_b_70hr','m12icr_700hr'],
    'wanted':'crdenr',
    'startno':600,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hr1kpc',
    'withoutr':0.01,
    'withinr':15.0,
    'maxlength':1.0,
    'nogrid':10
    }


inputplotdict['figCRpr_volave_m12i_z4kpc']={
    '_plotfunction':CRerzdis_volave_testinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'crpr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hr_z4kpc',
    'withoutr':0.01,
    'withinr':15.0,
    'maxlength':4.0,
    'nogrid':10,
    'addsurden':1
    }

inputplotdict['figCRpr_volave_m12i']={
    '_plotfunction':CRerzdis_volave_testinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'crpr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hr',
    'withoutr':0.01,
    'withinr':15.0,
    'maxlength':1.0,
    'nogrid':10,
    'addsurden':1
    }

inputplotdict['figCRpr_volave_m12m']={
    '_plotfunction':CRerzdis_volave_testinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'crpr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mcr700hr',
    'withoutr':0.01,
    'withinr':15.0,
    'maxlength':1.0,
    'nogrid':10,
    'addsurden':1
    }

inputplotdict['figCRpr_volave_m12i']={
    '_plotfunction':CRerzdis_volave_testinput,
    'dirneed':['m12imhdcvhr','m12icr_700hr'],
    'wanted':'crpr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hr_turz',
    'withoutr':0.01,
    'withinr':15.0,
    'maxlength':1.0,
    'nogrid':10,
    'addsurden':0,
    'title':'m12ires7000'
    }


inputplotdict['figCRpr_volave_m12f']={
    '_plotfunction':CRerzdis_volave_testinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'crpr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fcr700hr_turz',
    'withoutr':0.01,
    'withinr':15.0,
    'maxlength':1.0,
    'nogrid':10,
    'addsurden':0,
    'title':'m12fres7000'
    }


inputplotdict['figCRpr_volave_m12m']={
    '_plotfunction':CRerzdis_volave_testinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'crpr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mcr700hr_turz',
    'withoutr':0.01,
    'withinr':15.0,
    'maxlength':1.0,
    'nogrid':10,
    'addsurden':0,
    'title':'m12mres7000'
    }


inputplotdict['figCRpr_volave_m12m_HI']={
    '_plotfunction':CRerzdis_volave_testinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'crpr',
    'startno':600,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mcr700hr_HI',
    'withoutr':0.01,
    'withinr':15.0,
    'maxlength':1.0,
    'nogrid':10,
    'addsurden':0,
    'title':'m12mres7000',
    'HIonly': 1
    }

inputplotdict['figCRvr_volave_m12m_HI']={
    '_plotfunction':CRerzdis_volave_testinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'crvr',
    'startno':600,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mcr700hr_HI',
    'withoutr':0.01,
    'withinr':15.0,
    'maxlength':1.0,
    'nogrid':10,
    'addsurden':0,
    'title':'m12mres7000',
    'HIonly': 1
    }


inputplotdict['figCRpr_single_m12']={
    '_plotfunction':CRerzdis_single_testinput,
    'dirneed':[['m12imhdcvhr','m12icr_700hr'],['m12fmhdcvhr','m12fcr_700hr'],['m12mmhdcvhr','m12mcr_700hr']],
    'wanted':'crpr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12hr',
    'withoutr':0.01,
    'withinr':15.0,
    'maxlength':1.0,
    'nogrid':10,
    'addsurden':0,
    'title':'m12mres7000',
    'egyneed':'eth'
    }


inputplotdict['figCRvr_single_m12']={
    '_plotfunction':CRerzdis_single_testinput,
    'dirneed':[['m12imhdcvhr','m12icr_700hr'],['m12fmhdcvhr','m12fcr_700hr'],['m12mmhdcvhr','m12mcr_700hr']],
    'wanted':'crvr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12hr',
    'withoutr':0.01,
    'withinr':15.0,
    'maxlength':1.0,
    'nogrid':10,
    'addsurden':0,
    'title':'m12mres7000',
    'egyneed':'eth'
    }


inputplotdict['figCRpz_single_m12']={
    '_plotfunction':CRerzdis_single_testinput,
    'dirneed':[['m12imhdcvhr','m12icr_700hr'],['m12fmhdcvhr','m12fcr_700hr'],['m12mmhdcvhr','m12mcr_700hr']],
    'wanted':'crpz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12hr',
    'withoutr':0.01,
    'withinr':10.0,
    'maxlength':10.0,
    'nogrid':10,
    'addsurden':0,
    'title':'m12mres7000',
    'egyneed':'eth',
    'usehalfz':1
    }

inputplotdict['figCRpturz_single_m12']={
    '_plotfunction':CRerzdis_single_testinput,
    'dirneed':[['m12imhdcvhr','m12icr_700hr'],['m12fmhdcvhr','m12fcr_700hr'],['m12mmhdcvhr','m12mcr_700hr']],
    'wanted':'crpz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12hr',
    'withoutr':0.01,
    'withinr':10.0,
    'maxlength':10.0,
    'nogrid':10,
    'addsurden':0,
    'title':'m12mres7000',
    'egyneed':'etur',
    'usehalfz':1
    }

inputplotdict['figCRptur_single_m12']={
    '_plotfunction':CRerzdis_single_testinput,
    'dirneed':[['m12imhdcvhr','m12icr_700hr'],['m12fmhdcvhr','m12fcr_700hr'],['m12mmhdcvhr','m12mcr_700hr']],
    'wanted':'crpr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12hr',
    'withoutr':0.01,
    'withinr':15.0,
    'maxlength':1.0,
    'nogrid':10,
    'addsurden':0,
    'title':'m12mres7000',
    'egyneed':'etur',
    'usehalfz':1
    }

inputplotdict['figCRvtur_single_m12']={
    '_plotfunction':CRerzdis_single_testinput,
    'dirneed':[['m12imhdcvhr','m12icr_700hr'],['m12fmhdcvhr','m12fcr_700hr'],['m12mmhdcvhr','m12mcr_700hr']],
    'wanted':'crvr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12hr',
    'withoutr':0.01,
    'withinr':15.0,
    'maxlength':1.0,
    'nogrid':10,
    'addsurden':0,
    'title':'m12mres7000',
    'egyneed':'etur',
    'usehalfz':1
    }



inputplotdict['figCRpcrz_single_m12']={
    '_plotfunction':CRerzdis_single_testinput,
    'dirneed':[['m12imhdcvhr','m12icr_700hr'],['m12fmhdcvhr','m12fcr_700hr'],['m12mmhdcvhr','m12mcr_700hr']],
    'wanted':'crpz',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12hr',
    'withoutr':0.01,
    'withinr':10.0,
    'maxlength':10.0,
    'nogrid':10,
    'addsurden':0,
    'title':'',
    'egyneed':'ecr',
    'usehalfz':1
    }


inputplotdict['figCRpcr_single_m12']={
    '_plotfunction':CRerzdis_single_testinput,
    'dirneed':[['m12imhdcvhr','m12icr_700hr'],['m12fmhdcvhr','m12fcr_700hr'],['m12mmhdcvhr','m12mcr_700hr']],
    'wanted':'crpr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12hr',
    'withoutr':0.01,
    'withinr':15.0,
    'maxlength':1.0,
    'nogrid':10,
    'addsurden':0,
    'title':'m12mres7000',
    'egyneed':'ecr',
    'usehalfz':1
    }


inputplotdict['figCRvcr_single_m12']={
    '_plotfunction':CRerzdis_single_testinput,
    'dirneed':[['m12imhdcvhr','m12icr_700hr'],['m12fmhdcvhr','m12fcr_700hr'],['m12mmhdcvhr','m12mcr_700hr']],
    'wanted':'crvr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12hr',
    'withoutr':0.01,
    'withinr':15.0,
    'maxlength':1.0,
    'nogrid':10,
    'addsurden':0,
    'title':'m12mres7000',
    'egyneed':'ecr',
    'usehalfz':1
    }


inputplotdict['figCRpB_single_m12']={
    '_plotfunction':CRerzdis_single_testinput,
    'dirneed':[['m12imhdcvhr','m12icr_700hr'],['m12fmhdcvhr','m12fcr_700hr'],['m12mmhdcvhr','m12mcr_700hr']],
    'wanted':'crpr',
    'startno':580,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12hr',
    'withoutr':0.01,
    'withinr':15.0,
    'maxlength':1.0,
    'nogrid':10,
    'addsurden':0,
    'title':'m12mres7000',
    'egyneed':'eB',
    'usehalfz':1
    }


inputplotdict['figCRpr_volave_m11b']={
    '_plotfunction':CRerzdis_volave_testinput,
    'dirneed':['m11bmhdcv','m11bcr_700'],
    'wanted':'crpr',
    'startno':600,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m11bcr700',
    'withoutr':0.01,
    'withinr':15.0,
    'maxlength':0.25,
    'nogrid':10
    }


inputplotdict['figCRpr_volave_m11f']={
    '_plotfunction':CRerzdis_volave_testinput,
    'dirneed':['m11fmhdcv','m11fcr_700'],
    'wanted':'crpr',
    'startno':600,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m11fcr700',
    'withoutr':0.01,
    'withinr':15.0,
    'maxlength':0.25,
    'nogrid':10
    }


inputplotdict['figCRpr_volave_m11g']={
    '_plotfunction':CRerzdis_volave_testinput,
    'dirneed':['m11gmhdcv','m11gcr_700'],
    'wanted':'crpr',
    'startno':600,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m11gcr700',
    'withoutr':0.01,
    'withinr':15.0,
    'maxlength':0.25,
    'nogrid':10
    }


inputplotdict['figCRpr_m12i']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['m12imhdcvhr','m12icr_b_70hr','m12icr_700hr'],
    'wanted':'crpr',
    'startno':600,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hr',
    'withoutr':0.01,
    'withinr':8.0,
    'maxlength':0.2,
    'nogrid':10
    }


inputplotdict['figCRvi_m12i']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['m12imhdcvhr','m12icr_b_70hr','m12icr_700hr'],
    'wanted':'crvr',
    'startno':600,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hr',
    'withoutr':0.01,
    'withinr':8.0,
    'maxlength':0.2,
    'nogrid':10
    }

inputplotdict['figCRvi_m12m']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'crvr',
    'startno':600,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mcr700hr',
    'withoutr':0.01,
    'withinr':8.0,
    'maxlength':0.2,
    'nogrid':10
    }

inputplotdict['figCRvi_m12f']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['m12fmhdcvhr','m12fcr_700hr'],
    'wanted':'crvr',
    'startno':600,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12bcr700hr',
    'withoutr':0.01,
    'withinr':8.0,
    'maxlength':0.2,
    'nogrid':10
    }

inputplotdict['figCRezdist_m12ihr']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['fm12i','m12imhdcvhr','m12icr_b_70hr','m12icr_700hr'],
    'wanted':'crdenz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hr',
    'withoutr':3.0,
    'withinr':8.0,
    'maxlength':2.0,
    'nogrid':10
    }

inputplotdict['figCRezdist_m12ihr_7_10']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['fm12i','m12imhdcvhr','m12icr_b_70hr','m12icr_700hr'],
    'wanted':'crdenz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hr7_10',
    'withoutr':7.0,
    'withinr':10.0,
    'maxlength':5.0,
    'nogrid':10
    }


inputplotdict['figCRezdist_m12mhrin']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'crdenz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12mcr700hrin',
    'withoutr':0.01,
    'withinr':3.0,
    'maxlength':2.0,
    'nogrid':10
    }


inputplotdict['figCRezdist_m12ihrin']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['m12imhdcvhr','m12icr_b_70hr','m12icr_700hr'],
    'wanted':'crdenz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hrin',
    'withoutr':0.01,
    'withinr':3.0,
    'maxlength':2.0,
    'nogrid':10
    }

inputplotdict['figCRezdist_mwmrin']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['bwmwmrmhd','bwmwmrdc28mhd'],
    'wanted':'crdenz',
    'startno':490,
    'Nsnap':500,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'mwmrin',
    'withoutr':0.01,
    'withinr':3.0,
    'maxlength':10.0,
    'nogrid':10
    }


inputplotdict['figCRpzdist_mwmrin']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['bwmwmrmhd','bwmwmrdc28mhd'],
    'wanted':'crpz',
    'startno':490,
    'Nsnap':500,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'mwmrin',
    'withoutr':0.01,
    'withinr':3.0,
    'maxlength':10.0,
    'nogrid':10
    }


inputplotdict['figCRdpzdist_m12ihr']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['m12icr_700hr'],
    'wanted':'crdpz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hr',
    'withoutr':3.0,
    'withinr':8.0,
    'maxlength':10.0,
    'nogrid':10
    }

inputplotdict['figCRpzdist_m12itest']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['m12icr_700hr'],
    'wanted':'crpz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700test',
    'withoutr':3.0,
    'withinr':8.0,
    'maxlength':10.0,
    'nogrid':10
    }


inputplotdict['figrhogdist_m12ihr']={
    '_plotfunction':rhogzdistestinput,
    'dirneed':['m12imhdcvhr','m12icr_b_70hr','m12icr_700hr'],
    'wanted':'crdpz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hr',
    'withoutr':0.01,
    'withinr':2.0,
    'maxlength':10.0,
    'nogrid':10,
    'usehalfz':1
    }


inputplotdict['figrhogdist_m12i70hr']={
    '_plotfunction':rhogzdistestinput,
    'dirneed':['m12icr_b_70hr'],
    'wanted':'crdpz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr70hr',
    'withoutr':3.0,
    'withinr':8.0,
    'maxlength':10.0,
    'nogrid':10
    }



inputplotdict['figrhogcoxdist_m12i700hr']={
    '_plotfunction':rhogzdistestinput,
    'dirneed':['m12imhdcvhr','m12icr_b_70hr','m12icr_700hr'],
    'wanted':'crdpzcox',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hr',
    'withoutr':5.0,
    'withinr':10.0,
    'maxlength':12.0,
    'nogrid':6,
    'usehalfz':1
    }


inputplotdict['figrhogcoxdist_m12i700hrin']={
    '_plotfunction':rhogzdistestinput,
    'dirneed':['m12imhdcvhr','m12icr_b_70hr','m12icr_700hr'],
    'wanted':'crdpzcox',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hrin',
    'withoutr':0.1,
    'withinr':3.0,
    'maxlength':12.0,
    'nogrid':6,
    'usehalfz':1
    }


inputplotdict['figrhogcoxdist_m12i700hrin1kpc']={
    '_plotfunction':rhogzdistestinput,
    'dirneed':['m12icr_700hr'],
    'wanted':'crdpzcox',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hrin1kpc',
    'withoutr':0.1,
    'withinr':3.0,
    'maxlength':12.0,
    'nogrid':6,
    'usehalfz':1
    }

inputplotdict['figrhogdist_m12fhrtest']={
    '_plotfunction':rhogzdistestinput,
    'dirneed':['m12fcr_700hr'],
    'wanted':'crdpz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fcr700hrtest',
    'withoutr':6.0,
    'withinr':8.0,
    'maxlength':10.0,
    'nogrid':10,
    'usehalfz':1
    }


inputplotdict['figrhogdist_m12fhr']={
    '_plotfunction':rhogzdistestinput,
    'dirneed':['m12fcr_700hr'],
    'wanted':'crdpz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12fcr700hr',
    'withoutr':6.0,
    'withinr':8.0,
    'maxlength':10.0,
    'nogrid':10,
    'usehalfz':1
    }


inputplotdict['figCRerdist_mwmr']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['bwmwmr','bwmwmrmhd','bwmwmrdc28','bwmwmrdc29'],
    'wanted':'crdenr',
    'startno':490,
    'Nsnap':500,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'mwmr',
    'withoutr':0.01,
    'withinr':8.0,
    'maxlength':0.2,
    'nogrid':10
    }

inputplotdict['figCRerdist_smclr']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['bwsmclr','bwsmclrmhd','bwsmclrdc28','bwsmclrdc29'],
    'wanted':'crdenr',
    'startno':490,
    'Nsnap':500,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'smclr',
    'withoutr':0.01,
    'withinr':8.0,
    'maxlength':0.2,
    'nogrid':10
    }

inputplotdict['figSFRsur_m12imhdhr']={
    '_plotfunction':star_midplane_outtestinput,
    'dirneed':['m12imhdcvhr'],
    'wanted':'SFRsurden',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'normalizedsm':0,
    'M1labelneed':0,
    'M1runlabelneed':0,
    'resoneed':0,
    'diffusionsolverneed':0,
    'refereelabelneed':0,
    'newlabelneed':1,
    'strlabelneed':0,
    'showstarburst':0,
    'legendneed':1,
    'correctIa':0,
    'runlabelneed':0,
    'withincr':0
    }

inputplotdict['figSFRsur_m12ihr']={
    '_plotfunction':star_midplane_outtestinput,
    'dirneed':['m12imhdcvhr','m12icr_b_70hr','m12icr_700hr'],
    'wanted':'SFRsurden',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'normalizedsm':0,
    'M1labelneed':0,
    'M1runlabelneed':0,
    'resoneed':0,
    'diffusionsolverneed':0,
    'refereelabelneed':0,
    'newlabelneed':1,
    'strlabelneed':0,
    'showstarburst':0,
    'legendneed':1,
    'correctIa':0,
    'runlabelneed':0,
    'withincr':0
    }


inputplotdict['figSFRsur_m12mhr']={
    '_plotfunction':star_midplane_outtestinput,
    'dirneed':['m12mmhdcvhr','m12mcr_700hr'],
    'wanted':'SFRsurden',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'normalizedsm':0,
    'M1labelneed':0,
    'M1runlabelneed':0,
    'resoneed':0,
    'diffusionsolverneed':0,
    'refereelabelneed':0,
    'newlabelneed':1,
    'strlabelneed':0,
    'showstarburst':0,
    'legendneed':1,
    'correctIa':0,
    'runlabelneed':0,
    'withincr':0
    }


inputplotdict['figkslaw_m12mhr']={
    '_plotfunction':kslawtestinput,
    #'dirneed':['m12fmhdcvhr','m12fcr_700hr',\
    #          'm12imhdcvhr','m12icr_700hr','m12mmhdcvhr','m12mcr_700hr'],
    'dirneed':['fm12m','m12mcr_700hr','m12mmhdcvhr'],
    'wanted':'kslaw',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'legendneed':1,
    'withinr':10,
    'maxlength':10,
    'nogrid':10
    }


inputplotdict['figkslaw_m12ihr']={
    '_plotfunction':kslawtestinput,
    #'dirneed':['m12fmhdcvhr','m12fcr_700hr',\
    #          'm12imhdcvhr','m12icr_700hr','m12mmhdcvhr','m12mcr_700hr'],
    'dirneed':['fm12i','m12icr_700hr','m12imhdcvhr'],
    'wanted':'kslaw',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'legendneed':1,
    'withinr':10,
    'maxlength':10,
    'nogrid':10
    }


inputplotdict['figkslaw_m12ihrdc28']={
    '_plotfunction':kslawtestinput,
    #'dirneed':['m12fmhdcvhr','m12fcr_700hr',\
    #          'm12imhdcvhr','m12icr_700hr','m12mmhdcvhr','m12mcr_700hr'],
    'dirneed':['fm12i','m12icr_b_70hr','m12icr_700hr','m12imhdcvhr'],
    'wanted':'kslaw',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'dc28',
    'legendneed':1,
    'withinr':10,
    'maxlength':10,
    'nogrid':10
    }


inputplotdict['figkslaw_m12fhr']={
    '_plotfunction':kslawtestinput,
    #'dirneed':['m12fmhdcvhr','m12fcr_700hr',\
    #          'm12imhdcvhr','m12icr_700hr','m12mmhdcvhr','m12mcr_700hr'],
    'dirneed':['fm12f','m12fcr_700hr','m12fmhdcvhr'],
    'wanted':'kslaw',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'legendneed':1,
    'withinr':10,
    'maxlength':10,
    'nogrid':10
    }

inputplotdict['figkslaw_m11_m12']={
    '_plotfunction':kslawtestinput,
    'dirneed':['m11gcr_700','m11gmhdcv','m11dmhdcv','m11hmhdcv','m12bcr_700hr','m12wcr_700hr','m12rcr_700hr','m11bmhdcv',\
               'm11bcr_700','m11fmhdcv','m11fcr_700','fm11q',\
               'f476','f61','f46','fm12f','fm12i','fm12m','m12fmhdcvhr',\
               'm12fcr_700hr','m12mmhdcvhr','m12mcr_700hr',\
               'm12imhdcvhr','m12icr_700hr','f573'],
    'wanted':'kslaw',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'legendneed':1,
    'withinr':15,
    'withoutr':10,
    'maxlength':10,
    'nogrid':10
    }

inputplotdict['figkslaw_m11_m12']={
    '_plotfunction':kslawtestinput,
    'dirneed':['m11gcr_700','m11gmhdcv','m11dmhdcv','m11hmhdcv','m12bcr_700hr','m12wcr_700hr','m12rcr_700hr','m11bmhdcv',\
               'm11bcr_700','m11fmhdcv','m11fcr_700','fm11q',\
               'f476','f61','f46','fm12f','fm12i','fm12m','m12fmhdcvhr',\
               'm12fcr_700hr','m12mmhdcvhr','m12mcr_700hr',\
               'm12imhdcvhr','m12icr_700hr','f573'],
    'wanted':'kslaw',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'legendneed':1,
    'withinr':10,
    'withoutr':-10,
    'maxlength':10,
    'nogrid':10
    }




inputplotdict['figkslaw_m11_m12_out']={
    '_plotfunction':kslawtestinput,
    'dirneed':['m11gcr_700','m11gmhdcv','m11dmhdcv','m11hmhdcv','m12bcr_700hr','m12wcr_700hr','m12rcr_700hr','m11bmhdcv',\
               'm11bcr_700','m11fmhdcv','m11fcr_700','fm11q',\
               'f476','f61','f46','fm12f','fm12i','fm12m','m12fmhdcvhr',\
               'm12fcr_700hr','m12mmhdcvhr','m12mcr_700hr',\
               'm12imhdcvhr','m12icr_700hr','f573'],
    'wanted':'kslaw',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'gtr8kpc',
    'legendneed':1,
    'withinr':18,
    'withoutr':8,
    'maxlength':10,
    'nogrid':10
    }


inputplotdict['figsyn_m12i']={
    '_plotfunction':synchrotrontestinput,
    'dirneed':['m12imhdcvhr','m12icr_b_70hr','m12icr_700hr'],
    'wanted':'syncr',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':0.2,
    'withinr': 20.0,
    'withoutr': 0.1
    }

inputplotdict['figBmin_m12i']={
    '_plotfunction':synchrotrontestinput,
    'dirneed':['m12icr_700hr'],
    'wanted':'Bmin',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':0.2,
    'withinr': 20.0,
    'withoutr': 0.1
    }


inputplotdict['figfrm_m12i']={
    '_plotfunction':Faradayrotationtestinput,
    'dirneed':['m12imhdcvhr','m12icr_b_70hr','m12icr_700hr'],
    'wanted':'Faradayrotation',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'',
    'maxlength':10.0,
    'withinr': 20.0,
    'withoutr': 0.1
    }


inputplotdict['figCRcox_m12ihrin']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['m12imhdcvhr','m12icr_b_70hr','m12icr_700hr'],
    'wanted':'crpzcox',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hrin',
    'withoutr':0.01,
    'withinr':3.0,
    'maxlength':2.0,
    'nogrid':10
    }

inputplotdict['figCRcox_m12ihrout']={
    '_plotfunction':CRerzdistestinput,
    'dirneed':['m12imhdcvhr','m12icr_b_70hr','m12icr_700hr'],
    'wanted':'crpzcox',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':'m12icr700hrout',
    'withoutr':3.0,
    'withinr':6.0,
    'maxlength':2.0,
    'nogrid':10
    }


presetkeylist = [['withoutrnew',-1],['the_prefix','snapshot'],['the_suffix','.hdf5'],['usehalfz',0],['title',''],['addsurden',0],['egyneed',''],['addtitleinbox',0],['minsfr',1],['maxsfr',20],['extendedrange',0],['usepredata',0],['griddir','grid1kpc'],['cutcold',0],['parallel',0],['usekez',0],['normalized',1],['Tcut',-1],['withoutr',-1.0],['HIonly',0],['outHI',0],['usesnipshot',0],['extendedrange',0],['outz',0],['zup',1.0],['zdown',0.1],['phicut',0],['usehz',1],['useobs',0],['useverplot',0],['hotwarmmode',0],['trackstart',0]]

for key in list(inputplotdict):
    for pkey in presetkeylist:
        if not checkkey(inputplotdict[key],pkey[0]): inputplotdict[key][pkey[0]]=pkey[1]
    if inputplotdict[key]['fmeat']=='':
        inputplotdict[key]['fmeat']=inputplotdict[key]['dirneed'][-1]

