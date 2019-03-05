import collections
from importplotfunctionsISM import *

nested_dict = lambda: collections.defaultdict(nested_dict)
inputplotdict = nested_dict()


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
    '_plotfunction':CRezdistestinput,
    'dirneed':['m12mmhdcv','m12mcr_b_70','m12mcr_700'],
    'wanted':'crdenz',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':''
    }


inputplotdict['figanHmid_m12m']={
    '_plotfunction':nHmidplanetestinput,
    'dirneed':['m12mmhdcv','m12mcr_b_70','m12mcr_700'],
    'wanted':'nHmidplane',
    'startno':590,
    'Nsnap':600,
    'snapsep':10,
    'the_prefix':'snapshot',
    'the_suffix':'.hdf5',
    'fmeat':''
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



for key in inputplotdict:
    if inputplotdict[key]['fmeat']=='':
        inputplotdict[key]['fmeat']=inputplotdict[key]['dirneed'][-1]

