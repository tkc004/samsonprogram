import os

#plotnamelist=['gf12fmhdcv']

plotnamelist=['gf12fmhdcv','gf12fcr_700','gpm12fmhdcv','gpm12fcr_700',\
              'gf12mmhdcv','gf12mcr_700','gpm12mmhdcv','gpm12mcr_700']
for plotname in plotnamelist:
    os.system('qsub -v PLOTNAME={'+plotname+'} job.pbs')