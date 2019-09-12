import sys
import os

#runlist=[
#    'outssm12m580','outssm12m590','outssm12i580','outssm12i590','outssm12f580', 'outssm12f590',
#]

#runlist=['outssm12f580', 'outssm12f590']

#runlist=['gf12fmhdcvgrid0_25kpc590','gf12fmhdcvgrid0_25kpc595','gf12fcr_700grid0_25kpc595','gf12mmhdcvgrid0_25kpc593','gf12mmhdcvgrid0_25kpc596','gf12mcr_700grid0_25kpc590','gf12mcr_700grid0_25kpc595']
    
    
    
#runlist=[\
#         'outssm11g580',\
#         'outssm11g590',\
#         'outssm11d580',\
#         'outssm11d590',\
#         'outssm11h580',\
#         'outssm11h590',\
#         'outssm11f580',\
#         'outssm11f590',\
#         'outssm11b580',\
#         'outssm11b590',\
#        ]


#runlist=[\
#         'gpm11bgrid0_5kpcHI',\
#         'gpm11fgrid0_5kpcHI',\
#         'gpm11hgrid0_5kpcHI',\
#         'gpm11dgrid0_5kpcHI',\
#         'gpm11ggrid0_5kpcHI',\
#        ]
 
    
#runlist=[\
#         'gf11bgrid0_5kpc',\
#         'gf11fgrid0_5kpc',\
#         'gf11hgrid0_5kpc',\
#         'gf11dgrid0_5kpc',\
#         'gf11ggrid0_5kpc',\
#        ]


    
#runlist=[\
#         'gpm12imhdcvgrid0_25kpcHI',\
#         'gpm12icr_700grid0_25kpcHI',\

#         'gf12imhdcvgrid0_25kpc',\
#         'gf12icr_700grid0_25kpc',\
#         'gf12fmhdcvgrid0_25kpc',\
#         'gf12fcr_700grid0_25kpc',\
#         'gf12mmhdcvgrid0_25kpc',\
#         'gf12mcr_700grid0_25kpc',\
         
#         ]         
         

#runlist=[\
#         'outssm12i580',\
#         'outssm12i590',\
#         'outssm12f580',\
#         'outssm12f590',\
#         'outssm12m580',\
#         'outssm12m590',\
#        ]

runlist=[\
         'gpm12mmhdcvgrid0_25kpc580HI',\
         'gpm12mmhdcvgrid0_25kpc585HI',\
         'gpm12mmhdcvgrid0_25kpc590HI',\
         'gpm12mmhdcvgrid0_25kpc595HI',\
         'gpm12mcr_700grid0_25kpc580HI',\
         'gpm12mcr_700grid0_25kpc585HI',\
         'gpm12mcr_700grid0_25kpc590HI',\
         'gpm12mcr_700grid0_25kpc595HI',\
         'gpm12fmhdcvgrid0_25kpc580HI',\
         'gpm12fmhdcvgrid0_25kpc585HI',\
         'gpm12fmhdcvgrid0_25kpc590HI',\
         'gpm12fmhdcvgrid0_25kpc595HI',\
         'gpm12fcr_700grid0_25kpc580HI',\
         'gpm12fcr_700grid0_25kpc585HI',\
         'gpm12fcr_700grid0_25kpc590HI',\
         'gpm12fcr_700grid0_25kpc595HI',\
         #'gpm12imhdcvgrid0_25kpcHI',\
         #'gpm12imhdcvgrid0_25kpc590HI',\
         #'gpm12icr_700grid0_25kpcHI',\
         'gpm12icr_700grid0_25kpc590HI',\
         #'gpm12fmhdcvgrid0_25kpcHI',\
         'gpm12fmhdcvgrid0_25kpc590HI',\
         #'gpm12fcr_700grid0_25kpcHI',\
         #'gpm12fcr_700grid0_25kpc590HI',\
         #'gpm12mmhdcvgrid0_25kpcHI',\
#         'gpm12mmhdcvgrid0_25kpc590HI',\
#         'gpm12mcr_700grid0_25kpcHI',\
         #'gpm12mcr_700grid0_25kpc590HI'
        ]

#runlist=[\
#         'gpm12imhdcvgrid0_5kpcHI',\
#         'gpm12icr_700grid0_5kpcHI',\
#         'gpm12fmhdcvgrid0_5kpcHI',\
#         'gpm12fcr_700grid0_5kpcHI',\
#         'gpm12mmhdcvgrid0_5kpcHI',\
#         'gpm12mcr_700grid0_5kpcHI',\
#        ]


#runlist = [\
#         'gpm12fmhdcvHI',\
#          'gpm12fcr_700HI',\
#           'gpm12mmhdcvHI',\
#           'gpm12mcr_700HI',\
#           'gpm12imhdcvHI',\
#           'gpm12fcr_700HI',\
#          ]

#runlist = [\
#           'gf12mmhdcvgrid0_125kpc600'
#          ]

#runlist = [\
#           'gpm12fcr_700grid0_125kpc600',\
#           'gpm12fmhdcvgrid0_125kpc600',\
#           'gpm12fcr_700grid0_125kpc600',\
#          ]


#runlist = [\
#           'gpm12fmhdcvgrid1kpccc',\
#           'gpm12fcr_700grid1kpccc',\
#           'gpm12mmhdcvgrid1kpccc',\
#           'gpm12mcr_700grid1kpccc',\
#           'gpm12imhdcvgrid1kpccc',\
#           'gpm12fmhdcvgrid1kpccc',\
#          ]


#         'gf12imhdcv',\
#         'gf12icr_700',\
#         'gpm12imhdcv',\
#         'gpm12icr_700',\
#         'gf12fmhdcv',\
#         'gf12fcr_700',\
#         'gpm12fmhdcv',\
#         'gpm12fcr_700',\
#         'gf12mmhdcv',\
#         'gf12mcr_700',\
#         'gpm12mmhdcv',\
#         'gpm12mcr_700',\
           
           

#runlist=[\
#         'gf12imhdcvgrid2kpc',\
#         'gf12icr_700grid2kpc',\
#         'gpm12imhdcvgrid2kpc',\
#         'gpm12icr_700grid2kpc',\
#         'gf12fmhdcvgrid2kpc',\
#         'gf12fcr_700grid2kpc',\
#         'gpm12fmhdcvgrid2kpc',\
#         'gpm12fcr_700grid2kpc',\
#         'gf12mmhdcvgrid2kpc',\
#         'gf12mcr_700grid2kpc',\
#         'gpm12mmhdcvgrid2kpc',\
#         'gpm12mcr_700grid2kpc',\
#        ]


#runlist=[\
#         'gf12imhdcvgrid0_5kpc',\
#         'gf12icr_700grid0_5kpc',\
#         'gpm12imhdcvgrid0_5kpc',\
#         'gpm12icr_700grid0_5kpc',\
#         'gf12fmhdcvgrid0_5kpc',\
#         'gf12fcr_700grid0_5kpc',\
#         'gpm12fmhdcvgrid0_5kpc',\
#         'gpm12fcr_700grid0_5kpc',\
#         'gf12mmhdcvgrid0_5kpc',\
#         'gf12mcr_700grid0_5kpc',\
#         'gpm12mmhdcvgrid0_5kpc',\
#         'gpm12mcr_700grid0_5kpc',\
#        ]

#runlist=[\
#         'gf12imhdcvgrid0_25kpc',\
#         'gf12icr_700grid0_25kpc',\
#         'gpm12imhdcvgrid0_25kpc',\
#         'gpm12icr_700grid0_25kpc',\
#         'gf12fmhdcvgrid0_25kpc',\
#         'gf12fcr_700grid0_25kpc',\
#         'gpm12fmhdcvgrid0_25kpc',\
#         'gpm12fcr_700grid0_25kpc',\
#         'gf12mmhdcvgrid0_25kpc',\
#         'gf12mcr_700grid0_25kpc',\
#         'gpm12mmhdcvgrid0_25kpc',\
#         'gpm12mcr_700grid0_25kpc',\
         
#         ]

#runlist=[\
#        'gpm12mcr_700grid0_125kpc600',\
#        'gpm12mmhdcvgrid0_125kpc600',\
#        'gf12mcr_700grid0_125kpc600',\
#         'gf12mmhdcvgrid0_125kpc600',\
#         'gpm12fcr_700grid0_125kpc600',\
#         'gpm12fmhdcvgrid0_125kpc600',\
#         'gf12fcr_700grid0_125kpc600',\
#         'gf12fmhdcvgrid0_125kpc600',\
#         'gpm12icr_700grid0_125kpc600',\
#         'gpm12imhdcvgrid0_125kpc600',\
#         'gf12icr_700grid0_125kpc600',\
#         'gf12imhdcvgrid0_125kpc600',\
#        ]     

#runlist=[\
#        'gpm12mcr_700grid0_125kpc600HI',\
#        'gpm12mmhdcvgrid0_125kpc600HI',\
#        'gf12mcr_700grid0_125kpc600',\
#         'gf12mmhdcvgrid0_125kpc600',\
#         'gpm12fcr_700grid0_125kpc600HI',\
#         'gpm12fmhdcvgrid0_125kpc600HI',\
#         'gf12fcr_700grid0_125kpc600',\
#         'gf12fmhdcvgrid0_125kpc600',\
#         'gpm12icr_700grid0_125kpc600HI',\
#         'gpm12imhdcvgrid0_125kpc600HI',\
#         'gf12icr_700grid0_125kpc600',\
#         'gf12imhdcvgrid0_125kpc600',\
#        ] 

#runlist=[\
#        'gpm12mcr_700grid0_125kpc590HI',\
#        'gpm12mmhdcvgrid0_125kpc590HI',\
#        'gf12mcr_700grid0_125kpc590',\
#         'gf12mmhdcvgrid0_125kpc590',\
#         'gpm12fcr_700grid0_125kpc590HI',\
#         'gpm12fmhdcvgrid0_125kpc590HI',\
#         'gf12fcr_700grid0_125kpc590',\
#         'gf12fmhdcvgrid0_125kpc590',\
#         'gpm12icr_700grid0_125kpc590HI',\
#         'gpm12imhdcvgrid0_125kpc6590HI',\
#         'gf12icr_700grid0_125kpc590',\
#         'gf12imhdcvgrid0_125kpc590',\
#        ] 


#runlist=[\
#        'gpm12mcr_700grid0_125kpc580HI',\
#        'gpm12mmhdcvgrid0_125kpc580HI',\
#        'gf12mcr_700grid0_125kpc580',\
#         'gf12mmhdcvgrid0_125kpc580',\
#         'gpm12fcr_700grid0_125kpc580HI',\
#         'gpm12fmhdcvgrid0_125kpc580HI',\
#         'gf12fcr_700grid0_125kpc580',\
#         'gf12fmhdcvgrid0_125kpc580',\
#         'gpm12icr_700grid0_125kpc580HI',\
#         'gpm12imhdcvgrid0_125kpc580HI',\
#         'gf12icr_700grid0_125kpc580',\
#         'gf12imhdcvgrid0_125kpc580',\
#        ] 


#runlist = [\
#           'gpm12fmhdcvgrid0_125kpc590',\
#          ]

#runlist=[\
#        'gpm12mcr_700grid0_125kpc580',\
#        'gpm12mmhdcvgrid0_125kpc580',\
#        'gf12mcr_700grid0_125kpc580',\
#         'gf12mmhdcvgrid0_125kpc580',\
#         'gpm12fcr_700grid0_125kpc580',\
#         'gpm12fmhdcvgrid0_125kpc580',\
#         'gf12fcr_700grid0_125kpc580',\
#         'gf12fmhdcvgrid0_125kpc580',\
#         'gpm12icr_700grid0_125kpc580',\
#         'gpm12imhdcvgrid0_125kpc580',\
#         'gf12icr_700grid0_125kpc580',\
#         'gf12imhdcvgrid0_125kpc580',\
#        ]  
    
 

#runlist = ['gpm12fmhdcvgrid0_125kpc580']


for runname in runlist:
    finname = 'job.pbs'
    f=open(finname)
    dars = f.readlines()
    f.close()
    count = 0

    foutname = 'tempjob.pbs'
    g = open(foutname, 'w')
    for line in dars:
        g.write(line+'\n')
    g.write('#PBS -N '+runname+'\n')
    g.write('cd /home/tkc004/samsonprogram/tools/ISMmainplotfunction '+'\n')    
    g.write('pwd'+'\n')
    g.write('python testinputplotdataISMforjob.py -f '+runname+'\n')
    g.close()
    
    os.system('qsub tempjob.pbs')


