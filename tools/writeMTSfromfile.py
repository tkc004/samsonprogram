import os
import numpy as np
filelist=np.array(os.listdir('Pep/'))
filelist=sorted(filelist, key=str.lower)
particlelist=[]
print 'filelist', filelist
for filename in filelist:
	if 'particle' in filename:
		particlelist.append(filename)
print 'particlelist', particlelist
outname = 'MTSnaps.txt'
MTS= open(outname, "w")
for filename in particlelist:
        MTS.write('Pep/'+filename[0:22]+'\n')
MTS.close()
