import sys
import os
import os.path
Nsnap=601
count=0
single=1
while(count < Nsnap):
	tailno=0
	while (tailno < 4):
		strno = str(count)
		if count<10:
			strno='00'+str(count)
		elif count<100:
			strno='0'+str(count)
		if single==0:
			fname = 'snapdir_'+strno+'/snapshot_'+strno+'.'+str(tailno)+'.hdf5'
		else:
			fname = 'snapshot_'+strno+'.hdf5'
		if not os.path.exists(fname):
			print fname
		tailno+=1
	count = count +1
