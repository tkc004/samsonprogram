import sys
import os
import os.path
Nsnap=601
count=0
while(count < Nsnap):
	tailno=0
	while (tailno < 4):
		strno = str(count)
		if count<10:
			strno='00'+str(count)
		elif count<100:
			strno='0'+str(count)
		fname = 'snapdir_'+strno+'/snap_convertedshot_'+strno+'.'+str(tailno)
		if not os.path.exists(fname):
			print fname
		tailno+=1
	count = count +1
