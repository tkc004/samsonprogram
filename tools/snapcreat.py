import sys
import os
Nsnap=601
count=0
while(count < Nsnap):
        strno = str(count)
        if count<10:
                strno='00'+str(count)
        elif count<100:
                strno='0'+str(count)
        print 'doing ',strno
        mycommand = 'mkdir snapdir_'+strno
	mycommand2 = 'mv ../output/snapdir_'+strno+'/snap_convertedshot_'+strno+'.*  snapdir_'+strno
        os.system(mycommand)
	os.system(mycommand2)
        count = count +1
