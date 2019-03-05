import numpy as np

spmname='snapshot_times.txt'
spmfile=open(spmname,"r")
spmfile.readline()
spmfile.readline()
spmfile.readline()
dars = spmfile.readlines()

snapno=[]
redshift=[]

for line in dars:
	xsd = line.split()
	snapno = np.append(snapno, int(xsd[0]))
	redshift = np.append(redshift, float(xsd[2]))
	print 'xsd', xsd
spmfile.close()
snapno=np.array(snapno)
redshift=np.array(redshift)
print 'snapno', snapno
outname = 'MTSnaps.txt'
MTS= open(outname, "w")
for icount in range(len(snapno)):
	if snapno[icount]>-1:
		strno=str(int(snapno[icount]))
		if icount <10:
			strno =  '00' + strno
		elif icount<100:
			strno = '0'+strno
		if strno=='250':
			strz = '1.187'
		elif strno=='000':
			strz = '99.013'
		elif strno=='024':
			strz = '9.313'
		else:
			strz="%0.3f" % redshift[icount]
		MTS.write('Pep/snap'+strno+'RPep.z'+strz+'.AHF'+'\n')
MTS.close()
