import sys
import os

finname = 'iclist.txt'
f = open(finname)
dars = f.readlines()
f.close()
alist = []
Nlist = []

Nhalo = 0
for line in dars:
	xsd = line.split()
	Nlist.append(str(xsd[0]))
	Nhalo = Nhalo + 1

finext  = '.conf'
finpre  = 'ics_l8l_zbare'
finname = finpre+finext
f=open(finname)
dars = f.readlines()
f.close()
count = 0

inputlist = []

while (count < len(Nlist)):
	print '---------------------------- \n'*4 
	print 'DOING N ',Nlist[count]
	foutname = 'zooms_3rv/'+finpre + Nlist[count] + finext
	g = open(foutname, 'w')
	linecount = 0
	for line in dars:
		linecount = linecount + 1
		if (linecount == 21):
			print 'region_point_file       = zooms_3rv/parlist_halo_'+Nlist[count]+'.dat'
			g.write('region_point_file       = zooms_3rv/parlist_halo_'+Nlist[count]+'.dat'+'\n')
		elif (linecount == 57):
			print 'filename                = zooms_3rv/ics_halo_'+Nlist[count]+'.dat'
			g.write('filename                = zooms_3rv/ics_halo_'+Nlist[count]+'.dat'+'\n')
		else:
			print line.strip()
			g.write(line)
	g.close()
	count = count + 1
	inputlist.append(foutname)

count = 0
while(count < len(inputlist)):
	print 'doing ',inputlist[count]
	mycommand = './MUSIC zooms_3rv/ics_l8l_zbare'+Nlist[count]+'.conf'
	os.system(mycommand)
	count = count +1
