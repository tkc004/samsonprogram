#This python script is to generate a much finer output time for cosmological run
from gadget_lib.cosmo import *
#z=1.0
#t=quick_lookback_time(z)
#a=quick_scale_factor(t)
#print 't,z', t,1./a-1.
z=0.1694600
t=quick_lookback_time(z)
a=quick_scale_factor(t)
print 'a', a
tlist=t+np.linspace(0,0.2,num=200)
print 'tlist', tlist
alist=quick_scale_factor(tlist)
print 'alist', alist

finname = 'snapshot_scale-factors.txt'
f=open(finname)
dars = f.readlines()
f.close()
count = 0

inputlist = []

Nsnap=500
g = open('snapshot_scale-factors_smallt.txt', 'w')
linecount = 0
for line in dars:
	linecount = linecount + 1
	g.write(line)
	if linecount>Nsnap-1:
		break
for ncount in range(len(alist)):
	g.write(format(alist[ncount],'.7f')+'\n')
g.close()

