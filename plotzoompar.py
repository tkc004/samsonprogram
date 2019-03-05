import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

scatterplot=0
#spmname='/home/tkc004/scratch/m11b_res260_dmo_6rv/lagregion/parlist_halo_dm.dat'
spmname='/home/tkc004/scratch/m11b_res260_dmo_f6rv/zooms_8rv/parlist_halo.dat'
#spmname='/home/tkc004/scratch/m11b_res260_dmo_f6rv/zooms_7rv/parlist_halo.dat'
#spmname='/home/tkc004/scratch/m11b_res260_dmo_6rv/zooms_6rv/parlist_halo.dat'
#spmname='/home/tkc004/scratch/bare_l8l_10lv/zooms_6rv/parlist_halo_476.dat'
spmfile=open(spmname,"r")
spmfile.readline()
dars = spmfile.readlines()
xcor=[]
ycor=[]
zcor=[]
for line in dars:
	xsd =  line.split()
	xcor = np.append(xcor, float(xsd[0]))
	ycor = np.append(ycor, float(xsd[1]))
	zcor = np.append(zcor, float(xsd[2]))
spmfile.close()

if scatterplot==1:

	plt.scatter(xcor,ycor)
	#plt.imshow(a, cmap='hot', interpolation='nearest')
	plt.xlabel('x')
	plt.ylabel('y')
	plt.savefig('figures/zooms_xy.png')
	plt.clf()

	plt.scatter(xcor,zcor)
	plt.xlabel('x')
	plt.ylabel('z')
	plt.savefig('figures/zooms_xz.png')
	plt.clf()

	plt.scatter(ycor,zcor)
	plt.xlabel('y')
	plt.ylabel('z')
	plt.savefig('figures/zooms_yz.png')
	plt.clf()

else:
	xedges=np.linspace(0.,1.,num=1000)
	yedges=np.linspace(0.,1.,num=1000)
	#print 'xedges', xedges
	#print 'yedges', yedges
	#histogram reverse x y coordinate
	Hxy, xedges, yedges = np.histogram2d(ycor, xcor, bins=(xedges, yedges))
	plt.imshow(Hxy, interpolation='nearest', origin='low',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
	plt.xlabel('x')
	plt.ylabel('y')
	plt.savefig('figures/zooms_xy_his.png')
	plt.clf()

	Hxz, xedges, yedges = np.histogram2d(zcor, xcor, bins=(xedges, yedges))
	plt.imshow(Hxz, interpolation='nearest', origin='low',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
	plt.xlabel('x')
	plt.ylabel('z')
	plt.savefig('figures/zooms_xz_his.png')
	plt.clf()

	Hyz, xedges, yedges = np.histogram2d(zcor, ycor, bins=(xedges, yedges))
	plt.imshow(Hyz, interpolation='nearest', origin='low',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
	plt.xlabel('y')
	plt.ylabel('z')
	plt.savefig('figures/zooms_yz_his.png')
	plt.clf()
