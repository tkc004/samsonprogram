import numpy as np
def readfitdat(haloname, beginno, finalno, xmax_of_box, xno, theta='045',movdir='hdf5',silent=0,fixn=0,linsep=13):
	haloname=str(haloname)
	icount=0
	tnowl=[]
	gil=[]
	sagel=[]
	smassl=[]
	FeHl=[]
	xmax_of_boxl=[]
	xnol=[]
	xgalkpcl=[]
	ygalkpcl=[]
	surbril=[]
	reffl=[]
	sinl=[]
	axisrl=[]
	angl=[]
	gmagl=[]
	xgalkpc=500
	ygalkpc=500
	surbri=25.
	reff=2.
	ang=30.
	axisr=0.8
	sin=1.
	for i in range(beginno,finalno,1):
		flog = open('/home/tkc004/scratch/galfit-example/EXAMPLE/'+haloname+'/fit.log', 'r')
		dars = flog.readlines()
		flog.close()
		linecount = 0
		for line in dars:
			linecount = linecount + 1
			if (linecount == 8+linsep*icount):
				xsd = line.split()
				if silent==0:
					print 'linecount', linecount
					print 'xsd', xsd
		filemeat = movdir+'_s0'+str(i)+'_t'+theta+'_star'
		filenametxt='/home/tkc004/scratch/test_movie_maker/'+haloname+'r/'+movdir+'/'+filemeat+'.txt'
		ftxt = open(filenametxt, 'r')
		ftxt.readline()
		linetxt = ftxt.readline()
		ftxt.close()
		xsdtxt = linetxt.split()
		tnow = float(xsdtxt[0])
		gi = float(xsdtxt[1])
		sage = float(xsdtxt[2])
		smass = float(xsdtxt[3])
		FeH = float(xsdtxt[4])
		gmag =float(xsdtxt[5])
		try:
			reff = float(xsd[6])*xmax_of_box/xno
			surbri = float(xsd[5])
			if fixn==0:
				sin = float(xsd[7])
			else:
				sin =fixn
			axisr = float(xsd[8])
			ang = float(xsd[9])
			xgals = xsd[3]
			ygals = xsd[4]
			xgal = float(xgals[:-1])
			ygal = float(ygals[:-1])
			xgalkpc = float(xgal)*xmax_of_box/xno
			ygalkpc = float(ygal)*xmax_of_box/xno
		except (ValueError,IndexError):
			if silent==0:
				print 'some values are not fitted'
		surbriguess=surbri
		reffguess=reff
		angguess=ang
		aratguess=axisr
		nguess=sin
		icount = icount +1
		tnowl=np.append(tnowl,tnow)
		gil=np.append(gil,gi)
		sagel=np.append(sagel,sage)
		smassl=np.append(smassl,smass)
		FeHl=np.append(FeHl,FeH)
		xmax_of_boxl=np.append(xmax_of_boxl,xmax_of_box)
		xnol=np.append(xnol,xno)
		xgalkpcl=np.append(xgalkpcl,xgalkpc)
		ygalkpcl=np.append(ygalkpcl,ygalkpc)
		surbril=np.append(surbril,surbri)
		reffl=np.append(reffl,reff)
		sinl=np.append(sinl,sin)
		axisrl=np.append(axisrl,axisr)
		angl=np.append(angl,ang)
		gmagl=np.append(gmagl,gmag)
	return tnowl, gil, sagel, smassl, FeHl, xmax_of_boxl, xnol, xgalkpcl, ygalkpcl, surbril, reffl, sinl, axisrl, angl, gmagl


