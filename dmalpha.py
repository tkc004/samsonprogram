from readsnap_samson import *
from Sasha_functions import *
import sys
import matplotlib as mpl
from pathloc import *
mpl.use('Agg')
import matplotlib.pyplot as plt
from gadget_lib.cosmo import *
from samson_functions import *
mpl.use('Agg')
from matplotlib import rcParams
from dmalpha_tools import *
rcParams['axes.linewidth'] = 2
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True
#rcParams['text.usetex'] = True
rcParams['axes.unicode_minus']=False

rcParams.update({'figure.autolayout': True})
#dirneed=['fm11']
#tlist=[6.0]
usesnapno=0
#tlist=[400]
tlist=[13.9]
dirnr=[]
dirhr=[]
#dirneed=['m09','m10','1297','1146','573','553','476','m11','383']
#dirneed=['f1146']
#dirnr=['fm10q','fm11q','fm11v','f1146','f573','f553','f476','fm11','f383','f61']
#dirnr = ['fm10q','fm11q','fm11v','f1146','f573','f553','f476','fm11','f383','f61','fm12c','fm12f','fm12i','fm12m','fm12q']
#dirhr = ['fm11qCRK100','fm11qCRK1000']
#dirnr = ['fm11d','m11dmhdcv','m11dcr_b_70','m11hmhdcv','m11hcr_b_70','m11gcr_b_70','m11gmhdcv','m11fmhdcv', 'm11fcr_b_70','fm10q','m10qmhdcv','m10qcr_b_70','fm11q','f476','m11bmhdcv', 'm11bcr_b_70',\
#'m12imhdcv','m12icr_b_70','fm11v','f1146','f573','f553','fm11','f383','f61']
#dirnr = ['m11fmhdcv','m11fcr_b_70','fm11q','f476','fm11v','f1146','f573','f553','fm11','f383','f61']
#dirnr = ['fm11v','fm11d','f383','f553','fm11','m10vcr_700', 'm11bcr_700','m11dcr_700','m11fcr_700','m11gcr_700','m11hcr_700','m12icr_700','fm11q','fm11v','f1146','f61']
#dirnr = ['m11hcr_700','m12icr_700']
#dirnr = ['f553','f476']
#dirnr = ['fm10qmd','fm10vmd','fm11dmd','fm11emd', 'fm11hmd','fm11imd', 'fm11qmd','fm12fmd','fm12imd'\
#,'fm10q','f476','f383','f553','fm11q','fm11v','f1146','f61','fm12b']
#dirnr = [
#'fm10qmd','fm10vmd','fm11dmd','fm11emd', 'fm11hmd','fm11imd', 'fm11qmd','fm12fmd','fm12imd',\
#'m10qmhdcv','m12imhdcv','m10vmhdcv', 'm11bmhdcv','m11cmhdcv','m11dmhdcv','m11fmhdcv',\
#'m11hmhdcv','m11gmhdcv','m11v1mhdcv','m11v2mhdcv',\
#'m10qcr_b_70','m10vcr_b_70', 'm11bcr_b_70','m11dcr_b_70','m11vcr_b_70','m11v1cr_b_70',\
#'m12icr_b_70','m11fcr_b_70','m11hcr_b_70','m11gcr_b_70',\
#'m10vcr_700', 'm11bcr_700','m11dcr_700','m11fcr_700','m11gcr_700','m11hcr_700','m12icr_700',\
#'fm10q','f553','f573','f476','f383','fm11q','fm11v','fm11v1','fm11v2','f1146','f46','f61'
#,'f476','f383'
#]
dirnr = ['fm12c', 'fm12f', 'fm12i', 'fm12m', 'fm12q']
#dirnr = ['m10vcr_b_70', 'm11bcr_b_70','m11dcr_b_70','m12icr_b_70','m11fcr_b_70','m11hcr_b_70','m11gcr_b_70',\
#'m10qmhdcv','m12imhdcv','m10vmhdcv', 'm11bmhdcv','m11dmhdcv','m12imhdcv','m11fmhdcv','m11hmhdcv','m11gmhdcv',\
#'m10vcr_700', 'm11bcr_700','m11dcr_700','m11fcr_700','m11gcr_700','m11hcr_700','m12icr_700',\
#'fm11v','f383','f553','fm11','fm11q','fm11v','f1146','f61','fm12b']
#dirnr = ['m10vcr_b_70', 'm11bcr_b_70','m11dcr_b_70','m12icr_b_70','m11fcr_b_70','m11hcr_b_70','m11gcr_b_70','m10qmhdcv','m12imhdcv','m10vmhdcv', 'm11bmhdcv','m11dmhdcv','m12imhdcv','m11fmhdcv','m11hmhdcv','m11gmhdcv','fm11v','fm11d','f383','f553','fm11','m10vcr_700', 'm11bcr_700','m11dcr_700','m11fcr_700','m11gcr_700','m11hcr_700','m12icr_700','fm11q','fm11v','f1146','f61']
#dirnr = ['fm11v','fm11d','f383','f553','fm11','m10vmhdcv', 'm11bmhdcv','m11dmhdcv','m12imhdcv','m11fmhdcv','m11hmhdcv','m11gmhdcv','fm11q','fm11v','f1146','f61']
#dirnr = ['fm11v','fm11d','f383','f553','fm11','m10vcr_b_70', 'm11bcr_b_70','m11dcr_b_70','m12icr_b_70','m11fcr_b_70','m11hcr_b_70','m11gcr_b_70','fm11q','fm11v','f1146','f61']
#dirnr = ['fm11','m10vcr_700', 'm11bcr_700','m11dcr_700','m12icr_700','fm11q','fm11v','f1146','f61']
#dirnr = ['fm11q', 'fm11qCRK100','fm11qCRK1000']
#dirhr = ['f1146_hv','f573_hv','f383_hv']
dirneed = np.concatenate((dirnr,dirhr),axis=0)
#dirneed=['f1146','f1146_hv','f573_hv','f573','f553','f476','fm11','f383','f383_hv']
#dirneed=['f1146','f573','f553','f476','fm11','f383']
#tlist=[[13.4,8.3,13.5],[10.4,5.8,13.5],[10.7,2.3,13.5],[2.7,2.0,6.0],[2.0,2.0,6.5]]
#plotname='Oh15_fire2'
#plotname='Msalpha1251_fire1'
#todo='alpha37'
todo='Malpha1251'
#xaxis='mv'
#xaxis='ms'
xaxis='ms_mv'
outputlist=1
alpha12only=1
needanno=1

#plotname=todo+'_'+xaxis+'_'+'fire2sn400'
#plotname=todo+'_'+xaxis+'_'+'test'
#plotname=todo+'_'+xaxis+'_'+'fire2'
#plotname=todo+'_'+xaxis+'_'+'nocr70'
#plotname=todo+'_'+xaxis+'_'+'nomhd'
#plotname=todo+'_'+xaxis+'_'+'nomd'
plotname=todo+'_'+xaxis+'_'+'nocr700'
#plotname=todo+'_'+xaxis+'_'+'nocr700sn300'
#plotname=todo+'_'+xaxis+'_'+'nocr700noname'
#plotname=todo+'_'+xaxis+'_'+'mhdmd'
#plotname=todo+'_'+xaxis+'_'+'mhdsn400'
#plotname=todo+'_'+xaxis+'_'+'crtot'

if todo=='alpha37':
	alpha37list=[]
	mvirlist=[]
	mslist=[]
        alpha37list_hv=[]
        mvirlist_hv=[]
        mslist_hv=[]
        alpha37list_hv2=[]
        mvirlist_hv2=[]
        mslist_hv2=[]
        alpha37list_hv3=[]
        mvirlist_hv3=[]
        mslist_hv3=[]
        alpha37list_hv4=[]
        mvirlist_hv4=[]
        mslist_hv4=[]
	for ncount in range(len(dirneed)):
		try:
			runtodo = dirneed[ncount]
			haloinfo=cosmichalo(runtodo)
			highres=haloinfo['highres']
			if highres==1:
				for icount in range(len(tlist)):
					time=tlist[icount]
					print 'time', time
					Rvir37alpha, xdata, ydata, mvir, ms=dmalpha(runtodo,time,Rrange=[0.3,0.7],timesRv=0,usesnapno=usesnapno)
					alpha37list_hv.append(Rvir37alpha)
					mvirlist_hv.append(mvir)
					mslist_hv.append(ms)
			elif highres==2:
				for icount in range(len(tlist)):
					time=tlist[icount]
					print 'time', time
					Rvir37alpha, xdata, ydata, mvir, ms=dmalpha(runtodo,time,Rrange=[0.3,0.7],timesRv=0,usesnapno=usesnapno)
					alpha37list_hv2.append(Rvir37alpha)
					mvirlist_hv2.append(mvir)
					mslist_hv2.append(ms)
			elif highres==3:
				for icount in range(len(tlist)):
					time=tlist[icount]
					print 'time', time
					Rvir37alpha, xdata, ydata, mvir, ms=dmalpha(runtodo,time,Rrange=[0.3,0.7],timesRv=0,usesnapno=usesnapno)
                                alpha37list_hv3.append(Rvir37alpha)
                                mvirlist_hv3.append(mvir)
                                mslist_hv3.append(ms)
			else:
				for icount in range(len(tlist)):
					time=tlist[icount]
					print 'time', time
					Rvir37alpha, xdata, ydata, mvir, ms=dmalpha(runtodo,time,Rrange=[0.3,0.7],timesRv=0,usesnapno=usesnapno)
					alpha37list.append(Rvir37alpha)
					mvirlist.append(mvir)
					mslist.append(ms)
		except (KeyError,ValueError,IOError):
			continue
	print 'dirneed', dirneed	
	print 'Mvir', mvirlist
	print 'Mstar', mslist
	print 'log10(Mstar)', np.log10(mslist)
	print 'alpha37list', alpha37list
	print 'mslist_hv2', mslist_hv2
	print 'mslist_hv3', mslist_hv3
	plotvsTHINGS()
	if len(alpha37list)>0:
		xneed=np.log10(mslist)
		yneed=alpha37list
		plotalpha37(xneed,yneed)
	if len(alpha37list_hv)>0:
		xneed=np.log10(mslist_hv)
		yneed=alpha37list_hv
		plotalpha37(xneed,yneed,marker='s',needlabel=0)
        finname = 'figures/'+plotname+'.pdf'
        plt.savefig(finname)
        plt.clf()


if todo=='Malpha1251':
	alpha51list=[]
        alpha12list=[]
        mvirlist=[]
        mslist=[]
	labl = []
        alpha51list_hv=[]
        alpha12list_hv=[]
        mvirlist_hv=[]
        mslist_hv=[]
	labl_hv = []
        alpha51list_hv2=[]
        alpha12list_hv2=[]
        mvirlist_hv2=[]
        mslist_hv2=[]
	labl_hv2 = []
        alpha51list_hv3=[]
        alpha12list_hv3=[]
        mvirlist_hv3=[]
        mslist_hv3=[]
	labl_hv3 = []
        alpha51list_hv4=[]
        alpha12list_hv4=[]
        mvirlist_hv4=[]
        mslist_hv4=[]
	labl_hv4 = []
        for ncount in range(len(dirneed)):
		try:
			runtodo = dirneed[ncount]
			haloinfo=cosmichalo(runtodo)
			highres=haloinfo['highres']
			labelname=haloinfo['labelname']
			if highres==1:
				for icount in range(len(tlist)):
					time=tlist[icount]
					Rvir12alpha, xdata, ydata, mvir, ms=dmalpha(runtodo,time,Rrange=[0.01,0.02],timesRv=1,usesnapno=usesnapno)
					Rvir51alpha, xdata, ydata, mvir, ms=dmalpha(runtodo,time,Rrange=[0.005,0.01],timesRv=1,usesnapno=usesnapno)
					alpha12list_hv.append(Rvir12alpha)
					alpha51list_hv.append(Rvir51alpha)
					mvirlist_hv.append(mvir)
					mslist_hv.append(ms)
					labl_hv.append(labelname)
			elif highres==2:
				for icount in range(len(tlist)):
					time=tlist[icount]
					Rvir12alpha, xdata, ydata, mvir, ms=dmalpha(runtodo,time,Rrange=[0.01,0.02],timesRv=1,usesnapno=usesnapno)
					Rvir51alpha, xdata, ydata, mvir, ms=dmalpha(runtodo,time,Rrange=[0.005,0.01],timesRv=1,usesnapno=usesnapno)
					alpha12list_hv2.append(Rvir12alpha)
					alpha51list_hv2.append(Rvir51alpha)
					mvirlist_hv2.append(mvir)
					mslist_hv2.append(ms)
					labl_hv2.append(labelname)
			elif highres==3:
				for icount in range(len(tlist)):
					time=tlist[icount]
					Rvir12alpha, xdata, ydata, mvir, ms=dmalpha(runtodo,time,Rrange=[0.01,0.02],timesRv=1,usesnapno=usesnapno)
					Rvir51alpha, xdata, ydata, mvir, ms=dmalpha(runtodo,time,Rrange=[0.005,0.01],timesRv=1,usesnapno=usesnapno)
					alpha12list_hv3.append(Rvir12alpha)
					alpha51list_hv3.append(Rvir51alpha)
					mvirlist_hv3.append(mvir)
					mslist_hv3.append(ms)
					labl_hv3.append(labelname)
			elif highres==4:
				for icount in range(len(tlist)):
					time=tlist[icount]
					Rvir12alpha, xdata, ydata, mvir, ms=dmalpha(runtodo,time,Rrange=[0.01,0.02],timesRv=1,usesnapno=usesnapno)
					Rvir51alpha, xdata, ydata, mvir, ms=dmalpha(runtodo,time,Rrange=[0.005,0.01],timesRv=1,usesnapno=usesnapno)
					alpha12list_hv4.append(Rvir12alpha)
					alpha51list_hv4.append(Rvir51alpha)
					mvirlist_hv4.append(mvir)
					mslist_hv4.append(ms)
					labl_hv4.append(labelname)
			else:
				for icount in range(len(tlist)):
					time=tlist[icount]
					Rvir12alpha, xdata, ydata, mvir, ms=dmalpha(runtodo,time,Rrange=[0.01,0.02],timesRv=1,usesnapno=usesnapno)
					Rvir51alpha, xdata, ydata, mvir, ms=dmalpha(runtodo,time,Rrange=[0.005,0.01],timesRv=1,usesnapno=usesnapno)
					alpha12list.append(Rvir12alpha)
					alpha51list.append(Rvir51alpha)
					mvirlist.append(mvir)
					mslist.append(ms)
					labl.append(labelname)
		except (KeyError,ValueError,IOError):
			continue
	mslist=np.array(mslist)
	mvirlist=np.array(mvirlist)
        mslist_hv=np.array(mslist_hv)
        mvirlist_hv=np.array(mvirlist_hv)
        mslist_hv2=np.array(mslist_hv2)
        mvirlist_hv2=np.array(mvirlist_hv2)
        mslist_hv3=np.array(mslist_hv3)
        mvirlist_hv3=np.array(mvirlist_hv3)
        mslist_hv4=np.array(mslist_hv4)
        mvirlist_hv4=np.array(mvirlist_hv4)
        print 'dirneed', dirneed
        print 'Mvir', mvirlist
        print 'Mstar', mslist
        print 'log10(Mstar)', np.log10(mslist)
        print 'alpha12list', alpha12list
	print 'alpha51list', alpha51list
        print 'mslist_hv2', mslist_hv2
        print 'mslist_hv3', mslist_hv3
	print 'alpha12list_hv3', alpha12list_hv3
	if xaxis=='mv':
		if len(alpha12list)>0:
			xneed=mvirlist
			yneed=alpha12list
			y1need=alpha51list
			plotMvMsRvalpha(xneed,yneed,y1need,zr=0.,xaxis=xaxis,labl=labl,label='Hydro',needlabel=1,alpha12only=alpha12only,color='y')
		if len(alpha12list_hv)>0:
			xneed=mvirlist_hv
			yneed=alpha12list_hv
			y1need=alpha51list_hv
			plotMvMsRvalpha(xneed,yneed,y1need,zr=0.,xaxis=xaxis,labl=labl_hv,marker='s',needlabel=1,label='MHD',alpha12only=alpha12only,color='b')
                if len(alpha12list_hv2)>0:
                        xneed=mvirlist_hv2
                        yneed=alpha12list_hv2
                        y1need=alpha51list_hv2
                        plotMvMsRvalpha(xneed,yneed,y1need,zr=0.,xaxis=xaxis,labl=labl_hv2,marker='^',needlabel=1,label=r'$\kappa$=3e28',alpha12only=alpha12only,color='g')
                if len(alpha12list_hv3)>0:
                        xneed=mvirlist_hv3
                        yneed=alpha12list_hv3
                        y1need=alpha51list_hv3
                        plotMvMsRvalpha(xneed,yneed,y1need,zr=0.,xaxis=xaxis,labl=labl_hv3,marker='d',needlabel=1,label=r'$\kappa$=3e29',alpha12only=alpha12only,color='k')
	elif xaxis=='ms':
                if len(alpha12list)>0:
			xneed=mslist
			yneed=alpha12list
			y1need = alpha51list
			plotMvMsRvalpha(xneed,yneed,y1need,zr=0.,xaxis=xaxis,labl=labl,label='Hydro',needlabel=1,alpha12only=alpha12only,color='y')
                if len(alpha12list)>0:
                        xneed=mslist_hv
                        yneed=alpha12list_hv
                        y1need = alpha51list_hv
			plotMvMsRvalpha(xneed,yneed,y1need,zr=0.,xaxis=xaxis,marker='s',needlabel=1,labl=labl,label='MHD',alpha12only=alpha12only,color='b')
                if len(alpha12list_hv2)>0:
                        xneed=mslist_hv2
                        yneed=alpha12list_hv2
                        y1need=alpha51list_hv2
                        plotMvMsRvalpha(xneed,yneed,y1need,zr=0.,xaxis=xaxis,marker='^',needlabel=1,labl=labl,label=r'$\kappa$=3e28',alpha12only=alpha12only,color='g')
                if len(alpha12list_hv3)>0:
                        xneed=mslist_hv3
                        yneed=alpha12list_hv3
                        y1need=alpha51list_hv3
                        plotMvMsRvalpha(xneed,yneed,y1need,zr=0.,xaxis=xaxis,marker='d',needlabel=1,labl=labl,label=r'$\kappa$=3e29',alpha12only=alpha12only,color='k')
	elif xaxis=='ms_mv':
                if len(alpha12list)>0:
			xneed=mslist/mvirlist
			yneed=alpha12list
			y1need=alpha51list
			plotMvMsRvalpha(xneed,yneed,y1need,zr=0.,xaxis=xaxis,label='Hydro',needlabel=1,alpha12only=alpha12only,color='y')

                if len(alpha12list_hv)>0:
                        xneed=mslist_hv/mvirlist_hv
                        yneed=alpha12list_hv
                        y1need=alpha51list_hv
			plotMvMsRvalpha(xneed,yneed,y1need,zr=0.,xaxis=xaxis,marker='s',needlabel=1,label='MHD',alpha12only=alpha12only,color='b')

                if len(alpha12list_hv2)>0:
                        xneed=mslist_hv2/mvirlist_hv2
                        yneed=alpha12list_hv2
                        y1need=alpha51list_hv2
                        plotMvMsRvalpha(xneed,yneed,y1need,zr=0.,xaxis=xaxis,marker='^',needlabel=1,label=r'$\kappa$=3e28',alpha12only=alpha12only,color='g')
                if len(alpha12list_hv3)>0:
                        xneed=mslist_hv3/mvirlist_hv3
                        yneed=alpha12list_hv3
                        y1need=alpha51list_hv3
			print 'xneed,yneed,y1need', xneed,yneed,y1need
                        plotMvMsRvalpha(xneed,yneed,y1need,zr=0.,xaxis=xaxis,marker='d',needlabel=1,\
			label=r'$\kappa$=3e29',alpha12only=alpha12only,color='c')
                if len(alpha12list_hv4)>0:
                        xneed=mslist_hv4/mvirlist_hv4
                        yneed=alpha12list_hv4
                        y1need=alpha51list_hv4
                        print 'xneed,yneed,y1need', xneed,yneed,y1need
                        plotMvMsRvalpha(xneed,yneed,y1need,zr=0.,xaxis=xaxis,marker='*',needlabel=1,label='metal diffusion',alpha12only=alpha12only,color='brown')
		if needanno==1:
			text=np.concatenate((labl,labl_hv,labl_hv2,labl_hv3,labl_hv4), axis=None)
			x_data = np.concatenate((np.log10(mslist/mvirlist),np.log10(mslist_hv/mvirlist_hv),\
				np.log10(mslist_hv2/mvirlist_hv2),np.log10(mslist_hv3/mvirlist_hv3),\
				np.log10(mslist_hv4/mvirlist_hv4)), axis=None)
			y_data = np.concatenate((alpha12list,alpha12list_hv,alpha12list_hv2,alpha12list_hv3,
				alpha12list_hv4), axis=None)
			print 'text', text
			print 'x_data', x_data
			print 'y_data', y_data
			txt_width = 0.2; txt_height = 0.08;
			text_positions = get_text_positions(text, x_data, y_data, txt_width, txt_height)
			text_plotter(text, x_data, y_data, text_positions, txt_width, txt_height)
	plt.legend(fontsize=8)
        finname = plotloc+'figures/'+plotname+'.pdf'
        plt.savefig(finname)
        plt.clf()


if outputlist==1:
	f = open('data/dmalpha_info.txt', 'w')
	f.write('halo name    '+ 'Mvir    '+   'Mstar    '+ 'Mstar/Mvir   '+ 'alpha (0.01-0.02Rvir) '+ 'alpha (0.005-0.01) '+'\n')
	for i in range(len(mslist)):
		f.write(str(dirnr[i])+'      '+str(mvirlist[i])+'      '+str(mslist[i])+'      '+str(mslist[i]/mvirlist[i])+'      '+str(alpha12list[i])+'      '+str(alpha51list[i])+'  \n')
	for i in range(len(mslist_hv)):
                f.write(str(dirhr[i])+'      '+str(mvirlist_hv[i])+'      '+str(mslist_hv[i])+'      '+str(mslist_hv[i]/mvirlist_hv[i])+'      '+str(alpha12list_hv[i])+'      '+str(alpha51list_hv[i])+ '  \n')
