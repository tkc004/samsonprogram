from stdmodandoption import *
from mpl_toolkits.axes_grid1 import AxesGrid
import plot_setup as PS
from crgridfunction import crgridfunction

def rhoTgridtestinput(subdict):
        dirneed=subdict['dirneed']
        wanted=subdict['wanted']
        print 'wanted', wanted
        startno=subdict['startno']
        Nsnap=subdict['Nsnap']
        snapsep=subdict['snapsep']
        the_prefix=subdict['the_prefix']
        the_suffix=subdict['the_suffix']
        fmeat=subdict['fmeat']
        mtitle=subdict['mtitle'] 
        withinr=subdict['withinr']
        withoutr=subdict['withoutr']
        zup=subdict['zup']
        zdown=subdict['zdown']
        Tcut=subdict['Tcut']
        addtitleinbox=subdict['addtitleinbox']
        extendedrange=subdict['extendedrange']
        print 'dirneed', dirneed
        print 'extendedrange', extendedrange
        normalized=subdict['normalized']
        useverplot=subdict['useverplot']
        userad=subdict['userad']
        if useverplot==1:
            ncols=1
            nrows=len(dirneed) 
        else:
            nrows=1
            ncols=len(dirneed)       
        if len(dirneed)>3:
                nrows = 2
                ncols = 3
        print 'addtitleinbox', addtitleinbox
        print 'nrows, ncols', nrows, ncols

        fig, gridag = PS.setupaxesgrid(nrows=nrows, ncols=ncols)

        plotupt=[]
        plotit=[]

        for idir, runtodo in enumerate(dirneed):
            print 'dirneed', dirneed
            print 'idir, runtodo', idir, runtodo
            im, plotupt,plotit,cblabel,totalname\
            = crgridfunction(gridag,runtodo,wanted,idir,plotupt,plotit,startno,Nsnap,snapsep,cbcolor='YlOrRd',
                           rotface=1,
                           Tcut=Tcut, #K
                           highTcut=1e10,
                           vcut=-1e10, #outflow velocity cut
                           vhcut=1e10, #upper outflow velocity cut
                           withinr=withinr, 
                           withoutr=withoutr,
                           zup=zup,
                           zdown=zdown,
                           userad=userad,
                           trackgas=0,
                           normalized=normalized,
                           fmeat=fmeat,
                           extendedrange=extendedrange)
        #plotupt=['',mtitle,'']
        print 'plotupt', plotupt
        PS.finishaxesgrid(im,cblabel,gridag,plotit,plotupt,addtitleinbox=addtitleinbox)

        print 'totalname', totalname
        plt.savefig(totalname,bbox_inches='tight',pad_inches = 0.1)
        plt.clf()
        return None