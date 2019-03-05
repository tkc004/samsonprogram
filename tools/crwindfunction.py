

crwindfunction(wanted,i,runtodo,title='',)
if wanted=='vzwindm' or wanted=='Twindm':
    rcParams['figure.figsize'] = 5,4
    vcut=0.
    zup=30.
    zdown=20.
    withinr=20.
        info=outdirname(runtodo, i)
        rundir=info['rundir']
        runtitle=info['runtitle']
        slabel=info['slabel']
        snlabel=info['snlabel']
        dclabel=info['dclabel']
        resolabel=info['resolabel']
        the_snapdir=info['the_snapdir']
        Nsnapstring=info['Nsnapstring']
        havecr=info['havecr']
        Fcal=info['Fcal']
        iavesfr=info['iavesfr']
        timestep=info['timestep'] 
        cosmo=info['cosmo']
        haveB=info['haveB']
        color=info['color']
        ptitle=title
        if runtitle=='SMC':
            ptitle='Dwarf'
        elif runtitle=='SBC':
            ptitle='Starburst'
        elif runtitle=='MW':
            ptitle=r'$L\star$ Galaxy'
        lsn='solid'
        if haveB>0:
            lsn='dashed'
                    print 'the_snapdir', the_snapdir
                    print 'Nsnapstring', Nsnapstring
                    print 'havecr', havecr
                    if cosmo==1:
                            G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr,h0=1, cosmological=1)
                    else:
                            G = readsnapcr(the_snapdir, Nsnapstring, 0, snapshot_name=the_prefix, extension=the_suffix, havecr=havecr)
                    Gm = G['m']
        Gp = G['p']
        Gv = G['v']
        Tb = G['u']
        rho = G['rho']
        Neb = G['ne']
                    Gm_in_sun = Gm*1e10
        Gx = Gp[:,0]
        Gy = Gp[:,1]
        Gz = Gp[:,2]
        Gvz = Gv[:,2]
        cutv = Gvz*Gz/np.absolute(Gz)>vcut
        cutz = (np.absolute(Gz)>zdown) & (np.absolute(Gz)<zup)
        cutr = np.sqrt(Gx*Gx+Gy*Gy)<withinr
        cut = cutv*cutz*cutr
        if wanted=='Twindm':
            TrueTemp, converted_rho  = SF.convertTemp(Tb, Neb, rho)
            Logvar = np.log10(TrueTemp[cut])
                            xmin=3;xmax=7.0;
            if runtitle=='SMC':
                xmax=6.0;
        elif wanted=='vzwindm':
            Logvar = np.log10(Gvz[cut])
            xmin=0.5;xmax=3.0;
            if runtitle=='SMC':
                xmax=2.3;
        Gm_in_sun=Gm_in_sun[cut]
                    Logxaxis = np.linspace(xmin,xmax,num=nogrid)
                    Gml = []
                    for ib in range(nogrid-1):
                            cutg = (Logvar > Logxaxis[ib]) & (Logvar < Logxaxis[ib+1])
                            Gmcut = Gm_in_sun[cutg]
                            Gml = np.append(Gml, np.sum(Gmcut))
                    plt.plot((Logxaxis[1:]+Logxaxis[:-1])/2.,Gml, label=dclabel,color=color,ls=lsn,lw=2)
if runtitle=='SMC':
    plt.legend(loc='best',fontsize=9,ncol=3)
    plt.yscale('log')
if wanted=='Twindm':
            plt.xlabel(r'$\log (T[{\rm K}])$',fontsize=18)
            plt.ylabel(r'${\rm d}M_{\rm wind}/{\rm d}\log T\;[M_\odot]$',fontsize=18)
    figname='CRplot/Twindm/Twindm_'+runtodo+'_sn'+str(Nsnap)+'.pdf'
elif wanted=='vzwindm':
    plt.xlabel(r'$\log (v_z[{\rm km/s}])$',fontsize=18)
    plt.ylabel(r'${\rm d}M_{\rm wind}/{\rm d}\log v_z\;[M_\odot]$',fontsize=18)
    figname='CRplot/vzwindm/vzwindm_'+runtodo+'_sn'+str(Nsnap)+'.pdf'
plt.title(ptitle,fontsize=18)
plt.tight_layout()
print 'figname', figname
    plt.savefig(figname)
    plt.clf()