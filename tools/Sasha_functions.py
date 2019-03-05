import numpy as np
from tasz import tfora
import sys
from samson_const import *
from pathloc import *
#from samson_functions import mvee, ellipse

def rotateL_to_z(x_all,y_all,z_all,Lx,Ly,Lz,outangle=0):
        L2 = Lx*Lx+Ly*Ly+Lz*Lz
        theta = np.arccos(Lz/np.sqrt(L2))
        #if Lx<0.: theta*=-1.
        Lxy2 = Lx*Lx+Ly*Ly
        phi = np.arccos(Lx/np.sqrt(Lxy2))
        if Ly<0.: phi*=-1.
        if outangle==0:
                x, y, z = coordinates_rotate_Sasha(x_all, y_all, z_all, -theta, -phi)
                return x,y,z
        if outangle==1:
                return theta,phi

def coordinates_rotate_Sasha(x_all, y_all, z_all, theta, phi):
    #rotate the cordinate around z axis (phi is the angle between the x and the vector) and then around y (theta is the angle between z and the vector)
    x = np.cos(theta)*np.cos(phi)*x_all - np.cos(theta)*np.sin(phi)*y_all + np.sin(theta)*z_all
    y = np.sin(phi)*x_all + np.cos(phi)*y_all + 0.*z_all
    z = -np.sin(theta)*np.cos(phi)*x_all + np.sin(theta)*np.sin(phi)*y_all + np.cos(theta)*z_all

    #y = np.sin(phi)*x_all + np.cos(phi)*y_all + 0.*z_all
    #z =  -np.cos(phi)*np.sin(theta)*x_all + np.sin(phi)*np.sin(theta)*y_all + np.cos(theta)*z_all
    return x,y,z



def convertTemp(Tb, Neb, rho):
        #converts temperature and density to physical units (direct from readsnap.py, so no h factor)
        
        #print 'it is Sasha in tools'
        UnitLength_in_cm=kpc_in_cm
        UnitMass_in_g=Msun_in_g*1e10
        UnitVelocity_in_cm_per_s=km_in_cm
        UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
        UnitEnergy_in_cgs=UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2)

        MeanWeights = 4.0/(3*H_MASSFRAC+1+4*H_MASSFRAC* Neb) * PROTONMASS
        MeanWeights_con = 4.0/(3*H_MASSFRAC+1+4*H_MASSFRAC* Neb) 

        converted_rho = rho * UnitMass_in_g * Avogadro_const / MeanWeights_con
        converted_rho /= (UnitLength_in_cm * UnitLength_in_cm *UnitLength_in_cm)
        
        TrueTemp = Tb *  MeanWeights/BOLTZMANN * (GAMMA-1.) * UnitEnergy_in_cgs/ UnitMass_in_g
        return (TrueTemp, converted_rho)

def cut_stellar_age(a, N, age, omeganot, littleh):
        #returns the time since the previous snapshot, as well as an array that gives the actual age of all stellar particles in Myr (as opposed to formation time, as given by code)
        finname = 'output_times.txt'
        f = open(finname)
        SdarsA = np.loadtxt(f)
        prev_epoch = SdarsA[max((N-1), 0)]
        SFR_time_range = (tfora(a, omeganot, littleh) - tfora(prev_epoch, omeganot, littleh)) * 1e9
        print 'using exact SFR', a, prev_epoch, (SFR_time_range)/1e7, ' * 1e7 yrs'
        f.close()
        #print age
        time_in_Gyr = tfora(age, omeganot, littleh)
        age_in_Gyr = tfora(a, omeganot, littleh) - time_in_Gyr
        print ' i interpolated your ages'
        #print age_in_Gyr
        #print ' the times of formation are: ',time_in_Gyr
        return (SFR_time_range, age_in_Gyr)


def read_halo_catalog(rundir, Nsnapstring, redshiftstring, usembp='n',maindir='scratch',Sheaform=0):
        finname = '/home/tkc004/'+maindir+'/'+rundir+'/Pep/snap'+ Nsnapstring+ 'RPep.z'+redshiftstring+'.AHF_halos' #or whatever prefix
        f = open(finname)
        C = np.loadtxt(f)
        halID = C[:,0]
        x = C[:,5]
        y = C[:,6]
        z = C[:,7]
        vx = C[:,8]
        vy = C[:,9]
        vz = C[:,10]
        M = C[:,3]
        Vmax = C[:,16]
        Vsig = C[:,18]
        lambdano = C[:,19] #Bullock+01 definition
        Rvir = C[:,11]
        Lx = C[:,21]
        Ly = C[:,22]
        Lz = C[:,23]
        Mgas = C[:,53]
        if Sheaform==1:
                Mstar = C[:,64]
        else:
                Mstar = C[:,73]
        fMhires = C[:,37]
        Lxgas = C[:,56]
        Lygas = C[:,57]
        Lzgas = C[:,58]
        Lxstar = C[:,76]
        Lystar = C[:,77]
        Lzstar = C[:,78]
        #if (use_mbp=='y'): not even going to bother writing this for now - it isnt very accurate for halo center in Amiga. May want to use COM instead though
                
        return {'Lxgas':Lxgas, 'Lygas':Lygas, 'Lzgas':Lzgas,'Lxstar':Lxstar, 'Lystar':Lystar, 'Lzstar':Lzstar, 'Lx':Lx,'Ly':Ly,'Lz':Lz,'fMhires':fMhires,'ID':halID, 'x':x, 'y':y, 'z':z,\
 'vx':vx, 'vy':vy, 'vz':vz, 'M':M, 'Vmax':Vmax, 'Vsig':Vsig, 'lambdano':lambdano,\
'Rvir':Rvir, 'Mgas':Mgas, 'Mstar':Mstar}

        
def shell_check(dists, Rvir, R_bins, a, h, use_physical_bins='n', physical_extent=50):

        TheCuts = []
        bincount = 0
        ShellList = np.zeros(len(dists)) - 1 
        while (bincount < R_bins):
                DistMin = Rvir * bincount/R_bins
                DistMax = Rvir * (bincount+1)/R_bins
        
                #by default, the virial radius is cut into 10 bins. Here i'm also allowing for a fixed physical scale.
                if (use_physical_bins == 'y'):
                        print 'using physical bins ',physical_extent, 'kpc/h'
                        #gotta divide by a to go to comoving coordinates
                        Distmin = bincount/(Rbins) * physical_extent / a
                        Distmax = (bincount+1)/(Rbins) * physical_extent / a 

                ShellThickness = (DistMax - DistMin)*a/h

                Cut1 = dists > DistMin
                Cut2 = dists < DistMax
                Cuts = Cut1 * Cut2
                ShellList[Cuts] = bincount
                TheCuts.append(Cuts)
                bincount+=1
        TheCuts = np.array(TheCuts)
        return (TheCuts, ShellList, ShellThickness)

def calcdist2(x1, y1, z1, posx, posy, posz, boxsize):
#calculates distances of an array of points (particles) from a single given point (halo center)
#for your reference: x1, y1, z1 are floats for halo center or hwatever, posx, posy posz are arrays for positions of particles or whatever
         xdist = abs(posx - x1)
         ydist = abs(posy - y1)
         zdist = abs(posz - z1)
         
         #adjust for periodic boundary conditions
         x_straggler = xdist/boxsize > 0.5
         if (len(xdist[x_straggler])) > 0:
                print 'adjusting x stragglers ', len(xdist[x_straggler])
                xdist[x_straggler] = boxsize - xdist[x_straggler]
         y_straggler = ydist/boxsize > 0.5
         if (len(ydist[y_straggler])) > 0:
                ydist[y_straggler] = boxsize - ydist[y_straggler]               
                print 'adjusting y stragglers ', len(ydist[y_straggler])
         z_straggler = zdist/boxsize > 0.5
         if (len(zdist[z_straggler])) > 0:
                zdist[z_straggler] = boxsize - zdist[z_straggler]
                print 'adjusting z stragglers ', len(zdist[z_straggler])

         
         dists = pow(((xdist*xdist + ydist*ydist + zdist*zdist)), 0.5)
         return dists
         
def line_to_string(line, newline = 'y'):
        linestring = ''
        for q in line:
                linestring += "{0:.6g}".format(q) + '  '
        if (newline == 'y'):
                linestring+=' \n'
        return linestring
        
def sasha_max(a):
        if (len(a) > 0):
                return np.max(a)
        else:
                print 'warning there was an empty array max'
                return -99999999

def sasha_min(a):
        if (len(a) > 0):
                return np.min(a)
        else:
                print 'warning there was an empty array min'
                return 99999999

def check_Rvir_growth(halocount, a, Rvir, Vsig, Mass):
        #check to make sure that the snapshot you are using is not affected by a glitch in the halo catalog that causes a temporary dip in Rvir, Mvir, Vsig 
        H = read_halo_history(int(halocount))
        history_Rvirs = H['Rvir']
        history_M = H['M']
        history_Vsig = H['Vsig']
        history_Zs = H['redshift']
        history_As = 1.0 / (1.0 + history_Zs)
        history_Rvirs_phys = history_Rvirs * history_As
        RvirPhys = Rvir * a
        cut = history_As < a
        historical_max = sasha_max(history_Rvirs_phys[cut])  
        if (RvirPhys >= historical_max):
                print 'Rvir is already the biggest'
                return (Rvir, Vsig, Mass)
        maxindex = np.where(history_Rvirs_phys[cut]==historical_max)[0][0]
        print 'adjusting back to z = ',history_Zs[maxindex]
        print 'Rvir was ',history_Rvirs[maxindex]
        print 'instead of meager ',Rvir
        return (history_Rvirs[maxindex], history_Vsig[maxindex], history_M[maxindex])

def read_halo_history(rundir, halonostr='00', multifile='n',suffixadd='', hubble=1,comoving=1, maindir='scratch', singlesnap=0, atime=1, snumadd=0): #to output comoving distance comoving = 1
        redlist = []
        halolist = []
        xlist = []
        ylist = []
        zlist = []
        xvlist = []
        yvlist = []
        zvlist = []
        mvirlist = []
        mstarlist = []
        mgaslist = []
        rvirlist = []
        fMhireslist = []
        lambdalist = []
        Lxgaslist = []
        Lygaslist = []
        Lzgaslist = []
        Lxstarlist = []
        Lystarlist = []
        Lzstarlist = []
        dirname='/home/tkc004/'+maindir+'/'+rundir
        halofile=open(dirname+'/halos/halo_000'+halonostr+suffixadd+'.dat','r')
        halofile.readline()
        halofile.readline()
        dars = halofile.readlines()
        halofile.close()
        mvirno = 4
        xcenno = 6
        ycenno = 7
        zcenno = 8
        xvcenno = 9
        yvcenno = 10
        zvcenno = 11
        zredno = 0
        halono = 1
        mstarno = 74
        rvirno = 12
        fMhiresno = 38
        mgasno = 54
        lambdano=20
        Lxgasno = 57
        Lygasno = 58
        Lzgasno = 59
        Lxstarno = 77
        Lystarno = 78
        Lzstarno = 79
        if snumadd==1:
                mvirno+= 1
                xcenno+= 1
                ycenno+= 1
                zcenno+= 1
                xvcenno+= 1
                yvcenno+= 1
                zvcenno+= 1
                zredno+= 1
                halono+= 1
                mstarno+= 1
                rvirno+= 1
                fMhiresno+= 1
                mgasno+= 1
                lambdano+= 1
                Lxgasno+= 1
                Lygasno+= 1
                Lzgasno+= 1
                Lxstarno+= 1
                Lystarno+= 1
                Lzstarno+= 1
                
        for line in dars:
                xsd = line.split()
                mvir = float(xsd[mvirno])
                xcen = float(xsd[xcenno])
                ycen = float(xsd[ycenno])
                zcen = float(xsd[zcenno])
                xvcen = float(xsd[xvcenno])
                yvcen = float(xsd[yvcenno])
                zvcen = float(xsd[zvcenno])
                Lxgas = float(xsd[Lxgasno])
                Lygas = float(xsd[Lygasno])
                Lzgas = float(xsd[Lzgasno])
                Lxstar = float(xsd[Lxstarno])
                Lystar = float(xsd[Lystarno])
                Lzstar = float(xsd[Lzstarno])
                if comoving==1:
                        zred=0.0
                else:
                        zred=float(xsd[zredno])
                redlist=np.append(redlist,float(xsd[zredno]))
                halolist=np.append(halolist,int(xsd[halono]))
                mvirlist=np.append(mvirlist,mvir/hubble)
                mstarlist=np.append(mstarlist,float(xsd[mstarno])/hubble)
                xlist=np.append(xlist,xcen/hubble/(1.0+zred)) # xcen originally in comoving unit (kpc/h)
                ylist=np.append(ylist,ycen/hubble/(1.0+zred))
                zlist=np.append(zlist,zcen/hubble/(1.0+zred))
                xvlist=np.append(xvlist,xvcen) #originally in pecular velocity
                yvlist=np.append(yvlist,yvcen)
                zvlist=np.append(zvlist,zvcen)
                rvirlist=np.append(rvirlist,float(xsd[rvirno])/hubble/(1.0+zred)) # originally in comoving unit (kpc/h)
                fMhireslist= np.append(fMhireslist, float(xsd[fMhiresno]))
                mgaslist = np.append(mgaslist, float(xsd[mgasno])/hubble)
                lambdalist = np.append(lambdalist, float(xsd[lambdano])) #spin parameter defined in Bullock+2001
                Lxgaslist = np.append(Lxgaslist, float(xsd[Lxgasno]))
                Lygaslist = np.append(Lygaslist, float(xsd[Lygasno]))
                Lzgaslist = np.append(Lzgaslist, float(xsd[Lzgasno]))
                Lxstarlist = np.append(Lxstarlist, float(xsd[Lxstarno]))
                Lystarlist = np.append(Lystarlist, float(xsd[Lystarno]))
                Lzstarlist = np.append(Lzstarlist, float(xsd[Lzstarno]))
        redlist=np.array(redlist)
        halolist=np.array(halolist)
        mvirlist=np.array(mvirlist)
        mstarlist=np.array(mstarlist)
        mgaslist=np.array(mgaslist)
        xlist=np.array(xlist)
        ylist=np.array(ylist)
        zlist=np.array(zlist)
        xvlist=np.array(xvlist)
        yvlist=np.array(yvlist)
        zvlist=np.array(zvlist)
        rvirlist=np.array(rvirlist)
        fMhireslist=np.array(fMhireslist)
        lambdalist=np.array(lambdalist)
        Lxgaslist=np.array(Lxgaslist)
        Lygaslist=np.array(Lygaslist)
        Lzgaslist=np.array(Lzgaslist)
        Lxstarlist=np.array(Lxstarlist)
        Lystarlist=np.array(Lystarlist)
        Lzstarlist=np.array(Lzstarlist)
        if snumadd==1:
                redlist=redlist[::-1]
                halolist=halolist[::-1]
                mvirlist=mvirlist[::-1]
                mstarlist=mstarlist[::-1]
                mgaslist=mgaslist[::-1]
                xlist=xlist[::-1]
                ylist=ylist[::-1]
                zlist=zlist[::-1]
                xvlist=xvlist[::-1]
                yvlist=yvlist[::-1]
                zvlist=zvlist[::-1]
                rvirlist=rvirlist[::-1]
                fMhireslist=fMhireslist[::-1]
                lambdalist=lambdalist[::-1]
                Lxgaslist=Lxgaslist[::-1]
                Lygaslist=Lygaslist[::-1]
                Lzgaslist=Lzgaslist[::-1]
                Lxstarlist=Lxstarlist[::-1]
                Lystarlist=Lystarlist[::-1]
                Lzstarlist=Lzstarlist[::-1]
        if singlesnap==0:
                        return {'lambdalist':lambdalist, 'redshift':redlist, 'ID':halolist, 'x':xlist, 'y':ylist, 'z':zlist,\
                        'xv':xvlist, 'yv':yvlist, 'zv':zvlist, 'M':mvirlist, 'Ms':mstarlist, 'Mg':mgaslist, 'R':rvirlist,\
                        'fMhires':fMhireslist, 'Lxstar':Lxstarlist, 'Lystar':Lystarlist, 'Lzstar':Lzstarlist,\
                        'Lxgas':Lxgaslist, 'Lygas':Lygaslist, 'Lzgas':Lzgaslist}
        else:
                aAHF = 1.0/(1.0+redlist)
                rednow = 1.0/atime-1.0
                print 'atime', atime
                xcens = np.interp(atime,aAHF,xlist)
                ycens = np.interp(atime,aAHF,ylist)
                zcens = np.interp(atime,aAHF,zlist)
                xvcens = np.interp(atime,aAHF,xvlist)
                yvcens = np.interp(atime,aAHF,yvlist)
                zvcens = np.interp(atime,aAHF,zvlist)
                MAHF = np.interp(atime,aAHF,mvirlist)
                SmAHF = np.interp(atime,aAHF,mstarlist)
                GmAHF = np.interp(atime,aAHF,mgaslist)
                RAHF = np.interp(atime,aAHF,rvirlist)
                fAHF = np.interp(atime,aAHF,fMhireslist)
                lambdaspin = np.interp(atime,aAHF,lambdalist)
                ids = np.interp(atime,aAHF,halolist)    
                Lxgas = np.interp(atime,aAHF,Lxgaslist)
                Lygas = np.interp(atime,aAHF,Lygaslist)
                Lzgas = np.interp(atime,aAHF,Lzgaslist)
                Lxstars = np.interp(atime,aAHF,Lxstarlist)
                Lystars = np.interp(atime,aAHF,Lystarlist)
                Lzstars = np.interp(atime,aAHF,Lzstarlist)
                return {'lambdalist':lambdaspin, 'redshift':rednow, 'ID':ids, 'x':xcens, 'y':ycens, 'z':zcens,\
                'xv':xvcens, 'yv':yvcens, 'zv':zvcens, 'M':MAHF, 'Ms':SmAHF, 'Mg':GmAHF, 'R':RAHF, 'fMhires':fAHF,\
                'Lxstar':Lxstars,'Lystar':Lystars,'Lzstar':Lzstars,'Lxgas':Lxgas,'Lygas':Lygas,'Lzgas':Lzgas}


def read_halo_history_pep(rundir, finalno, Sheaform=0, snapsep=1, singlesnap=1, beginno=100,halonostr='00', multifile='n', hubble=1,comoving=1, maindir='scratch',firever=1, outputallhalos=0): #to output comoving distance comoving = 1
        redlist = []
        idlist = []
        xlist = []
        ylist = []
        zlist = []
        xvlist = []
        yvlist = []
        zvlist = []
        Mlist = []
        Mslist = []
        Mglist = []
        Rvlist = []
        fMlist = []
        lambdalist = []
        Lxgaslist=[]
        Lygaslist=[]
        Lzgaslist=[]
        Lxstarlist=[]
        Lystarlist=[]
        Lzstarlist=[]
        halono=int(halonostr)
        strnolist, zstrlist = snapshotzstr(firever=firever)
        if (singlesnap==1): beginno=finalno
        print 'beginno', beginno
        for Nsnap in range(beginno,finalno+1,snapsep):
                redshiftstring = str(zstrlist[Nsnap])
                #Nsnapstring = str(strnolist[Nsnap]) 
                Nsnapstring = str(Nsnap)
                hcat=read_halo_catalog(rundir, Nsnapstring, redshiftstring, usembp='n',maindir=maindir,Sheaform=Sheaform)
                idl = hcat['ID']
                mvirl = hcat['M']
                xl = hcat['x']
                yl = hcat['y']
                zl = hcat['z']
                xvl = hcat['vx']
                yvl = hcat['vy']
                zvl = hcat['vz']
                Msl = hcat['Mstar']
                Mgl = hcat['Mgas']
                Rl = hcat['Rvir']
                fMl = hcat['fMhires']
                lambdal =hcat['lambdano']
                Lxgasl = hcat['Lxgas']
                Lygasl = hcat['Lygas']
                Lzgasl = hcat['Lzgas']
                Lxstarl = hcat['Lxstar']
                Lystarl = hcat['Lystar']
                Lzstarl = hcat['Lzstar']
                redshift = float(redshiftstring)
                if comoving==1:
                        zred=0.0
                else:
                        zred=redshift
                haloid = halono
                mvirl = mvirl/hubble
                Msl = Msl/hubble
                Mgl = Mgl/hubble
                xl = xl/hubble/(1.0+zred)
                yl = yl/hubble/(1.0+zred)
                zl = zl/hubble/(1.0+zred)
                Rl = Rl/hubble/(1.0+zred)
                if outputallhalos==1 and singlesnap==1:
                        print 'output all halos'
                        return {'fMhires':fMl, 'redshift':redshift, 'ID':idl, 'x':xl, 'y':yl, 'z':zl,\
                        'xv':xvl, 'yv':yvl, 'zv':zvl, 'M':mvirl, 'Ms':Msl, 'Mg':Mgl, 'R':Rl,\
                        'Lxgas':Lxgasl, 'Lygas':Lygasl, 'Lzgas':Lzgasl,\
                        'lambdalist':lambdal,'Lxstar':Lxstarl, 'Lystar':Lystarl, 'Lzstar':Lzstarl}
                mvir = float(mvirl[halono])
                mstar = float(Msl[halono])
                mgas = float(Mgl[halono])
                xcen = float(xl[halono]) # xcen originally in comoving unit (kpc/h)
                ycen = float(yl[halono])
                zcen = float(zl[halono])
                xvcen = float(xvl[halono])
                yvcen = float(yvl[halono])
                zvcen = float(zvl[halono])
                rvir = float(Rl[halono])  # originally in comoving unit (kpc/h)
                lambdaB = float(lambdal[halono])
                Lxgas = float(Lxgasl[halono])
                Lygas = float(Lygasl[halono])
                Lzgas = float(Lzgasl[halono])
                Lxstar = float(Lxstarl[halono])
                Lystar = float(Lystarl[halono])
                Lzstar = float(Lzstarl[halono])
                xlist = np.append(xlist,xcen)
                ylist = np.append(ylist,ycen)
                zlist = np.append(zlist,zcen)
                xvlist = np.append(xvlist,xvcen)
                yvlist = np.append(yvlist,yvcen)
                zvlist = np.append(zvlist,zvcen)
                idlist = np.append(idlist,haloid)
                redlist = np.append(redlist,float(zstrlist[Nsnap]))
                Mlist = np.append(Mlist,mvir)
                Mslist = np.append(Mslist,mstar)
                Mglist = np.append(Mglist,mgas)
                Rvlist = np.append(Rvlist,rvir)
                fMlist = np.append(fMlist,fMl)
                lambdalist = np.append(lambdalist,lambdaB)
                Lxgaslist = np.append(Lxgaslist,Lxgas)
                Lygaslist = np.append(Lygaslist,Lygas)
                Lzgaslist = np.append(Lzgaslist,Lzgas)
                Lxstarlist = np.append(Lxstarlist,Lxstar)
                Lystarlist = np.append(Lystarlist,Lystar)
                Lzstarlist = np.append(Lzstarlist,Lzstar)
        xlist = np.array(xlist)
        ylist = np.array(ylist)
        zlist = np.array(zlist)
        xvlist = np.array(xvlist)
        yvlist = np.array(yvlist)
        zvlist = np.array(zvlist)
        idlist = np.array(idlist)
        redlist = np.array(redlist)
        Mlist = np.array(Mlist)
        Mslist = np.array(Mslist)
        Mglist = np.array(Mglist)
        Rvlist = np.array(Rvlist)
        fMlist = np.array(fMlist)
        Lxgaslist = np.array(Lxgaslist)
        Lygaslist = np.array(Lygaslist)
        Lzgaslist = np.array(Lzgaslist)
        Lxstarlist = np.array(Lxstarlist)
        Lystarlist = np.array(Lystarlist)
        Lzstarlist = np.array(Lzstarlist)
        lambdalist = np.array(lambdalist)
        if (singlesnap==1):
                return {'lambdalist':lambdaB, 'fMhires':fMl, 'redshift':redshift, 'ID':haloid,\
                'x':xcen, 'y':ycen, 'z':zcen, 'xv':xvcen, 'yv':yvcen, 'zv':zvcen,\
                'M':mvir, 'Ms':mstar, 'Mg':mgas, 'R':rvir, 'Lxstar':Lxstar, 'Lystar':Lystar,'Lzstar':Lzstar,\
                'Lxgas':Lxgas, 'Lygas':Lygas,'Lzgas':Lzgas}
        else:
                return {'lambdalist':lambdalist, 'fMhires':fMlist, 'redshift':redlist, 'ID':idlist,\
                'x':xlist, 'y':ylist, 'z':zlist, 'xv':xvlist, 'yv':yvlist, 'zv':zvlist,\
                'M':Mlist, 'Ms':Mslist, 'Mg':Mglist, 'R':Rvlist, 'Lxstar':Lxstarlist, 'Lystar':Lystarlist,'Lzstar':Lzstarlist,\
                'Lxgas':Lxgaslist, 'Lygas':Lygaslist,'Lzgas':Lzgaslist}



def snapshotzstr(firever=1):
        if firever==1:
                spmname=programdir+'/data/snapshottime_FIRE1.txt'
                spmfile=open(spmname,"r")
                dars = spmfile.readlines()
                snapno=[]
                redshift=[]
                for line in dars:
                        xsd = line.split()
                        snapno = np.append(snapno, int(xsd[0]))
                        redshift = np.append(redshift, float(xsd[2]))
                strzlist=np.array(redshift)
                strnolist=np.array(snapno)
        elif firever==2:
                spmname=programdir+'/data/snapshot_times.txt'
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
                spmfile.close()
                snapno=np.array(snapno)
                redshift=np.array(redshift)
                strzlist=[]
                strnolist=[]
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
                        strzlist = np.append(strzlist,strz)
                        strnolist = np.append(strnolist,strno)
        return strnolist, strzlist
