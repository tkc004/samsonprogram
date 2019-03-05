from Sasha_functions import *
import numpy as np
r =1.; x =0.5;
#print np.arccos(x/r)

#x_all=[1.,-1.];y_all=[1.,-1.];z_all=[1.,-1.]
#Lx=1.;Ly=1.;Lz=1#.;
#Lx=-2.44371633;Ly=-8.65157373;Lz=7.63155222;
Lx=1e-6;Ly=1e-6;Lz=1.0;
#Lx=-1.81442897;Ly=1.26445131;Lz=-0.698398290
x,y,z= rotateL_to_z(1.0,0.0,0.0,Lx,Ly,Lz)
print 'x,y,z', x,y,z
