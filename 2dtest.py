from __future__ import division
import numpy as np
from math import*

##two dimensional FDTD test
##TMz 

imp0 = 377.

sizex = 101
sizey = 81
maxtime = 300
cdtds = 1 / sqrt(2.0)

##
ceze = np.ones((sizey,sizex))
cezh = np.ones((sizey,sizex))
cezh[:]=cdtds / imp0

chxh = np.ones((sizey-1,sizex))
chxe = np.ones((sizey-1,sizex))
chxe[:] = cdtds / imp0

chyh = np.ones((sizey,sizex-1))
chye = np.ones((sizey,sizex-1))

chye[:] = cdtds / imp0

##
ez = np.zeros((sizey,sizex))
hx = np.zeros((sizey-1,sizex))
hy = np.zeros((sizey,sizex-1))

hx[:,:]= chxh[:,:]*hx[:,:]-chxe[:,:]*(ez[1:sizey,:]-ez[0:sizey-1,:])
hy[:,:]= chyh[:,:]*hy[:,:]+chye[:,:]*(ez[:,1:sizex]-ez[:,0:sizex-1])

ez[1:sizey-1,1:sizex-1]=ceze[1:sizey,1:sizex-1]*ez[1:sizey-1,1:sizex-1]+cezh[1:sizey-1,1:sizex-1]*((hy[1:sizey-1,1:sizex-1]-hy[1:sizey-1,0:sizex-2])-(hx[1:sizey-1,1:sizex-1]-hy[0:sizey-2,1:sizex-1]))  
