from __future__ import division
import numpy as np
from math import*

##one dimensional FDTD test
##harmonic wave

size = 200
N=100

ez=np.zeros(size)
hy=np.zeros(size)
imp0=377.0

#change this one
maxtime=250
#maxtime=1000

file=open('cos.txt','w')
for i in range(maxtime):
        hy[0:size-2]=hy[0:size-2]+(ez[1:size-1]-ez[0:size-2])/imp0

        ez[1:size-1]=ez[1:size-1]+(hy[1:size-1]-hy[0:size-2])*imp0

        ez[0]=cos(2*pi/N*i)

        file.write("%d %f %f\n"%(i,ez[50],hy[50]))
file.close()

print 'end'
