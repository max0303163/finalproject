from __future__ import division
import numpy as np
import math

##one dimensional FDTD test

size = 200

ez=np.zeros(size)
hy=np.zeros(size)
imp0=377.0

#change this one
#maxtime=250
maxtime=1000

file=open('oned1.txt','w')
for i in range(maxtime):
        hy[0:size-2]=hy[0:size-2]+(ez[1:size-1]-ez[0:size-2])/imp0

        ez[1:size-1]=ez[1:size-1]+(hy[1:size-1]-hy[0:size-2])*imp0

        ez[0]=math.exp(-(i-30)*(i-30)/100)

        file.write("%d %f\n"%(i,ez[50]))
file.close()

print 'end'

