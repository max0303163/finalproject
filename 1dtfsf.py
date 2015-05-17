from __future__ import division
import numpy as np
import math

################################################################################
##final project for general physics in NTU by 顏立峯 and 吳達懿
##dim : 1d
##type : gaussian
##snap : dynamic
##others : tfsf boundary
##command : do for [i=0:45]{load 'tfsf.'.i.'.plt';pause 0.2}
################################################################################

size = 200

ez=np.zeros(size)
hy=np.zeros(size)
imp0=377.0

#change this one
maxtime=1000

for i in range(maxtime):
        hy[0:size-2]=hy[0:size-2]+(ez[1:size-1]-ez[0:size-2])/imp0

        hy[49] -= math.exp(-(i-30)*(i-30)/100)/imp0

        ez[0]=ez[1]

        ez[size-1]=ez[size-2]=0

        ez[1:size-1]=ez[1:size-1]+(hy[1:size-1]-hy[0:size-2])*imp0


        ez[50]+=math.exp(-(i+0.5+0.5-30)*(i+0.5+0.5-30)/100)



        if i%10==0:
                file=open('tfsf50.%d.txt'%int(i/10),'w')
                for j in range(size):
                    file.write("%d %f\n"%(j,ez[j]))
                file.close()
print 'end'
