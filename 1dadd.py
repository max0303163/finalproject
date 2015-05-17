from __future__ import division
import numpy as np
import math

################################################################################
##final project for general physics in NTU by 顏立峯 and 吳達懿
##dim : 1d
##type : gaussian
##snap : dynamic
##others : add source
##command : set yrange[-0.3:1] ; set xrange[0:200] ; do for [i=0:99]{plot 'add50.'.i.'.txt' w l notitle ; pause 0.2} 
################################################################################
size = 200

ez=np.zeros(size)
hy=np.zeros(size)
imp0=377.0

#change this one
maxtime=1000

for i in range(maxtime):
        hy[size-1]=hy[size-2]

        hy[0:size-1]=hy[0:size-1]+(ez[1:size]-ez[0:size-1])/imp0

        ez[0]=ez[1]

        ez[1:size]=ez[1:size]+(hy[1:size]-hy[0:size-1])*imp0

        ez[50]+=math.exp(-(i-30)*(i-30)/100)



        if i%10==0:
                file=open('add50.%d.txt'%int(i/10),'w')
                for j in range(size):
                    file.write("%d %f\n"%(j,ez[j]))
                file.close()
print 'end'
