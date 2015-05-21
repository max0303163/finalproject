from __future__ import division
import numpy as np
from math import *

################################################################################
##final project for general physics in NTU by 顏立峯 and 吳達懿
##dim : 2d
##type : sin / gaussian
##snap : dynamic,from side
##others : none
##command : do for [i=0:45]{load '2d.'.i.'.plt'; pause 0.2 }
################################################################################
##parameter1
c = 299792458
u0 = pi*4E-7
e0 = 1/(u0*c**2)

##

dx = 1e-3
dt = dx/c/sqrt(2)

size = 5001

##
hx = np.zeros((size,size-1))
hy = np.zeros((size-1,size))
ez = np.zeros((size,size))

##
f = 278e+6

for  n in range(1000):

    hx[:,:] -= dt/(u0*dx)*(ez[:,1:]-ez[:,:size-1])

    hy[:,:] += dt/(u0*dx)*(ez[1:,:]-ez[:size-1,:])

    ez[1:size-1,1:size-1] += dt/(e0*dx)*(
        hy[1:,1:size-1]-hy[:size-2,1:size-1])-(
            hx[1:size-1,1:]-hx[1:size-1,:size-2]))

    #ez[100,100] = exp(-(n-8)**2/16)

    ez [100,100] = sin(2*pi*f*n*dt)

    if n %10 ==0:
        file = open ('2d.%d.txt'%int(n/10),'w')

        for i in range(size):

            for j in range(size):

                file.write("%f "%ez[i,j])

            file.write("\n")

        file.close()

    if n == 150:
        file = open ('2dside.txt','w')
        for i in range(size):
            file.write("%f %f\n"%((i-100)*dx,ez[100,i]))
        file.close()

print "end"
