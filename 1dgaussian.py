                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   from __future__ import division
import numpy as np
from math import *

################################################################################
##final project for general physics in NTU by 顏立峯 and 吳達懿
##dim : 1d
##type : gaussian
##snap : dynamic
##command : do for [i=0:19]{load '1dgaussian.'.i.'.plt'; pause 0.2 }
################################################################################

##parameter1
c = 299792458
u0 = pi*4E-7
e0 = 1/(u0*c**2)

##

dt = 0.18e-9
dx = c*dt

size = 401

##
hy = np.zeros(size)
ez = np.zeros(size)

##
f = 278e+6

for  n in range(200):

    hy[size-1] = hy[size-2]

    hy[:size-1] += dt/(u0*dx)*(ez[1:]-ez[:size-1])

    ez[0] = ez[1]

    ez[1:] += dt/(e0*dx)*(hy[1:]-hy[:size-1])

    ez[200] = exp(-(n-8)**2/16)

    for n%10 == 0:

        file = open ('1dgaussian.%d.txt'%int(n/10),'w')

        for i in range(size):

            file.write("%d %f\n"%(i,ez[i]))

        file.close()

print "end"
