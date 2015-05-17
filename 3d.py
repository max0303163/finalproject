from __future__ import division
import numpy as np
from math import *

################################################################################
##final project for general physics in NTU by 顏立峯 and 吳達懿
##dim : 3d
##type : sin
##snap : dynamic / all information
##command : do for [i=0:24]{load '3d.'.i.'.plt' ; pause 0.5}
################################################################################

##parameter1
c = 299792458
u0 = pi*4E-7
e0 = 1/(u0*c**2)

##

dx = 10.8e-2
dt = dx/c/sqrt(3)

size = 101

imp1 = dt/(u0*dx)

imp2 = dt/(e0*dx)

##
ex = np.zeros((size-1,size,size))
ey = np.zeros((size,size-1,size))
ez = np.zeros((size,size,size-1))

hx = np.zeros((size,size-1,size-1))
hy = np.zeros((size-1,size,size-1))
hz = np.zeros((size-1,size-1,size))

##
f = 278e+6



for  n in range(250):

    hx[:,:,:] += imp1*((ey[:,:,1:]-ey[:,:,:size-1])-(ez[:,1:,:]-ez[:,:size-1,:]))

    hy[:,:,:] += imp1*((ez[1:,:,:]-ez[:size-1,:,:])-(ex[:,:,1:]-ex[:,:,:size-1]))

    hz[:,:,:] += imp1*((ex[:,1:,:]-ex[:,:size-1,:])-(ey[1:,:,:]-ey[:size-1,:,:]))

    ex[:,1:size-1,1:size-1] += imp2*((hz[:,1:,1:size-1]-hz[:,:size-2,1:size-1])-(hy[:,1:size-1,1:]-hy[:,1:size-1,:size-2]))

    ey[1:size-1,:,1:size-1] += imp2*((hx[1:size-1,:,1:]-hx[1:size-1,:,:size-2])-(hz[1:,:,1:size-1]-hz[:size-2,:,1:size-1]))

    ez[1:size-1,1:size-1,:] += imp2*((hy[1:,1:size-1,:]-hy[:size-2,1:size-1,:])-(hx[1:size-1,1:,:]-hx[1:size-1,:size-2,:]))

    ez [50,50,50] = sin(2*pi*f*n*dt)

    if n%10 == 0:

        file = open('3d.%d.txt'%int(n/10),'w')

        for i in range(5):

            for j in range(5):

                for k in range(5):

                    file.write("%d %d %d %f %f %f %f %f %f\n "%(20*i,20*j,20*k,ex[20*i,20*j,20*k],ey[20*i,20*j,20*k],ez[10*i,10*j,10*k],
                                                                hx[20*i,20*j,20*k],hy[20*i,20*j,20*k],hz[20*i,20*j,20*k]))

        file.close()

print "end"
