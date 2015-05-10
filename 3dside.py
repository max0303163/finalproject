from __future__ import division
import numpy as np
from math import *

##1d test
#sin wave at the middle

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
        
    if n%100 == 0:
        file = open('3dside.txt','w')
        for i in range(size-1):
            file.write("%f %f\n"%((i-50)*dx,ez[i,50,50]))
        file.close()


print "end"

    
