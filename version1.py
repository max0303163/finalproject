from __future__ import division
import numpy as np
from math import *

##parameter1
c = 299792458  
u0 = pi*4E-7
e0 = 1/(u0*c**2)

##variable
sc= 1 / sqrt(2)
dist = 40E-10
size = 101
dx = dist/size

nodes = 10    ##at least 10 
lda = 490E-9  ##wave length

if lda/10 < dx :
    dx = lad/10
    size = int(dist/dx)

dt = sc*dx/c

maxtime=250


##matrix
ex = np.zeros((size,size-1))
ey = np.zeros((size-1,size))
hz = np.zeros((size-1,size-1))


file = open ('2dtest.txt','w')
##

ex[0,0]=1
for t in range(maxtime):

    hz[:,:] -= dt/(u0*dx)*((ex[1:,:]-ex[:size-1,:])
                           -(ey[:,1:]-ey[:,:size-1]))
                            
    ex[1:size-1,:] += dt/(e0*dx)*((hz[1:,:]-hz[0:size-2,:]))

    ey[:,1:size-1] -= dt/(e0*dx)*((hz[:,1:]-hz[:,0:size-2]))

    
    if t == 100 :
        for j in range(size):
                    file.write("%d %f\n"%(j,ey[40,j]))

##    if t %10 == 0:
##        file = open('2dtest.%d.txt'%int(t/10),'w')
file.close()

print "end"
                                

    















