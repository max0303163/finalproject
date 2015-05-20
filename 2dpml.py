from __future__ import division
import numpy as np
from math import *

################################################################################
##final project for general physics in NTU by 顏立峯 and 吳達懿
##dim : 2d
##type : sin / gaussian
##snap : dynamic,from side
##others : none
##command : do for [i=0:45]{load '2dpml.'.i.'.plt'; pause 0.2 }
################################################################################
##parameter1
c = 299792458
u0 = pi*4E-7
e0 = 1/(u0*c**2)

##

dx = 5.4e-2
dt = dx/c/sqrt(2)

size = 201

##
hx = np.zeros((size,size-1))
hy = np.zeros((size-1,size))
ez = np.zeros((size,size))

##
Ezx = np.zeros((size,size))
Ezy = np.zeros((size,size))
caEzx = np.ones((size,size))
caEzy = np.ones((size,size))
cbEzx = np.ones((size,size)) * dt/(e0*dx)
cbEzy = np.ones((size,size)) * dt/(e0*dx)
daHx = np.ones((size,size-1))
daHy = np.ones((size-1,size))
dbHx = np.ones((size,size-1)) * dt/(u0*dx)
dbHy = np.ones((size-1,size)) * dt/(u0*dx)

##
f = 278e+6

##
npmls = 10
rsize = size - npmls
sigmax = -3*e0*c*log(1e-5)/(2*dx*npmls)
rhomax = sigmax*(u0/e0)
sig = np.zeros((npmls))
rho = np.zeros((npmls))
ca = np.zeros((npmls))
cb = np.zeros((npmls))
da = np.zeros((npmls))
db = np.zeros((npmls))

for m in range(npmls):
    sig[m] = sigmax*((m+0.5)/(npmls+0.5))**2
    rho[m] = rhomax*((m+1)/(npmls+0.5))**2
    re = sig[m]*dt/e0
    rm = rho[m]*dt/u0
    ca[m] = exp(-re)
    cb[m] = -(exp(-re)-1)/sig[m]/dx
    da[m] = exp(-rm)
    db[m] = -(exp(-rm)-1)/rho[m]/dx


##
for i in range(1,size-1):
    for j in range(1,npmls+1):
        m = npmls+2-j-2
        caEzy[i,j] = ca[m]
        cbEzy[i,j] = cb[m]
        daHx[i,j-1] = da[m-1]
        dbHx[i,j-1] = db[m-1]

    for j in range(rsize,size-1):
        m = j+1-rsize-1
        caEzy[i,j] = ca[m]
        cbEzy[i,j] = cb[m]
        daHx[i,j] = da[m]
        dbHx[i,j] = db[m]

for j in range(1,size-1):
    for i in range(1,npmls+1):
        m = npmls+2-i-2
        caEzx[i,j] = ca[m]
        cbEzx[i,j] = cb[m]
        daHy[i-1,j] = da[m-1]
        dbHy[i-1,j] = db[m-1]

    for i in range(rsize,size-1):
        m = i+1-rsize-1
        caEzx[i,j] = ca[m]
        cbEzx[i,j] = cb[m]
        daHy[i,j] = da[m]
        dbHy[i,j] = db[m]



for  n in range(1000):

    hx[:,:] = daHx[:,:] * hx[:,:] - dbHx[:,:] * (ez[:,1:]-ez[:,:size-1])

    hy[:,:] = daHy[:,:] * hy[:,:] + dbHy[:,:] * (ez[1:,:]-ez[:size-1,:])

    Ezx[1:size-1,1:size-1] = caEzx[1:size-1,1:size-1] * Ezx[1:size-1,1:size-1] + cbEzx[1:size-1,1:size-1] * (hy[1:,1:size-1]-hy[:size-2,1:size-1])
    Ezy[1:size-1,1:size-1] = caEzy[1:size-1,1:size-1] * Ezy[1:size-1,1:size-1] - cbEzy[1:size-1,1:size-1] * (hx[1:size-1,1:]-hx[1:size-1,:size-2])
    ez[1:size-1,1:size-1] = Ezx[1:size-1,1:size-1] + Ezy[1:size-1,1:size-1]

    #ez[100,100] = exp(-(n-8)**2/16)

    ez [100,100] = sin(2*pi*f*n*dt)

    if n %10 ==0:
        file = open ('2dpml.%d.txt'%int(n/10),'w')

        for i in range(size):

            for j in range(size):

                file.write("%f "%ez[i,j])

            file.write("\n")

        file.close()

    if n == 150:
        file = open ('2dpmlside.txt','w')
        for i in range(size):
            file.write("%f %f\n"%((i-100)*dx,ez[100,i]))
        file.close()

print "end"
