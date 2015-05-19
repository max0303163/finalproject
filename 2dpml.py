from __future__ import division
import numpy as np
from math import *

################################################################################
##final project for general physics in NTU by 顏立峯 and 吳達懿
##dim : 2d
##type : sin / gaussian
##snap : dynamic,from side
##others : pml
##command : do for [i=0:45]{load '2d.'.i.'.plt'; pause 0.2 }
################################################################################
##parameter1
c = 299792458
u0 = pi*4E-7
e0 = 1/(u0*c**2)

##

dx = dy = 1
dt = dx/c/sqrt(2)

size = 201

##
hx = np.zeros((size,size-1))
hy = np.zeros((size-1,size))
ez = np.zeros((size,size))

##
f = 15e+6
sig = 0
sigm = 0

##
sym = sxm = sy = sx = np.zeros((size))
aex = bex = ahx = bhx = np.zeros((size,1))
aey = bey = ahy = bhy = np.zeros((1,size))
Phy = Phx = np.zeros((size-1,size-1))
Pex = Pey = np.zeros((size-2,size-2))

##
m = 4
sxmax = -(m+1)*log(1e-6)/2/377/(10*dx)

##
for mm in range(10):
    sy[mm+1] = sxmax*((11-mm-0.5)/10)**4
    sym[mm] = sxmax*((11-mm)/10)**4
    sy[size-1-mm] = sxmax*((11-mm-0.5)/10)**4
    sym[size-1-mm] = sxmax*((11-mm)/10)**4
    sx[mm+1] = sxmax*((11-mm-0.5)/10)**4
    sxm[mm] = sxmax*((11-mm)/10)**4
    sx[size-1-mm] = sxmax*((11-mm-0.5)/10)**4
    sxm[size-1-mm] = sxmax*((11-mm)/10)**4
    
##
aex[:,0] = np.exp(-sx[:]*dt/e0) - 1
bex[:,0] = aex[:,0] + 1
aey[0,:] = np.exp(-sy[:]*dt/e0) - 1
bey[0,:] = aey[0,:] + 1
ahx[:,0] = np.exp(-sxm[:]*dt/e0) - 1
bhx[:,0] = ahx[:,0] + 1
ahy[0,:] = np.exp(-sym[:]*dt/e0) - 1
bhy[0,:] = ahy[0,:] + 1

##
bhy = np.tile(bhy,(size-1,1))
ahy = np.tile(ahy,(size-1,1))
bey = np.tile(bey,(size-2,1))
aey = np.tile(aey,(size-2,1))
bhx = np.tile(bhx,(1,size-1))
ahx = np.tile(ahx,(1,size-1))
bex = np.tile(bex,(1,size-2))
aex = np.tile(aex,(1,size-2))


##

for  n in range(10):
    
    Phy[:,:] = bhy[:,:size-1]*Phy[:,:] + ahy[:,:size-1]*(ez[:size-1,1:]-ez[:size-1,:size-1])/dy

    Phx[:,:] = bhx[:size-1,:]*Phx[:,:] + ahx[:size-1,:]*(ez[1:,:size-1]-ez[:size-1,:size-1])/dx

    hx[:size-1,:] = hx[:size-1,:] - dt/(u0*dx)*(ez[:size-1,1:]-ez[:size-1,:size-1]) - dt/u0*Phy[:,:]

    hy[:,:size-1] = hy[:,:size-1] + dt/(u0*dx)*(ez[1:,:size-1]-ez[:size-1,:size-1]) + dt/u0*Phx[:,:]

    print Phy[100,:]

    Pex[:,:] = bex[1:size-1,:]*Pex[:,:] + aex[1:size-1,:]*(hy[:size-2,1:size-1]-hy[1:,1:size-1])/dx

    Pey[:,:] = bey[:,1:size-1]*Pey[:,:] + aey[:,1:size-1]*(hx[1:size-1,:size-2]-hx[1:size-1,1:])/dy

    

    ez[1:size-1,1:size-1] = ez[1:size-1,1:size-1] + dt/(e0*dx)*(
        (hy[1:,1:size-1]-hy[:size-2,1:size-1])-(
            hx[1:size-1,1:]-hx[1:size-1,:size-2])) + dt/u0*(Pex[:,:] - Pey[:,:])

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
