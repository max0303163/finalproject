# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
from math import *


################################################################################
##final project for general physics in NTU by 顏立峯 and 吳達懿
##dim : 1d
##type : gaussian
##snap : one time/fft
##other : ABC boundary/debye material
################################################################################
##parameter
c = 299792458
u0 = pi*4E-7
e0 = 1/(u0*c**2)
loe=0
lom=0
##
t = 0
dx = 20
dt = 66.7E-9

size = 401

##
edc = 10
einf = 1
tau = 1

##
hy = np.zeros(size)
ez = np.zeros(size)
dz = np.zeros(size)
dzp = np.zeros(size)
##
f = 278e+6
maxtime = 200*dt
record = np.zeros((int(maxtime/dt)))

while t <= maxtime:

    hy[size-1] = hy[size-2]

    hy[:size-1] += dt/(u0*dx)*(ez[1:]-ez[:size-1])

    dzp[:] = dz[:]

    dz[1:] += dt/dx*(hy[1:]-hy[:size-1])

    ez[0] = ez[1]

    ez[1:] = ((dt+2*tau)/(2*tau*e0*einf+e0*edc*dt))*dz[1:] + ((dt-2*tau)/(2*tau*e0*einf+e0*edc*dt)) *dzp[1:] + ((2*tau*e0*einf-e0*edc*dt)/(2*tau*e0*einf+e0*edc*dt)) * ez[1:]

    ez[200] = exp(-(t/dt-8)**2/16)

    record[t/dt] = ez[50]

    t += dt

yfft = np.fft.fft(record)
xfft = np.linspace(0.0,1/(2*dt),(int(maxtime/dt/2)))

file = open ('1ddebye.txt','w')

for i in range(size):

    file.write("%d %f\n"%(i,ez[i]))

file.close()

file = open ('1dfftdebye.txt','w')

for i in range((int(maxtime/dt/2))):

    file.write("%f %f\n"%(xfft[i],2/(int(maxtime/dt))*abs(yfft[i])))

file.close()
print "end"
