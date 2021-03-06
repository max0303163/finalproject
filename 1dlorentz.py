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

##
t = 0
dx = 20
dt = 66.7E-9

size = 401

##
edc = np.array([10])
einf = np.array([1])
f0 = np.array([70e+3])
w0 = f0*2*np.pi
fhi = w0/100

##
A1 = np.zeros(f0.size)
A2 = np.zeros(f0.size)
A3 = np.zeros(f0.size)
A1 = (2-(w0*dt)**2)/(1+fhi*dt)
A2 = (fhi*dt-1)/(fhi*dt+1)
A3 = e0*(edc-einf)*(w0*dt)**2/(1+fhi*dt)

##
C1 = 2*e0*einf/(2*e0*einf+0.5*np.sum(A3,dtype=np.float64))
C2 = 0.5*np.sum(A3,dtype=np.float64)/(2*e0*einf+0.5*np.sum(A3,dtype=np.float64))
C3 = 2*dt/(2*e0*einf+0.5*np.sum(A3,dtype=np.float64))

##

hy = np.zeros(size)
ez = np.zeros(size)
ezp = np.zeros(size)
ezpp = np.zeros(size)
Jz = np.zeros((size,f0.size))
Jzp = np.zeros((size,f0.size))
Jzpp = np.zeros((size,f0.size))
##
f = 300e+3
maxtime = 200*dt
record = np.zeros((int(maxtime/dt)))

while t <= maxtime:

    hy[size-1] = hy[size-2]

    hy[:size-1] += dt/(u0*dx)*(ez[1:]-ez[:size-1])

    ezpp[:] = ezp[:]

    ezp[:] = ez[:]

    ez[0] = ez[1]

    ez[1:] = C1[0]*ez[1:] + C2[0]*ezpp[1:] + C3[0]*(hy[1:]-hy[:size-1])/dx

    for m in range(f0.size):

        ez[1:] -= C3*0.5*((1+A1[m])*Jz[1:,m]+A2[m]*Jzp[1:,m])

    #ez[200] = exp(-(t/dt-8)**2/16)

    ez[200] = exp(-(t-3*20*dt)**2/(20*dt)**2)*sin(2*pi*f*(t-3*20*dt))

    Jzpp[:,:] = Jzp[:,:]

    Jzp[:,:] = Jz[:]

    for m in range(f0.size):

        Jz[:,m] = A1[m]*Jz[:,m] + A2[m]*Jzpp[:,m] + A3[m]*(ez[:]-ezp[:])/(2*dt)

    record[t/dt] = ez[50]

    t += dt

yfft = np.fft.fft(record)
xfft = np.linspace(0.0,1/(2*dt),(int(maxtime/dt/2)))

file = open ('1dlorentz.txt','w')

for i in range(size):

    file.write("%d %f\n"%(i,ez[i]))

file.close()

file = open ('1dfftlorentx.txt','w')

for i in range((int(maxtime/dt/2))):

    file.write("%f %f\n"%(xfft[i],2/(int(maxtime/dt))*abs(yfft[i])))

file.close()
print "end"
