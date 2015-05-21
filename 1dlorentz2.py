# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
from math import *
import cmath


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

sigma = 0
##
t = 0
dx = 20
dt = 66.7E-9

size = 401

##
edc = np.array([10])
einf = np.array([1])
f0 = np.array([1000e+3])
w0 = f0*2*np.pi
phi = w0/100

##
ap = -phi - cmath.sqrt(-1)*sqrt(w0**2-phi**2)
cp = cmath.sqrt(-1)*(edc-einf)*w0**2/(2*sqrt(w0**2-phi**2))
kp = (1+ap*dt/2)/(1-ap*dt/2)
bp = e0*cp*dt/(1-ap*dt/2)

##
sbeta = 0
for p in range(f0.size):
    sbeta += 2*bp[p].real
##

hy = np.zeros(size)
ez = np.zeros(size)
ezp = np.zeros(size)
Jp = np.zeros((size,f0.size),dtype = complex)
sJ = np.zeros(size)
##
f = 300e+3
maxtime = 200*dt
record = np.zeros((int(maxtime/dt)))

while t <= maxtime:

    ezp[:] = ez[:]

    ez[0] = ez[1]

    sJ[:] = 0
    for p in range(f0.size):

        sJ += ((1+kp[p]) * Jp[:,p]).real

    ez[1:] = ((2*e0*einf*sbeta-sigma*dt)/(2*e0*einf*sbeta+sigma*dt))*ez[1:] + 2*dt*(((hy[1:]-hy[:size-1])/dx)-sJ[1:])/(2*e0*einf+sbeta+sigma*dt)

    ez[200] = exp(-(t-3*20*dt)**2/(20*dt)**2)*sin(2*pi*f*(t-3*20*dt))

    #ez[200] = sin(2*pi*f*(t-3*20*dt))

    for p in range(f0.size):

        Jp[:,p] = kp[p] * Jp[:,p] + bp[p]*(ez[:]-ezp[:])/dt

    hy[size-1] = hy[size-2]

    hy[:size-1] += dt/(u0*dx)*(ez[1:]-ez[:size-1])

    record[t/dt] = ez[50]

    t += dt

yfft = np.fft.fft(record)
xfft = np.linspace(0.0,1/(2*dt),(int(maxtime/dt/2)))

file = open ('1dlorentz1000.txt','w')

for i in range(size):

    file.write("%d %f\n"%(i,ez[i]))

file.close()

file = open ('1dfftlorentz1000.txt','w')

for i in range((int(maxtime/dt/2))):

    file.write("%f %f\n"%(xfft[i],2/(int(maxtime/dt))*abs(yfft[i])))

file.close()
print "end"
