﻿from __future__ import division
import numpy as np
from math import *
import cmath
################################################################################
##final project for general physics in NTU by 顏立峯 and 吳達懿
##dim : 2d
##type : sin / gaussian
##snap : dynamic,from side
##others : none
##command : do for [i=0:45]{load 'ball.'.i.'.plt'; pause 0.2 }
################################################################################
##constants
c = 299792458
u0 = pi*4E-7
e0 = 1/(u0*c**2)
sigma2 = 3.50E+7
sigma1 = 8e-15

##parameter
dx = 5.4e-3
dt = dx/c/sqrt(2)
size = 201

f = 278e+7

##lorentz parameter
edc = np.array([2.52])
einf = np.array([1])
f0 = np.array([29.54E+14/2/np.pi])
w0 = f0*2*np.pi
phi = w0*0

##lorentz or debye parameter
ap = -phi - cmath.sqrt(-1)*sqrt(w0**2-phi**2)
cp = cmath.sqrt(-1)*(edc-einf)*w0**2/(2*sqrt(w0**2-phi**2))
kp = (1+ap*dt/2)/(1-ap*dt/2)
bp = e0*cp*dt/(1-ap*dt/2)

sbeta = 0
for p in range(f0.size):
    sbeta += 2*bp[p].real

C1 = (2*e0*einf*sbeta-sigma2*dt)/(2*e0*einf*sbeta+sigma2*dt)
C2 = (2*e0*einf+sbeta+sigma2*dt)


##plasma
me = 9.109e-31
qe = 1.602e-19
w_pla = 1.78e+13
v_pla = 4.4e+11
sigma_pla = 1j*e0*w_pla**2/((2*pi*f)**2-1j*(2*pi*f)*v_pla)

c_pla = -(w_pla**2)/2
a_pla = -v_pla
k_pla = (1+a_pla*dt/2)/(1-a_pla*dt/2)
b_pla = e0*c_pla*dt/(1-a_pla*dt/2)

kp_pla = np.zeros((size,size))
bp_pla = np.zeros((size,size))
Jp_pla = np.zeros((size,size))
C1_pla = np.zeros((size,size))
C2_pla = np.zeros((size,size))

kp_pla[:,:] = -1

C1_pla[:,:] = (2*e0-sigma1*dt)/(2*e0+sigma1*dt)
C2_pla[:,:] = (2*e0+sigma1*dt)

##TM mode
hx = np.zeros((size,size-1))
hy = np.zeros((size-1,size))
ez = np.zeros((size,size))

##lorentz field
ezp = np.zeros((size,size))
Jp = np.zeros((size,size,f0.size),dtype = complex)
sJ = np.zeros((size,size))

##PML field
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



##range for air
cutl = 10
cutr = 190
cutu = 10
cutd = 190

##range for pokemon
pokl = 50
pokr = 150
poku = 50
pokd = 150
##PML parameter
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

##PML setting
for m in range(npmls):
    sig[m] = sigmax*((m+0.5)/(npmls+0.5))**2
    rho[m] = rhomax*((m+1)/(npmls+0.5))**2
    re = sig[m]*dt/e0
    rm = rho[m]*dt/u0
    ca[m] = exp(-re)
    cb[m] = -(exp(-re)-1)/sig[m]/dx
    da[m] = exp(-rm)
    db[m] = -(exp(-rm)-1)/rho[m]/dx

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


## main part

for  n in range(500):

    ## magnetic field update
    hx[:,:] = daHx[:,:] * hx[:,:] - dbHx[:,:] * (ez[:,1:]-ez[:,:size-1])

    hy[:,:] = daHy[:,:] * hy[:,:] + dbHy[:,:] * (ez[1:,:]-ez[:size-1,:])

    ## PML region update
    Ezx[1:cutl,1:size-1] = caEzx[1:cutl,1:size-1] * Ezx[1:cutl,1:size-1] + cbEzx[1:cutl,1:size-1] * (hy[1:cutl,1:size-1]-hy[:cutl-1,1:size-1])
    Ezy[1:cutl,1:size-1] = caEzy[1:cutl,1:size-1] * Ezy[1:cutl,1:size-1] - cbEzy[1:cutl,1:size-1] * (hx[1:cutl,1:]-hx[1:cutl,:size-2])
    ez[1:cutl,1:size-1] = Ezx[1:cutl,1:size-1] + Ezy[1:cutl,1:size-1]

    Ezx[cutr:size-1,1:size-1] = caEzx[cutr:size-1,1:size-1] * Ezx[cutr:size-1,1:size-1] + cbEzx[cutr:size-1,1:size-1] * (hy[cutr:,1:size-1]-hy[cutr-1:size-2,1:size-1])
    Ezy[cutr:size-1,1:size-1] = caEzy[cutr:size-1,1:size-1] * Ezy[cutr:size-1,1:size-1] - cbEzy[cutr:size-1,1:size-1] * (hx[cutr:size-1,1:]-hx[cutr:size-1,:size-2])
    ez[cutr:size-1,1:size-1] = Ezx[cutr:size-1,1:size-1] + Ezy[cutr:size-1,1:size-1]

    Ezx[cutl:cutr,1:cutu] = caEzx[cutl:cutr,1:cutu] * Ezx[cutl:cutr,1:cutu] + cbEzx[cutl:cutr,1:cutu] * (hy[cutl:cutr,1:cutu]-hy[cutl-1:cutr-1,1:cutu])
    Ezy[cutl:cutr,1:cutu] = caEzy[cutl:cutr,1:cutu] * Ezy[cutl:cutr,1:cutu] - cbEzy[cutl:cutr,1:cutu] * (hx[cutl:cutr,1:cutu]-hx[cutl:cutr,:cutu-1])
    ez[cutl:cutr,1:cutu] = Ezx[cutl:cutr,1:cutu] + Ezy[cutl:cutr,1:cutu]

    Ezx[cutl:cutr,cutd:size-1] = caEzx[cutl:cutr,cutd:size-1] * Ezx[cutl:cutr,cutd:size-1] + cbEzx[cutl:cutr,cutd:size-1] * (hy[cutl:cutr,cutd:size-1]-hy[cutl-1:cutr-1,cutd:size-1])
    Ezy[cutl:cutr,cutd:size-1] = caEzy[cutl:cutr,cutd:size-1] * Ezy[cutl:cutr,cutd:size-1] - cbEzy[cutl:cutr,cutd:size-1] * (hx[cutl:cutr,cutd:size-1]-hx[cutl:cutr,cutd-1:size-2])
    ez[cutl:cutr,cutd:size-1] = Ezx[cutl:cutr,cutd:size-1] + Ezy[cutl:cutr,cutd:size-1]

    ##air update
    ezp[cutl:cutr,cutu:cutd] = ez[cutl:cutr,cutu:cutd]

    sJ[cutl:cutr,cutu:cutd] = 0
    sJ[cutl:cutr,cutu:cutd] += ((1+kp_pla[cutl:cutr,cutu:cutd]) * Jp_pla[cutl:cutr,cutu:cutd]).real

    ez[cutl:pokl,cutu:cutd] = C1_pla[cutl:pokl,cutu:cutd]*ez[cutl:pokl,cutu:cutd] + 2*dt*((((hy[cutl:pokl,cutu:cutd]-hy[cutl-1:pokl-1,cutu:cutd])-(hx[cutl:pokl,cutu:cutd]-hx[cutl:pokl,cutu-1:cutd-1]))/dx)-sJ[cutl:pokl,cutu:cutd])/C2_pla[cutl:pokl,cutu:cutd]
    ez[pokr:cutr,cutu:cutd] = C1_pla[pokr:cutr,cutu:cutd]*ez[pokr:cutr,cutu:cutd] + 2*dt*((((hy[pokr:cutr,cutu:cutd]-hy[pokr-1:cutr-1,cutu:cutd])-(hx[pokr:cutr,cutu:cutd]-hx[pokr:cutr,cutu-1:cutd-1]))/dx)-sJ[pokr:cutr,cutu:cutd])/C2_pla[pokr:cutr,cutu:cutd]
    ez[pokl:pokr,cutu:poku] = C1_pla[pokl:pokr,cutu:poku]*ez[pokl:pokr,cutu:poku] + 2*dt*((((hy[pokl:pokr,cutu:poku]-hy[pokl-1:pokr-1,cutu:poku])-(hx[pokl:pokr,cutu:poku]-hx[pokl:pokr,cutu-1:poku-1]))/dx)-sJ[pokl:pokr,cutu:poku])/C2_pla[pokl:pokr,cutu:poku]
    ez[pokl:pokr,pokd:cutd] = C1_pla[pokl:pokr,pokd:cutd]*ez[pokl:pokr,pokd:cutd] + 2*dt*((((hy[pokl:pokr,pokd:cutd]-hy[pokl-1:pokr-1,pokd:cutd])-(hx[pokl:pokr,pokd:cutd]-hx[pokl:pokr,pokd-1:cutd-1]))/dx)-sJ[pokl:pokr,pokd:cutd])/C2_pla[pokl:pokr,pokd:cutd]

    Jp_pla[cutl:cutr,cutu:cutd] = kp_pla[cutl:cutr,cutu:cutd] * Jp_pla[cutl:cutr,cutu:cutd] + bp_pla[cutl:cutr,cutu:cutd]*(ez[cutl:cutr,cutu:cutd]-ezp[cutl:cutr,cutu:cutd])/dt

    ##pokemon update
    ezp[pokl:pokr,poku:pokd] = ez[pokl:pokr,poku:pokd]

    sJ[pokl:pokr,poku:pokd] = 0
    for p in range(f0.size):

        sJ[pokl:pokr,poku:pokd] += ((1+kp[p]) * Jp[pokl:pokr,poku:pokd,p]).real

    ez[pokl:pokr,poku:pokd] = C1*ez[pokl:pokr,poku:pokd] + 2*dt*((((hy[pokl:pokr,poku:pokd]-hy[pokl-1:pokr-1,poku:pokd])-(hx[pokl:pokr,poku:pokd]-hx[pokl:pokr,poku-1:pokd-1]))/dx)-sJ[pokl:pokr,poku:pokd])/C2

    for p in range(f0.size):

        Jp[pokl:pokr,poku:pokd,p] = kp[p] * Jp[pokl:pokr,poku:pokd,p] + bp[p]*(ez[pokl:pokr,poku:pokd]-ezp[pokl:pokr,poku:pokd])/dt


    ##thunderbolt
    ez [50,100] = 1E+6/dx*(1-2*(np.pi*f*(n*dt-10*dt))**2)*exp(-(np.pi*f*(n*dt-10*dt))**2)

    ez [150,100] = 1E+6/dx*(1-2*(np.pi*f*(n*dt-10*dt))**2)*exp(-(np.pi*f*(n*dt-10*dt))**2)

    ez [100,50] = 1E+6/dx*(1-2*(np.pi*f*(n*dt-10*dt))**2)*exp(-(np.pi*f*(n*dt-10*dt))**2)

    ez [100,150] = 1E+6/dx*(1-2*(np.pi*f*(n*dt-10*dt))**2)*exp(-(np.pi*f*(n*dt-10*dt))**2)


    ##record data
    if n %10 ==0:
        file = open ('ball.%d.txt'%int(n/10),'w')

        for i in range(size):

            for j in range(size):

                file.write("%f "%ez[i,j])

            file.write("\n")

        file.close()

    if n == 150:
        file = open ('ballside.txt','w')
        for i in range(size):
            file.write("%f %f\n"%((i-100)*dx,ez[100,i]))
        file.close()

print "end"
