from __future__ import division
import numpy as np


# Number of samplepoints
N = 600
# sample spacing
T = 1.0 / 800.0
x = np.linspace(0.0, N*T, N)
y = np.sin(50.0 * 2.0*np.pi*x) + 0.5*np.sin(80.0 * 2.0*np.pi*x) + 0.8*np.sin(150.0 * 2.0*np.pi*x)
yf = np.fft.fft(y)
xf = np.linspace(0.0, 1.0/(2.0*T), N/2)

file = open('ffttest.txt','w')

for i in range(int(N/2)):
    file.write("%f %f\n" %(xf[i],2/N*abs(yf[i])))

file.close()

print "end"
