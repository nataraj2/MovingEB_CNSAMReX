import numpy as np
from math import *
#from pycse import bvp
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from scipy import *

data=loadtxt('DensityVsTime.txt')
datasize = size(data[:,0])
time = zeros(datasize);
volume = zeros(datasize);
density = zeros(datasize)
pressure = zeros(datasize)

time = data[:,0]
volume = data[:,1];
#density = data[:,2]
#pressure = data[:,3]

skipval=1

for i in arange(0,datasize,1):
	density[i]=1.226*volume[0]/volume[i]
	pressure[i] = 101325.0/1.226**1.4*density[i]**1.4

plt.figure(1)
plt.plot(data[1:datasize:skipval,0],data[1:datasize:skipval,2],marker='o',color='k',label='CNS')
plt.plot(time,density,linestyle='-',color='b',label='Exact')
plt.xlabel('Time (s)')
plt.ylabel('Density (kg/m$^3$)')
plt.gca().legend(loc=0)
imgfilename="./Images/DensityVsTime.png"
plt.savefig(imgfilename)

plt.figure(2)
plt.plot(data[1:datasize:skipval,0],data[1:datasize:skipval,3],marker='o',color='k',label='CNS')
plt.plot(time,pressure,linestyle='-',color='b',label='Exact')
plt.xlabel('Time (s)')
plt.ylabel('Pressure (N/m$^2$)')
plt.gca().legend(loc=0)
imgfilename="./Images/PressureVsTime.png"
plt.savefig(imgfilename)
plt.show()
