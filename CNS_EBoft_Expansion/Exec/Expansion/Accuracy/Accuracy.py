import numpy as np
from math import *
#from pycse import bvp
import matplotlib.pyplot as plt
from numpy import *

line = loadtxt('Expansion_512.curve')
x   = line[:,0]
s = line[:,1]
n = size(line[:,0])
errorvec = zeros(n)

error = 0.0
counter = 0

for i in arange(0,n,1):
	if(s[i]>0):	
		counter = counter+1
        	error = error + log(s[i])**2
		errorvec[i] = abs(log(s[i]))

error512 = (error/counter)**0.5
linfnorm512 = max(errorvec)
l1norm512 = sum(errorvec)/counter
print error512

line = loadtxt('Expansion_256.curve')
x   = line[:,0]
s = line[:,1]
n = size(line[:,0])
errorvec = zeros(n)

error = 0.0
counter = 0
for i in arange(0,n,1):
	if(s[i]>0):	
		counter = counter+1
        	error = error + log(s[i])**2
		errorvec[i] = abs(log(s[i]))

error256 = (error/counter)**0.5
linfnorm256 = max(errorvec)
l1norm256 = sum(errorvec)/counter
print error512

line = loadtxt('Expansion_128.curve')
x   = line[:,0]
s = line[:,1]
n = size(line[:,0])
errorvec = zeros(n)

error = 0.0
counter=0
for i in arange(0,n,1):
	if(s[i]>0):	
		counter = counter+1
        	error = error + log(s[i])**2
		errorvec[i] = abs(log(s[i]))

error128 = (error/counter)**0.5
linfnorm128 = max(errorvec)
l1norm128 = sum(errorvec)/counter
print error128

line = loadtxt('Expansion_064.curve')
x   = line[:,0]
s = line[:,1]
n = size(line[:,0])
errorvec = zeros(n)

error = 0.0
counter = 0
for i in arange(0,n,1):
	if(s[i]>0):	
		counter = counter+1
        	error = error + log(s[i])**2
		errorvec[i] = abs(log(s[i]))

error64 = (error/counter)**0.5
linfnorm64 = max(errorvec)
l1norm64 = sum(errorvec)/counter
print error64


line = loadtxt('Expansion_032.curve')
x   = line[:,0]
s = line[:,1]
n = size(line[:,0])
errorvec = zeros(n)

error = 0.0
counter = 0
for i in arange(0,n,1):
	if(s[i]>0):	
		counter = counter+1
        	error = error + log(s[i])**2
		errorvec[i] = abs(log(s[i]))

error32 = (error/counter)**0.5
linfnorm32 = max(errorvec)
l1norm32 = sum(errorvec)/counter
print error32

dx = zeros(5)
errorl1 = zeros(5)
errorl2 = zeros(5)

dx[0] =32
dx[1] = 64
dx[2] = 128
dx[3] = 256
dx[4] = 256

errorl2[0] = error32
errorl2[1] = error64
errorl2[2] = error128
errorl2[3] = error256
errorl2[4] = error256

#error[0] = linfnorm32
#error[1] = linfnorm64
#error[2] = linfnorm128
#error[3] = linfnorm256
#error[4] = linfnorm256

errorl1[0] = l1norm32
errorl1[1] = l1norm64
errorl1[2] = l1norm128
errorl1[3] = l1norm256
errorl1[4] = l1norm256


#plt.loglog(dx,errorl1,'-or',markersize=10,label='L1')
plt.loglog(dx,errorl2,'-sk',markersize=10,label='L2')


slope = zeros([2,2])

slope[0,0] = 10**1.5
slope[0,1] = 10**(-0.4)
slope[1,0] = 10**2.5
slope[1,1] = 10**(-1.4)

plt.plot(slope[:,0],slope[:,1],'--k',linewidth=2)
plt.xlim([20,300])
#plt.ylim([0.9*min(errorl2),1.1*max(errorl2)])
#plt.legend()
#plt.title("Original CNS, Accuracy=2.0",fontsize=30)
plt.ylabel(r'$L^2$ norm',fontsize=20)
plt.xlabel(r'Number of mesh points',fontsize=20)
plt.text(90, 0.15,'~$N^{-1}$', fontsize=20)
plt.savefig('AccuracyMovingEB.png')
plt.show()


