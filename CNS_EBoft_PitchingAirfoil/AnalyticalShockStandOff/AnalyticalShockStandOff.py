import numpy as np
from math import *
#from pycse import bvp
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from scipy import *

n = 100

ymin = -0.25
ymax =  0.25
dely = (ymax-ymin)/n
x = zeros(n)
y = zeros(n)

mach = 1.9411
theta = asin(1.0/mach)
R = 0.02
Delta = 0.386*R*exp(4.67/mach**2)
Rc = R*1.386*exp(1.8/(mach-1)**0.75)

for i in arange(0,n,1):
	y[i] = ymin + i*dely
	x[i] = R + Delta - Rc*1.0/tan(theta)**2*((1.0+(y[i]*tan(theta)/Rc)**2)**0.5-1)

plt.plot(-x/R,y/R,'o')
plt.xlim([-0.5/R,0.5/R])
plt.ylim([-0.25/R,0.25/R])
plt.show()



