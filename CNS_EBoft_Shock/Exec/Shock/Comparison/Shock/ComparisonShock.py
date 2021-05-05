import numpy as np
from math import *
#from pycse import bvp
import matplotlib.pyplot as plt

#The velocity of the piston is up=Ut/t0

gam=1.4
U=25.0
t0=2e-4
c0=340.15
xmin=-0.1
xmax=0.9
n=1000
delx=(xmax-xmin)/(n-1)

t=0.001

delt=5.0e-6

for filenum in np.arange(150,550,100):
	print filenum
	shock_num_filename="x_velocity%05d.txt"%filenum
        shock_num=np.loadtxt(shock_num_filename);
	n = np.size(shock_num[:,0])
	u = np.zeros(n)
	x = np.zeros(n)
	uexact=[]
	unum = []

        t=filenum*delt
        for i in np.arange(0,n,1):
		unum.append((shock_num[i,0],shock_num[i,1]))
                x[i]=shock_num[i,0]
		
                up=U*t/t0
		#if(x[i] > 0.5*U/t0*t**2 and 1/gam**2*(c0-(gam+1)/2*up)**2-2/gam*U/t0*(x[i]-c0*t) >= 0 and x[i]-c0*t<0.0):
		if(x[i] > 0.5*U/t0*t**2 and x[i]-c0*t<0.0):
			u[i]=-1.0/gam*(c0-(gam+1)/2*up)+(1/gam**2*(c0-(gam+1)/2*up)**2-2/gam*U/t0*(x[i]-c0*t))**0.5
		        #if(u[i]<0):
	               # 	u[i]=0.0
		else:
			u[i]=0.0
		uexact.append((x[i],u[i])) 

                #imgfilename="./Images/Shock%04d.png"%filenum
        imgfilename="./Images/ShockComparison.png"%filenum

	uexact=sorted(uexact, key = lambda x:float(x[0]))
	unum = sorted(unum, key = lambda x:float(x[0]))
		
	for i in np.arange(0,n,1):
		x[i] = uexact[i][0]
		u[i] = uexact[i][1]
		shock_num[i,0] = unum[i][0]
		shock_num[i,1] = unum[i][1]

        #plt.xlim([-0.45,0.45]);
        plt.ylim([-10.0,1.05*max(u)]);
        plt.plot(x,u,'k',linewidth=2)
        plt.plot(shock_num[0:-1:2,0],shock_num[0:-1:2,1],'ob',markersize=5)
        plt.xlabel('$x (m)$',fontsize=20)
        plt.ylabel('$u (m/s)$',fontsize=20)
        plt.gca().legend(('Exact','Present'),'upper left')
        plt.savefig(imgfilename)
	#plt.close()
	#plt.show()


