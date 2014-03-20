import numpy as np
import matplotlib.pyplot as plt
from math import *

Uinf = 1.0

# definition of the cylinder
R = 1.0
theta = np.linspace(0,2*pi,100)
xCylinder,yCylinder = R*np.cos(theta),R*np.sin(theta)

# plot cylinder
size = 4
plt.figure(num=None,figsize=(size,size))
plt.grid(True)
plt.xlim(-1.1,1.1)
plt.ylim(-1.1,1.1)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.plot(xCylinder,yCylinder,c='b',ls='-',lw=2);

class Panel:
    def __init__(self,xa,ya,xb,yb):
        self.xa,self.ya = xa,ya             # 1st end-point
        self.xb,self.yb = xb,yb             # 2nd end-point
        self.xc,self.yc = (xa+xb)/2,(ya+yb)/2       # control point
        self.length = sqrt((xb-xa)**2+(yb-ya)**2)   # length of the panel
        
        # orientation of the panel
        if (xb-xa<=0.): self.beta = acos((yb-ya)/self.length)
        elif (xb-xa>0.): self.beta = pi+acos(-(yb-ya)/self.length)
        
        self.sigma = 0.             # source strength
        self.Vt = 0.                # tangential velocity
        self.Cp = 0.                # pressure coefficient
        
Np = 10 # number of panels around cylinder

# defining the end points of the panels
xb = R*np.cos(np.linspace(0,2*pi,Np+1))
yb = R*np.sin(np.linspace(0,2*pi,Np+1))

# defining the panels
panel = np.empty(Np,dtype=object)
for i in range(Np):
    panel[i] = Panel(xb[i],yb[i],xb[i+1],yb[i+1])
    
# plotting the discretization
size = 6
plt.figure(num=None,figsize=(size,size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.plot(xCylinder,yCylinder,c='b',ls='-',lw=1)
plt.plot(xb,yb,c='#CD2305',ls='-',lw=2)
plt.scatter([p.xa for p in panel],[p.ya for p in panel],c='#CD2305',s=40)
plt.scatter([p.xc for p in panel],[p.yc for p in panel],c='k',s=40,zorder=3)
plt.legend(['cylinder','panels','end points','center points'],\
            loc='best',prop={'size':16})
plt.xlim(-1.1,1.1)
plt.ylim(-1.1,1.1);