import numpy as np
import matplotlib.pyplot as plt
from math import *

# creat mesh grid
N=50
xStart,xEnd = -2.0,2.0 
yStart,yEnd = -1.0,1.0
x = np.linspace(xStart,xEnd,N)
y = np.linspace(yStart,yEnd,N)
X,Y = np.meshgrid(x,y)

# Source
class Source:
    def __init__(self,strength,x,y):
        self.strength = strength
        self.x,self.y = x,y
    # get the velocity field
    def velocity(self,X,Y):
        self.u = self.strength/(2*pi)*(X-self.x)/((X-self.x)**2+(Y-self.y)**2)
        self.v = self.strength/(2*pi)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)
    # get the stream-function
    def streamFunction(self,X,Y):
        self.psi = self.strength/(2*pi)*np.arctan2((Y-self.y),(X-self.x))
        
strength = 1.0
xSource,ySource = 0.0,0.5


source1 = Source(strength,xSource,ySource)
source1.velocity(X,Y)
source1.streamFunction(X,Y)

source2 = Source(strength,xSource,-ySource)
source2.velocity(X,Y)
source2.streamFunction(X,Y)

# Superimposition
u = source1.u + source2.u
v = source1.v + source2.v
psi = source1.psi + source2.psi

# Plotting
size = 10
plt.figure(num=0, figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot( X,Y, u,v,\
               density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(source1.x, source1.y, c='#CD2305', s=80, marker='o')
plt.scatter(source2.x, source2.y, c='#CD2305', s=80, marker='D')
plt.axhline(0.0, color='k', linestyle='--', linewidth=4);


# Vortex
class Vortex:
    def __init__(self,strength,x,y):
        self.strength = strength
        self.x,self.y = x,y
        
    # get velocity field
    def velocity(self,X,Y):
        self.u = +self.strength/(2*pi)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)
        self.v = -self.strength/(2*pi)*(X-self.x)/((X-self.x)**2+(Y-self.y)**2)
        
    # get streamfunction
    def streamFunction(self,X,Y):
        self.psi = -self.strength/(4*pi)*np.log((X-self.x)**2+(Y-self.y)**2)

strengthVortex = 1.0
xVortex,yVortex = 0.0,0.5

Vortex1 = Vortex(strengthVortex,xVortex,yVortex)
Vortex1.velocity(X,Y)
Vortex1.streamFunction(X,Y)

Vortex2 = Vortex(-strengthVortex,xVortex,-yVortex)
Vortex2.velocity(X,Y)
Vortex2.streamFunction(X,Y)

# Superimposition
u = Vortex1.u + Vortex2.u
v = Vortex1.v + Vortex2.v
psi = Vortex1.psi + Vortex2.psi

# Plot
size = 10
plt.figure(num=0, figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot( X,Y, u,v,\
               density=2.0, linewidth=1, arrowsize=1, arrowstyle='->')
plt.scatter( Vortex1.x, Vortex1.y, c='#CD2305', s=80, marker='o')
plt.scatter( Vortex2.x, Vortex2.y, c='#CD2305', s=80, marker='D')
plt.axhline(0.0, color='k', linestyle='--', linewidth=4);

# Vortex pairs
# Vortex it self
strengthVortex = 1.0
xVortex3,yVortex3 = -0.1,0.5
xVortex4,yVortex4 = 0.1,0.5

vortex3 = Vortex(+strengthVortex,xVortex3,yVortex3)
vortex4 = Vortex(-strengthVortex,xVortex3,yVortex4)

vortex3.velocity(X,Y)
vortex3.streamFunction(X,Y)
vortex4.velocity(X,Y)
vortex4.streamFunction(X,Y)

# Vortex Reflection
vortexImage3 = Vortex(-strengthVortex,xVortex3,-yVortex3)
vortexImage4 = Vortex(+strengthVortex,xVortex4,-yVortex4)

vortexImage3.velocity(X,Y)
vortexImage3.streamFunction(X,Y)
vortexImage4.velocity(X,Y)
vortexImage4.streamFunction(X,Y)

# Superimposition
u = vortex3.u + vortex4.u + vortexImage3.u + vortexImage4.u
v = vortex3.v + vortex4.v + vortexImage3.v + vortexImage4.v
psi = vortex3.psi + vortex4.psi + vortexImage3.psi + vortexImage4.psi

# Plot
size = 10
plt.figure(num=0,figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,\
               density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(vortex3.x,vortex3.y,c='#CD2305',s=80,marker='o')
plt.scatter(vortex4.x,vortex4.y,c='g',s=80,marker='o')
plt.scatter(vortexImage3.x,vortexImage3.y,c='#CD2305',s=80,marker='D')
plt.scatter(vortexImage4.x,vortexImage4.y,c='g',s=80,marker='D')
plt.axhline(0.0,color='k',linestyle='--',linewidth=4);

# Doublet
Uinf = 1.0

uFreestream = Uinf*np.ones((N,N),dtype=float)
vFreestream = np.zeros((N,N),dtype=float)

psiFreestream = Uinf*Y

class Doublet:
    def __init__(self,strength,x,y):
        self.strength = strength
        self.x,self.y = x,y
        
    # get velocity field
    def velocity(self,X,Y):
        self.u = -self.strength/(2*pi)*\
                ((X-self.x)**2-(Y-self.y)**2)/((X-self.x)**2+(Y-self.y)**2)**2
        self.v = -self.strength/(2*pi)*\
                2*(X-self.x)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)**2
            
    # get stream function
    def streamFunction(self,X,Y):
        self.psi = -self.strength/(2*pi)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)
       
strengthDoublet = 1.0
xDoublet,yDoublet = 0.0,0.3

doublet = Doublet(strengthDoublet,xDoublet,yDoublet)
doublet.velocity(X,Y)
doublet.streamFunction(X,Y)

doubletImage = Doublet(strengthDoublet,xDoublet,-yDoublet)
doubletImage.velocity(X,Y)
doubletImage.streamFunction(X,Y)

# Superimposition 
u = uFreestream + doublet.u + doubletImage.u
v = vFreestream + doublet.v + doubletImage.v
psi = psiFreestream + doublet.psi + doubletImage.psi

# Plot
size = 10
plt.figure(num=0,figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,\
               density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(doublet.x,doublet.y,c='r',s=80,marker='o')
plt.scatter(doubletImage.x,doubletImage.y,c='r',s=80,marker='x')
plt.axhline(0.0,color='k',linestyle='--',linewidth=4);
plt.show()