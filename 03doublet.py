import numpy as np
import matplotlib.pyplot as plt
from math import*
from IPython.core.display import clear_output

N = 50                      # Number of points in each direction
xStart,xEnd = -2.0,2.0      # x-direction boundaries
yStart,yEnd = -1.0,1.0      # y-direction boundaries
x = np.linspace(xStart,xEnd,N)  # x 1D array
y = np.linspace(yStart,yEnd,N)  # y 1D array
X,Y = np.meshgrid(x,y)

class Doublet:
    def __init__(self,strength,x,y):
        self.strength = strength
        self.x,self.y = x,y
        
    # get velocity field
    def velocity(self,X,Y):
        self.u = -self.strength/(2*pi)*\
                ((X-self.x)**2-(Y- self.y)**2)/((X-self.x)**2+(Y-self.y)**2)**2
        self.v = -self.strength/(2*pi)*\
                2*(X-self.x)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)**2
                
    # get stream function
    def streamFunction(self,X,Y):
        self.psi = -self.strength/(2*pi)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)
        
kappa = 1.0  #strength of the doublet
xDoublet,yDoublet = 0.0,0.0  # location of the doublet

# function to compute the celocity components of a doublet
def getVelocityDoublet(strength,xd,yd,X,Y):
    u = - strength/(2*pi)*((X-xd)**2-(Y-yd)**2)/((X-xd)**2+(Y-yd)**2)**2
    v = - strength/(2*pi)*2*(X-xd)*(Y-yd)/((X-xd)**2+(Y-yd)**2)**2
    return u,v
    
# function to compute the stream-function of a doublet
def getStreamFunctionDoublet(strength,xd,yd,X,Y):
    psi = - strength/(2*pi)*(Y-yd)/((X-xd)**2+(Y-yd)**2)
    return psi
    
# computing the velocity components on the mesh grid
uDoublet,vDoublet = getVelocityDoublet(kappa,xDoublet,yDoublet,X,Y)

# computing the stream-function on the mesh grid
psiDoublet = getStreamFunctionDoublet(kappa,xDoublet,yDoublet,X,Y)

# plotting
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,uDoublet,vDoublet,\
               density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xDoublet,yDoublet,c='#CD2305',s=80,marker='o')


# Uniform flow past a doublet

Uinf = 1.0   # freestream speed
alpha = 0.0   # angle of attack (in degrees)

uFreestream = Uinf*cos(alpha*pi/180)*np.ones((N,N),dtype=float)
vFreestream = Uinf*sin(alpha*pi/180)*np.ones((N,N),dtype=float)

psiFreestream = Uinf*(cos(alpha*pi/180)*Y-sin(alpha*pi/180)*X)

# superimposition of the doublet on the freestream flow
u = uFreestream + uDoublet
v = vFreestream + vDoublet
psi = psiFreestream + psiDoublet

# plotting
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,\
               density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.contour(X,Y,psi,levels=[0.0],colors='#CD2305',linewidths=2,linestyles='solid')
plt.scatter(xDoublet,yDoublet,c='#CD2305',s=80,marker='o')

# compute the pressure coefficient
Cp = 1.0-(u**2+v**2)/Uinf**2

# plotting
size = 10
plt.figure(num=0,figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
contf = plt.contourf(X,Y,Cp,levels=np.linspace(-2.0,1.0,100),extend='both')
cbar = plt.colorbar(contf)
cbar.set_label(r'$C_p$',fontsize=16)
cbar.set_ticks([-2.0,-1.0,0.0,1.0])
plt.scatter(xDoublet,yDoublet,c='#CD2305',s=80,marker='o')
plt.contour(X,Y,psi,\
            levels=[0.0],\
            colors='#CD2305',linewidths=2,linestyles='solid')
            
            
plt.show()