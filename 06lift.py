import numpy as np
import matplotlib.pyplot as plt
from math import *

N = 50                            # Number of points in each direction
xStart,xEnd = -2.0,2.0            # x-direction boundaries
yStart,yEnd = -1.0,1.0            # y-direction boundaries
x = np.linspace(xStart,xEnd,N)    # x 1D-array
y = np.linspace(yStart,yEnd,N)    # y 1D-array
X,Y = np.meshgrid(x,y)            # generation of the mesh grid

kappa = 1.0                    # strength of the doublet
xDoublet,yDoublet = 0.0,0.0    # location of the doublet

Uinf = 1.0        # freestream speed

# function to compute the velocity components of a doublet
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

# freestream velocity components
uFreestream = Uinf*np.ones((N,N),dtype=float)
vFreestream = np.zeros((N,N),dtype=float)

# stream-function
psiFreestream = Uinf*Y

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
plt.scatter(xDoublet,yDoublet,c='#CD2305',s=80,marker='o')

# cylinder radius
R = sqrt(kappa/(2*pi*Uinf))
circle = plt.Circle((0,0),radius=R,color='#CD2305',alpha=0.5)
plt.gca().add_patch(circle)

# stagnation points
xStagn1,yStagn1 = +sqrt(kappa/(2*pi*Uinf)),0
xStagn2,yStagn2 = -sqrt(kappa/(2*pi*Uinf)),0
plt.scatter([xStagn1,xStagn2],[yStagn1,yStagn2],c='g',s=80,marker='o');

# add vortex
gamma = 4.5                  # strength of the vortex
xVortex0,yVortex = 0.0,0.0    # location of the vortex


# function to compute the velocity components of a vortex
def getVelocityVortex(strength,xv,yv,X,Y):
    u = + strength/(2*pi)*(Y-yv)/((X-xv)**2+(Y-yv)**2)
    v = - strength/(2*pi)*(X-xv)/((X-xv)**2+(Y-yv)**2)
    return u,v

# function to compute the stream-function of a vortex
def getStreamFunctionVortex(strength,xv,yv,X,Y):
    psi = strength/(4*pi)*np.log((X-xv)**2+(Y-yv)**2)
    return psi
    
# computing the velocity components on the mesh grid
uVortex,vVortex = getVelocityVortex(gamma,xVortex0,yVortex,X,Y)

# computing the stream-function on the mesh grid
psiVortex = getStreamFunctionVortex(gamma,xVortex0,yVortex,X,Y)

# superimposition of the doublet and the vortex on the freestream flow
u = uFreestream + uDoublet +uVortex
v = vFreestream + vDoublet +vVortex
psi = psiFreestream + psiDoublet +psiVortex

# cylinder radius
R0 = sqrt(kappa/(2*pi*Uinf))

# stagnation points
xStagn1,yStagn1 = +sqrt(R**2-(gamma/(4*pi*Uinf))**2),-gamma/(4*pi*Uinf)
xStagn2,yStagn2 = -sqrt(R**2-(gamma/(4*pi*Uinf))**2),-gamma/(4*pi*Uinf)

# plotting
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,\
               density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
circle = plt.Circle((0,0),radius=R,color='#CD2305',alpha=0.5)
plt.gca().add_patch(circle)
plt.scatter(xVortex0,yVortex,c='#CD2305',s=80,marker='o')
plt.scatter([xStagn1,xStagn2],[yStagn1,yStagn2],c='g',s=80,marker='o');

# change vortex location
xVortex1,yVortex1 = xStagn1,yStagn1

# function to compute the velocity components of a vortex
def getVelocityVortex(strength,xv,yv,X,Y):
    u = + strength/(2*pi)*(Y-yv)/((X-xv)**2+(Y-yv)**2)
    v = - strength/(2*pi)*(X-xv)/((X-xv)**2+(Y-yv)**2)
    return u,v

# function to compute the stream-function of a vortex
def getStreamFunctionVortex(strength,xv,yv,X,Y):
    psi = strength/(4*pi)*np.log((X-xv)**2+(Y-yv)**2)
    return psi
    
# computing the velocity components on the mesh grid
uVortex,vVortex = getVelocityVortex(gamma,xVortex1,yVortex1,X,Y)

# computing the stream-function on the mesh grid
psiVortex = getStreamFunctionVortex(gamma,xVortex1,yVortex1,X,Y)

# superimposition of the doublet and the vortex on the freestream flow
u = uFreestream + uDoublet +uVortex
v = vFreestream + vDoublet +vVortex
psi = psiFreestream + psiDoublet +psiVortex

# cylinder radius
R1 = sqrt(kappa/(2*pi*Uinf))

# stagnation points
xStagn1,yStagn1 = +sqrt(R**2-(gamma/(4*pi*Uinf))**2),-gamma/(4*pi*Uinf)
xStagn2,yStagn2 = -sqrt(R**2-(gamma/(4*pi*Uinf))**2),-gamma/(4*pi*Uinf)

# plotting
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,\
               density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
circle = plt.Circle((0,0),radius=R,color='#CD2305',alpha=0.5)
plt.gca().add_patch(circle)
plt.scatter(xVortex1,yVortex1,c='#CD2305',s=80,marker='o')
plt.scatter(xDoublet,yDoublet,c='#CD2305',s=80,marker='o')
plt.scatter([xStagn1,xStagn2],[yStagn1,yStagn2],c='g',s=80,marker='o');


# pressure and lift
theta = np.linspace(0,2*pi,100)
utheta = -2*Uinf*np.sin(theta)-gamma/(2*pi*R)

Cp = 1 - (utheta/Uinf)**2

# if there was no vortex
utheta_noVortex = -2*Uinf*np.sin(theta)
Cp_noVortex = 1 - (utheta_noVortex/Uinf)**2

# plotting
size = 6
plt.figure(figsize=(size,size))
plt.grid(True)
plt.xlabel(r'$\theta$',fontsize=18)
plt.ylabel(r'$C_p$',fontsize=18)
plt.xlim(theta.min(),theta.max())
plt.plot(theta,Cp,color='#CD2305',linewidth=2,linestyle='-')
plt.plot(theta,Cp_noVortex,color='g',linewidth=2,linestyle='-')
plt.legend(['with vortex','without vortex'],loc='best',prop={'size':16});

# lift and drag
theta = np.linspace(0,2*pi,100)
utheta = -2*Uinf*np.sin(theta)-gamma/(2*pi*R)

Pinf = 0
rho = 1
Uinf = 1
c = gamma/(2*pi*R)

# pressure with vortex
drag1 = -(Cp*0.5*Uinf**2)*np.sin(theta)*R
lift1 = (0.5*Cp*Uinf**2)*np.cos(theta)*R

# pressure without vortex
drag2 = 4/3*(np.sin(theta))**3+0.2*np.sin(theta)\
    -c*(np.cos(theta))**2*R*0.5
lift2 = (c*theta-.5*np.sin(2.*theta)\
    -3.2*np.cos(theta)+1/3*np.cos(3.0*theta))*R*0.5


# plot 
size = 6
plt.figure(figsize=(size,size))
plt.grid(True)
plt.xlabel(r'$\theta$',fontsize=18)
plt.ylabel(r'force',fontsize=18)
plt.xlim(theta.min(),theta.max())
#plt.ylim()
plt.plot(drag1,lift1,color='red',linewidth=2,linestyle='-')
plt.plot(drag2,lift2,color='g',linewidth=2,linestyle='-')

plt.show()