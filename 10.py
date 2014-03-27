import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from math import *

# read of the geometry from a data file
coords = np.loadtxt(fname='/Users/Eva/Documents/ip/code/resources/naca0012.dat')
xp,yp = coords[:,0],coords[:,1]

# plotting the geometry
valX,valY = 0.1,0.2
xmin,xmax = min(xp),max(xp)
ymin,ymax = min(yp),max(yp)
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.plot(xp,yp,'k-',linewidth=2);

# class Panel containing the info about one panel
class Panel:
    def __init__(self,xa,ya,xb,yb):
        self.xa,self.ya = xa,ya                     # 1st end-point
        self.xb,self.yb = xb,yb                     # 2nd end-point
        self.xc,self.yc = (xa+xb)/2,(ya+yb)/2       # control point
        self.length = sqrt((xb-xa)**2+(yb-ya)**2)   # length of the panel
        
        # orientation of the panel
        if (xb-xa<=0.): self.beta = acos((yb-ya)/self.length)
        elif (xb-xa>0.): self.beta = pi+acos(-(yb-ya)/self.length)
        
        # location of the panel
        if (self.beta<=pi): self.loc = 'extrados'
        else: self.loc = 'intrados'
        
        self.sigma = 0.                             # source strength
        self.vt = 0.                                # tangential velocity
        self.Cp = 0.                                # pressure coefficient

# function to descretize the geometry into panels
def definePanels(N,xp,yp):
    Np = len(xp)
    length = sum(np.sqrt((xp[1:Np]-xp[0:Np-1])**2+(yp[1:Np]-yp[0:Np-1])**2))
    length += sqrt((xp[0]-xp[Np-1])**2+(yp[0]-yp[Np-1])**2)
    lc = length/N
    
    x,y = np.empty(N,dtype=float),np.empty(N,dtype=float)
    x[0],y[0] = xp[0],yp[0]
    i = 0
    for j in range(N-1):
        xStart,yStart = x[j],y[j]
        xEnd,yEnd = xp[i+1],yp[i+1]
        d = sqrt((xStart-xEnd)**2+(yStart-yEnd)**2)
        if (d>=lc): d = 0
        elif (d<lc):
            while (i<Np-1):
                i += 1
                xStart,yStart = xEnd,yEnd
                if (i<Np-2): xEnd,yEnd = xp[i+1],yp[i+1]
                else: xEnd,yEnd = xp[0],yp[0]
                dd = sqrt((xStart-xEnd)**2+(yStart-yEnd)**2)
                if (d+dd<lc): d += dd        # +=
                else: break
        x[j+1] = xStart+(lc-d)/sqrt((xEnd-xStart)**2+(yEnd-yStart)**2)*(xEnd-xStart)
        y[j+1] = yStart+(lc-d)/sqrt((xEnd-xStart)**2+(yEnd-yStart)**2)*(yEnd-yStart)
                
    panel = np.empty(N,dtype=object)
    for i in range(N-1):
        panel[i] = Panel(x[i],y[i],x[i+1],y[i+1])
        panel[N-1] = Panel(x[N-1],y[N-1],x[0],y[0])
    
    return panel

N = 100                        # number of panels
panel = definePanels(N,xp,yp)  # discretization of the geometry into panels

# plotting the geometry with the panels
valX,valY = 0.1,0.2
xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
ymin,ymax = min([p.ya for p in panel]),max([p.ya for p in panel])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.plot(xp,yp,'k-',linewidth=2)
plt.plot(np.append([p.xa for p in panel],panel[0].xa),\
        np.append([p.ya for p in panel],panel[0].ya),\
        linestyle='-',linewidth=1,\
        marker='o',markersize=6,color='#CD2305');
        
# class Freestream containing the freestream conditions
class Freestream:
        def __init__(self,Uinf,alpha):
	   self.Uinf = Uinf                   # velocity magnitude
	   self.alpha = alpha*pi/180          # angle of attack (degrees --> radians)

# definition of the object freestream
Uinf = 1.0                               # freestream speed
alpha = 5.0                              # angle of attack (in degrees)
freestream = Freestream(Uinf,alpha)      # instantiation of the object freestream

# function to evaluate the integral Iij(zi)
def I(xci,yci,pj,dxdz,dydz):
    def func(s):
	return (+(xci-(pj.xa-sin(pj.beta)*s))*dxdz\
		+(yci-(pj.ya+cos(pj.beta)*s))*dydz)\
		/((xci-(pj.xa-sin(pj.beta)*s))**2\
	       + (yci-(pj.ya+cos(pj.beta)*s))**2)
    return integrate.quad(lambda s:func(s),0.,pj.length)[0]
	
# function to build the source matrix
def buildMatrix(p):
	N = len(p)
	A = np.empty((N,N),dtype=float)
	np.fill_diagonal(A,0.5)
	for i in range(N):
	   for j in range(N):
		if (i!=j):
		  A[i,j] = 0.5/pi*I(p[i].xc,p[i].yc,p[j],cos(p[i].beta),sin(p[i].beta))
	return A

# function to build the right hand-side of the linear system
def buildRHS(p,fs):
    N = len(p)
    B = np.zeros(N,dtype=float)
    for i in range(N):
	   B[i] = -fs.Uinf*cos(fs.alpha-p[i].beta)
    return B
	
A = buildMatrix(panel)					# calculate the singularity matrix
B = buildRHS(panel,freestream)			# calculate the freestream RHS

	
# solve the linear system
var = np.linalg.solve(A,B)
for i in range(len(panel)):
	panel[i].sigma = var[i]
	
# function to calculate the tangential velocity at each control point
def getTangentVelocity(p,fs,gamma):
	N = len(p)
	A = np.zeros((N,N),dtype=float)
	for i in range(N):
	   for j in range(N):
		if (i!=j):
		      A[i,j] = 0.5/pi*I(p[i].xc,p[i].yc,p[j],-sin(p[i].beta),cos(p[i].beta))
	B = fs.Uinf*np.sin([fs.alpha-pp.beta for pp in p])
	var = np.array([pp.sigma for pp in p])
	vt = np.dot(A,var)+B
	for i in range(N):
		p[i].vt = vt[i]
		
getTangentVelocity(panel,freestream,gamma)	# get tangential velocity

# function to calculate the pressure coefficient at each control point
def getPressureCoefficient(p,fs):
    for i in range(len(p)):
	p[i].Cp = 1-(p[i].vt/fs.Uinf)**2
	
getPressureCoefficient(panel,freestream)	# get pressure coefficient

# plotting the coefficient of pressure
valX,valY = 0.1,0.2
xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
Cpmin,Cpmax = min([p.Cp for p in panel]),max([p.Cp for p in panel])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = Cpmin-valY*(Cpmax-Cpmin),Cpmax+valY*(Cpmax-Cpmin)
plt.figure(figsize=(10,6))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('$C_p$',fontsize=16)
plt.plot([p.xc for p in panel if p.loc=='extrados'],\
		[p.Cp for p in panel if p.loc=='extrados'],\
		'ro-',linewidth=2)
plt.plot([p.xc for p in panel if p.loc=='intrados'],\
		[p.Cp for p in panel if p.loc=='intrados'],\
		'bo-',linewidth=1)
plt.legend(['extrados','intrados'],'best',prop={'size':14})
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.title('Number of panels : %d'%len(panel))
print min([p.Cp for p in panel])

# sum of all source/sink strengths
print '--> sum of source/sink strengths:',sum([p.sigma*p.length for p in panel])

# function to calculate the velocity field given a mesh grid
def getVelocityField(panel,freestream,gamma,X,Y):
    Nx,Ny = X.shape
    u,v = np.empty((Nx,Ny),dtype=float),np.empty((Nx,Ny),dtype=float)
    for i in range(Nx):
        for j in range(Ny):
            u[i,j] = freestream.Uinf*cos(freestream.alpha)\
				+ 0.5/pi*sum([p.sigma*I(X[i,j],Y[i,j],p,1,0) for p in panel])
            v[i,j] = freestream.Uinf*sin(freestream.alpha)\
				+ 0.5/pi*sum([p.sigma*I(X[i,j],Y[i,j],p,0,1) for p in panel])
    return u,v
    
# definition of the mesh grid
Nx,Ny = 20,20
valX,valY = 1.0,2.0
xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
ymin,ymax = min([p.ya for p in panel]),max([p.ya for p in panel])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
X,Y = np.meshgrid(np.linspace(xStart,xEnd,Nx),np.linspace(yStart,yEnd,Ny))

# get the velicity field on the mesh grid
u,v = getVelocityField(panel,freestream,gamma,X,Y)

# plotting the velocity field
size=12
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.streamplot(X,Y,u,v,density=1,linewidth=1,arrowsize=1,arrowstyle='->')
plt.fill([p.xa for p in panel],[p.ya for p in panel],'ko-',linewidth=2,zorder=2)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.title('Contour of velocity field');

# computing the pressure field
Cp = 1.0-(u**2+v**2)/freestream.Uinf**2

# plotting the pressure field
size=12
plt.figure(figsize=(1.1*size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
contf = plt.contourf(X,Y,Cp,levels=np.linspace(-2.0,1.0,100),extend='both')
cbar = plt.colorbar(contf)
cbar.set_label('$C_p$',fontsize=16)
cbar.set_ticks([-2.0,-1.0,0.0,1.0])
plt.fill([p.xc for p in panel],[p.yc for p in panel],'ko-',linewidth=2,zorder=2)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.title('Contour of pressure field');


plt.show()