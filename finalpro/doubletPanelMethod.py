import numpy as np
from scipy import integrate
from math import *
import matplotlib.pyplot as plt

# read geometry
coords = np.loadtxt(fname='/Users/Eva/Documents/python/code/finalpro/naca0012.dat')
xp,yp = coords[:,0],coords[:,1]

# classify each panel info
class Panel:
        def __init__(self,xa,ya,xb,yb):
		self.xa = xa
		self.xb = xb
		self.ya = ya
		self.yb = yb
		self.xc = (xa+xb)/2
		self.yc = (ya+yb)/2				# panel points 
		self.length = sqrt((xb-xa)**2+(yb-ya)**2)	# panel length
		self.beta = np.arctan((yb-ya)/(xb-xa))		# panel orientation

		self.mu = 0. 		# strength
		self.vt = 0.  			# tangential velocity
		self.cp = 0. 			# cp

# define panels
N=len(xp)
panel = np.empty(N,dtype = object)
for i in range(N):
        if(i==N-1):
                panel[i] = Panel(xp[-1],yp[-1],xp[i],yp[i])
        else:
	       panel[i] = Panel(xp[i],yp[i],xp[i+1],yp[i+1])


# freestream
class Freestream:
	def __init__(self, Uinf, alpha):
		self.Uinf = Uinf
		self.alpha = alpha*pi/180
		self.uinf=Uinf*cos(alpha*pi/180)
		self.vinf=Uinf*sin(alpha*pi/180)
		
Uinf  = 1.0
alpha = 5.0
freestream = Freestream(Uinf,alpha)


# plot air airfoil panels
plt.grid(True)
plt.plot(xp,yp,'k-',linewidth=2)
plt.title('Airfoil')
plt.show()


# doublet potential 
def PHICD(mu,xci,yci,pj):
	 xR=cos(panel.beta)*(panel.xb-panel.xa)-sin(panel.beta)*(panel.yb-panel.ya)
	 x=cos(panel.beta)*(xci-panel.xa)-sin(panel.beta)*(yci-panel.ya)
	 y=sin(panel.beta)*(xci-panel.xa)+cos(panel.beta)*(yci-panel.ya)
	 phi = -(mu/(2*pi))*((atan2(y,x-xR))-(atan2(y,x)))
	 return phi

# create LHS \N-1
A = np.empty((N,N),dtype=float)
for i in range(N-1):        
    for j in range(N):
        if (i!=j):
            A[i,j]=PHICD(1,panel[i].xc,panel[i].yc,panel[j])
        else:
            A[i,j]=0.5
            
# create RHS 
b=np.empty((N,1),dtype=float)
for i in range(N):
    b[i]=-freestream.uinf*panel[i].xc-freestream.vinf*panel[i].yc

# kutta condition rearrangement
A[N-1][0]=1
A[N-1][N-2]=-1
A[N-1][N-1]=1

# solve
mu=np.linalg.solve(A,b)

# store strength
for i in range(len(panel)):
    panel[i].mu=mu[i]
    
# solve for tangential velocity
B = np.zeros((N,N),dtype=float)
for i in range(N):
    for j in range(N):
        if (i!=j):
            B[i,j] = PHICD(1,panel[i].xc,panel[i].yc,panel[j])/2/pi
C = freestream.Uinf*np.sin([freestream.alpha-panel.beta])
Vt=mu+C
for i in range(N):
    panel[i].vt = Vt[i]
    
    
# cp
for i in range(len(panel)):
    panel[i].cp = 1-(panel[i].vt/freestream.Uinf)**2
    

# plot
figure()
plt.grid(True)
Plt.xlabel('x')
plt.ylabel('Cp')
plt.plot(panel[i].xc,panel[i].cp)
plt.show()

# velocity field
u = np.empty((N,N),dtype = float)
v = np.empty((N,N),dtype = float)
for i in range(N):
    for j in range(N):
        u[i,j] = freestream.uinf -panel.mu*0.5/pi*(PHICD(1,panel[i].xc,panel[i].yc,panel[j]))
        v[i,j] = freestream.vinf -panel.mu*0.5/pi*(PHICD(1,panel[i].xc,panel[i].yc,panel[j]))
        
# mesh grid
Nx,Ny = 20,20
valX,valY = 1.0,2.0
xmin,xmax = min(panel.xa),max(panel.xa)
ymin,ymax = min(panel.ya),max(panel.ya)

xStart,xEnd = xmin-valX*(xmax-xmin) ,xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)


X,Y = np.meshgrid(np.linspace(xStart,xEnd,Nx),np.linspace(yStart,yEnd,Ny))


# plot
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.streamplot(X,Y,u,v,density=1,linewidth=1,arrowsize=1,arrowstyle='->')
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)


