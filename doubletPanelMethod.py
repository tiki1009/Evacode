import numpy as np
from scipy import integrate
from math import *
import matplotlib.pyplot as plt

# read geometry
coords = np.loadtxt(fname='C:/Users/conans/Dropbox/selfshared/AERODYN/resources/naca0012.dat'
)
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

		self.sigma = 0. 		# strength
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
	def __inti__(self, Uinf, alpha):
		self.Uinf = Uinf
		self.alpha = alpha*pi/180
		self.uinf=Uinf*cos(alpha*pi/180)
		self.vinf=Uinf*sin(alpha*pi/180)
Uinf  = 1.0
alpha = 1.0
#freestream = Freestream(Uinf,alpha)


# plot air airfoil panels
plt.grid(True)
plt.plot([p.xb for p in panel],[p.yb for p in panel])
plt.plot([p.xa for p in panel],[p.ya for p in panel])
#plt.plot(xp,yp,'k-',linewidth=2)
plt.title('Airfoil')
plt.show()


# doublet potential 
def PHICD(mu,xci,yci,pj):
	 xR=cos(panel.beta)*(panel.xb-panel.xa)-sin(panel.beta)*(panel.yb-panel.ya)
	 x=cos(panel.beta)*(xci-panel.xa)-sin(panel.beta)*(yci-panel.ya)
	 y=sin(panel.beta)*(xci-panel.xa)+cos(panel.beta)*(yci-panel.ya)
	 phi = -(mu/(2*pi))*((atan2(y,x-xR))-(atan2(y,x)))
	 return phi

