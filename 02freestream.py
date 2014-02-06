import numpy as np
import matplotlib.pyplot as plt
from math import *

N = 200                          # Number of points in each direction
xStart,xEnd = -4.0,4.0           # x-direction boundaries
yStart,yEnd = -2.0,2.0           # y-direction boundaries
x = np.linspace(xStart,xEnd,N)   # x 1D-array
y = np.linspace(yStart,yEnd,N)   # y 1D-array
X,Y = np.meshgrid(x.y)           # Generation of the mesh grid

test