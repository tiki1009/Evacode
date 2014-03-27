import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from math import *

# function to read the coordinates to store into 2 1-D arrays
def readGeomtry():
    inFile = open('/Users/Eva/Documents/ip/code/resources/naca0012.dat','r')
    x,y = [],[]
    for line in inFile:
        data = line.split()
        x.append(float(data[0]))
        y.append(float(data[1]))
    x,y = np.array(x),np.array(y)
    inFile.close()
    return x,y
    
