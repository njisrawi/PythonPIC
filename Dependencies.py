from __future__ import division
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from science import *

L=1     #system length
NSP=2   #number of species
DT = 1e-2 #time step
NT = 10**4 #number of iterations
NG = 2**8 #number of grid points

Epsilon = 1 #1/epsilon0
A1 = 0 #compensation
A2 = 0 #smoothing factor
PlotIterStep=100

N=np.array([10,10.]) #number of particles in each species
WP=np.array([1,1.]) #plasma frequency - positive
WC=np.array([1,-1.]) #cyclotron frequency - signed
QM=np.array([1, -1/1000.])
