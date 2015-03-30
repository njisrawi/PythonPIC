from Dependencies import *
print "Initialization..."
L=1.     #system length
NSP=2   #number of species
DT = 1e-4 #time step
NT = 30 #number of iterations
NG = 31 #number of grid points
gridspacing=L/31.
N=100

##IterFinal=int(TFinal/DT)
##SamplingRate=IterFinal/1
##IterStep=IterFinal/SamplingRate
IPoints = np.arange(0, NT, 1)

Epsilon = 1 #1/epsilon0
A1 = 0 #compensation
A2 = 0 #smoothing factor
PlotIterStep=100

##N=np.array([10,10.]) #number of particles in each species
WP=np.array([1,1.]) #plasma frequency - positive
WC=np.array([1,-1.]) #cyclotron frequency - signed
QM=np.array([1, -1/1000.])
#

if __name__=="__main__":
    print "boo"
