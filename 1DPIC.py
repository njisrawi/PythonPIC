from Dependencies import *
N=100000
DeltaT=1e-5

def LinearInterpolation(x, Dx):
    result=1-np.abs(x)/Dx
    result[np.abs(x)>Dx]=0.
    return result

##def EBoundaryCondition(voltage, f, t):
##    return 

#grid setup
x=np.linspace(0,64,6400)
y=np.linspace(0,1,100)
z=np.linspace(0,1,100)
#hstack?

BVector=np.array([0,0,1])

#setup for one kind of particles; could be copied for ions and different kinds of particles could be put in list
r=np.zeros((N,3))
v=np.zeros((N,3))
a=np.zeros((N,3))
m=np.ones(N)
q=-np.ones(N)

#initial conditions
##r[:,0]=np.random.uniform(min(x), max(x), N)
##r[:,1:3]=np.random.uniform(0,1,(N,2))
#print (r)

#histogram showing
##plt.hist(r[:,2], bins=64, normed=True)
##plt.show()

#initial push back in preparation for leapfrog method
##v=v-a*DeltaT/2

#finish leapfrog integrator as 
