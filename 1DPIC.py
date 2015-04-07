from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from numpy import sin, cos, exp, pi
from numpy.fft import fft, ifft, fftfreq

from ScienceConstants import *
Dim=1
N=100
NG=10
Size=np.array([30])
DT=0.1
np.random.seed(1)
C=3e8
eps_0=1

def RelativisticCorrectionGamma(v):
    speed=np.linalg.norm(v, axis=1)
    denominator=np.repeat(1/np.sqrt(1-(speed/C)**2), Dim)
    return denominator

class grid(object):
    def __init__(self):
        self.L = Size[0]
        self.X = np.linspace(0,self.L, NG+1)
        self.dX=self.L/NG
        self.efield=np.zeros(NG+1)
        self.density=np.zeros(NG+1)
        self.freq=np.zeros(NG+1)
        self.pot=np.zeros(NG+1)
        self.init=True
    def update(self, list_of_species):
        self.density=np.zeros(NG+1)
        for species in list_of_species:
##            try:
            gridindex=(species.position[:,0]/self.L*NG).astype(int)
            interior = (gridindex+1)<NG
            print (gridindex+1)<NG
            
##            self.density[:NG][gridindex]+=(self.X[(gridindex+1)[interior]]-species.position[:,0])/self.dX
##            self.density[:NG][gridindex]+=(species.position[:,0]-self.X[gridindex[interior]])/self.dX
##            except:
##                #print species.position
##                print gridindex
##                print self.X
##            self.density[NG-1] += (self.X[NG] - species.position[:,0])/self.dX
##            self.density[0] += (species.position[:,0][-self.X[NG-1])/self.dX
            self.density *= species.charge
        self.density[NG]=self.density[0] #a może bez tego?

        if self.init:
            plt.plot(self.X, self.density)
            plt.ylabel("Density")
            plt.title("Density")
            plt.show()
            self.init=False


        #FOURIER TRANSFORM
        self.freq=self.L*fftfreq(NG, self.dX)
        self.freq[0]=0.01
        self.pot = np.real(ifft(fft(self.density)[0:NG]/self.freq[0:NG]**2/4./pi**2/eps_0))
        print self.pot
        self.efield = -np.gradient(self.pot)
class species(object):
    def __init__(self, mass, charge, position, velocity, number, name):
        self.name=name
        self.number=number
        self.mass=mass
        self.charge=charge
        self.position=position
        self.velocity=velocity
    
    def move(self, dt):
        self.position += self.velocity*dt
        for i in range(Dim):    #to na pewno można zrobić bez fora
            #out of positive range
            self.position[:,i][self.position[:,i]>Size[i]]-=Size[i]
            #out of negative range
            self.position[:,i][self.position[:,i]<0]+=Size[i]
    def accelerate(self, grid, dt):
        gridindex=(self.position/grid.L*NG).astype(int)
        EField=((grid.X[gridindex+1]-self.position)/grid.dX*grid.efield[gridindex] + (self.position-grid.X[gridindex])/grid.dX*grid.efield[gridindex+1])[:,0]
        acceleration=np.zeros((N,Dim))
        acceleration[:,0]=self.charge*EField/self.mass/(RelativisticCorrectionGamma(self.velocity))
        self.velocity+=acceleration*dt

        
        acceleration=self.charge
    def info(self):
        print "N=",self.number, self.name
        print "m=",self.mass, "q=", self.charge
        print "position\n",self.position
        print "velocity\n",self.velocity
    def step(self):
        self.accelerate(Grid, DT)
        self.move(DT)
Grid=grid()
electrons=species(1, -1, np.random.random((N,Dim))*Size, np.random.random((N,Dim))*Size/2, N, "electrons")
protons=species(1, 1, np.random.random((N,Dim))*Size, np.random.random((N,Dim))*Size/2, N, "protons")
Species=[protons,electrons]
Grid.update(Species)
for i in Species:
    i.accelerate(Grid, -0.5*DT)
for iterat in range(10):
    for i in Species:
        i.step()
    Grid.update(Species)


##def Time(i):
##    return TInitial+(TFinal-TInitial)/TFinal*i*DT
##
###Fields - continuous cases
##
####def EField(r):
####    return 0
####
####def BField(r):
####    result=np.zeros_like(r)
####    result[:,2] = r[:,1]
####    if DEBUG:
####        print "field", result
####    return result
##def BField(r):
##    return np.array([0,0,1])
####def Force(r, v, EField, BField): #Lorentz
####    f=q*(EField(r)+np.cross(v, BField(r)))
####    if DEBUG:
####        print "force", f
####    return f
##
###Particle movement
##
##def Force(r, v, q, EFieldArray, BFieldFunc):
####    print "vstack"
####    print np.vstack((q,q,q))
##    f1=np.vstack((q,q,q)).T*np.cross(v, BField(r))
##    f2=EFieldArray*q
##    f=f1
##    f[:,0]+=f2
##    for i in range(len(f1)):
##        print f1[i], f2[i]
####    print f
##    return f
##def InitialPush(r, v, q, m, EFieldArray, BFieldFunc):
##    return v+Force(r, v, q, EFieldArray, BFieldFunc)*DT/2./m
##
##def LeapfrogStep(r, v, q, m, EFieldArray, BFieldFunc):
##    return r + v*DT, v + Force(r, v, q, EFieldArray, BFieldFunc)*DT/m
##
###Diagnostics
##
##def KineticEnergy(v, m):
##    return np.sum(v*v*m)/2.
##
###Plotting
##
##def PlotTrajectory2D(rdata):
##    for j in range(len(rdata[0,0,:])):
##        plt.plot(rdata[j, 0, :], rdata[j, 1, :], "o-")
##    plt.xlabel("x")
##    plt.ylabel("y")
##    plt.grid()
##    plt.title("2D trajectory")
##    plt.show()
##    for j in range(len(rdata[0,0,:])):
##        plt.plot(rdata[j, 0, :], rdata[j, 2, :], "o-")
##    plt.xlabel("x")
##    plt.ylabel("z")
##    plt.grid()
##    plt.title("2D trajectory")
##    plt.show()
##    for j in range(len(rdata[0,0,:])):
##        plt.plot(rdata[j, 1, :], rdata[j, 2, :], "o-")
##    plt.xlabel("y")
##    plt.ylabel("z")
##    plt.grid()
##    plt.title("2D trajectory")
##    plt.show()
##
##def PlotTrajectory3D(rdata):
##    fig=plt.figure()
##    ax=fig.add_subplot(111,projection='3d')
##    for j in range(len(rdata[:,0,0])):
##        plt.plot(rdata[j, 0, :], rdata[j, 1, :], rdata[j,2,:], "o-")#Time(IPoints))
##    plt.grid()
##    plt.title("3D trajectory")
##    ax.set_xlabel('$x$', fontsize=16)
##    ax.set_ylabel('$y$', fontsize=16)
##    ax.set_zlabel('$z$', fontsize=16)
##    plt.show()    
##
##def EPlot(edata):
##    plt.plot(Time(IPoints),edata[0,0,:])
##    plt.xlabel("Time")
##    plt.ylabel("Energy")
##    plt.show()
##
##def EnforceBoundaryConditions(r):
##    result = r
##    result[:,0]=result[:,0]%L
##    result[:,1:]=result[:,1:]%W
##    return result
##L=60.     #system length
##W=1
##H=W
##NSP=2   #number of species
##DT = 1e-2 #time step
##NT = 100 #number of iterations
##NG = 31 #number of grid points
##gridspacing=L/31.
##N=100
##
####IterFinal=int(TFinal/DT)
####SamplingRate=IterFinal/1
####IterStep=IterFinal/SamplingRate
##IPoints = np.arange(0, NT, 1)
##
##Epsilon = 1 #1/epsilon0
##
##r=np.zeros((N,3))
##r[:,0]=np.linspace(L/10.,9*L/10.,N)
##r[:,1:]=W/2.
##print r
##v=np.zeros((N,3))
##m=np.ones((N,Dim))
##q=np.ones(N)
##q[N/2:]*=-1
##q*=1/10.
##gridpoints=np.linspace(0,L,NG)
##gridspacing=gridpoints[1]-gridpoints[0]
##
##rdata=r.copy()
##vdata=v.copy()
##edata=np.array(KineticEnergy(v, m))
##
##chargedensity=LinearInterpolateToGrid(r[:,0], q, gridpoints)
##phi=GaussSeidel1DSolver(chargedensity)
##field=EFieldSolver(phi)
##PlotFields(gridpoints, r[:,0], q, chargedensity, phi, field)
##particlefields=EFieldInterpolator(r[:,0], gridpoints, field)
##forces=Force(r,v,q,particlefields,BField)
##v=InitialPush(r,v,q,m,particlefields,BField)
##
##i=0
##while i<NT:
##    LastPercent=0
##    chargedensity=LinearInterpolateToGrid(r[:,0], q, gridpoints)
##    phi=GaussSeidel1DSolver(chargedensity)
##    field=EFieldSolver(phi)
##    particlefields=EFieldInterpolator(r[:,0], gridpoints, field)
##    forces=Force(r,v,q,particlefields,BField)
##    LastPercent=0
##    r, v = LeapfrogStep(r, v, q, m,particlefields,BField)
##    r = EnforceBoundaryConditions(r)
##    if i>0:
##        print i
##        rdata=np.dstack((rdata,r))
##        vdata=np.dstack((vdata,v))
##        edata=np.dstack((edata,KineticEnergy(v, m)))
##        if i%20==0:
##            PlotFields(gridpoints, r[:,0], q, chargedensity, phi, field)
##            PlotTrajectory3D(rdata)
##    i+=1
##print "Simulation finished"
##
##PlotTrajectory2D(rdata, N)
##PlotTrajectory3D(rdata, N)
