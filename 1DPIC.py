from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from numpy import sin, cos, exp, pi
from numpy.fft import fft, ifft, fftfreq
from datetime import datetime
#
from pylab import rcParams
rcParams['figure.figsize'] = 20, 10

#PHYSICAL CONSTANTS

C=299792458. #m/s
EPSILON0=8.854187817e-12 #farads/m
ELECTRONCHARGE=1.602176565e-19 #kg
ELECTRONMASS=9.10938215e-31 #kg
PROTONMASS=1.672621777e-27 #kg
BOLTZMANNKB=1.380648e-23 #J/K
##BOLTZMANNKBEV=8.6173324e-5 #eV/K

#SIMULATION PARAMETERS
Dim=1
N=100
NG=10
Size=np.array([1.])
DT=0.0001
np.random.seed(1)
RUNITERS=500
TIMEITERS=np.arange(RUNITERS+1)*DT

PARTICLESPERSUPER=1
SUPERPARTICLECHARGE=PARTICLESPERSUPER*ELECTRONCHARGE #for superparticles
SUPERPROTONMASS=PARTICLESPERSUPER*PROTONMASS
SUPERELECTRONMASS=PARTICLESPERSUPER*ELECTRONMASS
SnapshotEveryXIterations=100

ESTMAXCHARGE=N*SUPERPARTICLECHARGE/Size[0]/2/ELECTRONCHARGE/4
ESTMAXFIELD=N*SUPERPARTICLECHARGE/Size[0]/NG/Size[0]
INITSHOWDENS=True

time=datetime.now().time()
from time import gmtime, strftime
##RUNTIME="%.2d%.2d%.2d" %(time[0], time[1], time[2])
RUNTIME=strftime("%H.%M.%S")

def RelativisticCorrectionGamma(v):
    speed=np.linalg.norm(v, axis=1)
    denominator=np.repeat(1/np.sqrt(1-(speed/C)**2), Dim)
    return denominator

def zerofield(X):
    return 0



class grid(object):
    def __init__(self):
        self.L = Size[0]
        self.X, self.dX = np.linspace(0,self.L, NG, retstep=True)
        self.efield=np.zeros(NG)
        self.density=np.zeros(NG)
        self.freq=np.zeros(NG)
        self.pot=np.zeros(NG)        
        self.init=True
    def densityplot(self, filename, view=False):
        plt.title("Grid")
        plt.subplot(2,1,1)
        plt.plot(self.X, self.density/ELECTRONCHARGE, label="Density")
        plt.ylabel("Charge density [elementary charge]")
        plt.ylim(-ESTMAXCHARGE, ESTMAXCHARGE)
        plt.subplot(2,1,2)
        plt.plot(self.X,self.efield, label="eField")
        plt.ylabel("Field")
##        plt.ylim(-ESTMAXFIELD, ESTMAXFIELD)
        plt.xlabel("Position on the grid")
        plt.savefig(filename)
        if view:
            plt.show()
        plt.clf()

        
    def update(self, list_of_species, externalfield=zerofield):
        self.density=np.zeros(NG)
        for species in list_of_species:
            density=np.zeros_like(self.density)
            gridindex=(species.position/self.L*NG).astype(int)
            if ((gridindex>NG).any() or (gridindex<0).any()):
                print "Error: particles out of bounds"
                print "Positive indices"
                print gridindex[gridindex>NG], species.position[gridindex>NG]
                print "Negative indices"
                print gridindex[gridindex<0], species.position[gridindex<0]
            nonmaxcondition=gridindex<NG-1
            nonmaxgridindex=gridindex[nonmaxcondition]
            maxcondition=(gridindex==NG-1)
            maxgridindex=gridindex[maxcondition]
            maxzeros=np.zeros_like(maxgridindex)
            try:
                density[nonmaxgridindex]+=species.position[nonmaxcondition]-nonmaxgridindex
                density[nonmaxgridindex+1]+=nonmaxgridindex+1-species.position[nonmaxcondition]
            except IndexError, err:
                print ("Error: %s.\n" %str(err))
                pass
            density[maxgridindex]+=species.position[maxcondition]-maxgridindex
            density[maxzeros]+=maxgridindex+1-species.position[maxcondition]
            density *= species.charge
            self.density+=density

        self.density*=1./self.dX

        if self.init:
            self.densityplot(RUNTIME + "Initdensity.png", INITSHOWDENS)
            self.init=False

        #FOURIER TRANSFORM
        self.freq=self.L*fftfreq(NG, self.dX)
        self.freq[0]=0.01
        self.pot = np.real(ifft(fft(self.density)[0:NG]/self.freq[0:NG]**2/4./pi**2/EPSILON0))
        self.efield = -np.gradient(self.pot)
        self.efield+=externalfield(self.X)
class species(object):
    def __init__(self, mass, charge, position, velocity, number, name, color):
        self.name=name
        self.color=color
        self.number=number
        self.mass=mass
        self.charge=charge
        self.position=position
        self.velocity=velocity

        self.trajectories=np.copy(self.position)
        self.velocities=np.copy(self.velocity)
        self.temperatures=np.array(self.temperature())
    def move(self, dt):
        self.position += self.velocity*dt
        #wymuszanie warunków brzegowych
        for i in range(Dim):    #to na pewno można zrobić bez fora
            j=0
            while (self.position[:,i]>Size[i]).any() or (self.position[:,i]<0).any():
##                if j>1:
##                    print j
                self.position[:,i][self.position[:,i]>Size[i]]-=Size[i]
                self.position[:,i][self.position[:,i]<0]+=Size[i]
                j+=1
    def accelerate(self, grid, dt):
        gridindex=(self.position/grid.L*NG).astype(int)

        nonmaxcondition=gridindex<NG-1
        nonmaxgridindex = gridindex[nonmaxcondition]
        maxcondition=gridindex==NG-1
        maxgridindex=gridindex[maxcondition]
        maxzeros=np.zeros_like(maxgridindex)
        EField=np.zeros_like(self.position)
        EField[nonmaxcondition] = (grid.X[nonmaxgridindex+1]-self.position[nonmaxcondition])/grid.dX*grid.efield[nonmaxgridindex]+(self.position[nonmaxcondition]-grid.X[nonmaxgridindex])/grid.dX*grid.efield[nonmaxgridindex+1]
##        EField[maxcondition] = -((grid.X[maxgridindex-1]-self.position[maxcondition])/grid.dX*grid.efield[maxgridindex]+(self.position[maxcondition]-grid.X[maxgridindex])/grid.dX*grid.efield[maxgridindex-1]) #stara wersja bez periodycznych warunkow brzegowych z interpolacja w ujemna strone
        EField[maxcondition] = (grid.X[maxzeros]-self.position[maxcondition])/grid.dX*grid.efield[maxgridindex]+(maxgridindex+1-self.position[maxcondition])/grid.dX*grid.efield[maxzeros]

        acceleration=np.zeros((N,Dim))
        acceleration[:,0]=self.charge*EField[:,0]/self.mass/(RelativisticCorrectionGamma(self.velocity))

        self.velocity+=acceleration*dt
    def temperature(self):
        return 0.5*self.mass*np.sqrt(np.sum(self.velocity**2)/self.number)
    def info(self):
        print "N=",self.number, self.name
        print "m=",self.mass, "q=", self.charge
        print "position\n",self.position
        print "velocity\n",self.velocity
    def step(self):
        self.accelerate(Grid, DT)
        self.move(DT)
        self.trajectories=np.hstack((self.trajectories,self.position))
        self.velocities=np.hstack((self.velocities,self.velocity))
        self.temperatures=np.hstack((self.temperatures,self.temperature()))
def PlotAllTrajectories(ListOfSpecies):
    plt.title("Run history")
    for index, i in enumerate(ListOfSpecies):
        plt.subplot(len(ListOfSpecies)+2,1,index+1)
        plt.plot(TIMEITERS, i.velocities.T, color=i.color, alpha = 0.8)
        plt.ylabel('X velocity for ' + i.name + '[m/s]')

    plt.subplot(len(ListOfSpecies)+2,1,len(ListOfSpecies)+1)
    for i in ListOfSpecies:
        plt.plot(TIMEITERS, i.velocities.T, color=i.color, alpha = 0.5)
    plt.ylabel('X velocity [m/s]')

    plt.subplot(len(ListOfSpecies)+2,1,len(ListOfSpecies)+2)
    for i in ListOfSpecies:
        plt.plot(TIMEITERS, i.temperatures/BOLTZMANNKB, color=i.color)
    plt.ylabel("Temperature [K]")
    plt.xlabel("Time [s]")
    plt.savefig(RUNTIME + "Runhistory.png")
    plt.show()

Grid=grid()
#hot start
electrons=species(SUPERELECTRONMASS, -SUPERPARTICLECHARGE, np.random.random((N,Dim))*Size, np.random.random((N,Dim))*Size-Size/2, N, "electrons", "y")
protons=species(SUPERPROTONMASS, SUPERPARTICLECHARGE, np.random.random((N,Dim))*Size, np.random.random((N,Dim))*Size-Size/2, N, "protons", "b")
#cold start
##electrons=species(1., -SUPERPARTICLECHARGE, np.random.random((N,Dim))*Size, np.zeros((N,Dim)), N, "electrons", "y")
##protons=species(1., SUPERPARTICLECHARGE, np.random.random((N,Dim))*Size, np.zeros((N,Dim)), N, "protons", "b")

Species=[protons,electrons]
Grid.update(Species)
for i in Species:
    i.accelerate(Grid, -0.5*DT)
for iterat in range(RUNITERS):
    for i in Species:
        i.step()
##        print i.name, i.temperature()
    if iterat%SnapshotEveryXIterations==0:
        Grid.densityplot(RUNTIME + "DensityIter%d.png" %(iterat))
    Grid.update(Species)
PlotAllTrajectories(Species)
Grid.densityplot(RUNTIME + "Finaldensity.png", True)
