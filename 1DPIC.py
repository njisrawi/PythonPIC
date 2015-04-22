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
Dim=3
N=128
NG=32
Size=np.array([1e-3, 1e-4, 1e-4])
DT=1e-6
np.random.seed(1)
RUNITERS=5000
TIMEITERS=np.arange(RUNITERS+1)*DT
SnapshotEveryXIterations=1000
##MagneticFieldStrength=np.array([0,0,1.]) #teslas
##MagneticFieldStrength=1e-6

PARTICLESPERSUPER=100
SUPERPARTICLECHARGE=PARTICLESPERSUPER*ELECTRONCHARGE #for superparticles
SUPERPROTONMASS=PARTICLESPERSUPER*PROTONMASS
SUPERELECTRONMASS=PARTICLESPERSUPER*ELECTRONMASS

ESTMAXCHARGE=N/Size[0]*NG
ESTMAXFIELD=N*SUPERPARTICLECHARGE/Size[0]/NG/Size[0]
INITSHOWDENS=False

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
    def densityplot(self, listofspecies, filename, view=False):
##        plt.clf()
##        fig=plt.fig()
##        fig.title("Grid")
##        fig.subplot(2,1,1)
##        fig.plot(self.X, self.density/ELECTRONCHARGE, label="Density")
##        for i in listofspecies:
##            fig.scatter(i.position.T, np.zeros_like(i.position.T), c=i.color, alpha=0.3)
##        fig.ylabel("Charge density [elementary charge]")
####        plt.ylim(-ESTMAXCHARGE, ESTMAXCHARGE)
##        fig.subplot(2,1,2)
##        fig.plot(self.X,self.efield, label="eField")
##        fig.ylabel("Field")
####        plt.ylim(-ESTMAXFIELD, ESTMAXFIELD)
##        fig.xlabel("Position on the grid")
##        fig.savefig(filename)
##        if view:
##            fig.show()
##        return fig


        plt.title("Grid")
        plt.subplot(2,1,1)
        plt.plot(self.X, self.density/ELECTRONCHARGE, label="Density")
        for i in listofspecies:
            plt.scatter(i.position.T, np.zeros_like(i.position.T), c=i.color, alpha=0.3)
        plt.ylabel("Charge density [elementary charge]")
##        plt.ylim(-ESTMAXCHARGE, ESTMAXCHARGE)
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


##            gridindex=(species.position/self.L*NG).astype(int)
##            if ((gridindex>NG).any() or (gridindex<0).any()):
##                print "Error: particles out of bounds"
##                print "Positive indices"
##                print gridindex[gridindex>NG], species.position[gridindex>NG]
##                print "Negative indices"
##                print gridindex[gridindex<0], species.position[gridindex<0]
##            nonmaxcondition=gridindex<NG-1
##            nonmaxgridindex=gridindex[nonmaxcondition]
##            maxcondition=(gridindex==NG-1)
##            maxgridindex=gridindex[maxcondition]
##            maxzeros=np.zeros_like(maxgridindex)
##            try:
##                density[nonmaxgridindex]+=species.position[nonmaxcondition]-nonmaxgridindex
##                density[nonmaxgridindex+1]+=nonmaxgridindex+1-species.position[nonmaxcondition]
##            except IndexError, err:
##                print ("Error: %s.\n" %str(err))
##                pass
####            plt.plot(self.X, density, label=(species.name + "before"))
##            density[maxgridindex]+=species.position[maxcondition]-maxgridindex
####            density[maxzeros]+=maxgridindex+1-species.position[maxcondition]
####            density[maxgridindex]=0
####            density[maxzeros]=0
####            plt.plot(self.X, density, label=(species.name + "after"))
##            density *= species.charge
##            self.density+=density
####            plt.scatter(species.position, np.zeros_like(species.position), label=species.name)
####            plt.legend()
####            plt.show()
##        self.density*=1./self.dX
            for ppos in species.position[:,0]:
                index=int(ppos/self.dX)
                indexpos=self.X[index]
                w1=1-(ppos-indexpos)/self.dX
                w2=(ppos-indexpos)/self.dX
                if index<NG-1:
                    density[index]+=w1*species.charge
                    density[index+1]+=w2*species.charge
                else:
                    density[index]+=w1*species.charge
                    density[0]+=w2*species.charge              
            self.density+=density
        if (np.abs(self.density)>10000*ELECTRONCHARGE).any():
            print "DUN GOOFED"
            plt.plot(self.X, density)
            for species in list_of_species:
                plt.scatter(species.position, np.zeros_like(species.position), label=species.name, c=species.color, alpha=0.3)
            plt.legend()
            plt.show()
##            quit()
##        quit()
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
##                if j>50:
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
        acceleration[:,0]=EField[:,0]
        #MAGFIELD
##        acceleration[:,0]+=MagneticFieldStrength*self.velocity[:,1]
##        acceleration[:,1]-=MagneticFieldStrength*self.velocity[:,0]
        acceleration*=self.charge/self.mass
        self.velocity+=acceleration*dt
    def temperature(self):
##for i in range
        return self.mass*np.mean(np.square(self.velocity))/BOLTZMANNKB
    def info(self):
        print "N=",self.number, self.name
        print "m=",self.mass, "q=", self.charge
        print "position\n",self.position
        print "velocity\n",self.velocity
    def step(self):
        self.accelerate(Grid, DT)
        self.move(DT)
        self.trajectories=np.dstack((self.trajectories,self.position))
        self.velocities=np.dstack((self.velocities,self.velocity))
        self.temperatures=np.hstack((self.temperatures,self.temperature()))


#=========GENERAL DIAGNOSTICS===========

def PlotAllTrajectories(ListOfSpecies):
    plt.title("Run history")
    for index, i in enumerate(ListOfSpecies):
        plt.subplot(len(ListOfSpecies)+2,1,index+1)
        plt.plot(TIMEITERS, i.velocities[:,0,:].T, color=i.color, alpha = 0.8)
        plt.ylabel('X velocity for ' + i.name + '[m/s]')

    plt.subplot(len(ListOfSpecies)+2,1,len(ListOfSpecies)+1)
    for i in ListOfSpecies:
        plt.plot(TIMEITERS, i.velocities[:,2,:].T, color=i.color, alpha = 0.5)
    plt.ylabel('Y velocity [m/s]')

    plt.subplot(len(ListOfSpecies)+2,1,len(ListOfSpecies)+2)
    for i in ListOfSpecies:
        plt.plot(TIMEITERS, i.temperatures, color=i.color)
    plt.ylabel("Temperature [K]")
    plt.xlabel("Time [s]")
    plt.savefig(RUNTIME + "Runhistory.png")
    plt.show()

   
Grid=grid()
protondistribution=np.zeros((N,Dim))
electrondistribution=np.zeros((N,Dim))

###cosine start
##protondistribution[:,0]=np.linspace(Grid.L/10000, Grid.L*(1-1/10000), N)
##electrondistribution[:,0]=protondistribution[:,0]+0.1*Grid.L*np.sin(2*pi*protondistribution[:,0]/Grid.L)
##protons=species(SUPERPROTONMASS, SUPERPARTICLECHARGE, protondistribution, np.zeros((N,Dim)), N, "protons", "b")
##electrons=species(SUPERELECTRONMASS, -SUPERPARTICLECHARGE, electrondistribution, np.random.random((N,Dim))*Size-Size/2, N, "electrons", "y")

#hot start
##electrons=species(SUPERELECTRONMASS, -SUPERPARTICLECHARGE, np.random.random((N,Dim))*Size, np.random.random((N,Dim))*Size-Size/2, N, "electrons", "y")
##protons=species(SUPERPROTONMASS, SUPERPARTICLECHARGE, np.random.random((N,Dim))*Size, np.random.random((N,Dim))*Size-Size/2, N, "protons", "b")

#maxwellstart
INITTEMPERATURE=1
protons=species(SUPERPROTONMASS, SUPERPARTICLECHARGE, np.random.random((N,Dim))*Size, np.random.normal(0,np.sqrt(BOLTZMANNKB*INITTEMPERATURE/SUPERPROTONMASS), (N,Dim)), N, "protons", "b")
electrons=species(SUPERELECTRONMASS, -SUPERPARTICLECHARGE, np.random.random((N,Dim))*Size, np.random.normal(0,np.sqrt(BOLTZMANNKB*INITTEMPERATURE/SUPERELECTRONMASS), (N,Dim)), N, "electrons", "y")





#cold start
##electrons=species(1., -SUPERPARTICLECHARGE, np.random.random((N,Dim))*Size, np.zeros((N,Dim)), N, "electrons", "y")
##protons=species(1., SUPERPARTICLECHARGE, np.random.random((N,Dim))*Size, np.zeros((N,Dim)), N, "protons", "b")

Species=[protons,electrons]
Grid.update(Species)
for i in Species:
    i.accelerate(Grid, -0.5*DT)
for iterat in range(RUNITERS):
    print("Iteration", iterat)
    for i in Species:
        i.step()
##        print i.name, i.temperature()
    if iterat%SnapshotEveryXIterations==0:
        Grid.densityplot(Species, RUNTIME + "DensityIter%d.png" %(iterat))
    Grid.update(Species)
PlotAllTrajectories(Species)
Grid.densityplot(Species, RUNTIME + "Finaldensity.png", True)


##for i in range
