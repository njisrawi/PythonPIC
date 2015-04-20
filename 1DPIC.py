from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from numpy import sin, cos, exp, pi
from numpy.fft import fft, ifft, fftfreq

from ScienceConstants import *
Dim=1
N=100
NG=10
Size=np.array([30.])
DT=0.01
##np.random.seed(1)
C=3e8
eps_0=1
RUNITERS=3000
ELEMCHARGE=10.

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
        plt.plot(self.X, self.density)
        plt.ylabel("Charge density")
        plt.xlabel("Position on the grid")
        plt.title("Density plot")
        plt.ylim(-15,15)
        plt.savefig(filename)
        if view:
            plt.show()
        plt.clf()
        
    def update(self, list_of_species, externalfield=zerofield):
        self.density=np.zeros(NG)
        for species in list_of_species:
            density=np.zeros_like(self.density)
            gridindex=(species.position/self.L*NG).astype(int)
            nonmaxcondition=gridindex<NG-1
            nonmaxgridindex=gridindex[nonmaxcondition]
            maxcondition=(gridindex==NG-1)
            maxgridindex=gridindex[maxcondition]
            maxzeros=np.zeros_like(maxgridindex)
            try:
                density[nonmaxgridindex]+=species.position[nonmaxcondition]-nonmaxgridindex
                density[nonmaxgridindex+1]+=nonmaxgridindex+1-species.position[nonmaxcondition]
            except IndexError, err:
                print nonmaxgridindex
                print ("Error: %s.\n" %str(err))
                pass
            density[maxgridindex]+=species.position[maxcondition]-maxgridindex
            density[maxzeros]+=maxgridindex+1-species.position[maxcondition]
            density *= species.charge
            self.density+=density

        self.density*=1./self.dX

        if self.init:
            self.densityplot("Initdensity.png", True)
            self.init=False

        #FOURIER TRANSFORM
        self.freq=self.L*fftfreq(NG, self.dX)
        self.freq[0]=0.01
        self.pot = np.real(ifft(fft(self.density)[0:NG]/self.freq[0:NG]**2/4./pi**2/eps_0))
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
    def move(self, dt):
        self.position += self.velocity*dt
        #wymuszanie warunków brzegowych
        for i in range(Dim):    #to na pewno można zrobić bez fora
            self.position[:,i][self.position[:,i]>Size[i]]-=Size[i]
            self.position[:,i][self.position[:,i]<0]+=Size[i]
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
def PlotAllTrajectories(ListOfSpecies):
    plt.title("Run history")
    for index, i in enumerate(ListOfSpecies):
        plt.subplot(len(ListOfSpecies)+1,1,index+1)
        plt.plot(i.velocities.T, color=i.color)
        plt.ylabel('X velocity for ' + i.name)

    plt.subplot(len(ListOfSpecies)+1,1,len(ListOfSpecies)+1)
    for i in ListOfSpecies:
        plt.plot(i.velocities.T, color=i.color)
    plt.ylabel('X velocity')
    plt.xlabel("Iterations")
    plt.savefig("Runhistory.png")
    plt.show()

Grid=grid()
electrons=species(1., -ELEMCHARGE, np.random.random((N,Dim))*Size, np.random.random((N,Dim))*Size/2-Size/2, N, "electrons", "y")
protons=species(1., ELEMCHARGE, np.random.random((N,Dim))*Size, np.random.random((N,Dim))*Size/2-Size/2, N, "protons", "b")
Species=[protons,electrons]
Grid.update(Species)
for i in Species:
    i.accelerate(Grid, -0.5*DT)
for iterat in range(RUNITERS):
    for i in Species:
        i.step()
        if iterat%10==0:
            Grid.densityplot("DensityIter%d.png" %(iterat))
    Grid.update(Species)
PlotAllTrajectories(Species)
Grid.densityplot("Finaldensity.png", True)
