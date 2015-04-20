from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from numpy import sin, cos, exp, pi
from numpy.fft import fft, ifft, fftfreq

from ScienceConstants import *
Dim=1
N=1000
NG=10
Size=np.array([30.])
DT=0.1
np.random.seed(1)
C=3e8
eps_0=1
RUNITERS=3000
ELEMCHARGE=10.

def RelativisticCorrectionGamma(v):
    speed=np.linalg.norm(v, axis=1)
    denominator=np.repeat(1/np.sqrt(1-(speed/C)**2), Dim)
    return denominator

class grid(object):
    def __init__(self):
        self.L = Size[0]
        self.X, self.dX = np.linspace(0,self.L, NG, retstep=True)
        self.efield=np.zeros(NG)
        self.density=np.zeros(NG)
        self.freq=np.zeros(NG)
        self.pot=np.zeros(NG)
        self.init=True
    def update(self, list_of_species):
        self.density=np.zeros(NG)
        for species in list_of_species:
##            try:
            density=np.zeros_like(self.density)
            gridindex=(species.position/self.L*NG).astype(int)
            nonmaxcondition=gridindex<NG-1
            nonmaxgridindex=gridindex[nonmaxcondition]
            maxcondition=(gridindex==NG-1)
            maxgridindex=gridindex[maxcondition]
            maxzeros=np.zeros_like(maxgridindex)
##            print maxcondition
##            print maxgridindex
            density[nonmaxgridindex]+=species.position[nonmaxcondition]-nonmaxgridindex
            density[nonmaxgridindex+1]+=nonmaxgridindex+1-species.position[nonmaxcondition]
##            if maxcondition.any():
##                print maxcondition, maxgridindex, species.position[maxcondition]
##                print "FUCKUP"
####                break
            density[maxgridindex]+=species.position[maxcondition]-maxgridindex
            density[maxzeros]+=maxgridindex+1-species.position[maxcondition]
##            density[
##            print self.X
##            print species.position
##            self.density[:NG][nonmaxgridindex]+=(self.X[nonmaxgridindex]-species.position)/self.dX#LINE1
##            self.density[:NG][nonmaxgridindex]+=(species.position-self.X[gridindex[nonmaxgridindex]])/self.dX#LINE2
##            except:
##                #print species.position
##                print gridindex
##                print self.X
##            self.density[NG-1] += (self.X[NG] - species.position[:,0])/self.dX
##            self.density[0] += (species.position[:,0][-self.X[NG-1])/self.dX
            density *= species.charge
            self.density+=density
##        self.density[NG-1]=self.density[0] #a może bez tego?
        self.density*=1./self.dX

#FIELDS
##        gridindex=(self.position/grid.L*NG).astype(int)
##
##        nonmaxcondition=gridindex<NG-1
##        nonmaxgridindex = gridindex[nonmaxcondition]
##        maxcondition=gridindex==NG-1
##        maxgridindex=gridindex[maxcondition]
##
##        EField=np.zeros_like(self.position)
##        EField[nonmaxcondition] = (grid.X[nonmaxgridindex+1]-self.position[nonmaxcondition])/grid.dX*grid.efield[nonmaxgridindex]+(self.position[nonmaxcondition]-grid.X[nonmaxgridindex])/grid.dX*grid.efield[nonmaxgridindex+1]
##        EField[maxcondition] = -((grid.X[maxgridindex-1]-self.position[maxcondition])/grid.dX*grid.efield[maxgridindex]+(self.position[maxcondition]-grid.X[maxgridindex])/grid.dX*grid.efield[maxgridindex-1])


        if self.init:
            plt.plot(self.X, self.density)
            plt.ylabel("Density")
            plt.title("Density")
            plt.savefig("Initdensity.png")
            plt.clf()
            self.init=False


        #FOURIER TRANSFORM
##        print "FOURIER"
        self.freq=self.L*fftfreq(NG, self.dX)
##        print "freq", self.freq
        self.freq[0]=0.01
##        print "fixfreq", self.freq
##        print "dens", self.density
        self.pot = np.real(ifft(fft(self.density)[0:NG]/self.freq[0:NG]**2/4./pi**2/eps_0))
##        print "pot", self.pot
        self.efield = -np.gradient(self.pot)
##        print "efield", self.efield
##        self.efield=np.cos(self.X)
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

        EField=np.zeros_like(self.position)
        EField[nonmaxcondition] = (grid.X[nonmaxgridindex+1]-self.position[nonmaxcondition])/grid.dX*grid.efield[nonmaxgridindex]+(self.position[nonmaxcondition]-grid.X[nonmaxgridindex])/grid.dX*grid.efield[nonmaxgridindex+1]
        EField[maxcondition] = -((grid.X[maxgridindex-1]-self.position[maxcondition])/grid.dX*grid.efield[maxgridindex]+(self.position[maxcondition]-grid.X[maxgridindex])/grid.dX*grid.efield[maxgridindex-1])
##        maxgridrollzero=np.zeros_like(maxgridindex)
##        EField[maxcondition=(grid.X[maxgridrollzero]-self.position[nonmaxcondition
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
    plt.subplot(2,1,1)
    for i in ListOfSpecies:
        plt.plot(i.trajectories.T, label=i.name, color=i.color)
    plt.title("Run history")
    plt.ylabel('X position')

    plt.subplot(2,1,2)
    for i in ListOfSpecies:
        plt.plot(i.velocities.T, label=i.name, color=i.color)
    plt.title("Run history")
    plt.ylabel('X velocity')
    plt.xlabel("Iterations")
    
    plt.show()
Grid=grid()
electrons=species(1., -ELEMCHARGE, np.random.random((N,Dim))*Size, np.random.random((N,Dim))*Size/2, N, "electrons", "y")
protons=species(1., ELEMCHARGE, np.random.random((N,Dim))*Size, np.random.random((N,Dim))*Size/2, N, "protons", "b")
Species=[protons,electrons]
Grid.update(Species)
for i in Species:
    i.accelerate(Grid, -0.5*DT)
for iterat in range(RUNITERS):
    for i in Species:
        i.step()
    Grid.update(Species)
PlotAllTrajectories(Species)
