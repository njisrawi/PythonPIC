from Dependencies import *
from Initialization import *
print "Leapfrog systems online..."
DEBUG=False
Dim=3
def Time(i):
    return TInitial+(TFinal-TInitial)/TFinal*i*DT

#Fields - continuous cases

##def EField(r):
##    return 0
##
##def BField(r):
##    result=np.zeros_like(r)
##    result[:,2] = r[:,1]
##    if DEBUG:
##        print "field", result
##    return result
def BField(r):
    return np.array([0,0,1])
##def Force(r, v, EField, BField): #Lorentz
##    f=q*(EField(r)+np.cross(v, BField(r)))
##    if DEBUG:
##        print "force", f
##    return f

#Particle movement

def Force(r, v, q, EFieldArray, BFieldFunc):
##    print "vstack"
##    print np.vstack((q,q,q))
    f=np.vstack((q,q,q)).T*np.cross(v, BField(r))
    f[:,0]+=EFieldArray*q
##    print f
    return f
def InitialPush(r, v, q, m, EFieldArray, BFieldFunc):
    return v+Force(r, v, q, EFieldArray, BFieldFunc)*DT/2./m

def LeapfrogStep(r, v, q, m, EFieldArray, BFieldFunc):
    return r + v*DT, v + Force(r, v, q, EFieldArray, BFieldFunc)*DT/m

#Diagnostics

def KineticEnergy(v, m):
    return np.sum(v*v*m)/2.

#Plotting

def PlotTrajectory2D(rdata, N):
    for j in range(N):
        plt.plot(rdata[j, 0, :], rdata[j, 1, :], "-")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid()
    plt.title("2D trajectory")
    plt.show()
    for j in range(N):
        plt.plot(rdata[j, 0, :], rdata[j, 2, :], "-")
    plt.xlabel("x")
    plt.ylabel("z")
    plt.grid()
    plt.title("2D trajectory")
    plt.show()
    for j in range(N):
        plt.plot(rdata[j, 1, :], rdata[j, 2, :], "-")
    plt.xlabel("y")
    plt.ylabel("z")
    plt.grid()
    plt.title("2D trajectory")
    plt.show()

def PlotTrajectory3D(rdata, N):
    fig=plt.figure()
    ax=fig.add_subplot(111,projection='3d')
    for j in range(N):
        plt.plot(rdata[j, 0, :], rdata[j, 1, :], rdata[j,2,:])#Time(IPoints))
    plt.grid()
    plt.title("3D trajectory")
    plt.show()    

def EPlot(edata):
    print Time(IPoints).shape, edata.shape
    plt.plot(Time(IPoints),edata[0,0,:])
    plt.xlabel("Time")
    plt.ylabel("Energy")
    plt.show()

if __name__=="__main__":
    r=np.random.random((N,Dim))
    v=np.zeros((N,Dim))
    v[:, 0]=-5.
    m=np.ones((N, Dim))
    q=np.ones((N, Dim))
    v=InitialPush(r,v)

    #data acquisition
    rdata=r.copy()
    vdata=v.copy()
    edata=np.array(KineticEnergy(v))

    LastPercent=0
    while(i<=IterFinal):
        r, v = LeapfrogStep(r,v)
        if i%IterStep==0 and i>0:
            j=int((i/float(IterFinal)*100))
            if j>LastPercent:
                print ("%d%% done." %j)
                LastPercent=j
            rdata=np.dstack((rdata,r))
            vdata=np.dstack((vdata,v))
            edata=np.dstack((edata,KineticEnergy(v)))
        i+=1
    
    print "Simulation finished"

    PlotTrajectory2D(r)
    PlotTrajectory3D(r)
    EPlot(edata)
