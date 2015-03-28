from Dependencies import *
N=10 #number of bodies
Dim = 3 #number of dimensions
TInitial=0.
TFinal=100.
DT=1e-4
i=0
IterFinal=int(TFinal/DT)
SamplingRate=IterFinal/1000
IterStep=IterFinal/SamplingRate
IPoints = np.arange(0, IterFinal+1, IterStep)

def Time(i):
    return TInitial+(TFinal-TInitial)/TFinal*i*DT

def Force(r, v):
    return -r

def InitialPush(r, v):
    return v+Force(r,v)*DT/2./m

def LeapfrogStep(r, v):
    return r + v*DT, v + Force(r,v)*DT/m

if __name__=="__main__":
    r=np.ones((N,Dim))+np.random.random((N, Dim))
    v=np.random.random((N, Dim))
    m=np.ones((N, Dim))
    rdata=r.copy()
    v=InitialPush(r,v)
    
    while(i<=IterFinal):
        r, v = LeapfrogStep(r,v)
        if i%IterStep==0 and i>0:
            print ("%d%% done." %(i/float(IterFinal)*100))
            rdata=np.dstack((rdata,r))
        i+=1
    
    print "Simulation finished"
    for j in range(N):
        plt.plot(rdata[j, 0, :], rdata[j, 1, :], "-")
    plt.show()
    fig=plt.figure()
    ax=fig.add_subplot(111,projection='3d')
    for j in range(N):
        plt.plot(rdata[j, 0, :], rdata[j, 1, :], rdata[j,2,:])#Time(IPoints))
    plt.show()
