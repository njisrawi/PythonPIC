from Dependencies import *

DEBUG=False

def LinearInterpolateToGrid(position, weight, gridpoints):
    gridspacing=gridpoints[1]-gridpoints[0]
    interpolatedresult=np.zeros_like(gridpoints)
    i=0
    N=len(position)
    while i<N:
        odleglosc=np.abs(gridpoints-position[i])
        idmin=odleglosc.argmin()
        if idmin==len(gridpoints)-1:
            idmin1=idmin-1
            idmin2=idmin
            val1=np.abs(gridpoints[idmin2]-position[i])/gridspacing
            val2=np.abs(gridpoints[idmin1]-position[i])/gridspacing
            if DEBUG: print idmin, position[i], val1, val2
        elif idmin==0:
            idmin1=idmin
            idmin2=idmin+1
            val1=np.abs(gridpoints[idmin2]-position[i])/gridspacing
            val2=np.abs(gridpoints[idmin1]-position[i])/gridspacing
            if DEBUG: print idmin, position[i], val1, val2
        else:
            if odleglosc[idmin-1]>odleglosc[idmin+1]:
                idmin1=idmin
                idmin2=idmin+1
                val1=np.abs(gridpoints[idmin2]-position[i])/gridspacing
                val2=np.abs(gridpoints[idmin1]-position[i])/gridspacing
            else:
                idmin1=idmin-1
                idmin2=idmin
                val1=np.abs(gridpoints[idmin2]-position[i])/gridspacing
                val2=np.abs(gridpoints[idmin1]-position[i])/gridspacing
        interpolatedresult[idmin1]+=val1*weight[i]
        interpolatedresult[idmin2]+=val2*weight[i]
        i+=1
    return interpolatedresult/gridspacing

def GaussSeidel1DSolver(chargedensity):
    target=1e-6
    Delta=1.0
    phi=np.zeros(NG)
    phiprime=np.zeros(NG)
    while Delta>target:
        if DEBUG: print np.max(phiprime)
        for i in range(NG):
            if i==0 or i==NG-1:
                phiprime[i]=phi[i]
            else:
                phiprime[i]=((phi[i+1]+phi[i-1]) + chargedensity[i]/Epsilon)/2.
                if DEBUG: print phiprime[i]
        Delta = np.max(np.abs(phi-phiprime))
        phi,phiprime = phiprime,phi 
    return phi

def EFieldSolver(potential):
    i=0
    E=np.zeros_like(potential)
    print potential[1:NG-1].shape
    print potential[1:NG].shape
    print potential[0:NG-1].shape
    E[1:NG-1]=-(potential[2:NG]-potential[0:NG-2])/(2*gridspacing)
    E[0]=-(potential[1]-potential[0])/gridspacing
    E[NG-1]=-(potential[NG-1]-potential[NG-2])/gridspacing
    return E
##Continuous weighting function
##def FirstOrderWeighting(x, position, weight, dx):
##    wartosc=np.zeros_like(x)
##    i=0
##    N=len(position)
##    while i<N:
##        pwartosc=1-np.abs(x-position[i])/dx
##        pwartosc[pwartosc<0]=0
##        wartosc+=pwartosc*weight[i]
##        i+=1
##    return wartosc

if __name__=="__main__":
    np.random.seed(1)
    N=500
    L=10
    NG=101
    gridpoints=np.linspace(0,L,NG)
    gridspacing=gridpoints[1]-gridpoints[0]
    position=10*np.random.random(N)
    q=-1+2*np.random.random(N)
    chargedensity=LinearInterpolateToGrid(position, q, gridpoints)

    #FOURIERTRANSFORM
##    plt.plot(gridpoints, chargedensity, "ro-")
##    plt.scatter(gridpoints, np.zeros(NG))
##    plt.scatter(position, np.zeros(N), c='g')
####    plt.show()
##    
##    chargetransform=np.fft.rfft(chargedensity)
##    kvector=np.fft.rfftfreq(NG)*NG/
##    kvector=np.zeros_like(gridpoints)
##    kmin= 2*pi/L
##    for i in range(NG):
##        if i<NG/2:
##            kvector[i] = i*kmin
##        else:
##            kvector[i]=(i-NG)*kmin
##    print kvector, kvector2, kvector-kvector2
    
##    potentialtransform=chargetransform/kvector/kvector/Epsilon
    #if DEBUG:
##    print kvector, chargetransform
##    plt.plot(kvector, np.abs(chargetransform))
##    plt.plot(kvector, np.abs(potentialtransform))
##    plt.show()

##    if DEBUG: print "plotting potential"
##    potential=np.real(np.fft.irfft(potentialtransform))
##    if DEBUG: print potential
    plt.plot(gridpoints, chargedensity, "ro-", label="Charge")
##    plt.plot(gridpoints, potential, "y-", label="FFT")
    plt.scatter(gridpoints, np.zeros(NG), label="grid")
    plt.scatter(position, np.zeros(N), c='g', label="particles")

    phi=GaussSeidel1DSolver(chargedensity)
    plt.plot(gridpoints, phi, label="GaussSeidel solved")
    field=EFieldSolver(phi)
    plt.plot(gridpoints,field, label="Field")
    plt.legend()
    plt.grid()
    plt.show()
    
