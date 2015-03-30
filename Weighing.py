from Dependencies import *
from Initialization import *
print "Weighing and fields..."
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

def EFieldInterpolator(positions, gridpoints, field):
    gridspacing=gridpoints[1]-gridpoints[0]
    interpolatedresult=np.zeros_like(positions)
    i=0
    N=len(positions)
    while i<N:
        odleglosc=np.abs(gridpoints-positions[i])
        idmin=odleglosc.argmin()
        if idmin==len(gridpoints)-1:
            idmin1=idmin-1
            idmin2=idmin
            val1=np.abs(gridpoints[idmin2]-positions[i])/gridspacing
            val2=np.abs(gridpoints[idmin1]-positions[i])/gridspacing
        elif idmin==0:
            idmin1=idmin
            idmin2=idmin+1
            val1=np.abs(gridpoints[idmin2]-positions[i])/gridspacing
            val2=np.abs(gridpoints[idmin1]-positions[i])/gridspacing
        else:
            if odleglosc[idmin-1]>odleglosc[idmin+1]:
                idmin1=idmin
                idmin2=idmin+1
                val1=np.abs(gridpoints[idmin2]-positions[i])/gridspacing
                val2=np.abs(gridpoints[idmin1]-positions[i])/gridspacing
            else:
                idmin1=idmin-1
                idmin2=idmin
                val1=np.abs(gridpoints[idmin2]-positions[i])/gridspacing
                val2=np.abs(gridpoints[idmin1]-positions[i])/gridspacing
        if DEBUG: print idmin, positions[i], val1, val2
        interpolatedresult[i]=val1*field[idmin1]+val2*field[idmin2]
        i+=1
    return interpolatedresult

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
    E[1:NG-1]=-(potential[2:NG]-potential[0:NG-2])/(2*gridspacing)
    E[0]=-(potential[1]-potential[0])/gridspacing
    E[NG-1]=-(potential[NG-1]-potential[NG-2])/gridspacing
    return E

#Diagnostics

def PlotFields(gridpoints, position, q, chargedensity, phi, field):
    plt.plot(gridpoints, chargedensity, "ro-", label="Charge")
    plt.scatter(gridpoints, np.zeros(NG), label="grid")
    plt.scatter(position, q*3, c='g', label="particles")
    plt.plot(gridpoints, phi, label="GaussSeidel solved")
    plt.plot(gridpoints,field, label="Field")
    plt.legend()
    plt.grid()
    plt.show()

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
    N=100
    L=10
    NG=51
    gridpoints=np.linspace(0,L,NG)
    gridspacing=gridpoints[1]-gridpoints[0]
    position=L*np.random.random(N)
    q=-1+2*np.random.random(N)
    
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
##    potential=np.real(np.fft.irfft(potenti

##    altransform))
####    if DEBUG: print potential
####    plt.plot(gridpoints, potential, "y-", label="FFT")

    chargedensity=LinearInterpolateToGrid(position, q, gridpoints)
    phi=GaussSeidel1DSolver(chargedensity)
    field=EFieldSolver(phi)
    Forces=EFieldInterpolator(position, gridpoints, field)
    PlotFields(gridpoints, position, q, chargedensity, phi, field)
