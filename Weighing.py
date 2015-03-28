from Dependencies import *
DX=1e-1
np.random.seed(1)
Q=1.6e-19
if __name__=="__main__":
    x=np.random.random(40)*10 #particle locations
    x=np.sort(x)
    xgrid=np.arange(0,10,DX)
    y=Q/DX*np.interp(x-xgrid)
    yinterp=np.interp(xgrid, x, y)
    plt.plot(x,y)
    plt.plot(xgrid,yinterp, "o")
    plt.show()
