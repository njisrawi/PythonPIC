from Dependencies import *

if __name__=="__main__":
    x=np.linspace(0,10,1000)
    y=sin(x)
    xgrid=np.linspace(0,10,6)
    yinterp=np.interp(xgrid, x, y)
    plt.plot(x,y)
    plt.plot(xgrid,yinterp, "o")
    plt.show()
