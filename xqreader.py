import numpy as np
from scipy.interpolate import griddata

def xq2arr (file):
    if type(file)==str:
        file=open(file)
    [file.readline() for _ in (0,1,2)]
    out=np.loadtxt(file)
    file.close()
    return out

def xq2griddata (file):
    if type(file)==str:
        file=open(file)
    line=file.readline()
    (_,nx,ny)=[int(x)*1j for x in line.split()]
    print (_,nx,ny)
    [file.readline() for _ in (1,2)]
    dat=np.loadtxt(file)
    file.close()
    xmin,xmax=dat[:,0].min(),dat[:,0].max()
    ymin,ymax=dat[:,1].min(),dat[:,1].max()
    xgrid,ygrid=np.mgrid[xmin:xmax:nx,ymin:ymax:ny]
    print [xmin,xmax,nx,ymin,ymax,ny]
    print xgrid.shape
    print ygrid.shape
    print dat.shape
    return griddata(dat[:,0:2],dat[:,2], (xgrid,ygrid))
