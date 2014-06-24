import numpy as np
import matplotlib.pylab  as pl
import os
from scipy.interpolate import griddata
import utils as ut
import gzip

def xq2arr (file):
    try:
        f=gzip.open(file)
    except:
        f=open(file)
    out=np.loadtxt(file, skiprows=3)
    f.close()
    return out

def xq2griddata (file, factor=1):
    try:
        f=gzip.open(file)
        _,nx,ny=[int(x) for x in f.readline().split()]
    except:
        f=open(file)
        _,nx,ny=[int(x) for x in f.readline().split()]        
    dat=np.loadtxt(file, skiprows=3)
    xmin,xmax=dat[:,0].min(),dat[:,0].max()
    ymin,ymax=dat[:,1].min(),dat[:,1].max()
    xgrid,ygrid=np.mgrid[xmin:xmax:nx*1j/factor,ymin:ymax:ny*1j/factor]
    return griddata(dat[:,0:2],dat[:,2], (xgrid,ygrid))

def getdata (name, factor):
    try:
        rho=xq2griddata(name+'Rho.xq', factor)
    except:
        rho=xq2griddata(name+'Rho.xq.gz', factor)
    try:
        u1=xq2griddata(name+'U1.xq', factor)
    except:
        u1=xq2griddata(name+'U1.xq.gz', factor)
    try:
        u2=xq2griddata(name+'U2.xq', factor)
    except:
        u2=xq2griddata(name+'U2.xq.gz', factor)
    try:
        p=xq2griddata(name+'Pg.xq', factor)
    except:
        p=xq2griddata(name+'Pg.xq.gz', factor)
    try:
        L=xq2griddata(name+'L.xq', factor)
    except:
        L=xq2griddata(name+'L.xq.gz', factor)

    return (rho,u1,u2,p,L)    

def makeplot(name, vrate=30, drawPlot=True, factor=1, getshape=0):
    rho,u1,u2,p,L=getdata(name, factor)
    if drawPlot:
        pl.clf() 
        pl.imshow(np.log10(rho.T), origin='image', cmap=pl.cm.gist_heat_r)
        pl.colorbar()
        s=rho.shape
        pl.barbs(np.arange(0,s[0],vrate),np.arange(0,s[1],vrate),u1[::vrate,::vrate].T,u2[::vrate,::vrate].T)
    if getshape:
        try:
            f=gzip.open(name+'Rho.xq.gz')
            _,nx,ny=[int(x) for x in f.readline().split()]
        except:
            f=open(name+'Rho.xq')
            _,nx,ny=[int(x) for x in f.readline().split()] 
        dat=np.loadtxt(f, skiprows=2)
        xmin,xmax=dat[:,0].min(),dat[:,0].max()
        ymin,ymax=dat[:,1].min(),dat[:,1].max()
        return (rho,u1,u2,p,L ,((xmin,xmax),(ymin,ymax)))
    else:
        return (rho,u1,u2,p,L)

def plotfromdata((rho,u1,u2,p), vrate=30, cmap=pl.cm.gist_heat_r, interpolation='nearest', axes=False, **args):
    if axes:
        c=axes.imshow(np.log10(rho.T), origin='image', cmap=cmap, interpolation=interpolation, **args)
        
    else:
        pl.clf(); pl.imshow(np.log10(rho.T), origin='image', cmap=cmap, interpolation=interpolation, **args)
        c=pl.colorbar()
        s=rho.shape
        pl.quiver(np.arange(0,s[0],vrate),np.arange(0,s[1],vrate),u1[::vrate,::vrate].T,u2[::vrate,::vrate].T)
    return c

def makeStarwindHeader (name):
    rho=xq2arr(name+'Rho.xq')
    V=xq2arr(name+'V.xq')
    s=V.shape
    X=np.empty((s[0],3))
    X[:,0:2]=rho
    X[:,2]=V[:,1]
    ut.arr2h(X, 'star_arr', name+'.h')
    return 'done'
