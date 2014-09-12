import numpy as np
import matplotlib.pylab  as pl
import os
from scipy.interpolate import griddata
from scipy import constants as cns
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

def xq2griddata (file, fill_value=1e-99, shape=None, getshape=0):
    if not(shape):
        try:
            f=gzip.open(file)
            _,nx,ny=[int(x) for x in f.readline().split()]
        except:
            f=open(file)
            _,nx,ny=[int(x) for x in f.readline().split()]        
        dat=np.loadtxt(file, skiprows=2)
    else: 
        nx,ny=shape
        dat=np.loadtxt(file, skiprows=3)
    xmin,xmax=dat[:,0].min(),dat[:,0].max()
    ymin,ymax=dat[:,1].min(),dat[:,1].max()
    xgrid,ygrid=np.mgrid[xmin:xmax:nx*1j,ymin:ymax:ny*1j]
    if getshape:
        return griddata(dat[:,0:2],dat[:,2], (xgrid,ygrid), method='nearest'),((xmin,xmax),(ymin,ymax))
    else:
        return griddata(dat[:,0:2],dat[:,2], (xgrid,ygrid), method='nearest')

def getdata (name,  shape, fill_value=1e-99):
    try:
        rho=xq2griddata(name+'Rho.xq2',  fill_value, shape=shape)
    except:
        rho=xq2griddata(name+'Rho.xq2.gz',  fill_value, shape=shape)
    try:
        u1=xq2griddata(name+'U1.xq2',  fill_value, shape=shape)
    except:
        u1=xq2griddata(name+'U1.xq2.gz',  fill_value, shape=shape)
    try:
        u2=xq2griddata(name+'U2.xq2',  fill_value, shape=shape)
    except:
        u2=xq2griddata(name+'U2.xq2.gz',  fill_value, shape=shape)
    try:
        p=xq2griddata(name+'Pg.xq2',  fill_value, shape=shape)
    except:
        p=xq2griddata(name+'Pg.xq2.gz',  fill_value, shape=shape)
    try:
        L,s=xq2griddata(name+'L.xq2',  fill_value, shape=shape, getshape=1)
    except:
        L,s=xq2griddata(name+'L.xq2.gz',  fill_value, shape=shape, getshape=1)
    

    return (rho,u1,u2,p,L,s)    

def makeplot(name, vrate=30, drawPlot=True, factor=1, shape=None, numDensity=1):
    rho,u1,u2,p,L,s=getdata(name,  shape=shape)
    if drawPlot:
        pl.clf() 
        n=rho.T
        if numDensity: n=rho/(cns.m_p*(0.9*1+4*0.09+12*0.01))/1e3 #particles per cc
        pl.imshow(np.log10(n), origin='image', cmap=pl.cm.gist_heat_r)
        pl.colorbar()
        s=rho.shape
#        sqrtut=(u1[::vrate,::vrate]**2+u2[::vrate,::vrate]**2)**0.5
        pl.barbs(np.arange(0,s[0],vrate),np.arange(0,s[1],vrate),u1[::vrate,::vrate].T,u2[::vrate,::vrate].T)
    return (rho,u1,u2,p,L,s)

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
