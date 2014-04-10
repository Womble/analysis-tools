from __future__ import division
import utils as ut
import pyfits as P
import xqreader as xq
import numpy as np
from scipy import constants as cns
import pylab as pyl
import cPickle
import pprocess as pp
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib import rcParams
rcParams['text.usetex'] = True

from numpy import pi

m_sol=1.989e30
r_sol=6.955e8
L_sol=3.839e26

class dsave():
    def __init__(self, ts,rho,u1,u2,p,L):
        self.stamp=ts
        self.rho=rho
        self.u1=u1
        self.u2=u2
        self.p=p
        self.L=L

class disc ():

    def __init__ (self, inp, unitLength=5.5*r_sol, mu=1.3*cns.m_p, numberOfCellsInUnit=64, timestamp=False, factor=1, convert=0):
        if type(inp)==str : 
            rho,u1,u2,pg,L=xq.makeplot(inp, drawPlot=False, factor=factor)
            convert=1
        elif type(inp)==type(dsave(0,0,0,0,0,0)):
            timestamp,rho,u1,u2,pg,L=inp.stamp,inp.rho,inp.u1,inp.u2,inp.p,inp.L
            convert=0
        else : 
            try:
                ds=cPickle.load(inp)
                rho,u1,u2,pg,L=ds.rho,ds.u1,ds.u2,ds.p,ds.L
                convert=0
            except:
                rho,u1,u2,pg,L=inp
        self.S=rho.shape
        i=int(round(max(self.S)/1000))
        if i>1:
            rho=ut.degrade_arr(ut.degrade_arr(rho,0,i),1,i)
            u1=ut.degrade_arr(ut.degrade_arr(u1,0,i),1,i)
            u2=ut.degrade_arr(ut.degrade_arr(pg,0,i),1,i)
            pg=ut.degrade_arr(ut.degrade_arr(pg,0,i),1,i)
            L=ut.degrade_arr(ut.degrade_arr(L,0,i),1,i)
            
        if timestamp: self.timestamp=timestamp
        else        : self.timestamp=None

        R,Z=np.mgrid[0:self.S[0], 0:self.S[1]]
        R[R<=1e-6]=1e-6
        Z[Z<=1e-6]=1e-6
        r=np.sqrt(R*R+Z*Z)
        theta=np.arctan2(R,Z) #theta down from z axis

        self.R,self.Z=R*unitLength/numberOfCellsInUnit,Z*unitLength/numberOfCellsInUnit
        self.r=lambda :np.sqrt(self.R*2+self.Z**2)
        self.theta= lambda :np.arctan2(R,Z)
        self.unitLength=unitLength

        if convert:
            u1 *=self.unitLength
            u2 *=self.unitLength
            rho/=self.unitLength**3
            pg /=self.unitLength
            L  *=self.unitLength**2 #scaling to SI

#        p1,p2=u1*rho,u2*rho
        self.rho   =lambda :rho
        self.U1    =lambda :u1
        self.U2    =lambda :u2
        self.pg    =lambda :pg
        self.L     =lambda :L 
        self.p1    =lambda :u1*self.rho()
        self.p2    =lambda :u2*self.rho()
        self.N     =lambda :self.rho()/mu
        self.ke    =lambda :np.sqrt(self.p1()*self.p1()+self.p2()*self.p2())/self.rho()
        
        self.Uphi  =lambda :self.L()/(R*self.unitLength) #given by sqrt(R * MG)*R R has units of legnth and G has unit s of length^3
        self.T     =lambda :self.pg()/self.rho()/cns.Boltzmann*mu


    def pr (self):
        dtheta=self.theta-np.arctan2(self.U1(),self.U2())
        return np.sqrt(self.p1()**2+self.p2()**2)*np.cos(dtheta)

    def ptheta (self):
        dtheta=self.theta-np.arctan2(self.U1(),self.U2())
        return np.sqrt(self.p1()**2+self.p2()**2)*np.sin(dtheta)
      
    def makeplot(self, f=None):
        if not(f): f=pyl.figure()
        ax=f.add_subplot(111)
        xq.plotfromdata((self.rho(),self.U1(),self.U2(),self.pg()))
        if self.timestamp: f.axes[0].annotate('%.3e'%self.timestamp, (0,self.S[1]*.9), bbox=dict(fc='0.9'))
        f.show()
        return f

    def radialMassFlow(self, f=None):
        if not(f): f=pyl.figure()
        ax=f.add_subplot(111)
        M=sorted(abs(self.pr()).flat)
        mlow,mup=M[int(self.pr().size/200)],M[int(-self.pr().size/200)] #only use the middle 99% of data
        ut.imshowWslices(np.log10(self.pr()), cbar=True, vmin=np.log10(mlow), vmax=np.log10(mup), yslice=0.9, xslice=0.9)
        if self.timestamp: ax.annotate('%.3e'%self.timestamp, (0,self.S[1]*.9), bbox=dict(fc='0.9'))
        f.show()
        return f

    def kineticEnergy(self, f=None):
        if not(f): f=pyl.figure()
        ax=f.add_subplot(111)
        M=sorted(abs(self.ke()).flat)
        mlow,mup=M[int(self.ke().size/200)],M[int(-self.ke().size/200)] #only use the middle 99% of data
        ut.imshowWslices(np.log10(self.ke()), cbar=True, vmin=np.log10(mlow), vmax=np.log10(mup), yslice=0.9, xslice=0.9)
        if self.timestamp: ax.annotate('%.3e'%self.timestamp, (0,self.S[1]*0.9), bbox=dict(fc='0.9'))
        f.show()
        return f

    def degrade(self, factor):
        rho=ut.degrade_arr(ut.degrade_arr(self.rho(),0,factor),1,factor)
        u1 =ut.degrade_arr(ut.degrade_arr(self.p1(),0,factor),1,factor)/rho
        u2 =ut.degrade_arr(ut.degrade_arr(self.p2(),0,factor),1,factor)/rho
        pg =ut.degrade_arr(ut.degrade_arr(self.pg(),0,factor),1,factor)
        L  =ut.degrade_arr(ut.degrade_arr(self.L(),0,factor),1,factor)
        
        return disc((rho*self.unitLength**3,u1/self.unitLength,u2/self.unitLength,pg*self.unitLength,L/self.unitLength**2))

    def makeCRPs(self, **args):
        try:
            axes=args['axes']
        except KeyError:
            axes=pyl.gca()
        constRadiusPlot(self.rho(), **args )
        constRadiusPlot(self.pg(), **args )
        constRadiusPlot((self.U1()**2+self.U2()**2)*self.rho()/2, **args )
        constRadiusPlot((self.U1()*np.sin(self.theta())+self.U2()*np.cos(self.theta()))*self.rho(), **args )
        axes.legend((r'$\rho$',r'$P_g$',r'$KE$',r'$\rho v.\hat{r}$'), loc=6)
        axes.set_xlabel(r'$\frac{\theta}{2 \pi}$')
        pyl.show()

    def makeCTPs(self, **args):
        try:
            axes=args['axes']
        except KeyError:
            axes=pyl.gca()
        constThetaPlot(self.rho(), **args )
        constThetaPlot(self.pg(), **args )
        constThetaPlot((self.U1()**2+self.U2()**2)*self.rho()/2, **args )
        constThetaPlot((self.U1()*np.sin(self.theta())+self.U2()*np.cos(self.theta()))*self.rho(), **args )
        axes.legend((r'$\rho$',r'$P_g$',r'$KE$',r'$\rho v.\hat{r}$'), loc=1)
        axes.set_xlabel(r'$sin(\theta)$')
        pyl.show()    

    def massflux(self):
        mom_r_rq=ut.cartesian2polar((self.U1()*np.sin(self.theta())+self.U2()*np.cos(self.theta()))*self.rho())
        mom_r_rq[mom_r_rq<0]=0 #dont subtract material moving inward from the mas calculation
        arr=mom_r_rq[int(self.rho().shape[0]*0.9),:]*(0.9*self.unitLength)**2*np.sin(np.linspace(0,pi/2,mom_r_rq.shape[1]))
        return arr.sum()*3600*24*365.25/m_sol * 2


    def save(self, file):
        if type(file)==str:
            file=open(file, 'w')
        cPickle.dump(dsave(self.timestamp, self.rho(), self.U1(), self.U2(), self.pg(), self.L()),file,protocol=-1)
        return 0

def constRadiusPlot(arr, axes=None, radius=0.9, stlims=(0,1), normalise=False, logplot=False):
    r=int(min(arr.shape)*radius)
    rq_arr=ut.cartesian2polar(arr)
    if normalise : rq_arr/=abs(rq_arr[r,:].max())
    if axes:
        if logplot   : ret=axes.semilogy(np.linspace(0.5,0,rq_arr[r,:].size),rq_arr[r,:])
        else         : ret=axes.plot(np.linspace(0.5,0,rq_arr[r,:].size),rq_arr[r,:])
        axes.set_xlim(stlims)
    else:
        if logplot   : ret=pyl.semilogy(np.linspace(1,0,rq_arr[r,:].size),rq_arr[r,:])
        else         : ret=pyl.plot(np.linspace(1,0,rq_arr[r,:].size),rq_arr[r,:])
        pyl.gca().set_xlim(stlims)
    pyl.show()
    return ret

def constThetaPlot(arr, axes=None, Theta=0.95, normalise=False, logplot=False):
    rq_arr=ut.cartesian2polar(arr)
    q=rq_arr.shape[1]*(2/np.pi*Theta)
    if normalise : rq_arr/=abs(rq_arr[:,q].max())
    if not(axes) : axes=pyl.gca()
    if logplot : axes.semilogy(rq_arr[:,q])
    else       : axes.plot(rq_arr[:,q])
    pyl.show()
    return rq_arr

def convTests (discs):
    n=len(discs)*(len(discs)-1)/2
    nrow=int(np.sqrt(n))
    f=pyl.figure()
    grid=ImageGrid(f, 111, (nrow,int(n/nrow+0.5)), cbar_mode='single', cbar_location='top')
    k=0
    for i,d in enumerate(discs):
        if i<(len(discs)-1):
            for j in xrange(i+1, len(discs)):
                d2=discs[j]
                cbar=grid[k].imshow( np.abs(np.log10 (d.rho()/d2.rho())), origin='image', vmax=1)
                k+=1
    grid.cbar_axes[0].colorbar(cbar)
    pyl.show()
    _=raw_input('>>')
    k=0
    pyl.clf(); grid=ImageGrid(f, 111, (nrow,int(n/nrow+0.5)), cbar_mode='single', cbar_location='top')
    for i,d in enumerate(discs):
        if i<(len(discs)-1):
            for j in xrange(i+1, len(discs)):
                d2=discs[j]
                cbar=grid[k].imshow( np.abs(np.log10 (d.pr()/d2.pr() )), origin='image', vmax=1)
                k+=1
    grid.cbar_axes[0].colorbar(cbar)
    pyl.show()
    _=raw_input('>>')
    k=0
    pyl.clf(); grid=ImageGrid(f, 111, (nrow,int(n/nrow+0.5)), cbar_mode='single', cbar_location='top')
    for i,d in enumerate(discs):
        if i<(len(discs)-1):
            for j in xrange(i+1, len(discs)):
                d2=discs[j]
                cbar=grid[k].imshow( np.abs(np.log10 ( d.ke()*d2.rho()/d.rho()/d2.ke() )), origin='image', vmax=1)
                k+=1
    grid.cbar_axes[0].colorbar(cbar)
    pyl.show()

class timeSeries():
    def __init__(self, inp, factor=1):
        "inp wants to be a dictonary of form {timestamp:'xqRootName'} where timestamps are floats and xqRootNames are are strings containing the names of the xq files minus 'Rho.xq' at the end, or a filename of a saved timeseries"
        if type(inp)==type('string'):
            dlist=cPickle.load(open(inp))
            self.timeSeries=[disc(x) for x in dlist]
            for i in xrange(len(self.timeSeries)):
                self.timeSeries[i].timestamp=dlist[i].stamp
            self.timeSeries.sort(key=lambda x: x.timestamp)
        else:
            try:
                dlist=cPickle.load(inp)
                self.timeSeries=[disc(x) for x in dlist]
                for i in xrange(len(self.timeSeries)):
                    self.timeSeries[i].timestamp=dlist[i].stamp
                self.timeSeries.sort(key=lambda x: x.timestamp)
            except ValueError:
                l=[]
                #        tmp=[(key,inp[key]) for key in inp]
                #        def foo((k,v)):
                #            d=disc(v, timestamp=k)
                #            d.timestamp=k
                #            print k
                #            return d
                for key,val in inp.iteritems():
                    d=disc(val, timestamp=key, factor=factor)
                    d.timestamp=key
                    l.append(d)
                    print key
            #        l=pp.pmap(foo,tmp)
                self.timeSeries=sorted(l, key= lambda x: x.timestamp)

    def makeTimeSeriesPlots(self, save=False):
        "creates a series of plots, one for each member of the time series"
        for d in self.timeSeries:
            f=d.makeplot()
            if save: 
                pyl.savefig(save+'%.08d'%d.timestamp+'.png')
                pyl.close(f)

    def timeAverage(self):
        rhobar=np.array([x.rho() for x in self.timeSeries]).mean(0)
        U1bar=np.array([x.p1() for x in self.timeSeries]).mean(0)/rhobar
        U2bar=np.array([x.p2() for x in self.timeSeries]).mean(0)/rhobar
        pgbar=np.array([x.pg() for x in self.timeSeries]).mean(0)
        Lbar=np.array([x.L() for x in self.timeSeries]).mean(0)
        return disc([rhobar,U1bar,U2bar,pgbar,Lbar])

    def save(self, file):
        if type(file)==str:
            file=open(file, 'w')
        cPickle.dump([dsave(x.timestamp, x.rho(),x.U1(),x.U2(),x.pg(),x.L()) for x in self.timeSeries],file,protocol=-1)
        return 0
