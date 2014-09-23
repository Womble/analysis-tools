from __future__ import division
from math import *
import utils as ut
import pyfits as P
import xqreader as xq
import numpy as np
from scipy import constants as cns
from scipy.integrate import quad
import pylab as pyl
import dill
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib import rcParams

#rcParams['text.usetex'] = True

from numpy import pi

m_sol=1.989e30
r_sol=6.955e8
L_sol=3.839e26
h=cns.Planck
c=cns.speed_of_light
k=cns.Boltzmann
MMM=(0.909*1.008*cns.atomic_mass+0.088*4.002602*cns.atomic_mass) #mean molecular mass k

class dsave():
    def __init__(self, ts,rho,u1,u2,p,L):
        self.stamp=ts
        self.rho=rho
        self.u1=u1
        self.u2=u2
        self.p=p
        self.L=L

class disc ():

    def __init__ (self, inp, unitLength=5.5*r_sol, mu=1.3*cns.m_p, numberOfCellsInUnit=64, timestamp=False, factor=1, convert=0, lum=8500*L_sol, shape=(900,600)):
        if type(inp)==str : 
            rho,u1,u2,pg,L,s=xq.makeplot(inp, drawPlot=False, shape=shape) #expects to read in .xq files from mg named inp+ 'Rho' 'Pg' 'U1' 'U2' and  'L' for the density pressure velocities and spec ang mom of the region
            convert=1
        elif type(inp)==type(dsave(0,0,0,0,0,0)):
            timestamp,rho,u1,u2,pg,L=inp.stamp,inp.rho,inp.u1,inp.u2,inp.p,inp.L
            convert=0
        else : 
            try:
                ds=dill.load(inp)
                rho,u1,u2,pg,L,R,Z=ds.rho(),ds.U1(),ds.U2(),ds.pg(),ds.L(),ds.R,ds.Z
                self.R=R
                self.Z=Z
                convert=0
            except:
                rho,u1,u2,pg,L=inp
        self.S=rho.shape
#        i=int(round(max(self.S)/1000))
#        if i>1:
#            rho=ut.degrade_arr(ut.degrade_arr(rho,0,i),1,i)
#            u1=ut.degrade_arr(ut.degrade_arr(u1,0,i),1,i)
#            u2=ut.degrade_arr(ut.degrade_arr(u2,0,i),1,i)
#            pg=ut.degrade_arr(ut.degrade_arr(pg,0,i),1,i)
#            L=ut.degrade_arr(ut.degrade_arr(L,0,i),1,i)
            
        if timestamp: self.timestamp=timestamp
        else        : self.timestamp=None

        R,Z=np.mgrid[0:self.S[0]:rho.shape[0]*1j, 0:self.S[1]:rho.shape[1]*1j]
        R[R<=1e-6]=1e-6
        Z[Z<=1e-6]=1e-6
        r=np.sqrt(R*R+Z*Z)
        theta=np.arctan2(R,Z) #theta down from z axis

        if type(inp)==str:
            ((xmin,xmax),(ymin,ymax))=s
            self.R,self.Z=xmin+R/R.max()*(xmax-xmin)*unitLength,ymin+Z/Z.max()*(ymax-ymin)*unitLength
        #self.R,self.Z=R*unitLength/numberOfCellsInUnit,Z*unitLength/numberOfCellsInUnit
        self.r=lambda :np.sqrt(self.R**2+self.Z**2)
        self.theta= lambda :np.arctan2(Z,R)
        self.unitLength=unitLength

        if convert:
            u1 = u1*self.unitLength
            u2 = u2*self.unitLength
            rho= rho/self.unitLength**3
            pg = pg/self.unitLength
            L  = L*self.unitLength**2 #scaling to SI
            
            rho[np.isnan(rho)]=1e-30
            pg[np.isnan(pg)]=1e-30
            u1[np.isnan(u1)]=0
            u2[np.isnan(u2)]=0
            L[np.isnan(L)]=0

#        p1,p2=u1*rho,u2*rho
                  
        self._rho   =rho
        self._U1    =u1
        self._U2    =u2
        self._pg    =pg
        self._L     =L 

        self.rho   =lambda :self._rho
        self.U1    =lambda :self._U1
        self.U2    =lambda :self._U2
        self.pg    =lambda :self._pg
        self.L     =lambda :self._L 
        self.p1    =lambda :self.U1()*self.rho()
        self.p2    =lambda :self.U2()*self.rho()
        self.N     =lambda :self.rho()/mu
        self.ke    =lambda :(self.p1()*self.p1()+self.p2()*self.p2())/self.rho()/2
        
        self.Uphi  =lambda :(self.L()/((R+1/32.)*self.unitLength)) #given by sqrt(R * MG)*R R has units of legnth and G has unit s of length^3
        self.T     =lambda :self.pg()/self.rho()/cns.Boltzmann*mu
        self.Ut    =lambda :np.sqrt(self.U1()**2+self.U2()**2+self.Uphi()**2)
        self.Uplane=lambda :np.sqrt(self.U1()**2+self.U2()**2)

        self.lum=lum

    def pr (self):
        return self.p1()*np.cos(self.theta())+self.p2()*np.sin(self.theta())

    def ptheta (self):
        return self.p1()*np.sin(self.theta())+self.p2()*np.cos(self.theta())

    def makeplot(self, f=None, **args):
        if not(f): f=pyl.figure()
        ax=f.add_subplot(111)
        xq.plotfromdata((self.rho(),self.U1(),self.U2(),self.pg()), **args)
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

    def makeCRPs(self, rho=1, pg=1,ke=1,rhov=1,v=0, mach=0,**args):
        try:
            axes=args['axes']
        except KeyError:
            axes=pyl.gca()
        l=[]
        if rho: 
            constRadiusPlot(self.rho(), **args )
            l.append(r'$\rho$')
        if pg : 
            constRadiusPlot(self.pg(), **args )
            l.append(r'$p_g$')
        if ke : 
            constRadiusPlot((self.U1()**2+self.U2()**2)*self.rho()/2, **args )
            l.append(r'$E_k$')
        if mach:
            constRadiusPlot((self.U1()**2+self.U2()**2)*self.rho()/2/self.pg(), **args )
            l.append(r'$M$')
        if rhov:
            constRadiusPlot((self.U1()*np.sin(self.theta())+self.U2()*np.cos(self.theta()))*self.rho(), **args )
            l.append(r'$\rho \dot v$')
        if v:   
            constRadiusPlot((self.U1()*np.sin(self.theta())+self.U2()*np.cos(self.theta())), **args )
            l.append(r'$v$')
        axes.legend(l, loc=6)
        axes.set_xlabel(r'$\theta / \frac{\pi}{2}$')
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

    def massflux(self,sum=1, debug=0,stripSubSonic=1.0/3, frac=0.9):
        rho=self.rho().copy()
        if stripSubSonic : rho[self.ke()<(stripSubSonic*self.pg())]=1e-50
        mom_r_rq=ut.cartesian2polar((self.U1()*np.sin(self.theta())+self.U2()*np.cos(self.theta()))*rho)
        mom_r_rq[mom_r_rq<0]=0 #dont subtract material moving inward from the mass calculation
        arr=mom_r_rq[int(self.rho().shape[0]*frac),:]*np.sin(np.linspace(pi/2,0,mom_r_rq.shape[1])) #sin (phi) included here
        if debug: return(rho,arr.sum()*3600*24*365.25/m_sol * 4*pi*(frac*max(self.R.max(),self.Z.max()))**2 *(pi/2)/arr.size)
        elif sum :return arr.sum()*3600*24*365.25/m_sol * 4*pi*(frac*self.R.max())**2 *(pi/2)/arr.size
        else     :return arr*3600*24*365.25/m_sol       * 4*pi*(frac*self.R.max())**2 *(pi/2)/arr.size # r.d_phi * 2 pi r sin(phi)   * 2

    def discPolarRatio(self, rhoThresh=2,debug=0):
        phi=self.massflux(0)
        rhos=ut.cartesian2polar(self.rho())[self.rho().shape[0]*0.9,:]
        thresh=rhoThresh*sorted(rhos[0:int(rhos.size/2)])[int(rhos.size/4)]
        mask=rhos>thresh
        if debug : return [phi[mask],phi[np.logical_not(mask)]]
        return phi[mask].sum()/phi[np.logical_not(mask)].sum()

    def boundGas(self, centralMass):
        return (self.ke()+self.pg())>(centralMass*cns.gravitational_constant*self.rho()/self.r())

    def windPerfomance(self, sum=1, frac=0.9):
        "from sim2004 for spherical winds it is given by Phi.v_inf.c/L* where Phi is the mass loss rate, instead in do mass loss integral * v.c/L*"
        vr_rq=ut.cartesian2polar((self.U1()*np.cos(self.theta())+self.U2()*np.sin(self.theta())))
        vr_rq[vr_rq<0]=0
        massflux=self.massflux(sum=0, frac=frac)*m_sol/(3600*24*365.25)
        arr=massflux*vr_rq[vr_rq.shape[0]*frac,:]
        if sum :return arr.sum()  *3e8/self.lum
        else   :return arr   *3e8/self.lum

    def emmissionWeightedVelocity(self, mask=None, machLim=1):
        "returns velocity of supersonic gas weighted by R.rho^2 (rho^2 for recombination, R for amount of visible material"
        if mask==None: mask=self.ke()>(self.pg()*machLim)
        weights=(self.rho()*mask)**2*self.R
        return (self.Ut()*weights/weights.sum()).sum()
    
    def createHeader(self, name, factor=2):
        ut.arr2h(ut.degrade_arr(ut.degrade_arr(np.array([self.rho()*self.unitLength**3, self.U1()/self.unitLength,self.U2()/self.unitLength, self.pg()*self.unitLength, self.L()/self.unitLength**2]),1,factor),2,factor), 'innerDisc_arr', name)
        return 0

    def createPolarHeader(self, name, deconvert=1, convolve=0):
        s=ut.cartesian2polar(self.rho()).shape
        rho=self.rho()
        u1=self.U1()
        u2=self.U2()
        pg=self.pg()
        L=self.L()
        rho=ut.degrade_arr(ut.degrade_arr(ut.cartesian2polar(rho),0,s[0]/905.0),1,s[1]/905.0)
        u1=ut.degrade_arr(ut.degrade_arr(ut.cartesian2polar(u1),0,s[0]/905.0),1,s[1]/905.0)
        u2=ut.degrade_arr(ut.degrade_arr(ut.cartesian2polar(u2),0,s[0]/905.0),1,s[1]/905.0)
        pg=ut.degrade_arr(ut.degrade_arr(ut.cartesian2polar(pg),0,s[0]/905.0),1,s[1]/905.0)
        L=ut.degrade_arr(ut.degrade_arr(ut.cartesian2polar(L),0,s[0]/905.0),1,s[1]/905.0)
        s=rho.shape
        if convolve:
            rho=ut.beam_convolve(rho, convolve)
            pg=ut.beam_convolve(pg, convolve)
            u1=ut.beam_convolve(u1, convolve)
            u2=ut.beam_convolve(u2, convolve)
            L=ut.beam_convolve(L, convolve)
            s1=rho.shape[0]-s[0]
            s2=rho.shape[1]-s[1]
            rho=rho[s1/2:s1/2+s[0],s2/2:s2/2+s[1]]
            pg =pg [s1/2:s1/2+s[0],s2/2:s2/2+s[1]]
            u1 =u1 [s1/2:s1/2+s[0],s2/2:s2/2+s[1]]
            u2 =u2 [s1/2:s1/2+s[0],s2/2:s2/2+s[1]]
            L  = L [s1/2:s1/2+s[0],s2/2:s2/2+s[1]]

        if deconvert:
            u1 = u1/self.unitLength
            u2 = u2/self.unitLength
            rho= rho*self.unitLength**3
            pg = pg*self.unitLength
            L  = L/self.unitLength**2 #scaling to SI
            
        for a in [rho,pg,u1,u2,L]:
            pyl.semilogy(a[int(s[0]*0.9),:])
            raw_input('>')

        ut.arr2h(np.array([rho,pg,u1,u2,L])[:,int(s[0]*0.9),:], 'largeDisc', name)
        return np.array([rho,pg,u1,u2,L])

    def save(self, file):
        if type(file)==str:
            file=open(file, 'w')
        dill.dump(self,file,protocol=-1)
        return 0

    def n2r2(self):
        rho=self.rho().copy()
        rho[self.r()<(self.unitLength*1.25)]=rho.min()
        rho_pol=ut.cartesian2polar(rho*(1-2.0/3*np.exp(-13.6*cns.electron_volt/(k*self.T())))) #remove thermally ionised atoms from density
        n_pol=rho_pol/MMM
        r_pol=ut.cartesian2polar(self.r()/self.unitLength)
        return n_pol**2*r_pol**2

    def ionisation(self, T_eff, lyn=1, bal=1, skip=1):
        n2r2=self.n2r2()
        nphot_lyn=nphot(T_eff,3.3e15) #photons with energy for n=1->inf
        print nphot_lyn
        n2r2_sum=n2r2.copy()*2.076e-11*2.2/np.sqrt(ut.cartesian2polar(self.T()))
        ionise=np.zeros_like(n2r2)
        if lyn:
            for i in xrange(1,n2r2.shape[0]):
                if i>skip :n2r2_sum[i,:]+=n2r2_sum[i-1,:]
                else      :n2r2_sum[i,:]=0
#            pyl.clf();pyl.imshow(np.log10(n2r2_sum),vmax=np.log10(nphot_lyn), vmin=0);pyl.colorbar();pyl.show();raw_input()
            ionise=(n2r2_sum<(nphot_lyn))*1.0 #total ionisation when sum(n^2 r^2 dr) < number of >=13.6ev photons otherwise default to 1e-50
            for j in xrange(n2r2.shape[1]-1):
                for i in xrange(n2r2_sum.shape[0]-1):
                    if (n2r2_sum[i,j]>nphot_lyn):
                        ionise[i,j]=max((nphot_lyn-n2r2[i-1,j])/n2r2[i,j] , 1e-20) #partial ionisation
                        break
        if bal:
            i_B=ut.cartesian2polar(2.0/3*np.exp(-3.4*cns.electron_volt/(k*self.T()))) #boltzman distribution giving proportion of n=2 H atoms 
            ionise+=i_B #assume anything in the n=2 state is then ionised by radiation from the star, disc is thin the Ballmer photons
        ionise[ionise>1]=1
        ionise[ionise<1e-20]=1e-20
        return ut.polar2cartesian(ionise)[:self.rho().shape[0],:self.rho().shape[1]]

def emission(ne, T, nu):
    ne[ne<=0]=1e-50
    T [T<=0] =1e-50
    eps=6.8e-38*ne**2/np.sqrt(T)*np.exp(-cns.h*nu/cns.Boltzmann/T)
    return eps

def constRadiusPlot(arr, axes=None, radius=0.9, stlims=(0,1), normalise=False, logplot=False):
    r=int(min(arr.shape)*radius)
    rq_arr=ut.cartesian2polar(arr)
    if normalise : rq_arr/=abs(rq_arr[r,:].max())
    if axes:
        if logplot   : ret=axes.semilogy(np.linspace(0.5,0,rq_arr[r,:].size),rq_arr[r,:])
        else         : ret=axes.plot(np.linspace(0.5,0,rq_arr[r,:].size),rq_arr[r,:])
        axes.set_xlim(stlims)
        pyl.show()
    elif axes=='skip':
        print 'not drawing'
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
    if axes=='skip':
        print 'not drawing'
    elif logplot : 
        axes.semilogy(rq_arr[:,q])
        pyl.show()
    else            : 
        axes.plot(rq_arr[:,q])
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

def plank (T,v):
    return h*v**3/c**2 / (np.exp(h*v/(k*T))-1)

def nphot (T,vmin):
          "calculates the number of photons emitter per meter^2 per steradian with a freqency greater than vmin"
          tot=0
          m=max([np.log10(vmin)+20, 30])
          x=np.logspace(np.log10(vmin),m,100) 
          for i in xrange(99):
              tot+=quad(lambda v: plank(T,v)/(h*v) ,x[i],x[i+1])[0] #split the integral into 100 sections with equal spaced logs
          return tot

class timeSeries():
    def __init__(self, inp, factor=1):
        "inp wants to be a dictonary of form {timestamp:'xqRootName'} where timestamps are floats and xqRootNames are are strings containing the names of the xq files minus 'Rho.xq' at the end, or a filename of a saved timeseries"
        if type(inp)==type('string'):
            dlist=dill.load(open(inp))
            self.timeSeries=[disc(x) for x in dlist]
            for i in xrange(len(self.timeSeries)):
                self.timeSeries[i].timestamp=dlist[i].stamp
            self.timeSeries.sort(key=lambda x: x.timestamp)
        else:
            try:
                dlist=dill.load(inp)
                self.timeSeries=[disc(x) for x in dlist]
                for i in xrange(len(self.timeSeries)):
                    self.timeSeries[i].timestamp=dlist[i].stamp
                self.timeSeries.sort(key=lambda x: x.timestamp)
            except :
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
        kebar=np.array([x.p1()**2+x.p2()**2 for x in self.timeSeries]).mean(0)/rhobar/2
        pgbar=np.array([x.pg() for x in self.timeSeries]).mean(0)-(0.5*(p1bar**2+p2bar**2)/rhobar-kebar) #add/subtract missing ke from momentum averaging into pressure
        Lbar=np.array([x.L() for x in self.timeSeries]).mean(0)
        return disc([rhobar,U1bar,U2bar,pgbar,Lbar], convert=0, unitLength=self.timeSeries[0].unitLength, shape=self.timeSeries[0].rho().shape)

    def save(self, file):
        if type(file)==str:
            file=open(file, 'w')
        dill.dump([dsave(x.timestamp, x.rho(),x.U1(),x.U2(),x.pg(),x.L()) for x in self.timeSeries],file,protocol=-1)
        return 0
