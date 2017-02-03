from __future__ import division
from math import *
import utils as ut
from utils import lNot, lAnd, lOr, lXor
import xqreader as xq
import numpy as np
from scipy import constants as cns
from scipy.integrate import quad
from scipy.interpolate import UnivariateSpline,griddata
import pylab as pyl
import dill
import FreeFree as ff
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib import rcParams

#constants
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

    def __init__ (self, inp, unitLength=5.5*r_sol, mu=1.3*cns.m_p, timestamp=False, convert=1, lum=8500*L_sol, shape=(900,600), doLerp=0):

        if type(inp)==str:
            rhoarr=xq.xq2arr(inp+'Rho.xq2')
            pgarr =xq.xq2arr(inp+'Pg.xq2')
            u1arr =xq.xq2arr(inp+'U1.xq2')
            u2arr =xq.xq2arr(inp+'U2.xq2')
            Larr  =xq.xq2arr(inp+'L.xq2')

            Rmax=rhoarr[:,0].max()
            Zmax=rhoarr[:,1].max()
            
            rmax=max(Rmax,Zmax)
            rgrid,qgrid=np.mgrid[0:rmax:shape[0]*1j, 0:pi/2:shape[1]*1j]

            r=rgrid
            theta=qgrid
            R=rgrid*np.sin(qgrid)
            Z=rgrid*np.cos(qgrid)
            self.dtheta=pi/2/shape[1]

            rho=griddata(rhoarr[:,0:2],rhoarr[:,2],(R,Z), method='nearest')
            pg=griddata (pgarr[:,0:2], pgarr[:,2], (R,Z), method='nearest')
            u1=griddata (u1arr[:,0:2], u1arr[:,2], (R,Z), method='nearest')
            u2=griddata (u2arr[:,0:2], u2arr[:,2], (R,Z), method='nearest')
            L=griddata  (Larr[:,0:2],  Larr[:,2],  (R,Z), method='nearest')

            if doLerp: #doLerp flag redoes the interpolation as linear rather than nearest neighbour, then overwrites all values which arent returned nan
                rhoLin=griddata(rhoarr[:,0:2],rhoarr[:,2],(R,Z), method='linear')
                pgLin=griddata (pgarr[:,0:2], pgarr[:,2], (R,Z), method='linear')
                u1Lin=griddata (u1arr[:,0:2], u1arr[:,2], (R,Z), method='linear')
                u2Lin=griddata (u2arr[:,0:2], u2arr[:,2], (R,Z), method='linear')
                LLin=griddata  (Larr[:,0:2],  Larr[:,2],  (R,Z), method='linear')
                
                mask=lNot(reduce(lOr, [np.isnan(rhoLin),np.isnan(pgLin),np.isnan(u1Lin),np.isnan(u2Lin),np.isnan(LLin)]))
                #mask out any elements which have a nan for any variable 
                
                rho[mask]=rhoLin[mask]
                pg[mask] =pgLin[mask]
                u1[mask] =u1Lin[mask]
                u2[mask] =u2Lin[mask]
                L[mask]  =LLin[mask]

            self.S=rho.shape
    
        elif type(inp)==list or type(inp)==tuple:
            rho,pg,u1,u2,L,R,Z=inp
            r=np.sqrt(R*R+Z*Z)
            theta=np.arctan2(R,Z)

        elif type(inp)==type(nulldisc):
            rho=inp.rho()
            pg=inp.pg()
            u1=inp.U1()
            u2=inp.U2()
            L=inp.L()
            R=inp.R
            Z=inp.Z
            unitLength=inp.unitLength
            convert=0
            shape=rho.shape
            lum=inp.lum
            r=np.sqrt(R*R+Z*Z)
            theta=np.arctan2(R,Z)

            
        if timestamp: self.timestamp=timestamp
        else        : self.timestamp=None

        self.R,self.Z=R,Z
        self.unitLength=unitLength

        if convert:
            self.R*=unitLength
            self.Z*=unitLength
            r*=unitLength

            u1 = u1*self.unitLength
            u2 = u2*self.unitLength
            rho= rho/self.unitLength**3
            pg = pg/self.unitLength
            L  = L*self.unitLength**2 #scaling to SI
            
        self.dr=r[1,0]-r[0,0]
        self.r=lambda: r
        self.theta=lambda: theta
        self.dtheta=theta[1,1]-theta[0,0]
        self.unitLength=unitLength
        self.dphi=1e-9

        self.area=self.dr*r*self.dtheta         #area of a cell
        self.vol=self.area*r*np.sin(theta)*self.dphi #volume of a cell given dphi=1e-9 radians
                  
        self._rho   =rho
        self._rho[np.isnan(self._rho)]=np.nanmin(self._rho)
        self._U1    =u1
        self._rho[np.isnan(self._U1)]=np.nanmin(self._U1)
        self._U2    =u2
        self._rho[np.isnan(self._U2)]=np.nanmin(self._U2)
        self._pg    =pg
        self._rho[np.isnan(self._pg)]=np.nanmin(self._pg)
        self._L     =L 
        self._L[np.isnan(self._L)]=np.nanmin(self._L)
        #paranoid nan(droid) removal

        self.rho   =lambda :self._rho.copy()
        self.U1    =lambda :self._U1.copy()
        self.U2    =lambda :self._U2.copy()
        self.pg    =lambda :self._pg.copy()
        self.L     =lambda :self._L.copy() #only access copies of these via the function to prevent accidental modification
        self.p1    =lambda :self._U1*self._rho
        self.p2    =lambda :self._U2*self._rho
        self.N     =lambda :self._rho/mu
        
        self.Uphi  =lambda :(self._L/(self.R+1e-3*self.unitLength)) 
        self.Ur    =lambda :self._U1*np.sin(self.theta())+self._U2*np.cos(self.theta())
        self.Utheta=lambda :self._U1*np.cos(self.theta())-self._U2*np.sin(self.theta())  
        self.Ut    =lambda :np.sqrt(self._U1**2+self._U2**2+self.Uphi()**2)
        self.Uplane=lambda :np.sqrt(self._U1**2+self._U2**2)

        self.T     =lambda :self.pg()/self.rho()/cns.Boltzmann*mu
        self.lum=lum
        
        self._pg[self._pg<0.01*self.ke()]=0.01*self.ke()[self._pg<0.01*self.ke()]
        self._U1[self._pg<0.01*self.ke()]*=sqrt(0.99) # if pressure is < 0.01 ke floor it to that and remove it from the ke
        self._U2[self._pg<0.01*self.ke()]*=sqrt(0.99) # ie mach number of the gas is limited to 100
       
        self.ionisation=ut.memoize(self.ionisation) 

    def coerce2RZ (self, array, shape=(900,600)):
        scale=np.sqrt(shape[0]*shape[0]+shape[1]*shape[1])
        out=np.zeros(shape, dtype=np.float)
        lerp=ut.interpArray(array)
        
        for R in xrange(shape[0]):
            for z in xrange(shape[1]):
                out[R,z]=lerp[sqrt(R*R+z*z)/scale*array.shape[0]-1e-3,np.arctan2(R,z)*2/pi*array.shape[1]-1e-3]
        return out


    def ke (self, incRot=0):
        if incRot: ans=self.rho()*np.sqrt(self.U1()**2+self.U2()**2+self.Uphi()**2)/2
        else:      ans=(self.p1()*self.p1()+self.p2()*self.p2())/self.rho()/2
        return ans

    def mach (self, incRot=0):
        return self.ke(incRot)/self.pg()

    def pr (self):
        return self.p1()*np.cos(self.theta())+self.p2()*np.sin(self.theta())

    def ptheta (self):
        return self.p1()*np.sin(self.theta())+self.p2()*np.cos(self.theta())

    def mass (self):
        return self.rho()*self.vol

    def makeplot(self, f=None, **args):
        if not(f): f=pyl.figure()
        ax=f.add_subplot(111)
        xq.plotfromdata([ut.polar2cartesian(x).T for x in (self.rho(),self.U1(),self.U2(),self.pg())], **args)
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

    def massflux(self,sum=1, debug=0,stripSubSonic=False, unboundOnly=1 ,frac=0.9, mask=None, centralMass=10*m_sol):
        rho=self.rho()
        if np.any(mask) : rho=rho*mask
        if unboundOnly: rho[self.boundGas(centralMass)]=0
        if stripSubSonic : rho[self.ke(incRot=1)<(stripSubSonic*self.pg())]=0
        mom_r_rq=self.Ur()*rho
        mom_r_rq[mom_r_rq<0]=0 #dont subtract material moving inward from the mass calculation
        mom_r_rq*=self.vol
        mom_r_rq/=self.dr #rho.v *dv /dr =momflux*cellarea =massflux per cell
        mom_r_rq*=2*pi/self.dphi #mass flux per latitude
        arr=2*mom_r_rq[int(rho.shape[0]*frac),:] # double for both disc halves
#        arr*=np.sin(np.linspace(0,pi/2,arr.shape[0])) #sin (phi) included here
#        if debug: return(rho,arr.sum()*3600*24*365.25/m_sol * 4*pi*(frac*max(self.r().max(),self.Z.max()))**2 *(pi/2)/arr.size)
#        elif sum :return arr.sum()*3600*24*365.25/m_sol * 4*pi*(frac*self.r().max())**2 *(pi/2)/arr.size
#        else     :return arr*3600*24*365.25/m_sol       * 4*pi*(frac*self.r().max())**2 *(pi/2)/arr.size # r.d_phi * 2 pi r sin(phi)   * 2

        if debug: return(rho,arr.sum()*3600*24*365.25/m_sol) # * 4*pi*(frac*max(self.r().max(),self.Z.max()))**2 *(pi/2)/arr.size)
        elif sum :return arr.sum()*3600*24*365.25/m_sol #* 4*pi*(frac*self.r().max())**2 *(pi/2)/arr.size
        else     :return arr*3600*24*365.25/m_sol       #* 4*pi*(frac*self.r().max())**2 *(pi/2)/arr.size # r.d_phi * 2 pi r sin(phi)   * 2

    def discPolarRatio(self, centralMass,debug=0, frac=0.9):
        phi=self.massflux(sum=0, frac=frac)
        mask=self.discMaterial(centralMass)[self.rho().shape[0]*frac,:]
        if debug : return [phi[mask],phi[lNot(mask)]]
        return phi[mask].sum()/phi[lNot(mask)].sum()

    def boundGas(self, centralMass):
        if centralMass <1e6: centralMass*=m_sol #assume masses < 1e6 are in solar mass not kg
        return (self.ke(incRot=1)+self.pg())<(centralMass*cns.gravitational_constant*self.rho()/self.r())

    def discMaterial(self, centralMass):
        if centralMass <1e6: centralMass*=m_sol #assume masses < 1e6 are in solar mass not kg
        return (self.L()>(np.sqrt(self.R*centralMass*cns.G)/10)) # material is from the disc if it has at least 50% of its material with angular momentum > that of the disc at the stellar surface
    
    def avePolarVel(self, centralMass=10*m_sol, ionisedOnly=False, neutralOnly=False):
        mask = lAnd(lNot(self.discMaterial(centralMass)), lNot(self.boundGas(centralMass)))
        if ionisedOnly  :mask=mask*self.ionisation()
        elif neutralOnly:mask=mask*abs(1-self.ionisation())
        return (mask*self.Uplane()*self.vol).sum()/(mask*self.vol).sum()

    def aveDiscVel(self, centralMass=10*m_sol, ionisedOnly=False, neutralOnly=False):
        mask = lAnd(self.discMaterial(centralMass), lNot(self.boundGas(centralMass)))
        if ionisedOnly :mask=mask*self.ionisation()
        elif neutralOnly:mask=mask*abs(1-self.ionisation())
        return (mask*self.Uplane()*self.vol).sum()/(mask*self.vol).sum()

    def discIonFraction(self, centralMass):
        mask=lAnd(self.discMaterial(centralMass),lNot(self.boundGas(centralMass)))
        return (self.mass()*self.ionisation()*mask).sum()/(self.mass()*mask).sum()

    def windPerformance(self, sum=1, frac=0.9):
        "from sim2004 for spherical winds it is given by Phi.v_inf.c/L_* where Phi is the mass loss rate, instead in do mass loss integral * v.c/L*"
        vr_rq=self.Ur()
        vr_rq[vr_rq<0]=0
        massflux=self.massflux(sum=0, frac=frac)*m_sol/(3600*24*365.25)
        arr=massflux*vr_rq[int(vr_rq.shape[0]*frac),:]
        if sum :return arr.sum()*cns.c/self.lum
        else   :return arr      *cns.c/self.lum

    def massWeightedVelocity(self, mask=None, machLim=0.5, ioncut=True, incRot=True, sum=True):
        "returns velocity of supersonic gas weighted by mass"
        if type(mask)==np.ndarray: mask = mask*1.0 #implicit copy and coerce to float
        else: mask=np.ones_like(self.rho())
        
        if machLim: 
            mask*=(self.ke()>(self.pg()*machLim))
        if ioncut: 
            mask*=self.ionisation()
        
        weights=(self.vol*self.rho()*mask)
        if incRot : out= (self.Ut()*weights/weights.sum())
        else      : out= (self.Uplane()*weights/weights.sum())
        if sum : return out.sum()
        else   : return out

    def emissionWeightedVelocity(self, mask=None, machLim=0.5, ioncut=True, incRot=True, sum=True, retWeights=False):
        "returns velocity of supersonic gas weighted by R.rho^2 (rho^2 for recombination, R for amount of visible material"
        if type(mask)==np.ndarray: mask = mask*1.0 #implicit copy and coerce to float
        else: mask=np.ones_like(self.rho())
        
        if machLim: 
            mask*=(self.ke()>(self.pg()*machLim))
        if ioncut: 
            mask*=self.ionisation()
        
        weights=(self.rho()*mask)**2*self.vol
        if incRot : out= (self.Ut()*weights/weights.sum())
        else      : out= (self.Uplane()*weights/weights.sum())
        if retWeights :return (out,weights)
        if sum : return out.sum()
        else   : return out

    def emissionMeasure(self, sum=True, **kwargs):
        i=self.ionisation(**kwargs)
        out= ((self.rho()*i)**2*self.vol)
        if sum : return out.sum()
        else   : return out

    def n2r2(self, constRho=0):
        out=self.rho().copy()
        out*=self.rho()
        out/=(MMM**2)
        out[:constRho,:]=1e-99
        out*=self.r()
        out*=self.r()
        out/=(self.unitLength**2)
        return out

    def ionisation(self, T=None, debug=False, lyn=True, bal=True, zeros=10, inner=None, thermal=True, ionPot=13.5984, C_recomb=False):
        """calculate the ionisation state of the gas given a surface tmperature of the star
If t is not given then assume BB temp for 8500 L_sol and a radius of unitLength
if bal then remove atoms thermally excited into n=2 from the column of material
if thermal then remove thermally ionosed atoms from the column of material
zeros is the number or r cells to skip (to avoid trying to pass through the star itself)
if inner is an array then the last of the zero'd n2r2_sum cells from zeros is filled with an array interpolated from this
for the recombinations, we assume all ionised material is at 10^4K"""
        if T==None : T=T_eff(self.unitLength/r_sol, self.lum)
        if not(inner==None) :inner=UnivariateSpline(np.linspace(0,pi/2,len(inner)), inner, s=0)
        else                :inner=UnivariateSpline(np.arange(-10,10),np.arange(-10,10)*0)
        n2r2=self.n2r2(zeros)
        nphot_lyn=nphot(T,ionPot*cns.electron_volt/cns.h) #photons with energy for n=1->inf
        print nphot_lyn
        ionise=np.zeros_like(n2r2)
        if thermal:
            ionise+=np.exp(-ionPot*cns.electron_volt/(k*self.T())) #boltzman distribution giving thermally ionised atoms
            n2r2*=(1-ionise)**2 #remove thermally ionised from the column
        if bal:
            i_B=2.0/3*np.exp(-3.4*cns.electron_volt/(k*(self.T()))) #boltzman distribution giving proportion of n=2 H atoms 
            i_B[i_B>1]=1
            n2r2*=(1-i_B)**2 #remove atoms in thermally excited to n=2 from the column, assumes thin to lyman -> thin to balmer
        n2r2_sum=n2r2.copy()
        n2r2_sum[zeros,:]=inner(np.linspace(0,pi/2,n2r2_sum.shape[1]))
        if C_recomb:
            n2r2_sum*=4.27e-29*self.T()**-2.487*np.exp(12260/self.T())
        else:
            n2r2_sum*=2e-16*self.T()**-0.75
        n2r2_sum=(n2r2_sum*self.dr).cumsum(0)
        if debug:
            pyl.clf();pyl.imshow(np.log10(n2r2_sum), vmin=np.log10(nphot_lyn)-2);pyl.colorbar()
            pyl.show();raw_input('>')
        if lyn:
            mask=n2r2_sum<(nphot_lyn)
            ionise+=mask*1.0 #total ionisation when sum(n^2 r^2 dr) < number of >=13.6ev photons otherwise default to 1e-50
            borders=lXor(mask, np.roll(mask,1,0))
            borders[0,:]=False #gives a mask for the first element that is thick to the lyman column
            last=np.roll(borders,-1,0)           #gives a mask for the last  element that is thin  to the lyman column
            ionise[borders]+=(nphot_lyn-n2r2_sum[last])/(n2r2[borders]) #increase ionise by the ratio of the  contents of the border cell to the number of leftover photons
                              
        ionise[ionise>1]=1
        ionise[ionise<1e-20]=1e-20
        ionise[:zeros+1,:]=1e-20
        if debug : return ionise,n2r2_sum
        else     : return ionise

    def openingAngle(self,radius=2/3.0):
        s=int(max(self.pr().shape)*radius)
        x=getConstRadius(self.pr(), radius)
        q=np.linspace(0,pi/2,x.size)
        return (q * x*np.sin(q)).sum()/(x*np.sin(q)).sum()
    
    def makeCube(self, arr, size, zeros=3, verbose=0, scale=1):
        "interpolates the density of ionised material onto a cube of shape cubeshape, anti-alias by a factor of aa"
        length=self.R.max()
        nx=size
        ny=nz=size//2
        degrade=min(0, int(round(np.log10(min(arr.shape)/ny)-2)))
        for i in xrange(degrade):
            arr=ut.degrade_arr(ut.degrade_arr(arr,0),1)
        xgd,ygd,zgd=np.mgrid[-1:1:nx*1j,0:1:ny*1j, -0.5:0.5:nz*1j]
        zgd=np.abs(zgd)


        rgd=np.sqrt(xgd**2+ygd**2)

        qgd =(np.arctan2(rgd,zgd))
        qgd/=(0.5*pi)/arr.shape[1]

        rgd[...]=np.sqrt(xgd**2+ygd**2+zgd**2) #rgd converted to spherical r
        rgd*=arr.shape[0]

        rgd[rgd>=arr.shape[0]]=arr.shape[0]*0.999999
        qgd[qgd>=arr.shape[1]]=arr.shape[1]*0.999999


        lerp=ut.interpArray(arr)
        out=np.empty_like(rgd, dtype=np.float32)
        coords=iter(zip(rgd.flat,qgd.flat))
        for x in np.nditer(out, op_flags=['readwrite']):
            x[...]=lerp[coords.next()]

        if scale:
            scale=(arr*self.vol).sum()*2*pi/self.dphi*2 / (out*(2*length/nx)**3).sum()  #ensure sum(n_e**2) is constant after interpolation by scaling 
            if verbose:
                print "sum(rhoi**2*vol) %.3e"%((arr*self.vol).sum()*2*pi/self.dphi*2)
                print "sum(cube**2*vol) %.3e"%((out*(2*length/nx)**3).sum())
                print "scaling by %.3f"%scale
        else:
            scale=1
        return (out*scale)

    def ionisedRhoCube(self, size, zeros=3, verbose=0, surfTemp=None):
        "interpolates the density of ionised material onto a cube of shape cubeshape, anti-alias by a factor of aa"
        out=np.sqrt(self.makeCube((self.rho()*self.ionisation(surfTemp))**2, size, zeros,verbose,scale=1))
        out[-zeros:,...]=1e-99
        out[:zeros, ...]=1e-99
        out[:,-zeros:,:]=1e-99
        out[...,:zeros ]=1e-99
        out[...,-zeros:]=1e-99
        return out

    def neutRhoCube(self, size, zeros=3, verbose=0, surfTemp=None):
        "interpolates the temperature onto a cube of shape cubeshape, anti-alias by a factor of aa"
        out=self.makeCube(self.rho()*(1-self.ionisation(surfTemp)), size, zeros,verbose,scale=1)
        out[-zeros:,...]=1e-99
        out[:zeros, ...]=1e-99
        out[:,-zeros:,:]=1e-99
        out[...,:zeros ]=1e-99
        out[...,-zeros:]=1e-99
        return out

    def temperatureCube(self, size, zeros=3, verbose=0):
        "interpolates the temperature onto a cube of shape cubeshape, anti-alias by a factor of aa"
        T=self.T()
        return self.makeCube(self.T(), size, zeros,verbose,scale=0)

    def velocityCubes(self, size, zeros=3, verbose=0):
        V=np.array([self.makeCube(self.U1(), size, zeros,verbose,scale=0),  self.makeCube(self.U2(), size, zeros,verbose,scale=0),  self.makeCube(self.Uphi(), size, zeros,verbose,scale=0)])
        V_cart=np.zeros_like(V)
        x,y,z=np.mgrid[-size//2:size//2, 0:size//2, -size//4:size//4]
        theta=np.arctan2(x,y)
        V_cart[0,...]=V[0,...]*np.sin(theta)+V[2,...]*np.cos(theta)
        V_cart[1,...]=V[2,...]*np.sin(theta)+V[0,...]*np.cos(theta)
        V_cart[2,...]=V[1,...]
        V_cart[2,:,:,size//4:]*=-1 #invert v_z for z<0
        return V_cart

    def createFreeFreeRT(self, shape, zeros=5):
        return ff.freeFree(self.ionisedRhoCube(shape,zeros=zeros)/1000.0+1e-30,
                           self.neutRhoCube(shape,zeros=zeros)/1000.0+1e-30,
                           np.minimum(1.0e4*np.ones_like(self.neutRhoCube(shape,zeros=zeros)),
                                      self.temperatureCube(shape,zeros=zeros)),
                           self.velocityCubes(shape,zeros=zeros)/1000.0, 
                           self.R.max()*100*2/shape)
        
    def createHeader(self, name, factor=2):
        ut.arr2h(ut.degrade_arr(ut.degrade_arr(np.array([self.rho()*self.unitLength**3, self.U1()/self.unitLength,self.U2()/self.unitLength, self.pg()*self.unitLength, self.L()/self.unitLength**2]),1,factor),2,factor), 'innerDisc_arr', name)
        return 0

    def createPolarHeader(self, name, deconvert=1,frac=0.9):
        s=self.rho().shape
        if s[1]!=900: print "warning mg_g is expecting a 5x900 as the polar input to large models"
        rho=self.rho()[int(s[0]*frac),::-1]
        u1=self.U1()[int(s[0]*frac),::-1]
        u2=self.U2()[int(s[0]*frac),::-1]
        pg=self.pg()[int(s[0]*frac),::-1]
        L=self.L()[int(s[0]*frac),::-1]
        rmin=rho[:-100].min()
        u1min=rho[:-100].min()
        u2min=rho[:-100].min()
        pgmin=rho[:-100].min()
        lmin=rho[:-100].min()
        rho[rho<rmin]=rmin
        u1[u1<u1min]=u1min
        u2[u2<u2min]=u2min
        pg[pg<pgmin]=pgmin
        L[L<lmin]=lmin
#        s=rho.shape
#        if convolve:
#            rho=ut.beam_convolve(rho, convolve)
#            pg=ut.beam_convolve(pg, convolve)
#            u1=ut.beam_convolve(u1, convolve)
#            u2=ut.beam_convolve(u2, convolve)
#            L=ut.beam_convolve(L, convolve)
#            s1=rho.shape[0]-s[0]
#            s2=rho.shape[1]-s[1]
#            rho=rho[s1/2:s1/2+s[0],s2/2:s2/2+s[1]]
#            pg =pg [s1/2:s1/2+s[0],s2/2:s2/2+s[1]]
#            u1 =u1 [s1/2:s1/2+s[0],s2/2:s2/2+s[1]]
#            u2 =u2 [s1/2:s1/2+s[0],s2/2:s2/2+s[1]]
#            L  = L [s1/2:s1/2+s[0],s2/2:s2/2+s[1]]


        if deconvert:
            u1 = u1/self.unitLength
            u2 = u2/self.unitLength
            rho= rho*self.unitLength**3
            pg = pg*self.unitLength
            L  = L/self.unitLength**2 #scaling to SI
            
        ut.arr2h(np.array([rho,pg,u1,u2,L]), 'largeDisc', name)
        return np.array([rho,pg,u1,u2,L])

    def save(self, file):
        if type(file)==str:
            file=open(file, 'w')
        dill.dump(self,file,protocol=-1)
        return 0

    def disc2header(self,name):
        out=np.array([self.rho(),
                      self.ionisation(),
                      np.fmin(self.T(), np.ones_like(self.T())*1e4),
                      self.U1(),
                      self.U2(),
                      self.Uphi()])
        ut.arr2h(out, 'disc', name)
        return 0


def getConstRadius(arr, radius=0.9):
    polar=ut.cartesian2polar(arr)
    s=int(max(arr.shape)*radius)
    return polar[s,:]

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

def nphot (T,vmin, vmax=None):
          "calculates the number of photons emitted per meter^2 per steradian with a freqency greater than vmin"
          if not(vmax):
              vmax=100.0/(h/(k*T)) #only integrate up to (h nu)/kT =100 ~ photons per unit frequency down at 3e-38 of the peak, avoids exp overflows
          return quad(lambda v: plank(T,v)/(h*v) ,vmin,vmax)[0]

def T_eff (R, lum=L_sol):
    return pow(lum/cns.Stefan_Boltzmann/(4*pi*(R*r_sol)**2),0.25)


class timeSeries():
    def __init__(self, inp, shape=(900,600), unitLength=5.5*r_sol):
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
                    d=disc(val, timestamp=key, shape=shape, unitLength=UnitLength)
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

def mergeDiscs (d1, d2):
    "overwrites the smaller of the 2 discs over the correspoding region of the larger one, ensure they use the same unit length"
    if d1.r().max()>d2.r().max(): priority=1
    else                        : priority=0
    discs=(d1,d2)
    dbig,dsmall=discs[not(priority)],discs[priority]
    assert d1.R.shape[1] == d2.R.shape[1]
    shape=dbig.R.shape
    rhoL,u1L,u2L,pgL,LL=dbig.rho(),dbig.U1(),dbig.U2(),dbig.pg(),dbig.L()
    rhoS,u1S,u2S,pgS,LS=dsmall.rho(),dsmall.U1(),dsmall.U2(),dsmall.pg(),dsmall.L()
    rmax=dsmall.r().max()
    rsLarge=sorted([r for r in set(list(dbig.r().flat)) if r<=rmax])
    rsSmall=sorted([r for r in set(list(dsmall.r().flat)) if r<=rmax])
    n=len(rsLarge)
    for i,rL in enumerate(rsLarge):
        x='x'
        diff=+np.inf
        for j,rS in enumerate(rsSmall):
            if abs(rS-rL)<diff:
                x=j
                diff=abs(rS-rL)
        rhoL[i,int(shape[1]/sqrt(2)*(i/n)**2):]=rhoS[x,int(shape[1]/sqrt(2)*(i/n)**2):]
        pgL[i,int(shape[1]/sqrt(2)*(i/n)**2):] =pgS[x,int(shape[1]/sqrt(2)*(i/n)**2):]
        u1L[i,int(shape[1]/sqrt(2)*(i/n)**2):] =u1S[x,int(shape[1]/sqrt(2)*(i/n)**2):]
        u2L[i,int(shape[1]/sqrt(2)*(i/n)**2):] =u2S[x,int(shape[1]/sqrt(2)*(i/n)**2):]
        LL[i,int(shape[1]/sqrt(2)*(i/n)**2):]  =LS[x,int(shape[1]/sqrt(2)*(i/n)**2):]
        
    return disc((rhoL,pgL,u1L,u2L,LL, dbig.R, dbig.Z), unitLength=dbig.unitLength, shape=dbig.S, lum=dbig.lum, convert=0)


def createStarwindHeader(name):
    rho=xq.xq2arr(name+'Rho.xq2')
    pg=xq.xq2arr(name+'Pg.xq2')
    u=xq.xq2arr(name+'U1.xq2')
    rhos=rho[np.argsort(rho[:,0]),:]
    pgs=pg[np.argsort(pg[:,0]),:]
    us=u[np.argsort(u[:,0]),:]
    rhoSpline=UnivariateSpline(rhos[:,0],rhos[:,1],k=1,s=0) #k1 s0 spline is a first order spline that passes through eqach point, ie a lerp across all the points
    pgSpline=UnivariateSpline(pgs[:,0],pgs[:,1],k=1,s=0)
    uSpline=UnivariateSpline(us[:,0],us[:,1],k=1,s=0)
    out=np.array([np.linspace(0,10,2048),rhoSpline(np.linspace(0,10,2048)),uSpline(np.linspace(0,10,2048)),pgSpline(np.linspace(0,10,2048))]).T
    ut.arr2h(out,'star_arr',name+'.h')

