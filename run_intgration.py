import dumpParse as dp
from glob import glob
import numpy as np
from scipy.interpolate import Rbf
from scipy.integrate import odeint
from scipy.spatial import Delaunay
import utils as ut
from numpy import logical_and as lAnd
from numpy import logical_or as lOr
from collections import deque
import math

#dumps = glob('./fiducial/dump*.vtu')
#grids=[]
#for d in sorted(dumps):
#    n=d.split('/')[-1].split('.')[0].split('_')[-1]
#    grids.append(dp.Grid(*dp.parseDump(d), sourcefile='fiducial/source'+n+'.dat'))

class RBFstore():
    def __init__(self, array, vref=1e-5):
        self.RBFs=deque([None]*50, maxlen=50)
        self.array=array
        self.Xt=self.array[:,:3]
        self.U=self.array[:,3]
        self.W=self.array[:,4]
        self.minR=array[:,0].min()
        self.vref=vref
        self.constructFullHull()

    def clearRBFs(self):
        self.RBFs=deque([None]*50, maxlen=50)
    

    def constructFullHull(self):        
        zvals=set(self.array[:,1])
        tvals=set(self.array[:,2])
        newpoints=[]
        
        for t in tvals:
            arr=self.array[self.array[:,2]==t,:]
            for z in zvals:
                mask=arr[:,1]==z
                a=np.argpartition(arr[mask,:], 1, axis=0)[0,:]
                a[0]=0;a[3]=0      #append a new point on the axis with V_r=0 
                newpoints.append(a)#

        self.array=np.concatenate([self.array,np.array(newpoints)], axis=0)
                    
        Tmask=lOr(self.array[:,2]==self.array[:,2].max(),
                  self.array[:,2]==self.array[:,2].min())

        arr=self.array[Tmask,:]
        
        self.fullHull=Delaunay(arr[:,:3])        

    def addRBF(self, newXt, n=1000, qhull_opts="", output=True, makeHull=True):
        sampled=False; x=1.0;i=0
        while not(sampled):
            diff=self.Xt-newXt
            diff[:,2]*=self.vref*x
            diff*=diff
            diff=diff.sum(1)
            diff=np.sqrt(diff)
            index=np.argpartition(diff, n)[:n] #provide the indices of the n closest points
            points=self.array[index,:]

            rsampled=len(set(points[:,0]))>5
            zsampled=len(set(points[:,1]))>5
            tsampled=len(set(points[:,2]))>5
            
            if rsampled and zsampled and tsampled:
                sampled=True
            elif not(rsampled and zsampled):
                x*=2;i+=1
            else:
                x/=1.3;i+=1
            if i>=100:
                sampled=True
                print "max iterations reached, forcing convergence"
            

        if output: print("creating RBF at ", newXt)
        
        def dist_func(X1,X2):
            X=X1[:-1,...]-X2[:-1,...]
            T=(X1[-1,...]-X2[-1,...])*self.vref
            return np.sqrt((X**2).sum(0)+T**2)
        
        RBFu=Rbf(points[:,0], points[:,1], points[:,2], points[:,3],\
                 norm=dist_func, smooth=0, epsilon=5e2)

        RBFw=Rbf(points[:,0], points[:,1], points[:,2], points[:,4],\
                 norm=dist_func, smooth=0, epsilon=5e2)

        if makeHull:
            hull=Delaunay(RBFu.xi.T, qhull_options="QJ10000 "+qhull_opts) #creates a convex hull around the points used to construct the RBF
            RBFu.hull=hull           #a point cant be tested to be interior to the hull with:
            RBFw.hull=hull           #RBF.hull.find_simplex([r,z,t])>=0
            
            RBFu.centre=points.mean(0)
            RBFw.centre=points.mean(0)
        
        self.RBFs.appendleft([RBFu, RBFw])
        
class foo():
    def __init__(self, t=0, dt=3.15576e9, pos=[0,0], largeDx=7.0e7/2**5, store=False):
        self.t=t
        self.ti=int(math.floor(t/dt))
        self.dt=dt
        self.forewards=forewards
        self.pos=pos
        self.dx=largeDx
        self.store=store
        self.history=np.array(self.pos+[self.t]).reshape([1,3])
        self.getRBF(pos[0],pos[1],t)
        
    def getRBF(self,r,z,t):
        Xt=(r,z,t)
        found=False
        if all(v is None for v in self.store.RBFs):
            self.store.addRBF([r,z,t])
            u=self.store.RBFs[0][0];w=self.store.RBFs[0][1]
            found=True
        else:
            for rbf in self.store.RBFs:
                #look for an RBF containing the desired point
                if rbf and rbf[0].hull.find_simplex([r,z,t])>=0:
                    u=rbf[0];w=rbf[1]
                    found=True
                    break
        if not(found) and self.store.fullHull.find_simplex([r,z,t])>=0:
            #point is inside set of all points so create a new RBF for it
            self.store.addRBF([r,z,t])
            u=self.store.RBFs[0][0];w=self.store.RBFs[0][1]        
        elif not(found):
            #point is outside the set of all input points so extrapolate from closest RBF
            mindist=np.inf
            for rbf in self.store.RBFs:
                if rbf:
                    d=abs(rbf[0].hull.plane_distance(np.array([r,z,t]))).min()
                    if d<mindist :
                        mindist=d
                        u=rbf[0];w=rbf[1]
                        
        return u,w
            
    def return_derivs(self,Y,t):
      r,z=Y
      if np.sqrt(r**2+z**2)<1e3:
        theta=math.atan2(r,z)
        return [-1e-5*math.sin(theta),-1e-5*math.cos(theta)]
      else:
        sgnR=2*(r>=0)-1 #sign(r) with sign(0)=1
        r=abs(r)
        i=0;found=False

        u,w=self.getRBF(r,z,t)
        
        diffs=[u(r,z,t)*sgnR*1e-10,\
               w(r,z,t)*1e-10]

        if r<1e2 and abs(diffs[0])<1e-7 and abs(z)>1e3 : diffs[0]=sgnR*1e-7
        if abs(z) <1e3 and r<250e3:
            diffs[1]=diffs[1]*abs(z)/1e3
            diffs[0]=0-abs(diffs[0])
            
        return diffs
    
    def advance(self, n=1, nsteps=100):
        ts=np.linspace(self.t,self.t+self.dt*n,nsteps)
        integrated=odeint(self.return_derivs,self.pos,ts, hmax=self.dt/1000)
        self.pos=integrated[-1,:]
        self.t=ts[-1]
        self.history=np.concatenate([self.history,\
                                     np.concatenate([integrated[1:,:],\
                                                     ts[:-1].reshape([ts.size-1,1])],axis=1)]\
                                    ,axis=0)
        
