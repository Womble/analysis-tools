from random import random
import numpy as np
from numpy.random import randint
from scipy.optimize import fmin_bfgs

rndint=np.frompyfunc(randint,1,1)

_test_dat_=x=np.array([range(30),8-2*np.exp((-(np.arange(30)-12.)**2)/(2*3*3))])
#poolSize=500
#its=500

def test_cons():
    "sample constructor:constructs a new parameter list, constraints on params and weighting of random vars should go in here. All params in this list are optimised"
    flag=True
    while flag:
        ans=list((random()*100-50,random()*2,random()*30))
        if (ans[-1]-ans[0])>=0: flag=False
    return np.array(ans)

def _evaluate(f,lis):
    """takes a function f and a list of params lis and returns the sum of the squares of the difference between
f(_test_dat_[o,:],lis) and _test_dat_[1,:], ie sees how good a fit f(x,lis) is to the data"""
    ans=0
    for i in xrange(_test_dat_.shape[1]):
        ans+=(f(_test_dat_[0,i], *lis)-_test_dat_[1,i])**2
    return ans

def _foo(x,A,sigma,mu):
    "sample function with 3 mutate-able params A, sigma and mu. first param (here x) will be a 1d numpy array"
    c=8
    return A*np.exp((-(x-mu)**2)/(2*sigma*sigma))+c

def test_f(params):
    "sample fitness function, fits the gaussian funtion _foo to the data in _test_dat_"
    return _evaluate(_foo, params)

def test_cons2 ():
    return [random(), random(),random(), random()]

def test_f2 (X):
    v,x,y,w=X
    pi=np.pi
    r,theta=np.sqrt((x-0.5)**2+(y-0.5)**2),np.arctan2(y-0.25,x-0.25)
    return -1*np.cos(v-0.5)*np.cos(8*w-4)*np.cos(r*20*pi)*np.cos(theta*8)/(1+r+abs(v-0.5)+abs(w-0.5))**2

def _breed_(A,B,muterate=15):
    """takes two [fitness, params] lists and creates a new one by converting the params into uint8-arrays,
twisting them together at a random point and randomly mutating ~1/muterate of the uint16s"""
    origtype=type(A[0])
    aa,bb=A[1:].view(np.uint16),B[1:].view(np.uint16)
    mutations=(rndint(np.ones(len(aa))*muterate)==0).view(np.uint8) * np.array(rndint(np.ones(len(aa))*256), dtype=np.uint16)
    rr=randint(1,len(aa)-1)
    ans=np.concatenate((aa[:rr]+mutations[:rr], bb[rr:]+mutations[rr:])).view(origtype)
    return ans

def _detri_ (x): return int(np.sqrt(x)*np.sqrt(2))

def _hillWalk_(f,A,verbose=False, maxits=None):
    "takes a fitness function (of airity len(A-1)) and a list of [fitness, parameters] and moves it in parameter space up its local gradient"
    new=list(fmin_bfgs(f, A[1:],disp=False, maxiter=maxits))
    return np.concatenate(([f(new)], new))

class optimize ():
    """Class that uses a simple genetic optimiatoin routine to minimise n-dimensional functions avoiding local minima

    -f is a function which takes a list of n parameters and returns a single fitness value; this is the
value being minimised
    -cons is a function with no arguments which returns a new random guess for the parameters (i.e. the
output of this function can be given an f value
    -pool is either the length of the gene pool used in optimisation or the pool itself
    -its is the maximum number of major iterations for each run()
    -tolerance and satisfactory control early finishing of the optimization. If the fractional difference
in parameter space between the current top and the parameter list at poolSize/3 is less than tolerance or
the fitness of the top member is less than satisfactry the optimization ends
    -mutationRate is one in how many bytes will be mutated to a random byte
    -verbose causes output on the state of optimisation to be printed out at various intervals
"""
 #create
    def __init__ (self, f, cons, pool=200, its=200, tolerance=0.00001, satisfactory=1.0/5000, mutationRate=300, hillWalks=4, startWalk=False, finalWalk=True, verbose=False):
        self.verbose=verbose
        self.cons=cons
        if type(pool)==type([]): self.pool=pool
        else:
            self.pool=[]
            for _ in xrange(pool):
                c=self.cons()
                self.pool.append(np.concatenate(([f(c)], c)))
        self.poolSize=len(self.pool)
        self.pool.sort(key=lambda x: x[0])
        self.tri=self.poolSize*(self.poolSize+1)/2
        self.its=its
        self.tol=tolerance
        self.satisfactory=satisfactory
        self.muterate=mutationRate
        self.f=f
        self.walks=hillWalks
        self.sw=startWalk
        self.fw=finalWalk
    
    def tolTest(self):
        top=self.pool[0][1:]
        third=self.pool[self.poolSize/3][1:]
        return (np.sqrt(sum([(x-y)**2 for x,y in zip(top,third)]))/ np.sqrt(sum(top)))

    
    def convTest(self):
        "test if the top value is > satisfactory or fractional difference in parameter space between top answer and 1/3 of the way down the list is less than tol"
        s=self.pool[0][0]<self.satisfactory
        return s

    def run (self):
        "start optimizing"
        i=0
        if self.sw and not(self.convTest()):self.pool[0]=_hillWalk_(self.f, self.pool[0])
        while i<self.its and not(self.convTest()):
            if i%(self.its/50)==0: 
                if self.verbose: print [x[0] for x in self.pool[::self.poolSize/10]], str((i*100)/self.its)+'% '+str(self.tolTest())
            if self.walks and (i+1)%(self.its/self.walks)==0: 
                if self.verbose: print 'hill-walking'
                for j in xrange(self.poolSize/3): #only hill walk every 2nd element to (hopefully prevent convergence on one local maximum)
                    self.pool[j*2]=_hillWalk_(self.f,self.pool[j*2],self.verbose, maxits=10)
                if self.verbose: print 'done hill-walking'
     #breed
            else:
                for j in xrange(3*self.poolSize/4):
                    if self.tolTest<50*self.tol :m=self.muterate/100
                    else: m=self.muterate
                    r1=_detri_(randint(1,self.tri))-1
                    r2=(r1+randint(1,self.poolSize))%self.poolSize
                    x=randint(0,10)
                    if x==0:
                        new=self.cons()
                    else:
                        new=_breed_(self.pool[r1],self.pool[r2], m)
                    self.pool.append(np.concatenate(([self.f(new)], new)))
#            self.pool=self.pool[:3*self.poolSize/4]+self.pool[self.poolSize:] #cut bottom quarter from previous gen after giving them a chance to breed 
            self.pool.sort(key=lambda x: x[0]) #sort
            self.pool=self.pool[:self.poolSize]#cull to poolSize

            i+=1

        if self.verbose: print '100%'
        if self.fw: self.pool[0]=_hillWalk_(self.f,self.pool[0],self.verbose)
        return self.pool[0]
