import numpy as np
from scipy.weave import blitz
from random import random as rnd
from random import normalvariate as nv
import genetic as g
from complete import cube_convolve
import pyfits as P
import sys,os

def gauss2 (x, params):
     a1,s1,m1,a2,s2,m2,c=params
     func1=np.exp(-(x-m1)**2/(2*s1*s1))
     func2=np.exp(-(x-m2)**2/(2*s2*s2))
     func1*=a1
     func2*=a2
     func1+=func2
     func1+=c
     return func1	

def genetic_double_gaussian_fit(im) :
     shape=im.shape
     X=np.arange(float(shape[0]))
     res=np.empty((7,shape[1],shape[2]))
     for i in xrange(shape[1]):
#         f=open('log','a')
#         f.write(str(i)+'\n')
#         f.close()
         print i
         for j in xrange(shape[2]):
             c=sorted(im[:,i,j])[shape[0]/2]
             peak  =im[:,i,j].argmax()
             trough=im[:,i,j].argmin()

             def cons ():
                  if random.randint(0,1)==1:
                       return np.array((rnd()*10, rnd()*10, nv(peak,10),\
	                             rnd()*-10,rnd()*10, nv(peak,10),\
                                     c*nv(1, 0.1)))
                  else :
                       return np.array((rnd()*10, rnd()*10, nv(trough,10),\
                                     rnd()*-10,rnd()*10, nv(trough,10),\
                                     c*nv(1, 0.1)))

             def foo(params):
                d=im[:,i,j]
                return ((gauss2(X, params)-d)**2).sum()

             s=g.optimize(foo, cons, satisfactory=0.01, its=100)
             r=g.optimize(foo, cons, satisfactory=0.01, its=100)
             actualstdout = sys.stdout
             sys.stdout = open(os.devnull,'w')
             s.run()
             r.run()
             final=g.optimize(foo, cons, satisfactory=0.01, pool=s.pool+r.pool)
             sys.stdout.close()
             sys.stdout = actualstdout
             res[:,i,j]=final.pool[0][1:]
     return res

def disc_analysis (path):
    im=cube_convolve(P.getdata(path), 4)
    res=genetic_double_gaussian_fit(im)
    psplit=path.split('.')
    P.writeto(''.join(psplit+['fitted.fits']), res)
    ans=np.empty([4]+res.shape[1:])
    shape=res.shape
    for i in xrange(shape[1]):
        for j in xrange(shape[2]):
             ans[2:3,i,j]=res[2:6:3,i,j]
             if abs(res[0,i,j]) > abs(res[3,i,j])  and abs(res[0,i,j])/abs(res[6,i,j]) > 0.1:
                  ans[0,i,j]=res[2,i,j]
             if abs(res[0,i,j]) <= abs(res[3,i,j]) and abs(res[3,i,j])/abs(res[6,i,j]) > 0.1:
                  ans[1,i,j]=res[5,i,j]
    return ans
