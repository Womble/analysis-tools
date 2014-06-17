from numpy import arange, exp,array,empty
from scipy.weave import blitz
from random import random
import genetic as g

def gauss2 (x, params):
     a1,s1,m1,a2,s2,m2,c=params
     func1=exp(-(x-m1)**2/(2*s1*s1))
     func2=exp(-(x-m2)**2/(2*s2*s2))
     func1*=a1
     func2*=a2
     func1+=func2
     func1+=c
     return func1	

def mainfoo(im) :
     shape=im.shape
     X=arange(shape[0])
     res=empty((7,shape[1],shape[2]))
     for i in xrange(shape[1]):
         f=open('log','a')
         f.write(str(i)+'\n')
         f.close()
         print i
         for j in xrange(shape[2]):
             c=sorted(im[:,i,j])[shape[0]/2]
             peak=im[:,i,j].argmax()
             def cons ():
               return array((random()*10,random()*10,random()*10-20+peak, random()*-10,random()*10,random()*20-10+peak,c*(random()*0.2+0.9)))
             def foo(params):
                d=im[:,i,j]
                return ((gauss2(X, params)-d)**2).sum()
             s=g.optimize(foo, cons, satisfactory=0.01, verbose=False)
             s.run()
             res[:,i,j]=s.pool[0][1:]
     return res
