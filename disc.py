import numpy as np
from scipy import optimize
from random import random as rnd, randint as rndint
from random import normalvariate as nv
import genetic as g
from complete import cube_convolve
import pyfits as P
import sys,os
from scipy.weave import blitz
from matplotlib.cm import RdBu_r as cm
import matplotlib
from pylab import plot, imshow, clf, colorbar

cm.set_under((0.0,0.0,0.25))
cm.set_over ((0.25,0.0,0.0))
cm.set_bad  ((0.0,0.25,0.0))

last_plot=[]

def xyplot(lis): 
    x,y=apply(zip, lis)
    plot(x,y)

def gauss2 (x, params):
     """calculate the sum of 2 gaussians and a scalar
x is a scalar or a numpy array
params = amplitude1, sigma1, mean1, amplitude2, sigma2, mean2, baseline"""
     a1,s1,m1,a2,s2,m2,c=params
     func1=np.exp(-(x-m1)**2/(2*s1*s1))
     func2=np.exp(-(x-m2)**2/(2*s2*s2))
     func1*=a1
     func2*=a2
     func1+=func2
     func1+=c
     return func1

def gauss (x, params):
     a,s,m,c=params
     result=np.exp(-(x-m)**2/(2*s*s))
     result*=a
     result+=c
     return result

def almost_eq (a, b, diff=0.1):
     if (a <= b*(1+diff) and b*(1-diff) <= a) or\
        (b <= a*(1+diff) and a*(1-diff) <= b):
           return True
     else :return False

def genetic_double_gaussian_fit(im, verbose=False) :
     "fit double gaussian profiles to the each of the spectra in imcube im\nold, dont use"
     shape=im.shape
     X=np.arange(float(shape[0]))
     res=np.empty((7,shape[1],shape[2]))
     for i in xrange(shape[1]):
#         f=open('log','a')
#         f.write(str(i)+'\n')
#         f.close()
         print i
         if not(verbose):
              actualstdout = sys.stdout
              sys.stdout = open(os.devnull,'w')
         for j in xrange(shape[2]):
             d=im[:,i,j]
             c=sorted(d)[shape[0]/2]
             peak  =d.argmax()
             trough=d.argmin()
             mx,mn=d[peak]-c,d[trough]-c

             def cons ():
                  r=rndint(0,2)
                  if r==0:
                       return np.array((nv(mx,mx/4), rnd()*15, nv(peak,10),\
                                        nv(-mx,mx/4)*-10,rnd()*15, nv(peak,10),\
                                        c*nv(1, 0.1)))
                  elif r==1:
                       return np.array((nv(mx,mx/4)*10, rnd()*15, nv(peak,10),\
                                        nv(mn,mn/4)*-10,rnd()*15, nv(trough,10),\
                                        c*nv(1, 0.1)))
                  else :
                       return np.array((nv(mn,mn/4)*10, rnd()*15, nv(trough,10),\
                                        nv(mn,mn/4)*-10,rnd()*15, nv(trough,10),\
                                        c*nv(1, 0.1)))

             def foo(params):
                  #inlined gauss2 (this is the intermost part of the optimisation)
                  a1,s1,m1,a2,s2,m2,c=params
                  func1=np.exp(-(X-m1)**2/(2*s1*s1))
                  func2=np.exp(-(X-m2)**2/(2*s2*s2))
                  func1*=a1; func2*=a2; func1+=func2; func1+=c; func1-=d; func1*=func1
#                  blitz('func1=func1*a1+func2*a2+c')#not working for some reason
                  return func1.sum()

             s=g.optimize(foo, cons, satisfactory=0.01*251, tolerance=0.001, its=100, hillWalks=2, verbose=verbose)
             r=g.optimize(foo, cons, satisfactory=0.01*251, tolerance=0.001, its=100, hillWalks=2, verbose=verbose)
             s.run()
             r.run()
             final=g.optimize(foo, cons, satisfactory=1.0, tolerance=0.0001, pool=s.pool+r.pool)
             final.run()
             res[:,i,j]=final.pool[0][1:]
         if not(verbose):
              sys.stdout.close()
              sys.stdout = actualstdout
     return res

def double_gaussian_fit(im, verbose=False) :
     "fit double gaussian profiles to the each of the spectra in imcube im"
     shape=im.shape
     X=np.arange(float(shape[0]))
     res=np.empty((7,shape[1],shape[2]))
     for i in xrange(shape[1]):
#         f=open('log','a')
#         f.write(str(i)+'\n')
#         f.close()
         print i
         for j in xrange(shape[2]): 
             d=im[:,i,j]
             c=sorted(d)[shape[0]/2]
             peak  =d.argmax()
             trough=d.argmin()
             mx,mn=d[peak]-c,d[trough]-c

             def cons ():
                  r=rndint(0,2)
                  if r==0:
                       return np.array((nv(mx,mx/4), rnd()*15, nv(peak,10),\
                                        nv(-mx,mx/4)*-10,rnd()*15, nv(peak,10),\
                                        c*nv(1, 0.1)))
                  elif r==1:
                       return np.array((nv(mx,mx/4)*10, rnd()*15, nv(peak,10),\
                                        nv(mn,mn/4)*-10,rnd()*15, nv(trough,10),\
                                        c*nv(1, 0.1)))
                  else :
                       return np.array((nv(mn,mn/4)*10, rnd()*15, nv(trough,10),\
                                        nv(mn,mn/4)*-10,rnd()*15, nv(trough,10),\
                                        c*nv(1, 0.1)))

             def err(params,x,y):
                  a1,s1,m1,a2,s2,m2=params
                  func1=np.exp(-(x-m1)**2/(2*s1*s1))
                  func2=np.exp(-(x-m2)**2/(2*s2*s2))
                  func1*=a1; func2*=a2; func1+=func2; func1+=c; func1-=y;
                  return func1

             v,success=optimize.leastsq(err, (d[peak]-c, 10, peak , d[trough]-d[peak], 10, trough), args=(X,d))
             if success!=1:
                  def foo(params):
                      #inlined gauss2 (this is the intermost part of the optimisation)
                       a1,s1,m1,a2,s2,m2,c=params
                       func1=np.exp(-(X-m1)**2/(2*s1*s1))
                       func2=np.exp(-(X-m2)**2/(2*s2*s2))
                       func1*=a1; func2*=a2; func1+=func2; func1+=c; func1-=d; func1*=func1
                       return func1.sum()
    
                  if not(verbose):
                       actualstdout = sys.stdout
                       sys.stdout = open(os.devnull,'w')
                  #if leastsq doesnt give a goot fit then try to fit with a genetic alg
                  s=g.optimize(foo, cons, satisfactory=2.51, tolerance=0.001, its=100, hillWalks=2, verbose=verbose, startWalk=True)
                  r=g.optimize(foo, cons, satisfactory=2.51, tolerance=0.001, its=100, hillWalks=2, verbose=verbose, startWalk=True)
                  s.run()
                  r.run()
                  final=g.optimize(foo, cons, satisfactory=1.0, tolerance=0.0001, pool=s.pool+r.pool, finalWalk=True)
                  final.run()
                  if not(verbose):
                       sys.stdout.close()
                       sys.stdout = actualstdout
                  if final.pool[0][0]<=1.5:
                       res[:,i,j]=final.pool[0][1:]
                  else:
                       res[:,i,j]=np.nan
             else :
                  for k in xrange(len (v)): res[k,i,j]=v[k]
                  res[ -1,i,j]=c
     return res

def disc_analysis (path, threshold=0.1, sigma=4, verbose=False):
    """read in and peform a fitting on the fits file at path containing an imcube and write out
the results to a new fits file called Xfitted.fits of shape [7,x,y] where the 7 planes of the 
output are fitting parameters amplitude sigma and mean of the 2 gaussians and the baseline level"""
    im=cube_convolve(P.getdata(path), sigma)
    psplit=path.split('.')
    print(psplit)
#    data=[i for i in xrange(im.shape[1]) if im[:,i,:].sum()/im.shape[0]/im.shape[-1] > np.power(10, -1.5)]
    mn,mx=70,130#min(data),max(data)
#    res=scipy_double_gaussian_fit(im, verbose=verbose)
    res=double_gaussian_fit(im[:,mn:mx,:], verbose=verbose)
    ans=np.ones(res.shape, dtype=np.float)
    ans+=np.nan
    shape=res.shape
    for i in xrange(shape[1]):
        for j in xrange(shape[2]):
             ele=res[:,i,j]
             if abs(ele[0])>=abs(ele[3]) and not(1.6<abs(ele[1])<2.1 and 122<abs(ele[2])<128 and ele[0]<-4):
                  ans[:,i,j]=ele
             else:
                  ans[:,i,j]=(ele[3],ele[4],ele[5],ele[0],ele[1],ele[2],ele[6])
    try :P.writeto(''.join(psplit[:-1]+['fitted.fits']), ans)
    except IOError as msg: print msg 
    return ans



#displays
def pv (im, contSub=False, spatRes=0.625, velRes=0.075):
    global last_plot
#    fig=figure()
#    ax=matplotlib.image.NonUniformImage(fig, interpolation='nearest',  extent=[-velRes*p.shape[0]/2.0, velRes*p.shape[0]/2.0, -spatRes*p.shape[1]/2.0, spatRes*p.shape[1]/2.0])

    if (len (im.shape)==3):
        plane=float(raw_input("plane number for pv slice: "))
        axis=float(raw_input("axis number (y=1, x=2): "))
        global last_plot
        if   axis==1: p=im[:,plane,:].copy()
        elif axis==2: p=im[:,:,plane].copy()
        if contSub :
            for i in xrange(p.shape[1]):
                c=sorted(p[:,i])[p.shape[0]/2]
                p[:,i]-=c
                
        srt=sorted(p.flatten())
        l=len(srt)
        clf()
        imshow(p.transpose(), cmap=cm, vmin=srt[l/10], vmax=srt[49*l/50], interpolation='nearest')
        colorbar()
        last_plot=p.transpose()

    elif (len (im.shape)==2):
        p=im.copy()
        if contSub :
            for i in xrange(p.shape[1]):
                c=sorted(p[:,i])[p.shape[0]/2]
                p[:,i]-=c

        srt=sorted(p.flatten())
        l=len(srt)
        clf()
        imshow(p.transpose(), cmap=cm, vmin=srt[l/10], vmax=srt[99*l/100], interpolation='nearest')
        colorbar()
        last_plot=p.transpose()

def mom1map (fits, sensitivity=0.1):
     global last_plot
     f=fits[2,:,:].copy()
     mask=(abs(fits[0,...]/fits[-1,...]))<sensitivity
     f[mask]=np.nan
     clf()
     imshow(f, cmap=cm, vmin=0, vmax=250)
     colorbar()
     last_plot=f
