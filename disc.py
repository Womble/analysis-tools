import numpy as np
from scipy import optimize
from random import random as rnd, randint as rndint, normalvariate as nv
import genetic as g
from complete import cube_convolve
import pyfits as P
import sys,os
from complete import cube_convolve
from matplotlib.cm import RdBu_r as cm
import matplotlib
from pylab import plot, imshow, clf, colorbar, figure, draw ,savefig
import pprocess

cm.set_under((0.0,0.0,0.35))
cm.set_over ((0.35,0.0,0.0))
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

def almost_eq (a, b, diff=0.01):
    if type(a)==type(b)==np.ndarray:
        return (abs(a-b) <= diff) | (abs((a-b)/a) <= diff) 
    else:
        return abs(a-b) <= diff   or abs((a-b)/a) <= diff

def pvfit (arr, f, init_args_guess):
    """arr: a 2D array 
f: a function taking at least 2 parameters, the first of which is the x-cord eg f=lambda x,m,c: m*x+c for straight line fit
init_args_guess: a list of initial guess for the arfs in f eg [1,0] for f=1*x+0

fits the [1:] params of f to minimise the distance between function and pixel weighted by pixel value
"""
    pos=np.mgrid[0:arr.shape[0],0:arr.shape[1]]

    def func(args):
        f2=lambda x: f(x, *args)
        return ((f2(pos[0,...])-pos[1,...])**2 * arr).sum()        

    return optimize.fmin_bfgs(func, init_args_guess, maxiter=100)



def decorate_image(F, pixels=201,imres=0.005):
    F.axes[0].set_xticks([x*pixels for x in [0,0.25,0.5,0.75,1]])
    F.axes[0].set_xticklabels([round(10000*(x-0.5)*pixels*imres)/10000. for x in [0,0.25,0.5,0.75,1]])
    F.axes[0].set_yticks([x*pixels for x in [0,0.25,0.5,0.75,1]])
    F.axes[0].set_yticklabels([round(10000*(x-0.5)*pixels*imres)/10000. for x in [0,0.25,0.5,0.75,1]])
    F.axes[0].set_xlabel('x "')
    F.axes[0].set_ylabel('y "')

def decorate_pv (F, pixels=201, imres=0.005, nchan=251, velres=100):
    F.axes[0].set_xticks([x*nchan for x in [0,0.25,0.5,0.75,1]])
    F.axes[0].set_xticklabels([round(10000*(x-0.5)*nchan*velres/1000.)/10000. for x in [0,0.25,0.5,0.75,1]])
    F.axes[0].set_yticks([x*pixels for x in [0,0.25,0.5,0.75,1]])
    F.axes[0].set_yticklabels([round(10000*(x-0.5)*pixels*imres)/10000. for x in [0,0.25,0.5,0.75,1]])
    F.axes[0].set_ylabel('x "')
    F.axes[0].set_xlabel('v km/s')

def produce_figs (path, name=''):
    if name=='': name=''.join(path.split('.')[:-1])
    hdr=P.getheader(path)
    im=P.getdata(path)
    cube_convolve(im,1.0)
    nchan=hdr.get('NAXIS3')
    vcentre=nchan/2
    pixels=hdr.get('NAXIS2')
    pcentre=pixels/2
    imres=round(hdr.get('CDELT2')*3600*10000)/10000
    velres=hdr.get('CDELT3')
    im_cs=im-im[-1,...]
    
    mom0=im_cs[:round(vcentre-500/velres),...].sum(0)*velres/1000.+im_cs[round(vcentre+500/velres):,...].sum(0)*velres/1000. #messing about with centre is to avoid evelope emission/absorbtion, needs to talk to paola about whether or not this is kosher
    extent=sorted(abs(mom0).flat)[mom0.size*999/1000]
    F=figure();clf();imshow(mom0, interpolation='nearest', cmap=cm , vmin=-extent, vmax=extent, origin='image');c=colorbar();c.set_label('K km/s')
    decorate_image(F, pixels,imres)
    draw()
    savefig(name+'_contSub.png')
    
    cbar=pv(im_cs[:,pcentre,:],contSub=False, spatRes=imres, velRes=velres/1000., cutFrac=0.003)
    cbar.set_label('K')
    decorate_pv(F, pixels,imres, nchan,velres)
    draw()
    savefig(name+'_PV_centre.png')
    
    im_cs[round(vcentre-500/velres):round(vcentre+500/velres),...]=(im_cs[round(vcentre-500/velres),...]+im_cs[round(vcentre+500/velres),...])/2
    mask=abs(mom0)>abs(mom0).max()/1000
    cbar=mom1map(im_cs*mask, contSub=False, velwidth=velres/1000.)
    cbar.set_label('km/s')
    decorate_image(F, pixels,imres)
    draw()
    savefig(name+'_mom1.png')

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

def double_gaussian_fit_wCentral(im, verbose=False) :
     "fit double gaussian profiles to the each of the spectra in imcube im"
     shape=im.shape
     X=np.arange(float(shape[0]))
     res=np.empty((7,shape[1],shape[2]))
     temp=np.empty(X.shape)
     for i in xrange(shape[1]):
#         f=open('log','a')
#         f.write(str(i)+'\n')
#         f.close()
#     def blah (i):
         print i
         for j in xrange(shape[2]): 
             d=im[:,i,j]
             c=sorted(d)[shape[0]/2]
             temp=d.copy()
             centre=shape[0]/2-8,shape[0]/2+9
             temp[centre[0]:centre[1]]=[temp[centre[0]]*(1-x/16.)+temp[centre[1]]*x/16. for x in xrange (17)]
             peak  =temp.argmax()
             trough=temp.argmin()
             mx,mn=temp[peak]-c,temp[trough]-c



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

             v,success=optimize.leastsq(err, (d[peak]-c, 10, peak , d[trough]-d[peak], 10, trough), args=(X,temp))
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

             if verbose :
                 plot(d)
                 plot (temp)
                 plot (gauss2(np.arange(251.), res[:,i,j]))
                       
     return res

def disc_analysis (path, threshold=0.1, sigma=4, verbose=False):
    """read in and peform a fitting on the fits file at path containing an imcube and write out
the results to a new fits file called Xfitted.fits of shape [7,x,y] where the 7 planes of the 
output are fitting parameters amplitude sigma and mean of the 2 gaussians and the baseline level"""
    im=cube_convolve(P.getdata(path), sigma)
    print('starting: '+path)
#    data=[i for i in xrange(im.shape[1]) if im[:,i,:].sum()/im.shape[0]/im.shape[-1] > np.power(10, -1.5)]
    mn,mx=70,130#min(data),max(data)
#    res=scipy_double_gaussian_fit(im, verbose=verbose)
    res=double_gaussian_fit_wCentral(im[:,mn:mx,:], verbose=verbose)
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
    try :
        print('writting: '+path)
        P.writeto(path[:-5]+'_fitted.fits', ans)
    except IOError as msg: print msg 
    return ans

def continuumSubtract(spec, method=lambda x: sorted(x)[int(len(x)/2)]):
    return spec-method(spec)


#displays
def pv (im, contSub=True, spatRes=0.625, velRes=0.075, cutFrac=0.01, fractional=False, cmap=cm, **kwargs):
    global last_plot

    def contfind(spec):
        f=lambda x:x[0]
	return f(spec)
#        f=lambda x: np.sqrt(np.abs(spec-x)).sum()
#        return optimize.fmin(f, spec[0], disp=0)[0]

    if (len (im.shape)==3):
        plane=float(raw_input("plane number for pv slice: "))
        axis=float(raw_input("axis number (y=1, x=2): "))
        global last_plot
        if   axis==1: p=im[:,plane,:].copy()
        elif axis==2: p=im[:,:,plane].copy()
        if contSub :
            for i in xrange(p.shape[1]): p[:,i]=continuumSubtract(p[:,i], contfind)
        clf()
        if fractional:
            #assuming here for fractional that channel 0 is all cont emission
            if   axis==1: q=im[0,plane,:].reshape((p.shape[1],1))
            elif axis==2: q=im[0,:,plane].reshape((p.shape[1],1))
            extent=abs(p.transpose()/q).max()
            imshow(p.transpose()/q, cmap=cmap, vmax=extent, vmin=-extent, interpolation='nearest', origin='image')
            last_plot=p.transpose()/q
            return colorbar()

        else:
            srt=sorted(list(p[:120,:].flatten())+list(p[130:,:].flatten()))
            l=len(srt)
            extent=max(abs(srt[int(l*cutFrac)]),abs(srt[int((1-cutFrac)*l)]))
            imshow(p.transpose(),   cmap=cmap, vmax=extent, vmin=-extent, interpolation='nearest', origin='image')
            last_plot=p.transpose()
            return colorbar()


    elif (len (im.shape)==2):
        p=im.copy()
        if contSub :
            for i in xrange(p.shape[1]): p[:,i]=continuumSubtract(p[:,i], contfind)

        srt=sorted(list(p[:120,:].flatten())+list(p[130:,:].flatten()))
        l=len(srt)
        extent=max(abs(srt[int(l*cutFrac)]),abs(srt[int((1-cutFrac)*l)]))
        clf()
        imshow(p.transpose(), cmap=cmap, vmin=-extent, vmax=extent, interpolation='nearest', origin='image', **kwargs)
        last_plot=p.transpose()
        return colorbar()


def mom0map (im, velwidth=0.075, **kwargs):
     global last_plot
     p=im.sum(0)*velwidth
     clf();imshow(p, cmap=cm, vmax=p.max()*1.01, interpolation='nearest', origin='image', **kwargs);cbar=colorbar()
     last_plot=p
     cbar.set_label('K km/s')
     draw()
     return cbar

def mom1map (im, velwidth=0.075, contSub=True, **kwargs):
     global last_plot
     v=np.linspace(-im.shape[0]/2.0*velwidth, im.shape[0]/2.0*velwidth, im.shape[0]).reshape((im.shape[0],1,1))
     mcont=im.copy()
     if contSub:
         for i,j in [(x,y) for x in xrange(im.shape[1]) for y in xrange(im.shape[2])]:
             mcont[:,i,j]=continuumSubtract(mcont[:,i,j])
     x=mcont.copy()
     mcont=(abs(mcont)*v).sum(0)/abs(mcont).sum(0)
     mcont[np.isnan(mcont)]=0
     last_plot=mcont
     r=max([abs(x) for x in mcont.flatten() if x<1e99])
     clf();imshow(mcont, cmap=cm, vmax=r, vmin=-r, interpolation='nearest', origin='image', **kwargs);
     cbar=colorbar()
     cbar.set_label('km/s')
     draw()
     return cbar

def mom2map (im, velwidth=75, **kwargs):
     global last_plot
     v=np.linspace(-im.shape[0]/2.0*velwidth, im.shape[0]/2.0*velwidth, im.shape[0]).reshape((im.shape[0],1,1))
     mcont=im.copy()
     for i,j in [(x,y) for x in xrange(im.shape[1]) for y in xrange(im.shape[2])]:
         mcont[:,i,j]=continuumSubtract(mcont[:,i,j])
     mom1map(mcont, velwidth,False)
     m1=last_plot
     mcont=np.sqrt((abs(mcont)*(v-m1)**2).sum(0)/abs(mcont).sum(0))
     mcont[mcont==np.nan]=0
     last_plot=mcont
     clf();imshow(mcont, cmap=cm, vmin=0, interpolation='nearest', origin='image', **kwargs);cbar=colorbar()
     cbar.set_label('km/s')
     draw()
     return cbar
