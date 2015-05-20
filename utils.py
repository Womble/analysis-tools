import numpy as np
import pylab as pyl
from math import pi
import math
import os
from scipy import ndimage
from scipy import constants as cns
from scipy.interpolate import griddata
from numpy.fft import fft2,ifft2, irfft2, rfft2

class interpArray(np.ndarray):
    """interpolateably nd array: fractional indices give linear interpolations between points 
For an N dimensional array uses (1+N) times as much memory
Cannot use fractional indexing in slices, only for single points"""
    def __new__(cls, inpArray):
        obj = np.asarray(inpArray).view(cls)
        obj.grads = np.gradient(inpArray)
        return obj

    def __getitem__(self, index):
        try:
            if float in map(type,index):
                if slice in map(type,index): 
                    raise TypeError()
                Lind=tuple(map(int,index))
                lower=self.view(np.ndarray).__getitem__(Lind)
                deltas=map(lambda x: x%1, index)
                for i,d in enumerate(deltas):
                    lower+=deltas[i]*self.grads[i][Lind]
                return lower
            else:
                raise TypeError()
        except TypeError:
            return self.view(np.ndarray).__getitem__(index)

def div (V):
    "assumes vector field will be a rank n+1 array with shape [n,x_1,x_2,...,x_n]"
    gv=np.gradient(V)
    x=gv[0,0,...]
    for i in xrange(1,len(V.shape)):
        x+=gv[i,i,...]
    return x

def laplacian (X):
    gx=np.gradient(X)
    shape=X.shape
    if len(shape)==1:
        return np.gradient(gx)
    else:
        out=np.zeros_like(X)
        for i in xrange(len(shape)):
            out+=np.gradient(gx[i])[i]
        return out

def imshowWslices(data, cbar=False, yslice=0.5, xslice=0.5, **args):
    shape=data.shape
    if not(len(shape)==2): raise valueError
    xslice=int(xslice*shape[0])
    yslice=int(yslice*shape[1])
    # Define the locations for the axes
    left, width = 0.12, 0.55
    bottom, height = 0.12, 0.55
    bottom_h = left_h = left+width+0.02

    # Set up the geometry of the three plots
    rect_temperature = [left, bottom, width, height] # dimensions of temp plot
    rect_histx = [left, bottom_h, width, 0.25] # dimensions of x-histogram
    rect_histy = [left_h, bottom, 0.25, height] # dimensions of y-histogram

    # Set up the size of the figure
    fig = pyl.figure(1, figsize=(9.5,9))

    # Make the three plots
    axTemperature = pyl.axes(rect_temperature) # temperature plot
    cax = axTemperature.imshow(data, interpolation='nearest', origin='image', **args)

    axHistx = pyl.axes(rect_histx) # x histogram
    axHisty = pyl.axes(rect_histy) # y histogram
    if cbar : pyl.colorbar(cax, ax=axTemperature, orientation='horizontal')
    
    nullfmt = pyl.NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

#    axHisty.set_ylim([0,shape[1]])
#    axHistx.set_xlim([0,shape[0]])

    axHistx.plot(np.arange(shape[1]), data[xslice,:], linestyle='steps-mid')
    axHisty.plot(data[:,yslice], np.arange(shape[0]), linestyle='steps-mid')

    axHisty.set_yticklabels(['%.1e'%x for x in axHisty.get_xticks()])
    axHistx.set_xticklabels(['%.1e'%x for x in axHisty.get_yticks()])

    pyl.draw()
    pyl.show()

def arr2h (data, name, opath='temp.h'):
    S=np.array(data.shape)
    SS=np.array([S[i:].prod() for i in xrange(-S.size+1,0)]) #size of row, plane, cube, hypercube etc
    c=0
    of=open(opath,'w')
    print SS
    of.write('const static double '+name+''.join(tuple([str([x]) for x in data.shape]))+' = {')
    for x in data.flat:
        c+=1
        ns=c%SS
        for n in ns: #if we are are at the start of a block open a brace
            if n==1: of.write('{')
        of.write("%.3e"%x) #write the element
        if (ns!=0).all()or c==0: of.write(', ') #if we're not at the end of a block write a comma
        elif c!=0: #if we are at the end of a block close a brace and maybe write a new line
            for (i,n) in enumerate(ns[::-1]):
                if n==0: 
                    of.write('}')
                    if i==ns.size-1: 
                        if c!=data.size: of.write(',\n')
                        else: of.write('\n')
                    else:
                        if ns[::-1][i+1]==0: 
                            of.write('\n')
                        else: of.write(',\n')
    of.write('};\n')#close the array  

def getfits (path='./', *iname):
    return [x for x in os.walk(path).next()[2] if len(x)>=5 and (x[-5:]=='.fits' or x[-8:]=='.fits.gz') and all([y in x for y in ['']+list(iname)])]

def almost_eq (a, b, diff=1e-6):
    "for comparing floats, evaluates to 10 if the factional difference between a and b is less than diff"
    return np.logical_or(a==b, (np.logical_or((abs((a-b)/a) <= diff) , (abs((a-b)/b) <= diff))))

def next_pow2 (x):
    return np.power(2, int(np.ceil(np.log2(x))))

def JypPx2Temp (J,  freq, cellWidth):
    """X is intensity in janskys per cell
freq is the frequency of the light in Hz
cellWidth is the size of one cell in arcsec"""
    lamb=cns.speed_of_light/freq*1000 #lambda in mm
    return 13.6 * (lamb/(cellWidth*pi/4))**2 * J

    T = 13.6 * (lamb/(cellWidth*pi/4))**2 * J

def Temp2JypPx (T, freq, cellWidth):
    """T is brightness temperature in Kelvin
freq is the frequency of the light in Hz
cellWidth is the size of one cell in arcsec"""
    lamb=cns.speed_of_light/freq*1000
    return T / (13.6 * (lamb/(cellWidth*pi/4))**2)

def stripStokes (im):
    "casa annoyingly creates fits files in a shape [v,1,x,x] rather than [v,x,x] so this removes that"
    return im.reshape(list((im.shape[0],))+list(im.shape[2:]))

def FWHM2sigma (FWHM): return FWHM/(2*np.sqrt(2*np.log(2)))
    
def garray(shape, sigma, normalise=1):
    midi,midj=shape[0]/2,shape[1]/2
    X,Y=np.mgrid[-midi:midi:shape[0]*1j, -midj:midj:shape[1]*1j]
    r=np.sqrt(X*X+Y*Y)
    out= np.exp(-1*r**2 / (2*sigma*sigma))
    if normalise :out=out/out.sum()
    return out
#exp(-1*(i-midi)**2 / (2*sigma*sigma))*exp(-1*(j-midj)**2 / (2*sigma*sigma))/(np.sqrt(2*pi*sigma**2))  

def beam_convolve(arr, sigma, normalise=1):
    "convoles a 2D image with a gaussian profile with sigma in px"
    if len(arr.shape)!=2 or 3*sigma > max(arr.shape): raise ValueError ("arr is not 2d or beam is too wide")
    else: 
        shape=arr.shape
        gauss_mask=garray(shape,sigma,normalise)
        s=[y*2 for y in gauss_mask.shape]
        ftg=rfft2(gauss_mask, s)
        return irfft2(rfft2(arr,s)*ftg)  

def cube_convolve(imcube, sigma, inplace=0, normalise=1):
    "performs a convolution with a gaussian beam of width sigma on each yz plane of the cube"
    if not(inplace) : imcube=imcube.copy()
    shape=imcube.shape[1:]
    if len(shape)!=2:
        raise ValueError ("cube is not a cube")
    gauss_mask=garray(shape,sigma, normalise=normalise)
    s=[next_pow2(y*2+1) for y in gauss_mask.shape]
    ftg=fft2(gauss_mask, s).reshape(s)
    for i in xrange(imcube.shape[0]):
        imcube[i,...]=np.real(ifft2(fft2(imcube[i,...],s)*ftg)[shape[0]/2:3*shape[0]/2, shape[1]/2:3*shape[1]/2])
    return imcube

def convolve_3dg(arr1, sigma):
    s=arr1.shape
    arrs=np.mgrid[0:s[0],0:s[1],0:s[2]]
    arr2=np.exp( -((arrs**2).sum(0)/float(sigma)**2) )
    arr2/=arr2.sum()
    return np.fft.irfftn(np.fft.rfftn(arr1)*np.fft.rfftn(arr2))

def convolve (arr1, arr2):
    "convolves 2 arrays together with fft, arrays will be zero padded to equal size"
    if max(len(arr1.shape), len(arr2.shape)) > 2: raise ValueError("only dealing with 2d convolves here thankyou")
    s=(int(max(arr1.shape[0],arr2.shape[0])*1.5),int(max(arr1.shape[1],arr2.shape[1])*1.5))
    return irfft2(rfft2(arr1,s)*rfft2(arr2,s))

def strech_arr(arr, axis=0, factor=2, conserveSum=False):
    "increase the dimensions of array arr by factor along axis, if conserveSum is True then the array if divided by factor so that its sum remains the same"
    axes=[x for x in arr.shape[:axis]]+[int(arr.shape[axis]*factor)]+[x for x in arr.shape[axis+1:]]
    new=np.empty(axes)
    l=len(axes)
    for i in xrange(arr.shape[axis]):
        sl1=[slice(0,x) for x in arr.shape[:axis]]+[slice(i,i+1)]
        sl2=[slice(0,x) for x in arr.shape[:axis]]+[slice(i*factor,(i+1)*factor)]
        try:
            sl1+=[slice(0,axes[x]) for x in arr.shape[axis+1:]]
        except:
            None
        try:
            sl2+=[slice(0,axes[x]) for x in arr.shape[axis+1:]]
        except:
            None
        new[sl2]=arr[sl1]
    if conserveSum : new/=factor
    return new    

def degrade_arr (arr, axis=0, factor=2, conserveSum=False):
    "reduce the dimensions of array arr by factor along axis, if conserveSum is True then the array if multiplied by factor so that its sum remains the same"
    s=list(arr.shape)
    s[axis]=int(s[axis]//factor)
    new=np.zeros(s, dtype=arr.dtype)
    for i in xrange(s[axis]):
        sl1=[slice(_) for _ in s]
        sl1[axis]=i
        sl2=[slice(_) for _ in s]
        sl2[axis]=slice(int(i*factor),int((i+1)*factor))
        new[sl1]=arr[sl2].mean(axis)
    if conserveSum : new*=factor
    return new

def Rx(theta):
    "returns the matrix for a rotation of theta radians about x"
    ct=math.cos(theta)
    st=math.sqrt(1-ct*ct)
    return np.matrix([[1,0,0],[0,ct,-st],[0,st,ct]])

def Ry(theta):
    "returns the matrix for a rotation of theta radians about y"
    ct=math.cos(theta)
    st=math.sqrt(1-ct*ct)
    return np.matrix([[ct,0,st],[0,1,0],[-st,0,ct]])

def Rz(theta):
    "returns the matrix for a rotation of theta radians about z"
    ct=math.cos(theta)
    st=math.sqrt(1-ct*ct)
    return np.matrix([[ct,-st,0],[st,ct,0],[0,0,1]])

def cartesian2polar (grid):
    rmax=max(grid.shape)
    r,theta=np.mgrid[0:rmax:rmax*1j, 0:pi/2:rmax*1j]
    out=np.zeros_like(r)
    ndimage.map_coordinates(grid, [r*np.cos(theta),r*np.sin(theta)], output=out, mode='nearest', order=1)
    return out

def lin2loglog (grid, low=-3, shape=[], order=3):
    xmax,ymax=grid.shape
    if not(shape): shape=grid.shape
    oxmax,oymax=shape
    out=np.zeros(shape)
    x,y=np.mgrid[low:0:oxmax*1j, low:0:oymax*1j]
    ndimage.map_coordinates(grid, [10**x*xmax,10**y*ymax], output=out, mode='nearest', order=order)
    return out

def polar2cartesian (grid):
    "converts and interpolates a 2D polar (r,phi) grid to a cartesian  one"
    X,Y=np.mgrid[0.0:grid.shape[0], 0:grid.shape[0]]
    out=np.zeros_like(X)
    ndimage.map_coordinates(grid, [np.sqrt(X*X+Y*Y),np.arctan2(Y,X)/(pi/2)*grid.shape[1]], output=out, mode='nearest', order=1)
    return out

def cartesian2cylindical (grid, z=0):
    """converts and interpolates a 3D cartesian grid to a cylindrical one
z is the number of the z axis, defualts to 0 (first)"""
    if not(len(grid.shape)==3): raise IndexError 
    if z==0: return np.array([x for x in  pp.pmap(cartesian2polar, [grid[i,...] for i in xrange(grid.shape[0])], limit=6)]) 
    if z==1: return np.array([x for x in  pp.pmap(cartesian2polar, [grid[:,i,:] for i in xrange(grid.shape[1])], limit=6)]) 
    if z==2: return np.array([x for x in  pp.pmap(cartesian2polar, [grid[...,i] for i in xrange(grid.shape[2])], limit=6)]) 

def polarSlice2Cube (Slice, zaxis=1, size=None, mode='nearest'):
    s=max(Slice.shape)
    if not(size):
        size=s
    g=x,y,z=np.mgrid[-s:s:size*1j,-s:s:size*1j,-s:s:size*1j]
    R=np.sqrt(x*x+y*y)
    if zaxis==1: return ndimage.map_coordinates(Slice, (R,abs(z)), mode=mode, order=1)
    else       : return ndimage.map_coordinates(Slice, (abs(z),R), mode=mode, order=1)

_quadadd= lambda x: np.sqrt((x*x).sum())
def quad_add (array, axis=0):
    return np.apply_along_axis(_quadadd, axis, array)

def niceim (arr, cutFrac=0.001, log=False, centreZero=False, **kwargs):
    if log : arr=np.log10(arr)
    a=sorted(arr.flat)
    prod=np.product(arr.shape)
    kwargs['vmin']=a[int(prod*cutFrac)]
    kwargs['vmax']=a[int(prod*(1-cutFrac))]
    if centreZero:
        lim=max(abs(kwargs['vmin']), abs(kwargs['vmax']))
        kwargs['vmin']=-lim
        kwargs['vmax']= lim
    pyl.imshow(arr, **kwargs)

class memoize():
    "a class for creating memoized instances of functions"
    def __init__(self, function, dict={}):
        self.f=function
	if dict: self.d=d
        else   : self.d={}
    def __call__(self, *args, **kargs):
        key = (args, frozenset(sorted(kargs.items())))
        try:
            return self.d[key]
        except KeyError:
            ret=self.f(*args, **kargs)
            self.d[key]=ret
            return ret

def funcDiff (f, n=0):
    "provides the partial differential of f with respect to the nth parameter"
    def fprime (*args):
        eps=np.abs(np.finfo(float).eps)          #sqrt(machine epsilon)*x gives a good balance
        h=max(10*eps, np.sqrt(eps)*abs(args[n])) #between small delta and rounding errors
        l1=list(args[:n])+[args[n]-h]+list(args[n+1:])
        f1=f(*l1)
        l2=list(args[:n])+[args[n]+h]+list(args[n+1:])
        f2=f(*l2)
        return (f2-f1)/(2*h)
    return fprime

def compileClib(lib, Oflag='-O3'):
    if '/' in lib and not(lib[:2]=='./' and not('/' in lib[2:])):
        os.system('cp '+lib+' ./')
    if not(lib[-2:] == '.c' or lib[-4:]=='.cpp'):
        raise IndexError('"'+lib+'" doesnt look like a .c or .cpp file, rename it')
    if lib[-1:] == 'c': lname=lib[:-2]
    else              : lname=lib[:-4]
    os.system('gcc -fPIC -shared '+Oflag+' -o '+lname+'.so '+lib)
