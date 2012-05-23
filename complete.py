import numpy,pyfits,os,tarfile,Gnuplot,pprocess
from math import exp,log,sqrt,pi,ceil
from numpy.fft import rfft2,irfft2

g=Gnuplot.Gnuplot()
g('set style data histeps')

def getfits (path='./'):
    return [x for x in os.walk(path).next()[2] if len(x)>=5 and x[-5:]=='.fits']

def garray(shape, sigma):
    a=numpy.empty(shape, dtype=float)
    midi,midj=shape[0]/2,shape[1]/2
    for i in xrange(shape[0]):
        for j in xrange(shape[1]):
            a[i,j]=exp(-1*sqrt((i-midi)**2+(j-midj)**2)**2 / (2*sigma*sigma))/(sigma*sqrt(2*pi))**2
#exp(-1*(i-midi)**2 / (2*sigma*sigma))*exp(-1*(j-midj)**2 / (2*sigma*sigma))/(sqrt(2*pi*sigma**2))    
    return a

def FWHM2sigma (FWHM): return FWHM/(2*sqrt(2*log(2)))

def beam_convolve(arr, sigma):
    "convoles a 2D image with a gaussian profile with sigma in px"
    if len(arr.shape)!=2 or 3*sigma > max(arr.shape): raise ValueError ("arr is not 2d or beam is too wide")
    else: 
        shape=arr.shape
        gauss_mask=garray(shape,sigma)
        s=[y*2 for y in gauss_mask.shape]
        ftg=rfft2(gauss_mask, s)
        return irfft2(rfft2(arr,s)*ftg)

def cube_convolve(imcube, sigma):
    "performs a convolution with a gaussian beam of width sigma on each xy plane of the cube"
    shape=imcube.shape[1:]
    if len(shape)!=2 or 3*sigma > max(shape): 
        raise ValueError ("cube is not a cube or beam is too wide")
    gauss_mask=garray(shape,sigma)
    s=[y*2+1 for y in gauss_mask.shape]
    ftg=rfft2(gauss_mask, s)
    ref=(ceil(shape[0]/2.0),ceil(shape[0]/2.0)+shape[0],ceil(shape[1]/2.0),ceil(shape[1]/2.0)+shape[1])
#    cwrap=lambda x: convolve(x[0],x[1])
#    pprocess.pmap(cwrap, [(ftg,imcube[i,...]) for i in xrange(imcube.shape[0])])
    for i in xrange(imcube.shape[0]):
        imcube[i,:,:]=irfft2(rfft2(imcube[i,:,:],s)*ftg)[ref[0]:ref[1],ref[2]:ref[3]]
    return imcube

def spec_at(imcube, pos, chanwidth=10):
    x,y=pos
    shape=imcube.shape
    bandwidth=shape[0]*chanwidth
    if len(shape)!=3: raise ValueError("imcube must be a cube")
    dat=[(i*chanwidth-bandwidth/2,imcube[i,y,x]) for i in xrange(shape[0])]
    g.plot(dat)
    return dat

def plotfits (path, sigma, chanwidth, replot=False):
    if replot : g.replot(line(pyfits.getdata(path),sigma,chanwidth))
    else      : g.plot  (line(pyfits.getdata(path),sigma,chanwidth))

def lineold (imcube, sigma, chanwidth=10):
    "produces a spectrum by convolving each slice of imcube with a gaussian of width sigma and returning the value of the central pixel for each slice"
    shape=imcube.shape
    bandwidth=shape[0]*chanwidth
    if len(shape)!=3: raise ValueError("imcube must be a cube")
    gauss_mask=garray(shape[1:],sigma)
    s=[y*2 for y in gauss_mask.shape]
    ftg=rfft2(gauss_mask, s)
    
    return [(i*chanwidth-bandwidth/2,irfft2(rfft2(imcube[i,:,:],s)*ftg)[s[0]/2,s[1]/2])
            for i in xrange(shape[0])]

def spec_at(imcube, pos, chanwidth=10):
    "produces an (unconvolved) spectrum at the x,y pixel given"
    bandwidth=imcube.shape[0]*chanwidth
    x,y=pos
    return [(i*chanwidth-bandwidth/2, imcube[i,x,y]) for i in xrange(imcube.shape[0])]

def line (imcube, sigma, chanwidth=10):
    "produces a spectrum by convolving each slice of imcube with a gaussian of width sigma and returning the value of the central pixel for each slice, but cheats by just convolving the central pixel not the whole image"
    shape=imcube.shape
    bandwidth=shape[0]*chanwidth
    if len(shape)!=3: raise ValueError("imcube must be a cube")
    gauss_mask=garray(shape[1:],sigma)
    return [(i*chanwidth-bandwidth/2,(imcube[i,:,:]*gauss_mask).sum())
            for i in xrange(shape[0])]
    

def datafile(path):
    lines=[]
    f=open(path)
    for line in f: 
        if line.strip()[0] in '1234567890.-+': lines.append([numpy.float64(x) for x in line.split()])
    f.close()
    return numpy.array(lines)


def convolve (arr1, arr2):
    "convolves 2 arrays together with fft, arrays will be zero padded to equal size"
    if max(len(arr1.shape), len(arr2.shape)) > 2: raise ValueError("only dealing with 2d convolves here thankyou")
    s=(int(max(arr1.shape[0],arr2.shape[0])*1.5),int(max(arr1.shape[1],arr2.shape[1])*1.5))
    return irfft2(rfft2(arr1,s)*rfft2(arr2,s))
 
def single_line (path, sigma=27/(2*sqrt(2*log(2))), chanwidth=100, replot=False):
    im=pyfits.getdata(path)
    data=line(im,sigma,chanwidth)
    g('set style data histeps')
    if replot : g.replot(data)
    else : g.plot(data)
    

def linespectra (arr, freqs, sigma=4, channelWidth=20, kms=False, source_speed=0): #nb sigma is given in px (can be fractional)
    """arr should be an array of shape (x,:,>pix,>pix)
freqs an array or list of nums of length x"""
    shifts=[int(round((freqs[-1]-freqs[i])*299792458/(channelWidth*freqs[-1]))) for i in xrange(len(freqs))]
#    print shifts
    x=[[] for _ in xrange(arr.shape[0])]
    mid=arr.shape[2]/2.0-0.5

    gauss_mask=garray(arr.shape[-2:],sigma)
    s=[y*2 for y in gauss_mask.shape]
    ftg=rfft2(gauss_mask, s)

    for i in xrange(len(x)):
        for j in xrange(arr.shape[1]):
            convolved=irfft2(rfft2(arr[i,j,:,:],s)*ftg)
            x[i].append(convolved[s[0]/2,s[1]/2])

    padding=abs(max(shifts))
    padded=[0 for _ in xrange(arr.shape[1]+padding*2+2)]
    for i in xrange(len(x[0])):
        for j in xrange(len(x)):
            try:
                padded[i+shifts[j]+padding]+=x[j][i]
            except IndexError : 
                print j,i,len(x),len(x[j])
                None
    if kms: return [((i-150)*20/1000.0,x) for i,x in enumerate(padded)]
    else:   return [((i-150)*20,x)        for i,x in enumerate(padded)]

def main (tarobj):
    files=[tarobj.extractfile(x) for x in tarobj.getmembers() if x.name[-5:]=='.fits']
    files.sort(key=lambda x:x.name)
    print files[0].name
    tmp=pyfits.getdata(files[0])
    im=numpy.empty([len(files)]+list(tmp.shape), dtype=float)
    im[0,:,:,:]=tmp
    for i,f in enumerate(files[1:]):
        print f.name
        im[i+1,:,:,:]=pyfits.getdata(f)
        f.close()
    tarobj.close()

    freqs=[93.171616,93.171913,93.172050,93.173477,93.173775,93.173965,93.176254]
    x=numpy.zeros([7]+list(tmp.shape), dtype=float)
    x[0, :, :, :]=im[0 ,  :, :, :]
    x[1, :, :, :]=im[1:3, :, :, :].sum(0)
    x[2, :, :, :]=im[3:6, :, :, :].sum(0)
    x[3, :, :, :]=im[6:8, :, :, :].sum(0)
    x[4, :, :, :]=im[8:9, :, :, :].sum(0)
    x[5, :, :, :]=im[9:12,:, :, :].sum(0)
    x[6, :, :, :]=im[12:, :, :, :].sum(0)

    return linespectra(x,freqs, sigma=13/(2*sqrt(2*log(2))), kms=True)

def process_tar (tar):
    return [(x-15.75,y) for x,y in main(tar)]

def all_tars():
    tars=[]
    for  _,_,x in os.walk('./'):
        tars=sorted([tarfile.open(f, 'r:gz') for f in x if f[-7:]=='.tar.gz'], key=lambda x:x.name)
    lists={}
    for i,tar in enumerate(tars):
        print "****starting "+tar.name+"("+str(i+1)+" of "+str(len(tars))+")****"
        lists[tar.name]=(process_tar(tar))
    return lists

def diff (lis, data):
    tot=0
    for x1,y1 in lis:
        cont=True
        for x2,y2 in data:
            if cont and x2>x1: 
                tot+=(y1-y2)**2
                cont=False
    return numpy.sqrt(tot)

def get_data(infl):
    lines=[]
    for line in infl:
        l=line.split()
        lines.append((float(l[0]),float(l[1])))
    return lines

def dict2shiftedList (dic, shift, dat=False):
    ans=[]
    for key in dic:
        ans.append([dic[key],key])
    for i in xrange(len(ans)):
        ans[i][0]=[(x+shift,y) for x,y in ans[i][0]]
    if dat:
        ans.sort(key=lambda x:diff(x[0], dat))
    return ans

#old version w/o convolution
def oldlinespectra (arr, freqs, pix=1, channelWidth=20, kms=False, source_speed=0, sigma=0):
    """arr should be an array of shape (x,:,>pix,>pix)
freqs an array or list of nums of length x"""
    shifts=[int(round((freqs[-1]-freqs[i])*299792458/(channelWidth*freqs[-1]))) for i in xrange(len(freqs))]
    x=range(arr.shape[0])
    mid=arr.shape[2]/2.0-0.5
    if pix<1 or arr.shape[1]<1: raise ValueError
    elif pix%2==1 and arr.shape[2]%2==1:
        for i in xrange(len(x)):
            x[i]=[arr[i,j, int(mid+(pix-1)/2):int(mid+(pix-1)/2+1), int(mid-(pix-1)/2):int(mid+(pix-1)/2+1)].sum() for j in xrange(arr.shape[1])]
    elif pix%2==0 and arr.shape[1]%2==0:
        for i in xrange(len(x)):
            mid+=0.5
            x[i]=[arr[i,j, int(mid-pix/2):int(mid+pix/2), int(mid-pix/2):int(mid+pix/2)].sum() for j in xrange(arr.shape[1])]
    elif pix==1 and  arr.shape[1]%2==0:
        for i in xrange(len(x)):
            x[i]=[arr[i,j, int(mid-0.5):int(mid+1.5), int(mid-0.5):int(mid+1.5)].sum()/4 for j in xrange(arr.shape[1])]
    else: raise "not implemented"
    padding=abs(max(shifts))
    padded=[0 for _ in xrange(arr.shape[1]+padding*2+2)]
    for i in xrange(len(x[0])):
       for j in xrange(len(x)):
           try:
               padded[i+shifts[j]+padding]+=x[j][i]
           except IndexError : 
               print j,i,len(x),len(x[j])
               None
    if kms: return [((i-150)*20/1000.0,x) for i,x in enumerate(padded)]
    else:   return [((i-150)*20,x)        for i,x in enumerate(padded)]
