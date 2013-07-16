#!/usr/bin/python
from numpy import *
from scipy import interpolate
from pylab import plot, loglog, semilogx, semilogy, clf, imshow
import pylab as plt
from matplotlib.pyplot import gca


def datafile(path):
    lines=[]
    if isinstance(path, str): f=open(path)
    else: f=path
    for line in f: 
        l=line.strip()
        if len(l) > 0 and l[0] in '1234567890.-+': lines.append([float(x) for x in l.split()])
    f.close()
    return array(lines)

def _f(a,b):
    if   a[0]<b[0] : return -1
    elif a[0]>b[0] : return 1
    else           : return _f(a[1:], b[1:])

def xyplot(lis, **kwargs): 
    x,y=apply(zip, lis)
    plot(x,y, **kwargs)

def sortrows (arr, f):
    "sorts array arr by the ascending value of f applied to each row"
    x=list(arr)
    return array(sorted(x, key=f))

lookup={'d':3, 'D':3, 't':4, 'T':4, 'x':5, 'X':5}
for i in xrange(100):lookup[i]=i+6

def map2grid(dat, value, xpars,ypars,zpars, logScale=False, **kwargs):
    xmin,xmax,xstep=xpars
    ymin,ymax,ystep=ypars
    zmin,zmax,zstep=zpars
    one=ones((xstep,ystep,zstep))
    if 'X' in logScale.upper():
        X=(logspace(log10(xmin),log10(xmax),xstep)).reshape(xstep,1,1)*one
    else:
        X=(linspace(xmin,xmax,xstep)).reshape(xstep,1,1)*one

    if 'Y' in logScale.upper():
        Y=(logspace(log10(ymin),log10(ymax),ystep)).reshape(1,ystep,1)*one
    else:
        Y=(linspace(ymin,ymax,ystep)).reshape(1,ystep,1)*one

    if 'Z' in logScale.upper():
        Z=(logspace(log10(zmin),log10(zmax),zstep)).reshape(1,1,zstep)*one
    else:
        Z=(linspace(zmin,zmax,zstep)).reshape(1,1,zstep)*one

    points=dat[:,:3]
    vals=dat[:,lookup[value]]
    return interpolate.griddata(points,vals,(X,Y,Z))

def __main__():
    g1=g2=g3=g4=g5=0
    path=raw_input("enter path: ")

    try :dat=datafile (path)
    except IOError: 
        print "IOERROR"
        dat=datafile(raw_input("enter path: "))

    sdat=empty((dat.shape[0],dat.shape[1]-2))
    sdat[:,1:]=dat[:,3:]
    sdat[:,0]=sqrt((dat[:,:3]**2).sum(1))

    sdat=array(sorted(list(sdat), cmp=_f))

    clf()
    print 'grid points in rz'
#    plot(sqrt((dat[:,:2]**2).sum(1)), dat[:,2], ',')

    r,z = sqrt((dat[:,:2]**2).sum(1)), dat[:,2]
    H, xedges, yedges = histogram2d(r, z, bins=(20, 100))
    clf();plt.imshow(log10(H.T+0.1), vmin=0, interpolation='nearest',origin='image', cmap=plt.cm.gist_heat)
#    gca().set_xticks([x for x in gca().get_xticks() if x>=0 and x<len(xedges)])
#    gca().set_yticks([y for y in gca().get_yticks() if y>=0 and y<len(yedges)])
#    gca().set_xticklabels(['%.2f'%xedges[i]/2.5e13 for i in gca().get_xticks()])
#    gca().set_yticklabels(['%.2f'%yedges[i]/2.5e13 for i in gca().get_yticks()])
    plt.colorbar()
    plt.draw()

    raw_input("enter to continue")

    print 'grid points in log rz'
    H, xedges, yedges = histogram2d(log10(r), log10(abs(z)), bins=75)
    clf();plt.imshow(log10(H.T+0.1), vmin=0, interpolation='nearest',origin='image', cmap=plt.cm.gist_heat);plt.colorbar();
    gca().set_xticks([x for x in gca().get_xticks() if x>=0 and x<len(xedges)])
    gca().set_yticks([y for y in gca().get_yticks() if y>=0 and y<len(yedges)])
    gca().set_xticklabels(['%.2f'%(xedges[i]-log10(1.49598e11)) for i in gca().get_xticks()])
    gca().set_yticklabels(['%.2f'%(yedges[i]-log10(1.49598e11)) for i in gca().get_yticks()])
    plt.draw()

    raw_input("enter to continue")

    points=sorted([sqrt(x*x+y*y+z*z) for x,y,z in dat[:,0:3]])
    print 'cumulative distribution of points in r'
    clf()
    semilogx(points, xrange(len(points)), drawstyle='steps-mid')
    raw_input("enter to continue")

    bins=logspace(log10(min(points)),log10(max(points)), 51)
    histo=histogram(points,bins)
    print 'histogram (log spaced) grid points in r'
    clf();semilogx(histo[1][:-1],histo[0],drawstyle='steps-mid')
    raw_input("enter to continue")

    print "nh2 density"
    #h,m=dat[0,3],dat[0,5]
    #g.replot(Gnuplot.Data(sorted([(sqrt(x*x+y*y+z*z),nh2) for x,y,z,nh2        in dat[:,0:4] if nh2>0]), with_='dots'))
    clf();loglog (sdat[:,0],sdat[:,1],'.')
    gca().lines[-1].set_markeredgecolor(gca().lines[-1].get_markerfacecolor())
    raw_input("enter to continue")

    print "molecular abundance"
    loglog (sdat[:,0], sdat[:,3], '.')
    gca().lines[-1].set_markeredgecolor(gca().lines[-1].get_markerfacecolor())
    raw_input("enter to continue")

    print "molecular density"
    loglog (sdat[:,0], sdat[:,1]*sdat[:,3], '.')
    gca().lines[-1].set_markeredgecolor(gca().lines[-1].get_markerfacecolor())
    raw_input("enter to continue")

    try:
        statWeights=[int (x) for x in raw_input("input statistical weights of levels (or leave blank for simple rotor levels): ").split()]
        g1,g2,g3,g4,g5=statWeights
    except:
        g1,g2,g3,g4,g5=statWeights=1,3,5,7,9
    print statWeights

    try:
        print "population inversion points"
        inv10x,inv10y=apply(zip, [(sqrt(x*x+y*y),abs(z)) for x,y,z,nh2,t,nX,db,vx,vy,vz,n1,n2,n3 in dat[:,0:13] if not(1.01*n1/g1>=n2/g2)])
        print "{} inversions out of {} in 1->0".format(len(inv10x), dat.shape[0])
    except ValueError:
        print 'no inversions in 1-0'
        inv10x,inv10y=[],[]
    try:
        inv21x,inv21y=apply(zip, [(sqrt(x*x+y*y),abs(z)) for x,y,z,nh2,t,nX,db,vx,vy,vz,n1,n2,n3 in dat[:,0:13] if not(1.01*n2/g2>=n3/g3)])
        print "{} inversions out of {} in 2->1".format(len(inv21x), dat.shape[0])
    except ValueError:
        print 'no inversions in 2-1'
        inv21x,inv21y=[],[]
    try:
        inv32x,inv32y=apply(zip, [(sqrt(x*x+y*y),abs(z)) for x,y,z,nh2,t,nX,db,vx,vy,vz,n1,n2,n3,n4 in dat[:,0:14] if not(1.01*n3/g3>=n4/g4)])
        print "{} inversions out of {} in 3->2".format(len(inv32x), dat.shape[0])
    except ValueError:
        print 'no inversions in 3-2'
        inv32x,inv32y=[],[]

    if inv10x or inv21x or inv32x:
        clf()
        loglog(inv10x,inv10y,'x')
        loglog(inv21x,inv21y,'+')
        loglog(inv32x,inv32y,'*')
        raw_input("enter to continue")

        clf()
        plot(inv10x,inv10y,'x')
        plot(inv21x,inv21y,'+')
        plot(inv32x,inv32y,'*')
        raw_input("enter to continue")



    clf();plot([sqrt(x*x+y*y) for x,y in dat[:,:2]],[n1*g2/(n2*g1) for x,y,z,nh2,t,nX,db,vx,vy,vz,n1,n2,n3 in dat[:,0:13]],'.');semilogy()
    raw_input("enter to continue")

    print "checking for negative values"
    try: 
        negs=[(sqrt(x*x+y*y),abs(z)) for x,y,z,h,T,db,vx,vy,vz,m,n1,n2 in dat[:,0:8] if (h<=0 or T<=0 or m<=0 or n1<=0 or n2<=0)]
        if negs:
            clf()
            xyplot(negs, linestyle='.')
        print [(i, (h,T,m,n1,n2)) for (i,(x,y,z,h,T,db,vx,vy,vz,m,n1,n2)) in enumerate(dat[:,0:8]) if (h<=0 or T<=0 or m<=0 or n1<=0 or n2<=0)]
    except : print "ok"
    raw_input("enter to continue")

    print "kinetic and excitation temperature"
    kb=1.3806503e-23
    try : de=float(raw_input("enter n1-n2 energy gap (in Hz ev or J)"))
    except ValueError: de=1e11
    if de>1e6: de*=6.626068e-34
    elif de>1: de*=1.60217646e-19
    clf()
    semilogx(sqrt((dat[:,:3]**2).sum(1)),(log(dat[:,11]*g1/dat[:,10]/g2)*(-kb/de))**-1, 'r.')
    raw_input("enter to continue")
    semilogx(sqrt((dat[:,:3]**2).sum(1)), dat[:,4], 'g.')
    raw_input("enter to continue")

    print "level populations"
    clf()
    for i in xrange(min(5,dat.shape[1]-6)):
        loglog(sdat[:,0],sdat[:,8+i]/statWeights[i], '.')
        gca().lines[-1].set_markeredgecolor(gca().lines[-1].get_markerfacecolor())
#        loglog(sdat[:,0],sdat[:,4+i])
        raw_input("level {} (press enter to continue)".format(i))

    raw_input("enter to finish")
    return dat
