#!/usr/bin/python
from numpy import *
from pylab import plot, loglog, semilogx, clf

def datafile(path):
    lines=[]
    if isinstance(path, str): f=open(path)
    else: f=path
    for line in f: 
        if line.strip()[0] in '1234567890.-+': lines.append([float(x) for x in line.split()])
    f.close()
    return array(lines)

def _f(a,b):
    if   a[0]<b[0] : return -1
    elif a[0]>b[0] : return 1
    else           : return _f(a[1:], b[1:])

def xyplot(lis, **kwargs): 
    x,y=apply(zip, lis)
    plot(x,y, **kwargs)

def __main__():
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
    plot(sqrt((dat[:,:2]**2).sum(1)), abs(dat[:,2]), ',')

    raw_input("enter to continue")

    print 'grid points in log rz'
    clf();loglog(sqrt((dat[:,:2]**2).sum(1)), abs(dat[:,2]), ',')
    raw_input("enter to continue")

    points=sorted([sqrt(x*x+y*y+z*z) for x,y,z in dat[:,0:3]])
    print 'cumulative distribution of points in r'
    clf()
    semilogx(points, xrange(len(points)), drawstyle='steps-mid')
    raw_input("enter to continue")

    lrange=max(points)-min(points)
    bins=[min(points)+exp(i/50.)/exp(1)*lrange for i in xrange(51)]
    histo=list(histogram(points,bins)[0])+[0]
    print 'histogram (log spaced) grid points in r'
    clf();semilogx(bins,histo,drawstyle='steps-mid')
    raw_input("enter to continue")

    print "nh2 density"
    #h,m=dat[0,3],dat[0,5]
    #g.replot(Gnuplot.Data(sorted([(sqrt(x*x+y*y+z*z),nh2) for x,y,z,nh2        in dat[:,0:4] if nh2>0]), with_='dots'))
    clf();loglog (sdat[:,0],sdat[:,1],'.')
    raw_input("enter to continue")

    print "molecular abundance"
    loglog (sdat[:,0], sdat[:,3], '.')
    raw_input("enter to continue")

    print "molecular density"
    loglog (sdat[:,0], sdat[:,1]*sdat[:,3], '.')
    raw_input("enter to continue")

    try:
        statWeights=[int (x) for x in raw_input("input statistical weights of levels (or leave blank for simple rotor levels): ").split()]
        g1,g2,g3,g4,g5=statWeights
    except:
        g1,g2,g3,g4,g5=statWeights=1,3,5,7,9
    print statWeights

    try:
        print "population inversion points"
        inv10x,inv10y=apply(zip, [(sqrt(x*x+y*y),abs(z)) for x,y,z,_,_,_,n1,n2,n3 in dat[:,0:9] if not(n1/g1>=n2/g2)])
        print "{} inversions out of {} in 1->0".format(len(inv10x), dat.shape[0])
    except ValueError:
        print 'no inversions in 1-0'
        inv10x,inv10y=[],[]
    try:
        inv21x,inv21y=apply(zip, [(sqrt(x*x+y*y),abs(z)) for x,y,z,_,_,_,n1,n2,n3 in dat[:,0:9] if not(n2/g2>=n3/g3)])
        print "{} inversions out of {} in 1->0".format(len(inv21x), dat.shape[0])
    except ValueError:
        print 'no inversions in 2-1'
        inv21x,inv21y=[],[]
    try:
        inv32x,inv32y=apply(zip, [(sqrt(x*x+y*y),abs(z)) for x,y,z,_,_,_,n1,n2,n3,n4 in dat[:,0:10] if not(n3/g3>=n4/g4)])
        print "{} inversions out of {} in 1->0".format(len(inv32x), dat.shape[0])
    except ValueError:
        print 'no inversions in 3-2'
        inv32x,inv32y=[],[]

    clf()
    loglog(inv10x,inv10y,'x')
    loglog(inv21x,inv21y,'+')
    loglog(inv32x,inv32y,'*')

    raw_input("enter to continue")

    print "checking for negative values"
    try: 
        clf()
        xyplot([(sqrt(x*x+y*y),abs(z)) for x,y,z,h,T,m,n1,n2 in dat[:,0:8] if (h<=0 or T<=0 or m<=0 or n1<=0 or n2<=0)], linestyle='.')
        print [(i, (h,T,m,n1,n2)) for (i,(x,y,z,h,T,m,n1,n2)) in enumerate(dat[:,0:8]) if (h<=0 or T<=0 or m<=0 or n1<=0 or n2<=0)]
    except : print "ok"
    raw_input("enter to continue")

    print "kinetic and excitation temperature"
    kb=1.3806503e-23
    try : de=float(raw_input("enter n1-n2 energy gap (in Hz ev or J)"))
    except ValueError: de=1e11
    if de>1e6: de*=6.626068e-34
    elif de>1: de*=1.60217646e-19
    clf()
    semilogx(sqrt((dat[:,:3]**2).sum(1)),(log(dat[:,7]*g1/dat[:,6]/g2)*(kb/de)*-1)**-1, 'r.')
    semilogx(sqrt((dat[:,:3]**2).sum(1)), dat[:,4], 'g.')
    raw_input("enter to continue")

    print "level populations"
    clf()

    for i in xrange(min(5,dat.shape[1]-6)):
        loglog(sdat[:,0],sdat[:,4+i]/statWeights[i], '.')
        loglog(sdat[:,0],sdat[:,4+i])
        raw_input("level {} (press enter to continue)".format(i))

    raw_input("enter to finish")
    return dat
