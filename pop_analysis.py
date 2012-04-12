#!/usr/bin/python
import Gnuplot, numpy
from sys import argv
g=Gnuplot.Gnuplot()

def datafile(path):
    lines=[]
    if isinstance(path, str): f=open(path)
    else: f=path
    for line in f: 
        if line.strip()[0] in '1234567890.-+': lines.append([float(x) for x in line.split()])
    f.close()
    return numpy.array(lines)




try: path=argv[1]
except ValueError: path=raw_input("enter path: ")

try :
    dat=datafile (path)
except IOError: 
    dat=datafile(raw_input("enter path: "))
g('reset;set style data dots')
g('unset logscale')
print 'grid points in rz'
g.plot([(numpy.sqrt(x*x+y*y),abs(z)) for x,y,z in dat[:,0:3]])
raw_input("enter to continue")

g('set logscale xy')
print 'grid points in log rz'
g.replot()
raw_input("enter to continue")

points=sorted([numpy.sqrt(x*x+y*y+z*z) for x,y,z in dat[:,0:3]])
min,max=numpy.log(points[0]),numpy.log(points[-1])
g('set xrange ['+str(numpy.exp(min))+':'+str(numpy.exp(max))+']')
g('unset logscale;set logscale x')
print 'cumulative distribution of points in r'
g.plot([(r,i) for i,r in enumerate(points)])
raw_input("enter to continue")

g('set style data steps')
lrange=max-min
bins=[numpy.exp(min+lrange*i/100.) for i in xrange(101)]
histo=list(numpy.histogram(points,bins)[0])+[0]
print 'histogram (log spaced) grid points in r'
g.plot([(bins[i],histo[i]) for i in xrange(len(bins))])
raw_input("enter to continue")

print "nh2 density"
g('set logscale xy')
h,m=dat[0,3],dat[0,5]
g.replot(Gnuplot.Data(sorted([(numpy.sqrt(x*x+y*y+z*z),nh2) for x,y,z,nh2        in dat[:,0:4] if nh2>0]), with_='dots'))
raw_input("enter to continue")

print "molecular density"
g.replot(Gnuplot.Data(sorted([(numpy.sqrt(x*x+y*y+z*z),nh2*nm) for x,y,z,nh2,_,nm in dat[:,0:6]if nh2>0 and nm>0]), with_='dots'))
raw_input("enter to continue")

g.plot(Gnuplot.Data(sorted([(numpy.sqrt(x*x+y*y+z*z),nh2*nm) for x,y,z,nh2,_,nm in dat[:,0:6]if nh2>0 and nm>0]), with_='dots'))
raw_input("enter to continue")
g.replot(Gnuplot.Data(sorted([(numpy.sqrt(x*x+y*y+z*z),nm) for x,y,z,nh2,_,nm in dat[:,0:6]if nh2>0 and nm>0]), with_='dots'))
raw_input("enter to continue")
g.replot(Gnuplot.Data(sorted([(numpy.sqrt(x*x+y*y+z*z),nh2) for x,y,z,nh2,_,nm in dat[:,0:6]if nh2>0 and nm>0]), with_='dots'))
raw_input("enter to continue")

print "population inversion points"
g('reset;set logscale xy; set style data points; set pointsize 0.5')
try: g.plot  ([(numpy.sqrt(x*x+y*y),abs(z)) for x,y,z,_,_,_,n1,n2,n3 in dat[:,0:9] if not(n1>=n2)])
except AssertionError : print "no masing in 1-0"
try: g.replot([(numpy.sqrt(x*x+y*y),abs(z)) for x,y,z,_,_,_,n1,n2,n3 in dat[:,0:9] if not(n2>=n3)])
except AssertionError : print "no masing in 2-1"
g('set pointsize 1.0')
raw_input("enter to continue")

print "checking for negative values"
try: 
    g.plot  ([(numpy.sqrt(x*x+y*y),abs(z)) for x,y,z,h,T,m,n1,n2 in dat[:,0:8] if (h<=0 or T<=0 or m<=0 or n1<=0 or n2<=0)])
    print [(i, (h,T,m,n1,n2)) for (i,(x,y,z,h,T,m,n1,n2)) in enumerate(dat[:,0:8]) if (h<=0 or T<=0 or m<=0 or n1<=0 or n2<=0)]
except AssertionError : print "ok"
raw_input("enter to continue")

print "kinetic and excitation temperature"
kb=1.3806503e-23
try : de=float(raw_input("enter n1-n2 energy gap (in Hz ev or J)"))
except ValueError: de=1e11
if de>1e6: de*=6.626068e-34
elif de>1: de*=1.60217646e-19
g('unset logscale; set logscale x; set style data dots')
g.plot([(numpy.sqrt(x*x+y*y+z*z),(numpy.log(n2/n1)*(kb/de)*-1)**-1) for x,y,z,_,_,_,n1,n2 in dat[:,0:8]])
g.replot([(numpy.sqrt(x*x+y*y+z*z),T) for x,y,z,_,T in dat[:,0:5]])
raw_input("enter to continue")

g('set logscale xy; set style data dots; set yrange [1e-5:1.1]')
g.plot([(numpy.sqrt(x*x+y*y+z*z),n1) for x,y,z,_,_,_,n1 in dat[:,0:7]])
for i in xrange(5):
    g.replot([(numpy.sqrt(dat[j,0]*dat[j,0]+dat[j,1]*dat[j,1]+dat[j,2]*dat[j,2]),dat[j,7+i])for j in xrange(len(dat[:,0]))])
#    g.replot([(numpy.sqrt(dat[j,0]*dat[j,0]+dat[j,1]*dat[j,1]+dat[j,2]*dat[j,2]),dat[j,6:].sum())for j in xrange(len(dat[:,0]))])

raw_input("enter to finish")
