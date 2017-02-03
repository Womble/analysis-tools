import vtk
from vtk.util import numpy_support as VN
import numpy as np
from math import *
from matplotlib.mlab import griddata

lAnd=np.logical_and
lOr =np.logical_or
griddistancescale=1e10
msol_g=1.988435e33
year_s=365.25*24*3600
k_cgs=1.38065e-16
MMM_g=2.12597e-24
G_cgs=6.674e-8
au_cm=1.496e13
Lsol_cgs=3.848e33

def parseDump (filename):
    reader=vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    output=reader.GetOutput()
    celldata=output.GetCellData()
    cells=output.GetCells()
    points=output.GetPoints()
    return cells,celldata,points

def parseSource(filename):
    X=np.fromfile(filename).astype(float32)
    radius,lum,teff,mass,mdot=X[14],X[16]/Lsol_cgs,X[18],X[22]/msol_g,X[26]/msol_g*year_s
    return radius,lum,teff,mass,mdot# in R_sol,L_sol,K,msol,msol/s
    
def randomUnitVector2():
    flag=True

    while flag:
        vec=np.random.rand(2)
        if sqrt((vec**2).sum())<=1: flag=False

    return vec/sqrt((vec**2).sum())

def randomUnitVector2Cylindrical():
    flag=True

    while flag:
        vec=np.random.rand(3)
        if sqrt((vec**2).sum())<=1: flag=False

    vec[0]=sqrt(vec[0]**2+vec[2]**2)
    vec[2]=0
    if np.random.randint(2)==0: vec[1]*=-1
    return vec/sqrt((vec**2).sum())

def randomUnitVector3():
    flag=True

    while flag:
        vec=np.random.rand(3)
        if sqrt((vec**2).sum())<=1: flag=False

    return vec/sqrt((vec**2).sum())

class Grid():
    def __init__ (self, cells, celldata, points, sourcefile=None, starMass=1):
        cellpoints=np.delete(VN.vtk_to_numpy(cells.GetData()), slice(None,None,5)) # remove the first element and evey fith following that, cells is the the form [nPoints_1, point1_1, point2_1, ..., points_1, nPoints_2, ...] as all our cells have 4 points we can remove the 4 which occurs every 5 numbers
        loc=VN.vtk_to_numpy(points.GetData())
        pos=[]
        sizes=[]
        for i in xrange(cellpoints.size/4):
            j=i*4
            size=np.inf
            pos.append((loc[cellpoints[j+0],:]+ \
                             loc[cellpoints[j+1],:]+ \
                             loc[cellpoints[j+2],:]+ \
                             loc[cellpoints[j+3],:])/4.0)
            for k in xrange(3):
                for l in xrange(k,4):
                    d=loc[cellpoints[j+k]]-loc[cellpoints[j+l]]
                    s=sqrt((d**2).sum())
                    if s<size and s>0 : size=s
            sizes.append(size)
        self.cellsize=np.array(sizes)
        self.points=loc
        self.pos=np.array(pos)
        self.cellvol=self.cellsize**2*2*pi*self.pos[:,0]
        self.rho=VN.vtk_to_numpy(celldata.GetArray('rho'))
        self.pressure=VN.vtk_to_numpy(celldata.GetArray('pressure'))
        self.temp=VN.vtk_to_numpy(celldata.GetArray('temperature'))
        self.dustFrac=VN.vtk_to_numpy(celldata.GetArray('dust1'))
        self.rhou=VN.vtk_to_numpy(celldata.GetArray('rhou'))
        self.rhov=VN.vtk_to_numpy(celldata.GetArray('rhov'))
        self.rhow=VN.vtk_to_numpy(celldata.GetArray('rhow'))
        self.jnu=VN.vtk_to_numpy(celldata.GetArray('jnu'))
        self.HI=VN.vtk_to_numpy(celldata.GetArray('HI'))
        self.phi=VN.vtk_to_numpy(celldata.GetArray('phi'))
        self.KE_lin=0.5*(self.rhou**2+self.rhow**2)/self.rho**2
        self.KE_rot=0.5*(self.rhov/(griddistancescale*self.pos[:,0]))**2/self.rho**2
        self.KE=self.KE_rot+self.KE_lin
        self.boundness=abs(self.phi)/(self.pressure+self.KE)
        if sourcefile:
            self.starRadius,self.starLum,self.starTeff,self.starMass,self.starMdot=parseSource(sourcefile)
        else:
            self.mass=starMass*msol_g
        
        
    def massWithinRadius(self, radius):
        if radius>1e11 : radius/=1e10
        cellradius=np.sqrt((self.pos**2).sum(1))
        mass=0
        for i in xrange(self.rho.size):
            if cellradius[i]<(radius-sqrt(3/4)*self.cellsize[i]):
                mass+=self.rho[i]* (self.pos[i,0]*2*pi) *self.cellsize[i]**2 * griddistancescale**3
            elif cellradius[i]<(radius+sqrt(3/4)*self.cellsize[i]):
                frac=(radius-cellradius[i]+(cellsize[i]*sqrt(3/4))) / sqrt(3)*self.cellsize[i]
                frac=min(1,max(0,frac))
                mass+=self.rho[i]* (self.pos[i,0]*2*pi) *self.cellsize[i]**2 * griddistancescale**3\
                      *frac
            else:
                None
        mass/=msol_g
        return mass

    def returnNearestCellRef(self, pos):
        dmin=np.inf
        ref=np.nan
        for i in xrange(self.pos[:,0].size):
            dist=sqrt(((pos-self.pos[i,:])**2).sum())
            if dist<dmin:
                ref=i
                dmin=dist
        return ref

    def returnNearestCellRefs(self, positions):
        dmin=np.ones_like(positions[:,0])*np.inf
        ref=np.zeros(positions.shape[0], dtype=np.int)*np.nan
        distances=np.zeros_like(positions)
        for i in xrange(self.pos[:,0].size):
            distances*=0
            distances+=positions
            distances-=self.pos[i,:]
            distances*=distances
            dist=np.sqrt(distances.sum(1))
            mask=dist<dmin
            dmin[mask]=dist[mask]
            ref[mask]=i
        return ref

    def fluxThroughRadius(self, radius, n=1e3, inwardsOnly=False, outwardsOnly=False):
        n=int(n)
        da=4.0*pi*(radius*griddistancescale)**2 / n #cm^2
        totflux=0.0
        vecs=np.array([randomUnitVector2Cylindrical() for i in xrange(n)])
        refs=self.returnNearestCellRefs(vecs*radius)
        for i in xrange(n):
            vec=vecs[i]
            ref=refs[i]
            flux= -1*(self.rhou[ref]*vec[0]+\
                      self.rhow[ref]*vec[1])*da #g.cm^-3
            if inwardsOnly:  flux=max(flux,0)
            if outwardsOnly: flux=min(flux,0)
            totflux+=flux

        return totflux/msol_g*year_s

    def griddata(self, data, Rmax=7e7, Zmin=0, Zmax=3.5e7, n=200):
        R,z=np.mgrid[0:Rmax:n*1j,Zmin:Zmax:n*1j]
        return griddata(self.pos[:,0], self.pos[:,1],data, R,z, interp='linear')

    def discMask(self):
        return lAnd(self.KE_rot>self.KE_lin, self.KE_rot>abs(self.phi))
    
    def ToomreUnstable(self):
        discMask=self.discMask()
        Rmax=self.pos[discMask][:,0].max()
        gdRho=self.griddata(self.rho, Rmax=Rmax, Zmin=-Rmax/2, Zmax=Rmax/2)
        gdTemp=self.griddata(self.temp, Rmax=Rmax, Zmin=-Rmax/2, Zmax=Rmax/2)
        gdBoundness=self.griddata(np.maximum(self.boundness,1),\
                                  Rmax=Rmax, Zmin=-Rmax/2, Zmax=Rmax/2)
        gdBoundness[gdBoundness<0.5]=0.0 #material with ke+pressure >|2*phi| is excluded from the disc

        Rmax*=griddistancescale #Rmax now in cm rather than code units
        boundMassSigma=[]
        boundMassWeightedTemp=[]
        for i in xrange(200):
            boundMassSigma.append((gdRho[i,:]*gdBoundness[i,:]).sum()*Rmax/200)
            boundMassWeightedTemp.append((gdTemp[i,:]*gdRho[i,:]*gdBoundness[i,:]).sum()\
                                         *Rmax/200/boundMassSigma[-1])
        boundMassSigma=np.array(boundMassSigma)
        boundMassWeightedTemp=np.array(boundMassWeightedTemp)

        Rs=np.linspace(0,Rmax,200)
        Omega=np.sqrt(G_cgs*self.mass/Rs**3)
        cs=np.sqrt(5.0/3*k_cgs*boundMassWeightedTemp/MMM_g)
        Q=Omega*cs/(boundMassSigma*G_cgs*pi)
        return Rs/au_cm,Q

    def outflowMomentum(self, rmin=0, boundnessLim=1):
        totflux=0.0
        for i in xrange(self.rho.size):
            if self.boundness[i]<boundnessLim:
                vec=self.pos[i,:]
                modvec=np.sqrt((vec**2).sum())

                if modvec>=rmin:
                    vec/=modvec
                    flux= (self.rhou[i]*vec[0]+\
                           self.rhow[i]*vec[1])*self.cellvol[i]*griddistancescale**3 #g.cm^-3
                    flux=max(flux,0)
                    totflux+=flux

        return totflux/msol_g/1e5 #momentum in msol.km/s
    
    def discMass(self, boundnessLim=1):
        mask=self.discMask()
        mask=lAnd(mask,self.boundness>boundnessLim)
        mass=(self.rho[mask]*self.cellvol[mask]).sum()*griddistancescale**3
        return mass/msol_g
    
    def outflowMass(self, boundnessLim=1):
        totmass=0.0
        for i in xrange(self.rho.size):
            if self.boundness[i]<boundnessLim:
                vec=self.pos[i,:]
                vec/=np.sqrt((vec**2).sum())

                if (self.rhou[i]*vec[0]+self.rhow[i]*vec[1])>0:
                    totmass+=self.rho[i]*self.cellvol[i]*griddistancescale**3
            return totmass/msol_g

    def virial(self, rmin=0, rmax=np.inf):
        totmass=0
        totvir=0
        for i in xrange(self.rho.size):
            r=sqrt((self.pos[i,:]**2).sum())
            if r>rmin and r<rmax:
                mass=self.rho[i]*self.cellvol[i]
                totvir+=self.boundness[i]*mass
                totmass+=mass
        return totvir/totmass

    
