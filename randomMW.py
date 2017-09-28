import numpy as np
from scipy.special import lambertw
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from astropy import units as u
import scipy.integrate as integrate

# plotting
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class expdisk():
    def generate_disk(self):
        GCdist = coord.Galactocentric.galcen_distance.value*1000.

        h = self.h
        H = self.H
        rho0 = self.rho0

        theta = np.random.uniform(np.pi,2*np.pi,N)

        # the r component
        rho_r = lambda r : np.exp(-r/H)/np.exp(-8300./H) # normalised density
        pdf_r = lambda r : 2*np.pi * r * rho_r(r) # cylindrical coordiantes
        cdf_r = lambda l : integrate.quad(pdf_r,0.,l)
        rgrid = np.hstack(([0],np.logspace(0,7,10000)))
        totr = cdf_r(np.inf)[0]
        icdf_r = lambda n: np.interp(n,[cdf_r(_r)[0]/totr for _r in rgrid],rgrid)
        
        # the z compoment
        rho_z = lambda z : np.exp(-z/h)/np.exp(-27./h) # normalised density
        cdf_z = lambda l : integrate.quad(rho_z,0,l) # cdf
        rgrid = np.hstack(([0],np.logspace(-2,5,10000)))
        totz = cdf_z(np.inf)[0]
        icdf_z = lambda n: np.interp(n,[cdf_z(_r)[0]/totz for _r in rgrid],rgrid)
    
        # total number of points in disk
        N = rho0*cdf_r(np.inf)*2*cdf_z(np.inf) # factor 2 is from z and -z

        r = icdf_r[np.uniform(N)] # make r
        z = icdf_z[np.uniform(N)] # make z
        z *= (-1)**np.random.randint(0,2,N) # give z positive and negative values

        # make x,y
        x = r*np.sin(theta)
        y = r*np.cos(theta)


        # density around the sun:
        self.c = SkyCoord(x=x*u.pc,y=y*u.pc,z=z*u.pc, frame='galactocentric')
        self.scale = scale
        self.eq = self.c.transform_to(coord.ICRS)
        self.gc = self.c.transform_to(coord.Galactic)



    def __init__(self,rho,h,H,maxdist):
        GCdist = coord.Galactocentric.galcen_distance.value*1000.
        self.rho = rho
        self.h = h
        self.H = H
        self.maxdist = maxdist
        self.C = (np.exp(-GCdist/self.H)/self.H**2 * 
                    np.exp(-27./self.h)/(2*self.h*2*np.pi))
        self.calc_scale()
        self.Nfromdensity(rho)



class halo():
    def generate_halo(self):
        p = self.p
        k = self.k
        rc = self.rc # pc
        maxdist = self.maxdist
        rs2 = 8300**2+(27./k)**2 # distance from galactic center
        rho0 = self.rho0
        rho = lambda r: rho0*((rc**2+r**2)/(rc**2+rs2))**(-p/2)

        
        # calculate the integrand
        pdf = lambda r: 4*np.pi*rho(r)*r**2
        cdf = lambda l : integrate.quad(pdf,0,l)
        N = int(cdf(maxdist)[0])
        self.N = N
        r = np.linspace(0,maxdist,1000) # needed for interpolation
        icdf = lambda n: np.interp(n,[cdf(_r)[0]/N for _r in r],r)

        # radial distribution
        theta = np.random.uniform(0,2*np.pi,N)
        phi = np.random.uniform(0,np.pi,N)
        r = icdf(np.random.rand(N))

        # make x,y
        x = r*np.sin(theta)*np.sin(phi)
        y = r*np.cos(theta)*np.sin(phi)
        z = k*r*np.cos(phi) # is multiplied by c for oblateness

        # density around the sun:
        self.c = SkyCoord(x=x*u.pc,y=y*u.pc,z=z*u.pc, frame='galactocentric')
        self.eq = self.c.transform_to(coord.ICRS)
        self.gc = self.c.transform_to(coord.Galactic)

    def __init__(self,rho0,p,k,rc,maxdist):
        self.GCdist = coord.Galactocentric.galcen_distance.value*1000.
        self.rho0 = rho0 # halo stellar density around sun
        self.p = p #power
        self.k = k #oblateness <1
        self.rc = rc # core radius
        self.maxdist = maxdist # maxdist of halo

