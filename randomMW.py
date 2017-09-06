import numpy as np
from scipy.special import lambertw
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from astropy import units as u

# plotting
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class disk():
    def generate_disk(self):
        scale = 1
        GCdist = coord.Galactocentric.galcen_distance.value*1000.

        h = self.h
        H = self.H
        N = self.N
        maxdist = self.maxdist

        # generate a simple dist
        # generate distances from the core, according to exponential with scale dist H
        # the distribution of r is according to x*exp(-x/h) (because 2D)
        if self.maxdist < GCdist:
            g = lambda x: (H-np.exp(-x/H)*(H+x))/H
            minval = g(GCdist-self.maxdist)
            maxval = g(GCdist+self.maxdist)
            r = H*(-1*lambertw(((np.random.uniform(minval,maxval,N)-1)*np.e**-1),-1).real-1)
            print (maxval-minval)
            scale /= (maxval-minval)
        else:
            r = H*(-1*lambertw(((np.random.rand(N)-1)*np.e**-1),-1).real-1)

        # radial distribution
        if maxdist < GCdist:
            t = np.tan(maxdist/(GCdist))
            theta = np.random.uniform(1.5*np.pi-t,1.5*np.pi+t,N)
            scale *= (np.pi/t)
        else:
            theta = np.random.uniform(np.pi,2*np.pi,N)

        # make x,y
        x = r*np.sin(theta)
        y = r*np.cos(theta)

        # make random z distribution given scale height h
        z = h*np.log(np.random.rand(N))
        #z = h*2*np.arctan(np.tanh(0.5*np.random.rand(N)))/(0.5*np.pi) # if the distribution is sech(x/h) instead of exp(-x/h)
        z *= (-1)**np.random.randint(0,2,N) # above or below disk

        # density around the sun:
        self.c = SkyCoord(x=x*u.pc,y=y*u.pc,z=z*u.pc, frame='galactocentric')
        self.scale = scale
        self.eq = self.c.transform_to(coord.ICRS)
        self.gc = self.c.transform_to(coord.Galactic)
        self.density = (self.N*self.scale * np.exp(-8300./self.H)/self.H**2 * 
                            np.exp(-27./self.h)/(2*self.h*2*np.pi))
        self.dmask = self.eq.distance.pc < self.maxdist

    def calc_scale(self):
            scale = 1
            GCdist = coord.Galactocentric.galcen_distance.value*1000    
            if self.maxdist < GCdist:
                g = lambda x: (self.H-np.exp(-x/self.H)*(self.H+x))/self.H
                minval = g(GCdist-self.maxdist)
                maxval = g(GCdist+self.maxdist)
                scale /= (maxval-minval)
            else:
                pass

            # radial distribution
            if self.maxdist < GCdist:
                t = np.tan(self.maxdist/(GCdist))
                scale *= (np.pi/t)
            else:
                theta = np.random.uniform(np.pi,2*np.pi,N)
                
            self.scale = scale

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

    def plot_overview(self,*args):
        plt.plot(self.c.x[self.dmask],self.c.y[self.dmask],'r.',*args)
        plt.plot(self.c.x[~self.dmask],self.c.y[~self.dmask],'b.',*args)
        plt.axes().set_aspect('equal', 'datalim')

    def plot_view(self,*args):
        plt.plot(self.eq.ra[self.dmask],self.eq.dec[self.dmask],
            marker='.',*args,lw=0)

    def Nfromdensity(self,density):
        self.N = int(density/self.scale / (np.exp(-8300./self.H)/self.H**2 * 
                            np.exp(-27./self.h)/(2*self.h*2*np.pi) ))
