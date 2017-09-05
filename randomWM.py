import numpy as np
from scipy.special import lambertw
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from astropy import units as u

# plotting
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def generate_disk(N,h=200,H=3000,maxdist=10**9):
    scale = 1
    
    # generate a simple dist

    # generate distances from the core, according to exponential with scale dist H
    # the distribution of r is according to x*exp(-x/h) (because 2D)
    if maxdist < 8300:
        cendist = coord.Galactocentric.galcen_distance.value*1000-maxdist
        minval = (H-np.exp(-cendist/H)*(H+cendist))/H
        r = H*(-1*lambertw(((np.random.rand(N)-1)*np.e**-1),-1).real-1)
        scale /= minval
        print 1./minval
    else:
        r = H*(-1*lambertw(((np.random.rand(N)-1)*np.e**-1),-1).real-1)

    # radial distribution
    if maxdist < 8300:
        theta = np.random.uniform(np.pi,2*np.pi,N)
    else:
        t = np.tan(maxdist/(coord.Galactocentric.galcen_distance.value*1000))
        theta = np.random.uniform(0+t,2*np.pi-t,N)
        print np.pi/t
        scale *= (np.pi/t)

    # make x,y
    x = r*np.sin(theta)
    y = r*np.cos(theta)

    # make random z distribution given scale height h
    z = h*np.log(np.random.rand(N))
    #z = h*2*np.arctan(np.tanh(0.5*np.random.rand(N)))/(0.5*np.pi) # if the distribution is sech(x/h) instead of exp(-x/h)
    z *= (-1)**np.random.randint(0,2,N) # above or below disk

    return np.c_[x,y,z],scale



def test_galaxy(N=10**7,blim=0,distlim=10**9):
    thin = generate_disk(int(0.98*N),h=200,H=3000)
    c_thin = SkyCoord(x=thin[:,0]*u.pc,y=thin[:,1]*u.pc,z=thin[:,2]*u.pc, frame='galactocentric')
    eq_thin = c_thin.transform_to(coord.ICRS)
    gc_thin = c_thin.transform_to(coord.Galactic)
    
    thick = generate_disk(int(0.02*N),h=1000,H=3000)
    c_thick = SkyCoord(x=thick[:,0]*u.pc,y=thick[:,1]*u.pc,z=thick[:,2]*u.pc, frame='galactocentric')
    eq_thick = c_thick.transform_to(coord.ICRS)
    gc_thick = c_thick.transform_to(coord.Galactic)
    
    # make masks
    m_thin = (abs(gc_thin.b.degree) > blim)&(abs(gc_thin.distance) < distlim*u.pc)
    m_thick = (abs(gc_thick.b.degree) > blim)&(abs(gc_thick.distance) < distlim*u.pc)
    
    # ratio of thin ot thick disk
    print("Total %d" %(np.sum(m_thin)+np.sum(m_thick)))
    print("ratio %g" %(np.sum(m_thin)*1./np.sum(m_thick)))
    
    #3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot(thin[m_thin,0],thin[m_thin,1],thin[m_thin,2],'k,')
    ax.plot(thick[m_thick,0],thick[m_thick,1],thick[m_thick,2],'k,')
     
    ax.set_xlim(-10**4,10**4)
    ax.set_ylim(-10**4,10**4)
    ax.set_zlim(-10**4,10**4)
    
    
    fig = plt.figure()
    plt.plot(eq_thin[m_thin].ra.deg,eq_thin[m_thin].dec.deg,'k,')
    plt.plot(eq_thick[m_thick].ra.deg,eq_thick[m_thick].dec.deg,'r,')
    
    
    plt.show()
