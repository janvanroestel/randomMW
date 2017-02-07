import numpy as np

def generate_disk(N,h=200,H=3000):
    # generate a simple disk

    # generate distances from the core
    _r_min = 1-np.exp(r_min/-H)
    scale /= (1-_r_min)
    r = -H*np.log(1-np.random.uniform(_r_min,1,N))

    # radial distribution
    theta = np.random.uniform(0,2*np.pi,N)

    # make x,y
    x = r*np.cos(theta)
    y = r*np.sin(theta)

    # make random z distribution given scale height h
    z = -h*np.log(1-np.random.rand(N))


    return x,y,z
