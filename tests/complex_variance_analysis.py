import numpy as np
from matplotlib import pyplot as plt
from movie import Movie

## Analysis of motion-based contrast methods

### Sources of motion

from oct_simulator import *

dt = 1e-4

source = GaussianSpectrum(1e-3,740e-9,940e-9,1024,80e-9)
sample = Sample(xmax=0.0e-3,ymax=0.0e-3,zmax=1.0e-4,dt=dt)
s = Scatterer(x0=0.0e-3,y0=0.0e-3,z0=5.0e-5,vx=0.0,vy=0.0,vz=30.0e-3,coef=1.0e-3)
sample.add_scatterer(s)
oct = OCT(source,sample,xmax=0.0,ymax=0.0)

adu = oct.compute_spectrum()

#m = Movie(avifn='test.mp4')

for k in range(100):
#    print sample
    plt.subplot(1,2,1)
    sample.plot_sample('xz',False)
    plt.subplot(1,2,2)
    oct.plot_spectrum(False)
    plt.pause(.1)
#    m.add(plt.gcf())
    adu = oct.update()
    #plt.cla()
    #plt.plot(adu)
    #plt.pause(.1)
#m.make()
