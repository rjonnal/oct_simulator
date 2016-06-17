import numpy as np
from matplotlib import pyplot as plt
from movie import Movie

## Analysis of motion-based contrast methods

### Sources of motion

from oct_simulator import *

dt = 1e-6

source = GaussianSpectrum(1e-3,740e-9,940e-9,1024,80e-9)
sample = Sample(0.0,0.0,1.0e-4,dt=dt)
sample.add_scatterer(Scatterer(0.0,0.0,5.0e-5,0.0,0.0,30.0e-3,1.0e-3))
oct = OCT(source,sample)

adu = oct.compute_spectrum()

m = Movie(avifn='test.mp4')

for k in range(100):
    plt.subplot(1,2,1)
    sample.plot_sample('xz',False)
    plt.subplot(1,2,2)
    oct.plot_spectrum(False)
    plt.pause(.1)
    m.add(plt.gcf())
    oct.update()

m.make()
