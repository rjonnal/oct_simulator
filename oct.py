import numpy as np
from matplotlib import pyplot as plt
from oct_simulator import GaussianSpectrum,Sample

class OCT:

    def __init__(self,source,sample,reference_fraction=0.9,ref_z_position=0.0,xmax=1.0e-3,ymax=1.0e-3,dx=1e-5,dy=1e-5,dt=1e-5,pupil_diameter=6.75e-3):
        self.source = source
        self.sample = sample
        self.r_frac = reference_fraction
        self.s_frac = 1.0 - self.r_frac
        self.r_z_position = ref_z_position
        self.xmax = xmax
        self.ymax = ymax
        self.dx = dx
        self.dy = dy
        self.x = 0.0
        self.y = 0.0
        self.L = self.source.get_center_wavelength()
        self.Rx = 1.22*self.L*16.67e-3/pupil_diameter

    def update(self):
        print self.x,self.y
        self.sample.update()
        self.x = self.x + self.dx
        if self.x>self.xmax:
            self.x = 0.0
            self.y = self.y + self.dy
            if self.y>self.ymax:
                self.y = 0.0

    def plot(self):
        self.sample.plot('xy',do_pause=False)
        plt.plot(self.x,self.y,'ro')
        plt.pause(.0001)


    
        

if __name__=="__main__":

    dt = 1e-5
    xmax = 5e-4
    ymax = 1e-6
    zmax = 5e-4
    
    source = GaussianSpectrum(1e-3,740e-9,940e-9,1024,100e-9)
    sample = Sample(xmax,ymax,zmax,dt)
    for k in range(10):
        sample.add_random_scatterer()
    oct = OCT(source,sample,xmax=xmax,ymax=ymax,dt=dt)

    for k in range(10000):
        oct.update()
        oct.plot()
