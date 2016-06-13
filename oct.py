import numpy as np
from matplotlib import pyplot as plt
from oct_simulator import GaussianSpectrum,Sample

class OCT:

    def __init__(self,source,sample,reference_fraction=0.9,ref_z_position=0.0,xmax=1.0e-3,ymax=1.0e-3,npx=512,dx=1e-5,dy=1e-5,dz=1e-8,dt=1e-5,pupil_diameter=6.75e-3):
        self.source = source
        self.sample = sample
        self.r_frac = reference_fraction
        self.s_frac = 1.0 - self.r_frac
        self.r_z_position = ref_z_position
        self.xmax = xmax
        self.ymax = ymax
        self.npx = npx
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.x = 0.0
        self.y = 0.0
        self.L = self.source.get_center_wavelength()
        self.Rx = 1.22*self.L*16.67e-3/pupil_diameter
        self.psf_sigma = 0.42*self.L*16.67e-3/pupil_diameter
        max_rad = np.sqrt(xmax**2+ymax**2)
        q = np.linspace(0,max_rad,2**16)
        self.radial_psf = np.exp(-q**2/(2*(self.psf_sigma**2)))
        self.beam_radius = q[np.where(self.radial_psf<.05)[0][0]]
        # make an empty depth profile and axis vector; these are
        # filled with scattering coefficients at imaging time
        # first compute maximum depth:
        

    def psf_radial_intensity(self,q):
        return np.exp(-q**2/(2*(self.psf_sigma**2)))


    def find_scatterers(self):
        effective_scatterers = []
        for s in self.sample.scatterers:
            rad = np.sqrt((self.y-s.y)**2+(self.x-s.x)**2)
            if rad<=self.beam_radius:
                effective_scatterers.append(s)
        return effective_scatterers


    def compute_aline(self):
        scatterers = self.find_scatterers()
        

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
        for s in self.find_scatterers():
            plt.plot(s.x,s.y,'bs')
        plt.pause(.1)


    
        

if __name__=="__main__":

    dt = 1e-5
    xmax = 5e-5
    ymax = 5e-5
    zmax = 5e-4
    
    source = GaussianSpectrum(1e-3,740e-9,940e-9,1024,100e-9)
    sample = Sample(xmax,ymax,zmax,dt)
    for k in range(1000):
        sample.add_random_scatterer()
    oct = OCT(source,sample,xmax=xmax,ymax=ymax,dt=dt)

    for k in range(10000):
        oct.update()
        oct.plot()
