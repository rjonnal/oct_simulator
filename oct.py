import numpy as np
from matplotlib import pyplot as plt
from oct_simulator import GaussianSpectrum,Sample
import sys,os
from scipy import interpolate

h = 6.626068e-34 # J*s, or  m^2 kg / s
c = 299792458.0 # m/s
e = 1.602176487e-19 #  C


class OCT:

    def __init__(self,source,sample,reference_fraction=0.9,ref_z_position=0.0,xmax=1.0e-3,ymax=1.0e-3,npx=512,dx=1e-5,dy=1e-5,dz=1e-6,dt=1e-5,pupil_diameter=6.75e-3):
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
        L1 = source.L[0]
        LN = source.L[-1]
        n = 1.38
        Sz = 1.0/((2*n)/L1 - (2*n)/LN)
        self.zmax = Sz*npx

        self.z_axis = np.arange(0.0,self.zmax,dz)
        self.z = np.zeros(self.z_axis.shape)
        

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

        in_z_axis = []
        in_z = []

        self.z = self.z*0.0
        
        for s in scatterers:
            if s.z<=self.zmax:
                in_z_axis.append(s.z)
                rad = np.sqrt((self.y-s.y)**2+(self.x-s.x)**2)
                blur = np.exp(-rad**2/(2*(self.psf_sigma**2)))
                coef = s.coef*blur*self.s_frac
                in_z.append(coef)

                left = np.where(self.z_axis<s.z)[0][-1]
                right = np.where(self.z_axis>=s.z)[0][0]

                tot = self.z_axis[right]-self.z_axis[left]
                left_weight = 1.0 - (np.abs(self.z_axis[left]-s.z))/tot
                right_weight = 1.0 - (np.abs(self.z_axis[right]-s.z))/tot
                self.z[left] = left_weight*coef
                self.z[right] = right_weight*coef


        k,Pk = self.source.get_k()
        kmat = k[:,np.newaxis]
        xnmat = self.z_axis[np.newaxis,:]
        xnk2 = 2*kmat*xnmat
        cosxnk2 = np.cos(xnk2)
        sqrtrncosxnk2 = np.sqrt(self.z)*cosxnk2
        sqrtrrsqrtrncosxnk2 = 2*np.sqrt(self.r_frac)*np.sum(sqrtrncosxnk2,axis=1)
        sumrn = np.sum(self.z)
        sumterm = self.r_frac + sumrn + sqrtrrsqrtrncosxnk2

        photons = Pk*sumterm*(2*np.pi/k)/h/c
        
        

        
        plt.cla()
        plt.plot(photons)
        plt.pause(.1)
        


            

    def update(self):
        self.sample.update()
        self.x = self.x + self.dx
        if self.x>self.xmax:
            self.x = 0.0
            self.y = self.y + self.dy
            if self.y>self.ymax:
                self.y = 0.0

        self.compute_aline()

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
        #oct.plot()
