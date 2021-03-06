import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate

class Spectrum:

    def __init__(self,N=1024):
        self.N = N
        self.k_computed = False

    def set_L(self,L0=790e-9,dL=.1e-9):
        self.L = np.arange(self.N)*dL+L0

    def set_P(self,P):
        self.P = P

    def get_k(self):
        if not self.k_computed:
            self.k_interpolate()
        return self.k,self.Pk

    def k_interpolate(self):
        self.k = np.linspace(2*np.pi/self.L[0],2*np.pi/self.L[-1],self.N)
        k_rev = self.k[::-1]
        P_rev = self.P[::-1]

        k_in = (2*np.pi)/self.L
        k_in_rev = k_in[::-1]
        
        interpolator = interpolate.interp1d(k_in_rev,P_rev)
        Pk_rev = interpolator(k_rev)
        self.Pk = Pk = Pk_rev[::-1]
        self.k_computed = True
        

    def set_P_gaussian(self,total_power=1e-3,L_sigma=20e-9):
        self.L_sigma = L_sigma
        self.fwhm = self.sigma_to_fwhm(L_sigma)
        g = np.exp(-(self.L-np.mean(self.L))**2/(2*L_sigma**2))
        self.P = g/np.sum(g)*total_power

    def sigma_to_fwhm(self,sigma):
        return 2*np.sqrt(2*np.log(2))*sigma

    def fwhm_to_sigma(self,fwhm):
        return fwhm/2/np.sqrt(2*np.log(2))


    def get_fwhm(self):
        valid = np.where(self.P>np.max(self.P)/2.0)[0]
        return self.L[valid[-1]] - self.L[valid[0]]

    def plot(self,mode='L'):
        if mode=='L':
            plt.plot(self.L,self.P)
            plt.show()

    def get_center_wavelength(self):
        return np.sum(self.L*self.P)/np.sum(self.P)

class GaussianSpectrum(Spectrum):

    def __init__(self,total_power,L_start,L_end,N,L_fwhm):
        Spectrum.__init__(self,N)
        dL = (L_end-L_start)/N
        L_sigma = Spectrum.fwhm_to_sigma(self,L_fwhm)
        Spectrum.set_L(self,L0=L_start,dL=dL)
        Spectrum.set_P_gaussian(self,total_power,L_sigma)
        
    
if __name__=="__main__":

    g = GaussianSpectrum(1.0e-3,740e-9,940e-9,1024,100e-9)
    g.plot()
