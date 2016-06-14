import numpy as np
from matplotlib import pyplot as plt
from oct_simulator import GaussianSpectrum,Sample
import sys,os
from scipy import interpolate
import xml.etree.cElementTree as ET
from xml.etree import ElementTree
from xml.dom import minidom
from octopod import *

h = 6.626068e-34 # J*s, or  m^2 kg / s
c = 299792458.0 # m/s
e = 1.602176487e-19 #  C


class OCT:

    def __init__(self,source,sample,reference_fraction=0.9,ref_z_position=0.0,xmax=1.0e-3,ymax=1.0e-3,npx=512,dx=1e-5,dy=1e-5,dz=1e-6,nvol=1,nbm=1,dt=1e-5,pupil_diameter=6.75e-3,photons_per_adu=5e9):
        self.source = source
        self.sample = sample
        self.r_frac = reference_fraction
        self.s_frac = 1.0 - self.r_frac
        self.r_z_position = ref_z_position
        self.xmax = xmax
        self.ymax = ymax
        self.npx = int(npx)
        self.dx = dx
        self.nx = int(self.xmax/self.dx)
        self.dy = dy
        self.ny = int(self.ymax/self.dy)
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
        self.nvol = int(nvol)
        self.nbm = 1
        self.photons_per_adu = photons_per_adu

    def write_xml(self,fn):
        monsterlist = ET.Element("MonsterList")
        monsterlist.set("version","1.o")
        monster = ET.SubElement(monsterlist, "Monster")
        name = ET.SubElement(monster,"Name")
        name.text = "Goblin"
        
        attrib = {"Data_Acquired_at":"01/01/2000 0:00:00 PM"}
        time = ET.SubElement(monster,"Time",attrib)

        attrib = {"Width":"%d"%(self.npx),
                  "Height":"%d"%(self.nx),
                  "Number_of_Frames":"%d"%(self.ny),
                  "Number_of_Volumes":"%d"%(self.nvol)}
        volumesize = ET.SubElement(monster,"Volume_Size",attrib)

        xsr = 3168*self.xmax/300e-6
        ysr = 2112*self.ymax/300e-6
        xso = 0
        yso = 0
        attrib = {"X_Scan_Range":"%d"%(xsr),"Y_Scan_Range":"%d"%(ysr),
                  "X_Scan_Offset":"%d"%xso,"Y_Scan_Offset":"%d"%yso,
                  "Number_of_BM_scans":"%d"%(self.nbm)}
        scanningparameters = ET.SubElement(monster,"Scanning_Parameters",attrib)

        c2 = 0.0
        c3 = 0.0
        attrib = {"C2":"%0.2e"%c2,"C3":"%0.2e"%c3}
        dispersionparameters = ET.SubElement(monster,"Dispersion_Parameters",attrib)
        
        
        #tree = ET.ElementTree(monsterlist)
        #tree.write("filename.xml")

        def prettify(elem):
            """Return a pretty-printed XML string for the Element.
            """
            rough_string = ElementTree.tostring(elem, 'utf-8')
            reparsed = minidom.parseString(rough_string)
            return reparsed.toprettyxml(indent="  ")
    
        fid = open(fn,'w')
        fid.write(prettify(monster))
        fid.close()

        
    def psf_radial_intensity(self,q):
        return np.exp(-q**2/(2*(self.psf_sigma**2)))


    def find_scatterers(self):
        effective_scatterers = []
        for s in self.sample.scatterers:
            rad = np.sqrt((self.y-s.y)**2+(self.x-s.x)**2)
            if rad<=self.beam_radius:
                effective_scatterers.append(s)
        return effective_scatterers


    def compute_spectrum(self):
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
        self.adu = photons/self.photons_per_adu
        return self.adu

    def update(self):
        print self.x,self.y
        adu = self.compute_spectrum()
        self.sample.update()
        self.x = self.x + self.dx
        if self.x>=self.xmax:
            self.x = 0.0
            self.y = self.y + self.dy
            if self.y>=self.ymax:
                self.y = 0.0
        return adu

    def acquire(self,fn='test',do_plot=False):
        out = np.zeros((self.nvol,self.ny,self.nx,self.npx),dtype=np.uint16)
        for ivol in range(self.nvol):
            for iy in range(self.ny):
                for ix in range(self.nx):
                    adu = self.update()
                    if do_plot:
                        self.plot()
                    out[ivol,iy,ix,:] = adu

        fnroot = os.path.splitext(fn)[0]
        dfn = fnroot+'.unp'
        xfn = fnroot+'.xml'
        out.tofile(dfn)
        self.write_xml(xfn)
                    
    def plot(self):
        plt.subplot(1,2,1)
        self.sample.plot('xy',do_pause=False)
        plt.plot(self.x,self.y,'ro')
        for s in self.find_scatterers():
            plt.plot(s.x,s.y,'bs')
        plt.subplot(1,2,2)
        plt.cla()
        plt.plot(self.source.k,self.adu)
        plt.pause(.1)


    
        

if __name__=="__main__":

    dt = 1e-5
    xmax = 5e-5
    ymax = 5e-5
    zmax = 5e-4
    npx = 1024
    
    source = GaussianSpectrum(1e-3,740e-9,940e-9,npx,100e-9)
    sample = Sample(xmax,ymax,zmax,dt)
    for k in range(1000):
        sample.add_random_scatterer()
    oct = OCT(source,sample,xmax=xmax,ymax=ymax,dt=dt,npx=npx)
    oct.acquire('test')
    
