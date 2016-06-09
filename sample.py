import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Scatterer:

    def __init__(self,x0,y0,z0,vx=0.0,vy=0.0,vz=0.0,coef=1.0e-6):
        # coef is the fraction of photons returned to the pupil
        self.x = x0
        self.y = y0
        self.z = z0
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.coef = coef

    def update(self,dt,xmax,ymax,zmax):
        self.x = (self.x + self.vx*dt)%xmax
        self.y = (self.y + self.vy*dt)%ymax
        self.z = (self.z + self.vz*dt)%zmax

    def reset(self):
        self.x = self.x0
        self.y = self.y0
        self.z = self.z0


class Sample:

    def __init__(self,xmax,ymax,zmax,t0=0.0,dt=1e-3):
        self.xmin = 0.0
        self.xmax = xmax
        self.ymin = 0.0
        self.ymax = ymax
        self.zmin = 0.0
        self.zmax = zmax
        self.scatterers = []
        self.t = t0
        self.coef_max = -np.inf
        self.coef_min = np.inf
        self.dt = dt

    def add_scatterer(self,scatterer):
        if scatterer.coef<self.coef_min:
            self.coef_min = scatterer.coef
        if scatterer.coef>self.coef_max:
            self.coef_max = scatterer.coef
        self.scatterers.append(scatterer)

    def out_of_bounds(self,s):
        return not (s.x>self.xmin and s.x<self.xmax and s.y>self.ymin and
                s.x<self.ymax and s.z>self.zmin and s.z<self.zmax)

    def get_vec(self,dim):
        out = []
        for s in self.scatterers:
            out.append(eval('s.%s'%dim))
        return np.array(out)
    
    def plot3d(self,nsteps=1000):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlim((self.xmin,self.xmax))
        ax.set_ylim((self.ymin,self.ymax))
        ax.set_zlim((self.zmin,self.zmax))
        for k in range(nsteps):
            ax.cla()
            x = self.get_vec('x')
            y = self.get_vec('y')
            z = self.get_vec('z')
            self.update()
            ax.scatter(x, y, z)
            plt.pause(.1)
            
        plt.close()
            
    def plot0(self,proj='xy'):
        d1 = proj[0]
        d2 = proj[1]
        plt.cla()
        for s in self.scatterers:
            x = eval('s.%s'%d1)
            y = eval('s.%s'%d2)
            if self.coef_max==self.coef_min:
                alpha = 1.0
            else:
                alpha = (s.coef-self.coef_min)/(self.coef_max-self.coef_min)
            plt.plot(x,y,'go',alpha=alpha)
        plt.xlim((eval('self.%smin'%d1),eval('self.%smax'%d1)))
        plt.ylim((eval('self.%smin'%d2),eval('self.%smax'%d2)))
        plt.title('t=%0.3f'%self.t)
        plt.pause(.0001)

    def plot(self,proj='xy',do_pause=True):
        d1 = proj[0]
        d2 = proj[1]
        plt.cla()
        x = self.get_vec(d1)
        y = self.get_vec(d2)
        #alpha = self.get_vec('coef')
        #alpha = (alpha - self.coef_min)/(self.coef_max - self.coef_min)
        #plt.plot(x,y,'go',alpha=alpha)
        plt.plot(x,y,'go')
        plt.xlim((eval('self.%smin'%d1),eval('self.%smax'%d1)))
        plt.ylim((eval('self.%smin'%d2),eval('self.%smax'%d2)))
        plt.title('t=%0.3f'%self.t)
        if do_pause:
            plt.pause(.0001)


    def add_random_scatterer(self):
        x = np.random.rand()*(self.xmax-self.xmin)+self.xmin
        y = np.random.rand()*(self.ymax-self.ymin)+self.ymin
        z = np.random.rand()*(self.zmax-self.zmin)+self.zmin
        vx = np.random.randn()*1e-3
        vy = np.random.randn()*1e-3
        vz = np.random.randn()*1e-3
        coef = np.random.rand()*1e-6
        s = Scatterer(x,y,z,vx,vy,vz,coef=coef)
        self.add_scatterer(s)

    def update(self):
        self.t = self.t + self.dt
        for s in self.scatterers:
            s.update(self.dt,self.xmax,self.ymax,self.zmax)
    
if __name__=="__main__":

    s = Sample(1e-3,1e-3,1e-3,dt=1e-3)
    for k in range(100):
        s.add_random_scatterer()

    #s.plot3d()
    for k in range(1000):
        s.update()
        s.plot()
