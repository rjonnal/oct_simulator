import numpy as np
from matplotlib import pyplot as plt

class Phasor:

    colors = 'rgbkcmy'
    pid = 0
    
    def __init__(self,magnitude=1.0,angle=0.0):
        self.magnitude = magnitude
        self.angle = angle
        self.x = np.cos(self.angle)*self.magnitude
        self.y = np.sin(self.angle)*self.magnitude
        self.pid = Phasor.pid
        Phasor.pid = Phasor.pid + 1
        self.color = Phasor.colors[self.pid%len(Phasor.colors)]
        
    def plot(self,ax):
        hl = 0.1
        px = np.cos(self.angle)*(self.magnitude-hl)
        py = np.sin(self.angle)*(self.magnitude-hl)
        ax.arrow(0.0, 0.0, px, py, head_width=0.05, head_length=hl, fc=self.color, ec=self.color)
        xlim = list(ax.get_xlim())
        xlim[0] = min(xlim[0],-self.magnitude)
        xlim[1] = max(xlim[1],self.magnitude)
        ylim = list(ax.get_ylim())
        ylim[0] = min(ylim[0],-self.magnitude)
        ylim[1] = max(ylim[1],self.magnitude)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_aspect('equal')

    def plot_context(self,ax):
        ax.add_artist(plt.Circle((0,0),1.0,ec=[0.5,0.5,0.5],fill=False))
        for ang in np.arange(0,2,0.25):
            x = np.cos(ang*np.pi)
            y = np.sin(ang*np.pi)
            plt.text(x,y,'%0.2f $\pi$'%ang)
        
class PhasorGroup:

    colors = 'rgbkcmy'
    
    def __init__(self,complex_array=np.array([])):
        self.complex_array = complex_array
        self.magnitude = np.abs(complex_array)
        self.phase = np.angle(complex_array)%(2*np.pi)
        self.x = np.real(complex_array)
        self.y = np.imag(complex_array)

    def circular_variance(self):
        X = np.mean(np.cos(self.phase))
        Y = np.mean(np.sin(self.phase))
        return 1.0 - np.sqrt(X**2 + Y**2)

    def phase_variance(self):
        return np.var(self.phase)

    def complex_variance(self):
        return np.abs(np.var(self.complex_array))

    def perturb(self,magnitude_std=.02,phase_std=np.pi*.02):
        self.phase = self.phase + np.random.standard_normal(self.phase.shape)*phase_std
        self.magnitude = self.magnitude + np.random.standard_normal(self.magnitude.shape)*magnitude_std
        neg_idx = np.where(self.magnitude<0)[0]
        self.magnitude[neg_idx] = -self.magnitude[neg_idx]
        self.phase[neg_idx] = (self.phase[neg_idx]+np.pi)%(2*np.pi)
        self.complex_array = self.magnitude*np.exp(self.phase*1j)
        self.x = np.real(self.complex_array)
        self.y = np.imag(self.complex_array)
        
    def get_text_alignment(self,angle):
        if angle==0.0:
            ha = 'left'
            va = 'center'
        elif 0.0<angle<np.pi/2:
            ha = 'left'
            va = 'bottom'
        elif angle==np.pi/2:
            ha = 'center'
            va = 'bottom'
        elif np.pi/2<angle<np.pi:
            ha = 'right'
            va = 'bottom'
        elif angle==np.pi:
            ha = 'right'
            va = 'center'
        elif np.pi<angle<1.5*np.pi:
            ha = 'right'
            va = 'top'
        elif angle==1.5*np.pi:
            ha = 'center'
            va = 'top'
        else:
            ha = 'left'
            va = 'top'
        return ha,va
    
    def plot(self,ax,pad=1.0,markup_color=[0.8,0.8,0.8],print_coords=False):
        rad = np.max(self.magnitude)
        for k in np.arange(1.0,rad,1.0):
            ax.add_artist(plt.Circle((0,0),k,ec=markup_color,fill=False))
        for ang in np.arange(0,2,0.25):
            x = np.cos(ang*np.pi)*rad
            y = np.sin(ang*np.pi)*rad
            ha,va = self.get_text_alignment(ang*np.pi)
            ax.text(x,y,'%0.2f $\pi$'%ang,ha=ha,va=va,color=markup_color)
        xlim = [-rad-pad,rad+pad]
        ylim = [-rad-pad,rad+pad]
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_aspect('equal')

        hl = 0.1
        for angle,magnitude in zip(self.phase.ravel(),self.magnitude.ravel()):
            color = [angle/2.0/np.pi,(magnitude%2.0)/2.0,(angle+magnitude)%1.0]
            px = np.cos(angle)*(magnitude-hl)
            py = np.sin(angle)*(magnitude-hl)
            ax.arrow(0.0, 0.0, px, py, head_width=0.05, head_length=hl, fc=color, ec=color)
            if print_coords:
                ax.text(px,py,'(%0.1f,%0.1f)'%(px,py),color=color)


