import numpy as np
from matplotlib import pyplot as plt
from oct_simulator import *
import sys,os
import shutil

noisy = 'noisy' in sys.argv
try:
    multiplier = float(sys.argv[-1])
except:
    multiplier = 1e-3
    
dt = 1e-5
xmax = 1e-4
ymax = 0
#zmax = 1.3e-3
zmax = 1e-4
npx = 1024
N_scatterers = 1000
dx = xmax/float(N_scatterers)
N_frames = 100
outdir = '/home/rjonnal/data/Dropbox/Private/oct_simulations/flow_example_frames'

do_simulation = True
stack_frames = False or do_simulation

source = GaussianSpectrum(1e-3,740e-9,940e-9,npx,100e-9)
sample = Sample(xmax,ymax,zmax,dt=dt)
vessel_start = xmax/3.0
vessel_end = xmax*2.0/3.0


approx_z_sampling = 1.3e-6
vz_base = approx_z_sampling/dt*multiplier



for k in range(N_scatterers):
    #x0 = np.random.rand()*xmax
    x0 = dx/2.0 + dx*k
    y0 = 0.0
    z0 = np.random.rand()*zmax
    #z0 = 1e-5
    vx = 0.0
    vy = 0.0
    vz = 0.0
    if vessel_start<x0<vessel_end:
        #vz = np.random.randn()*5e-5+5e-4
        vz = vz_base + np.random.randn()*vz_base*0.05
    sample.add_scatterer(Scatterer(x0,y0,z0,vx,vy,vz,coef=1e-3))

oct = OCT(source,sample,xmax=xmax,ymax=ymax,dt=dt,npx=npx,dx=xmax/100.0,noisy=True)
cutoff = round(zmax/oct.Sz/1.38)

if noisy:
    noise_string = 'noisy'
else:
    noise_string = 'noiseless'

vz_base_string = 'vz-%0.3e'%vz_base

outdir = '%s_%s_%s'%(outdir,noise_string,vz_base_string)

if os.path.exists(outdir):
    shutil.rmtree(outdir)

os.mkdir(outdir)



if do_simulation:
    for k in range(N_frames):
        b = oct.acquire_bscan()
        outfn = os.path.join(outdir,'%03d.npy'%k)
        print outfn
        np.save(outfn,b)
        b = b[npx/2:npx/2+cutoff,:]
        amp = np.abs(b)
        phase = np.angle(b)
        plt.subplot(1,2,1)
        plt.cla()
        plt.imshow(amp,interpolation='none',aspect='auto')
        plt.subplot(1,2,2)
        plt.cla()
        plt.imshow(phase,interpolation='none',aspect='auto')
        plt.pause(.0000001)

sys.exit()
def SNR(im):
    noise = np.std(np.hstack((im[:,0:30],im[:,70:])))
    return im/noise

def plot_SNR(im,method):
    snr = SNR(im)
    plt.figure()
    plt.imshow(snr,interpolation='none',aspect='auto')
    plt.colorbar()
    plt.title(method)
    fn = 'SNR_%s.png'%method.replace('(','_').replace(')','').replace(' ','_').replace('$','').replace('\\','')
    plt.savefig(fn)
    #plt.pause(1)
    #plt.close()
    plt.show()

def angular_statistics(phase_stack):
    """Return the angluar (circular) statistics for a stack of phasors:
       X, Y, r, cosa, sina, theta:
       X and Y are rectangular coordinates of the mean angle
       r is the length of the mean unit phasor
       cosa and sina are cosine and sine of the mean angle
       theta is the mean angle"""
    X = np.mean(np.cos(phase_stack),axis=0)
    Y = np.mean(np.sin(phase_stack),axis=0)
    r = np.sqrt(X**2+Y**2)
    cosa = X/r
    sina = Y/r
    theta = np.arctan(sina/cosa)
    return X,Y,r,cosa,sina,theta
    
if stack_frames:
    for iframe in range(N_frames):
        print iframe
        fn = os.path.join(outdir,'%03d.npy'%iframe)
        frame = np.load(fn)[npx/2:npx/2+cutoff,:]
        sz,sx = frame.shape
        if iframe==0:
            stack = np.zeros((N_frames,sz,sx),dtype=np.complex64)
        stack[iframe,:,:] = frame

    np.save(os.path.join(outdir,'stack.npy'),stack)
else:
    stack = np.load(os.path.join(outdir,'stack.npy'))
    stack = stack[:,1:,:]
    X, Y, r, cosa, sina, theta = angular_statistics(np.angle(stack))



    plt.subplot(2,2,1)
    im = np.abs(np.var(stack,axis=0))
    plt.imshow(SNR(im),interpolation='none',aspect='auto')
    plt.title('complex variance SNR')
    plt.colorbar()
    plt.subplot(2,2,2)
    im = np.var(np.angle(stack),axis=0)
    plt.imshow(SNR(im),interpolation='none',aspect='auto')
    plt.title('phase variance SNR')
    plt.colorbar()
    plt.subplot(2,2,3)
    im = np.var(np.abs(stack),axis=0)
    plt.imshow(SNR(im),interpolation='none',aspect='auto')
    plt.title('amplitude variance SNR')
    plt.colorbar()
    plt.subplot(2,2,4)
    plt.imshow(SNR(1.0-r),interpolation='none',aspect='auto')
    plt.title('circular variance SNR')
    plt.colorbar()
    plt.show()
    # cvar = np.abs(np.var(stack,axis=0))
    # plot_SNR(cvar,'absolute complex variance')

    # cvar = np.abs(np.var(stack[::3],axis=0))
    # plot_SNR(cvar,'absolute complex variance 1 of 3')
    
    # cvar = np.var(np.angle(stack),axis=0)
    # plot_SNR(cvar,'phase variance')

    # cvar = np.var(np.abs(stack),axis=0)
    # plot_SNR(cvar,'amplitude variance')

    # cvar = np.var(np.angle(stack[::3]),axis=0)
    # plot_SNR(cvar,'phase variance 1 of 3')

    # cvar = np.var(np.abs(stack[::3]),axis=0)
    # plot_SNR(cvar,'amplitude variance 1 of 3')

    # plot_SNR(X,'X')
    # plot_SNR(Y,'Y')
    # plot_SNR(r,'r')
    # plot_SNR(cosa,'cosa')
    # plot_SNR(sina,'sina')
    # plot_SNR(theta,'theta')

    plot_SNR(2*np.pi-r,'$2\pi - r$')
