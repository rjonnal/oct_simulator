import numpy as np
from matplotlib import pyplot as plt
from oct_simulator import *
import sys,os

dt = 1e-5
xmax = 1e-4
ymax = 0
#zmax = 1.3e-3
zmax = 1e-4
npx = 1024
N_scatterers = 1000
dx = xmax/float(N_scatterers)
N_frames = 30
outdir = './flow_example_frames'

do_simulation = False
stack_frames = False or do_simulation

source = GaussianSpectrum(1e-3,740e-9,940e-9,npx,100e-9)
sample = Sample(xmax,ymax,zmax,dt=dt)
vessel_start = xmax/3.0
vessel_end = xmax*2.0/3.0

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
        vz = np.random.randn()*5e-5+5e-4
    sample.add_scatterer(Scatterer(x0,y0,z0,vx,vy,vz,coef=1e-3))

oct = OCT(source,sample,xmax=xmax,ymax=ymax,dt=dt,npx=npx,dx=xmax/100.0)
cutoff = round(zmax/oct.Sz/1.38)

if do_simulation:
    for k in range(N_frames):
        b = oct.acquire_bscan()
        outfn = os.path.join(outdir,'%03d.npy'%k)
        print outfn
        np.save(outfn,b)
        continue
        b = b[npx/2:npx/2+cutoff,:]
        amp = np.abs(b)
        phase = np.angle(b)
        plt.figure()
        plt.subplot(1,2,1)
        plt.imshow(amp,interpolation='none',aspect='auto')
        #plt.colorbar()
        plt.subplot(1,2,2)
        plt.imshow(phase,interpolation='none',aspect='auto')
        #plt.colorbar()
    #plt.show()


def SNR(im,method):
    noise = np.std(np.hstack((im[:,0:30],im[:,70:])))
    plt.figure()
    plt.imshow(im/noise,interpolation='none',aspect='auto')
    plt.colorbar()
    plt.title(method)
    fn = 'SNR_%s.png'%method.replace('(','_').replace(')','').replace(' ','_')
    plt.savefig(fn)
    #plt.pause(1)
    #plt.close()
    plt.show()
    
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

    cvar = np.abs(np.var(stack,axis=0))
    SNR(cvar,'abs(var(c))')

    cvar = np.abs(np.var(stack[::3],axis=0))
    SNR(cvar,'abs(var(c)) skip 2')
    
    cvar = np.var(np.angle(stack),axis=0)
    SNR(cvar,'var(phase)')

    cvar = np.var(np.abs(stack),axis=0)
    SNR(cvar,'var(amp)')

    cvar = np.var(np.angle(stack[::3]),axis=0)
    SNR(cvar,'var(phase) skip 2')

    cvar = np.var(np.abs(stack[::3]),axis=0)
    SNR(cvar,'var(amp) skip 2')

