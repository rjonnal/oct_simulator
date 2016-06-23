import numpy as np
from matplotlib import pyplot as plt
from oct_simulator import *
import sys,os
from movie import Movie

#x = (np.random.randn(30)*.01+1.0)*np.exp((.1*np.random.randn(30)*1j))
x = np.ones(30)*np.exp(np.pi/4.0*1j)

rootgroup = PhasorGroup(x)

pg_list = []

start_t = 50
n_steps = 1000

for k in range(n_steps):
    pg_list.append(PhasorGroup(rootgroup.complex_array))
    if k>=start_t:
        rootgroup.perturb()

circv_list = [pg.circular_variance() for pg in pg_list]
pv_list = [pg.phase_variance() for pg in pg_list]
compv_list = [pg.complex_variance() for pg in pg_list]

mags = [pg.magnitude for pg in pg_list]
maxmag = np.max(mags)

t = np.arange(n_steps) - start_t


#ax1 = plt.subplot(1,2,1)
#ax2 = plt.subplot(1,2,2)
l1 = .05
w1 = .5
b1 = .05
h1 = .85
l234 = .65
b2 = .65
b3 = .35
b4 = .05
h234 = .3
w234 = .3


t1 = t[0]
t2 = t[-1]

m = Movie('test.mp4')

def fixticks(ax):
    labels = ax.get_yticks()
    outlabels = []
    for lab in labels:
        try:
            outlabels.append('%0.1f'%lab)
        except:
            outlabels.append('')
    ax.set_yticklabels(outlabels)
        

for pg,now in zip(pg_list,t):
    fig = plt.figure(figsize=(9,6))
    ax1 = plt.axes([l1,b1,w1,h1])
    ax2 = plt.axes([l234,b2,w234,h234])
    ax3 = plt.axes([l234,b3,w234,h234])
    ax4 = plt.axes([l234,b4,w234,h234])

    ax_list = [ax1,ax2,ax3,ax4]

    pg.plot(ax1)
    ax1.set_xlim((-maxmag-.5,maxmag+.5))
    ax1.set_ylim((-maxmag-.5,maxmag+.5))
    
    ax2.plot(t,circv_list,'r')
    ax2.set_ylabel('circular variance')
    ax2.set_xticks([])
    ax2.set_xlim((t1,t2))

    ax3.plot(t,pv_list,'g')
    ax3.set_ylabel('phase variance')
    ax3.set_xticks([])
    ax3.set_yticks(ax3.get_yticks()[:-1])
    ax3.set_xlim((t1,t2))

    ax4.plot(t,compv_list,'b')
    ax4.set_xlim((t1,t2))
    ax4.set_yticks(ax4.get_yticks()[:-2])
    ax4.set_ylabel('complex variance')

    for ax in ax_list[1:]:
        ax.axvline(now,c=[.5,.5,.5])
        fixticks(ax)
        
    m.add(fig)

    plt.close()
    print now


m.make()
sys.exit()

p1 = Phasor(1.0,np.pi/4.0)
p2 = Phasor(1.0,np.pi/3.0)
ax = plt.axes()
ax.set_xlim((-1,1))
ax.set_ylim((-1,1))
p1.plot_context(ax)
p1.plot(ax)
p2.plot(ax)
plt.show()
sys.exit()
        
        
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
