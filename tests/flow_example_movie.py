import numpy as np
from matplotlib import pyplot as plt
from oct_simulator import *
import sys,os
from movie import Movie

frame_dir = './flow_example_slow_frames'
N_frames = 30
cutoff = 60
npx = 1024

try:
    stack = np.load(os.path.join(frame_dir,'stack.npy'))
except IOError:
    for iframe in range(N_frames):
        print iframe
        fn = os.path.join(frame_dir,'%03d.npy'%iframe)
        frame = np.load(fn)
        frame = frame[npx/2:npx/2+cutoff,:]
        sz,sx = frame.shape
        if iframe==0:
            stack = np.zeros((N_frames,sz,sx),dtype=np.complex64)
        stack[iframe,:,:] = frame
    np.save(os.path.join(frame_dir,'stack.npy'),stack)


astack = np.abs(stack)
smax = np.max(stack)
smin = np.min(stack)

mov = Movie('oct_slow_flow.avi')
fig = plt.figure()

for k in range(N_frames):
    plt.cla()
    plt.imshow(astack[k,:,:],interpolation='none',aspect='auto')
    plt.clim((smin,smax))
    if k==0:
        plt.colorbar()
    plt.pause(1)

sys.exit()
    


plt.figure()
plt.imshow(np.abs(stack[0,:,:]))
plt.figure()
plt.imshow(np.abs(stack[1,:,:]))
plt.show()
sys.exit()
sz,sy,sx = stack.shape

intervals = range(1,10)

variances = []


def cvar(mat,axis=0):
    a = np.angle(mat)
    X = np.cos(a)
    Y = np.sin(a)
    Xa = np.mean(X,axis=axis)
    Ya = np.mean(Y,axis=axis)
    mag = np.sqrt(Xa**2+Ya**2)
    return 1 - mag

for interval in intervals:
    temp = np.zeros((sz-interval,sy,sx),dtype=np.complex64)
    for idx2 in range(interval,sz):
        idx1 = idx2-interval
        f2 = stack[idx2,:,:]
        f1 = stack[idx1,:,:]
        df = f2 - f1
        temp[idx1,:,:] = df
    temp = np.abs(np.diff(temp,axis=0))
    v = np.mean(temp,axis=0)
    plt.imshow(v)
    plt.colorbar()
    plt.show()
        
        
        
