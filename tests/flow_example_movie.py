import numpy as np
from matplotlib import pyplot as plt
from oct_simulator import *
import sys,os
from movie import Movie
from glob import glob

root_dir = '/home/rjonnal/data/Dropbox/Private/oct_simulations/*'

dlist = glob(root_dir)

for frame_dir in dlist:
    junk,tail = os.path.split(frame_dir)
    tag = tail.replace('flow_example_frames_','')
    N_frames = 100
    cutoff = 60
    npx = 1024

    for iframe in range(N_frames):
        print iframe
        fn = os.path.join(frame_dir,'%03d.npy'%iframe)
        frame = np.load(fn)
        frame = frame[npx/2:npx/2+cutoff,:]
        sz,sx = frame.shape
        if iframe==0:
            stack = np.zeros((N_frames,sz,sx),dtype=np.complex64)
        stack[iframe,:,:] = frame


    astack = np.abs(stack)
    smax = np.max(stack)
    smin = np.min(stack)

    avifn = os.path.join(frame_dir,'frames.avi')
    make_movie = not os.path.exists(avifn)

    if make_movie:
        mov = Movie(avifn)
        fig = plt.figure()

        for k in range(N_frames):
            plt.cla()
            plt.imshow(astack[k,:,:],interpolation='none',aspect='auto')
            plt.clim((smin,smax))
            if k==0:
                plt.colorbar()
            mov.add(fig)
            plt.pause(.1)

        mov.make()


    sz,sy,sx = stack.shape

    intervals = range(1,80)

    variances = []
    control_variances = []

    mads = []
    control_mads = []

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
        mad = np.abs(np.diff(temp,axis=0))
        mad = np.mean(mad,axis=0)

        v = np.var(temp,axis=0)
        
        mads.append(np.mean(mad[:,34:66]))
        control_mads.append(np.mean(mad[:,:33]))
        
        variances.append(np.sqrt(np.mean(v[:,34:66])))
        control_variances.append(np.sqrt(np.mean(v[:,:33])))

    plt.figure(figsize=(6,8))
    
    plt.subplot(2,1,1)
    plt.plot(intervals,variances,'g-',label='flow')
    plt.plot(intervals,control_variances,'k:',label='no flow',alpha=0.8)
    plt.xlabel('interval (# of B-scans)')
    plt.ylabel('RMS')
    plt.title(tag)
    plt.legend()
    #plt.savefig('./flow_variance_plots/%s_vista_variance.png'%tag)
    #plt.close()
    
    plt.subplot(2,1,2)
    plt.plot(intervals,mads,'g-',label='flow')
    plt.plot(intervals,control_mads,'k:',label='no flow',alpha=0.8)
    plt.xlabel('interval (# of B-scans)')
    plt.ylabel('MAD')
    #plt.title(tag)
    plt.legend()
    
    plt.savefig('./flow_variance_plots/%s.png'%tag)
    plt.close()
