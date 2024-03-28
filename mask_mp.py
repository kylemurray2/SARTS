#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 15:10:44 2023

Vertical/horizontal decomposition for Oahu

Resamples descending velocities to the ascending grid

@author: km
"""

import numpy as np
import h5py 
from scipy.interpolate import griddata
from matplotlib import pyplot as plt
from SARTS import util, config, bridge, makeMap
import os


def get_mask(mintPy_dir='MintPy',spatial_coh_thresh  = 0.6,temporal_coh_thresh = 0.8,numtriplets_thresh_ratio  = 0.05,flipud=False):  
    '''
    Makes a final time series/velocity mask for MintPy workflow
    Masks based on:
        spatial coherence
        temporal coherence
        non-zero phase closure (ratio of ifgs)
        water mask
        shadow mask

    '''
    # conn comp files 
    filename = os.path.join(mintPy_dir, 'inputs/ifgramStack.h5')
    ds = h5py.File(filename,'r+')
    concomp_stack = np.asarray(ds['connectComponent'])
    concomp_stack = np.asarray(concomp_stack)
    concomp_stack[concomp_stack>1] = 0
    cc_frac = np.count_nonzero(concomp_stack,axis=0)/concomp_stack.shape[0]
    plt.figure();plt.imshow(cc_frac,vmin=.5,vmax=1)
    # a = np.zeros_like(cc_frac)
    # a [cc_frac > .9 ] = 1
    # plt.figure();plt.imshow(a)

    # shadow and water mask
    filename = os.path.join(mintPy_dir, 'inputs/geometryRadar.h5')
    ds = h5py.File(filename,'r+')
    waterMask_a = np.asarray(ds['waterMask'])
    shadowMask =  np.asarray(ds['shadowMask'])
    ds.close()
    
    # avgSpatialCoh
    filename = os.path.join(mintPy_dir, 'avgSpatialCoh.h5')
    ds = h5py.File(filename,'r+')   
    avgSpatialCoh = np.asarray(ds['coherence'])
    ds.close()
    
    # temporalCoherence
    filename = os.path.join(mintPy_dir, 'temporalCoherence.h5')
    ds = h5py.File(filename,'r+')   
    temporalCoherence = np.asarray(ds['temporalCoherence'])
    ds.close()
    
    # maskConnComp
    filename = os.path.join(mintPy_dir, 'maskConnComp.h5')
    ds = h5py.File(filename,'r+')   
    maskConnComp_a = np.asarray(ds['mask'])
    ds.close()

    # Triplet non zero
    filename = os.path.join(mintPy_dir, 'numTriNonzeroIntAmbiguity.h5')
    ds = h5py.File(filename,'r+')   
    numTriNonzeroIntAmbiguity_a = np.asarray(ds['mask']) 
    ds.close()
    plt.figure();plt.imshow(np.flipud(numTriNonzeroIntAmbiguity_a),vmin=None,vmax=200)
    
    # Make the other masks
    numtriplets_thresh = int(concomp_stack.shape[0]*numtriplets_thresh_ratio)
    msk_spatial_coh = avgSpatialCoh > spatial_coh_thresh
    msk_temporal_coh = temporalCoherence > temporal_coh_thresh 
    msk_numtriplets = numTriNonzeroIntAmbiguity_a < numtriplets_thresh
    
    # Combine into final full maks
    stacked_masks = np.array([msk_spatial_coh, msk_temporal_coh, msk_numtriplets,~shadowMask,waterMask_a])
    mask_full = np.all(stacked_masks, axis=0)
    
    
    if flipud:
        msk_spatial_coh = np.flipud(msk_spatial_coh)
        msk_temporal_coh = np.flipud(msk_temporal_coh)
        msk_numtriplets = np.flipud(msk_numtriplets)
        waterMask_a = np.flipud(waterMask_a)
        shadowMask = np.flipud(shadowMask)
        mask_full = np.flipud(mask_full)
    
    fig,ax = plt.subplots(3,2)
    ax[0,0].imshow(msk_spatial_coh,cmap='gray_r');ax[0,0].set_title('Spatial coherence')
    ax[0,1].imshow(msk_temporal_coh,cmap='gray_r');ax[0,1].set_title('Temporal coherence')
    ax[1,0].imshow(msk_numtriplets,cmap='gray_r');ax[1,0].set_title('Non-zero phase closure')
    ax[1,1].imshow(waterMask_a,cmap='gray_r');ax[1,1].set_title('Water')
    ax[2,0].imshow(~shadowMask,cmap='gray_r');ax[2,0].set_title('Shadows')
    ax[2,1].imshow(mask_full,cmap='gray_r');ax[2,1].set_title('Mask all')
    for i in range(3):
        for j in range(2):
            ax[i, j].set_xticks([])  # Remove x-axis tick marks
            ax[i, j].set_yticks([])  # Remove y-axis tick marks
            ax[i, j].set_xticklabels([])  # Remove x-axis tick labels
            ax[i, j].set_yticklabels([])  # Remove y-axis tick labels
    plt.tight_layout()
    plt.show()
    
    return mask_full

# Test
# mintPy_dir='MintPy_1_4_filtfull'
# spatial_coh_thresh  = 0.6
# temporal_coh_thresh = 0.8
# numtriplets_thresh_ratio  = 0.05
# flipud=True
# mask_full = get_mask(mintPy_dir=mintPy_dir,spatial_coh_thresh  = spatial_coh_thresh,temporal_coh_thresh =temporal_coh_thresh,numtriplets_thresh_ratio  =numtriplets_thresh_ratio,flipud=flipud)
