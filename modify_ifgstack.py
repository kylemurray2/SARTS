#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 14:56:55 2024

@author: km
"""

import numpy as np
from matplotlib import pyplot as plt
import os
import h5py



filename = 'MintPy/inputs/geometryRadar.h5'
ds = h5py.File(filename,'r+')
watermask = np.asarray(ds['waterMask'])
ds.close()

ifg_stack_file = 'MintPy/inputs/ifgramStack.h5'
ds = h5py.File(ifg_stack_file,'r+')
dropIfgram = np.asarray(ds['dropIfgram'])
unw_stack = np.asarray(ds['unwrapPhase'])
connectComponent = np.asarray(ds['connectComponent'])

# unw_stack_clean = unw_stack[dropIfgram,:,:]

# for ii in range(unw_stack_clean.shape[0]):
#     plt.figure()
#     plt.imshow(connectComponent[ii,:,:])
#     plt.title(ii)
    
mask = (connectComponent == 0) | (connectComponent == 1)
connectComponent[mask] = 1
connectComponent[:,~watermask] = 0

    

ds['connectComponent'][...] = connectComponent

ds.close()

# for ii in range(connectComponent.shape[0]):
#     plt.figure()
#     plt.imshow(connectComponent[ii,:,:])
#     plt.title(ii)