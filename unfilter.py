#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 12:12:02 2023

1. unwrap the downlooked/filtered pairs, 
2. "upsample" them, 
3. add the same 2pi wraps to the full res unwrapped data,

@author: km
"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('svg')
import os
from osgeo import gdal
import localParams

ps = localParams.getLocalParams()
plt.close('all')

# Load the original IFG
fn_ifg = 'Fringe/PS_DS/sequential1/20160204_20160416/fine_lk.int.vrt'
ds = gdal.Open(fn_ifg)
fine_lk = ds.GetVirtualMemArray()
plt.figure();plt.imshow(np.angle(fine_lk),cmap='rainbow')

# Lod the filtered IFG
fn_ifg_filt = 'Fringe/PS_DS/sequential1/20160204_20160416/filt_lk.int.vrt'
ds = gdal.Open(fn_ifg_filt)
filt_lk = ds.GetVirtualMemArray()
plt.figure();plt.imshow(np.angle(filt_lk),cmap='rainbow')

# Load the filtered UNW 
fn_unw_filt = 'Fringe/PS_DS/sequential1/20160204_20160416/filt_lk.unw.vrt'
ds = gdal.Open(fn_unw_filt)
filt_lk_unw = ds.GetVirtualMemArray()
plt.figure();plt.imshow(filt_lk_unw[1,:,:],cmap='rainbow')

ds = None

wraps =  filt_lk_unw[1,:,:] - np.angle(filt_lk)
