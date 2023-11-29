#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 10:34:47 2023

@author: km
"""

import numpy as np
from matplotlib import pyplot as plt
import os
from osgeo import gdal

#
ds = gdal.Open('./Fringe/coreg_stack/slcs_base.vrt')
stack = ds.GetVirtualMemArray()

ts = np.angle(stack[:,500,9000])
np.where(ts==0)
plt.figure()
plt.plot(ts,'.')


slc = np.angle(stack[200,:,:])
np.where(slc==0)

plt.figure()
plt.imshow(np.angle(stack[200,:,:]))