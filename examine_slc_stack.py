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
from SARTS import config

ps = config.getPS() 


ds = gdal.Open('./Fringe/coreg_stack/slcs_base.vrt')
stack = ds.GetVirtualMemArray()

px = int(ps.nx/2)
py = int(ps.ny/2)

ts = np.angle(stack[:,py,px])
np.where(ts==0)
plt.figure()
plt.plot(ts,'.')


awds_stack = []
for d in ps.dates:
    print(d)
    ds = gdal.Open('./Fringe/adjusted_wrapped_DS/'+ d + '.slc.vrt')
    awds_stack.append(ds.GetVirtualMemArray()[500,12323])
ts = np.asarray(awds_stack)
plt.figure()
plt.plot(ts,'.')

slc = np.angle(stack[200,:,:])
np.where(slc==0)

plt.figure()
plt.imshow(np.angle(stack[200,:,:]))