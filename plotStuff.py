#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 12:58:37 2022

Some examples/ starting point to Load and plot Fringe and  MintPy outputs
Mintpy also has many useful plotting functions. 

@author: km
"""

import numpy as np
import glob
from osgeo import gdal
from SARTS import config,util,makeMap
from matplotlib import pyplot as plt
import h5py
import cartopy.crs as ccrs
import argparse
# import mintpy.utils as mpu
# from mintpy.objects.coord import coordinate
# from mintpy.utils import utils as ut
# from Scripts.GPS import getGpsTs, getGpsRates
from astropy.convolution import Gaussian2DKernel,convolve
import os

mintDir = './MintPy/'
ps = config.getPS()

# Geometry files
filename = mintDir + 'inputs/geometryRadar.h5'
ds = h5py.File(filename,'r+')   
lons = np.asarray(ds['longitude'])
lats = np.asarray(ds['latitude'])
hgt = np.asarray(ds['height'])
waterMask = np.asarray(ds['waterMask'])
ds.close()


# maskConnComp
filename = mintDir + 'maskConnComp.h5'
ds = h5py.File(filename,'r+')   
maskConnComp = np.asarray(ds['mask'])
ds.close()

# avgSpatialCoh
filename = mintDir + 'avgSpatialCoh.h5'
ds = h5py.File(filename,'r+')   
avgSpatialCoh = np.asarray(ds['coherence'])
ds.close()
avgSpatialCoh[waterMask<1] = np.nan

# Temporal Coh
filename = mintDir + 'temporalCoherence.h5'
ds = h5py.File(filename,'r+')   
temporalCoherence = np.asarray(ds['temporalCoherence'])
ds.close()
temporalCoherence[waterMask<1] = np.nan

filename = mintDir + 'velocity.h5' 
ds = h5py.File(filename,'r+')   
velocity = np.asarray(ds['velocity']) *1000 #convert to mm
ds.close()
velocity[waterMask<1] = np.nan
plt.scatter(1600,3113,s=40,color='red')

#Mask
thresh = .7
msk_time = np.zeros(velocity.shape)
msk_space= np.zeros(velocity.shape)
msk_time[temporalCoherence>thresh] += 1
msk_space[avgSpatialCoh>thresh] += 1
msk_water=waterMask.copy()
msk_concom=maskConnComp.copy()
msk_sum = msk_time + msk_space + msk_water + msk_concom
fig,ax = plt.subplots(2,2,figsize=(7,7))
ax[0,0].imshow(msk_time,cmap='magma');ax[0,0].set_title('msk_time')
ax[0,1].imshow(msk_space,cmap='magma');ax[0,1].set_title('msk_space')
ax[1,0].imshow(msk_water,cmap='magma');ax[1,0].set_title('msk_water')
ax[1,1].imshow(msk_concom,cmap='magma');ax[1,1].set_title('msk_concom')

ps = argparse.Namespace()
ps.nyl = velocity.shape[0]
ps.nxl = velocity.shape[1]
ps.minlat = lats.min()
ps.maxlat = lats.max()
ps.minlon = lons.min()
ps.maxlon = lons.max()
ymin,ymax,xmin,xmax = 0,ps.nyl,0,ps.nxl

msk = np.ones(velocity.shape)
msk[waterMask==0] = 0
msk[maskConnComp<1] =0
msk[avgSpatialCoh<0.7] = 0

phs_model = util.phaseElev(velocity, hgt,msk, 0, ps.nyl, 0, ps.nxl,makePlot=True)
#img, hgt,msk, ymin, ymax, xmin, xmax,makePlot=False


velocity_c = velocity-phs_model
velocity_c[msk==0] = np.nan
velocity[msk==0] = np.nan

vmin,vmax = -5,5
bg = 'World_Shaded_Relief'
pad=0
zoomLevel = 12
decimate=2

title = 'Velocity'
makeMap.mapImg(velocity[::decimate,::decimate], lons[::decimate,::decimate], lats[::decimate,::decimate], vmin,vmax, pad,zoomLevel, title, bg, cm='RdBu_r',plotFaults= False,alpha=1)

title='Velocity corrected'
makeMap.mapImg(velocity_c[::decimate,::decimate], lons[::decimate,::decimate], lats[::decimate,::decimate], vmin,vmax, pad,zoomLevel, title, bg, cm='RdBu_r',plotFaults= False,alpha=1)

title = 'Phase-elevation dependence'
phs_model[waterMask==0] = np.nan
makeMap.mapImg(phs_model[::decimate,::decimate], lons[::decimate,::decimate], lats[::decimate,::decimate], vmin,vmax, pad,zoomLevel, title, bg, cm='RdBu_r',plotFaults= False,alpha=1)

title = 'Mean spatial coherence'
makeMap.mapImg(avgSpatialCoh[::decimate,::decimate], lons[::decimate,::decimate], lats[::decimate,::decimate], .25,.95, pad,zoomLevel, title, bg, cm='magma',plotFaults= False,alpha=1)

#Spatial Filter
kernel = Gaussian2DKernel(x_stddev=1)
velFilt = convolve(velocity, kernel)
velFilt[msk==0] = np.nan
title='velocity'
makeMap.mapImg(velFilt[::decimate,::decimate], lons[::decimate,::decimate], lats[::decimate,::decimate],vmin,vmax, pad,zoomLevel, title, bg, cm='RdBu_r',plotFaults= False,alpha=1)
    
# timeseries.h5
filename = mintDir + 'timeseries.h5'
dsts = h5py.File(filename,'r+')  

vmin,vmax = -20,20
plt.figure();plt.imshow(velocity,vmin=vmin,vmax=vmax)

# Pixel coordinates you want to plot
py,px =  200,200
plt.figure();plt.plot(ps.dec_year,dsts['timeseries'][:,py,px],'.');plt.title('title')
