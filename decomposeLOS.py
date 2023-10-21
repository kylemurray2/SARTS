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
from SARTS import util, config, bridge
# import cartopy.crs as ccrs
from osgeo import gdal
from astropy.convolution import Gaussian2DKernel,convolve
import os

island='Oahu'
doBridge=True

ascDir = os.path.join('/d/HI/S1/Asc', island)
desDir = os.path.join('/d/HI/S1/Des', island)

os.chdir('/d/HI/S1/verthorz/'+island)

ps_a = config.getPS(ascDir)
ps_d = config.getPS(desDir)

mpdir_a = os.path.join(ascDir, 'MintPy_delaunay_cohthresh')
mpdir_d = os.path.join(desDir, 'MintPy_delaunay_cohthresh')

filename = os.path.join(mpdir_a, 'inputs/geometryRadar.h5')
ds = h5py.File(filename,'r+')
lons_a = np.asarray(ds['longitude'])
lats_a = np.asarray(ds['latitude'])
az_a =  np.asarray(ds['azimuthAngle'])
inc_a =  np.asarray(ds['incidenceAngle'])
waterMask_a = np.asarray(ds['waterMask'])
ds.close()

#Land cover
filename = os.path.join(ascDir, 'merged/geom_reference/landCover_lk.rdr.vrt')
ds = gdal.Open(filename)   
lc = ds.GetVirtualMemArray()
plt.figure();plt.imshow(lc,cmap='jet')

# avgSpatialCoh
filename = os.path.join(mpdir_a, 'avgSpatialCoh.h5')
ds = h5py.File(filename,'r+')   
avgSpatialCoh_a = np.asarray(ds['coherence'])
ds.close()

# temporalCoherence
filename = os.path.join(mpdir_a, 'temporalCoherence.h5')
ds = h5py.File(filename,'r+')   
temporalCoherence_a = np.asarray(ds['temporalCoherence'])
ds.close()

# maskConnComp
filename = os.path.join(mpdir_a, 'maskConnComp.h5')
ds = h5py.File(filename,'r+')   
maskConnComp_a = np.asarray(ds['mask'])
ds.close()

# maskConnComp
filename = os.path.join(mpdir_d, 'maskConnComp.h5')
ds = h5py.File(filename,'r+')   
maskConnComp_d = np.asarray(ds['mask'])
ds.close()

filename = os.path.join(mpdir_a, 'velocity.h5')
ds = h5py.File(filename,'r+')   
velocity_a = np.asarray(ds['velocity']) *1000 #convert to mm
ds.close()





filename    = os.path.join(mpdir_d, 'inputs/geometryRadar.h5')
ds          = h5py.File(filename,'r+')   
lons_d      = np.asarray(ds['longitude'])
lats_d      = np.asarray(ds['latitude'])
az_d        = np.asarray(ds['azimuthAngle'])
inc_d       = np.asarray(ds['incidenceAngle'])
waterMask_d = np.asarray(ds['waterMask'])
ds.close()

filename    = os.path.join(mpdir_d, 'velocity.h5')
ds          = h5py.File(filename,'r+')   
velocity_d  = np.asarray(ds['velocity']) *1000 #convert to mm
ds.close()

# velocity_d[np.isnan(velocity_d)] = 0
# velocity_a[np.isnan(velocity_a)] = 0

vmin,vmax = -8,8
plt.figure();plt.imshow(velocity_a,vmin=vmin,vmax=vmax)
plt.figure();plt.imshow(velocity_d,vmin=vmin,vmax=vmax)

# velocity_d[maskConnComp_d==0] = np.nan
# velocity_a[maskConnComp_a==0] = np.nan
# plt.figure();plt.imshow(velocity_d,vmin=-10,vmax=10)
# plt.figure();plt.imshow(velocity_a,vmin=-10,vmax=10)
# plt.show()
# decimate=2;vmin=-10;vmax=10;pad=0;zoomLevel=12;title='veld';bg='World_Shaded_Relief'
# makeMap.mapImg(velocity_a[::decimate,::decimate], lons_a[::decimate,::decimate], lats_a[::decimate,::decimate], vmin,vmax, pad,zoomLevel, title, bg, cm='RdBu_r',plotFaults= False,alpha=1)

velocity_d_resamp = griddata((lons_d.ravel(),lats_d.ravel()), velocity_d.ravel(), (lons_a,lats_a), method='linear')
az_d_resamp = griddata((lons_d.ravel(),lats_d.ravel()), az_d.ravel(), (lons_a,lats_a), method='linear')
inc_d_resamp = griddata((lons_d.ravel(),lats_d.ravel()), inc_d.ravel(), (lons_a,lats_a), method='linear')

if doBridge:
    minPix = 4000
    velocity_a_trim,mask_a_trim = bridge.main(velocity_a, maskConnComp_a, minPix)
    velocity_d_resamp_trim,mask_a_trim = bridge.main(velocity_d_resamp, maskConnComp_a, minPix)

# mask with asc conncomp mask
velocity_d_resamp_trim_msk = velocity_d_resamp_trim.copy()
velocity_d_resamp_trim_msk[~mask_a_trim] = np.nan
velocity_a_trim_msk = velocity_a_trim.copy()
velocity_a_trim_msk[~mask_a_trim] = np.nan


plt.figure();plt.imshow(maskConnComp_d)
plt.figure();plt.imshow(mask_a_trim)
plt.show()

# Light spatial filter
kernel = Gaussian2DKernel(x_stddev=.2)
velocity_a_trim_msk_filt = convolve(velocity_a_trim_msk, kernel)
velocity_d_resamp_trim_msk_filt = convolve(velocity_d_resamp_trim_msk, kernel)

plt.figure();plt.imshow(velocity_a_trim_msk_filt,vmin=-6,vmax=6)
plt.figure();plt.imshow(velocity_d_resamp_trim_msk_filt,vmin=-6,vmax=6)
plt.figure();plt.imshow(velocity_d_resamp,vmin=-26,vmax=26)

npix = len(velocity_a_trim_msk_filt.ravel())

# smooth  = .4
# vertHorVec = np.zeros((npix,2))
# for ii in range(npix):
#     if ii%200000 == 0:
#         print("\r" + str(100*np.round((ii/npix),2)) + '%', end='', flush=True)
#     asc = velocity_a_trim_msk_filt.ravel()[ii]
#     des = velocity_d_resamp_trim_msk_filt.ravel()[ii]
#     psi_a = az_a.ravel()[ii]
#     psi_d = az_d_resamp.ravel()[ii]
#     theta_a = inc_a.ravel()[ii]
#     theta_d = inc_d_resamp.ravel()[ii]
#     vertHorVec[ii,:] = util.invertVertHor(asc,des,psi_a,theta_a,psi_d,theta_d,smooth)
    
# east = vertHorVec[:,0].reshape(velocity_a.shape)
# vert = vertHorVec[:,1].reshape(velocity_a.shape)

# if not os.path.isdir('Npy'):
#     os.mkdir('Npy')
# np.save(ps_d.workdir + '/Npy/vert2.npy',vert)
# np.save(ps_d.workdir + '/Npy/east2.npy',east)

vert = np.load(ps_d.workdir + '/Npy/vert2.npy')
east = np.load(ps_d.workdir + '/Npy/east2.npy')

eastm = east.copy()
vertm = vert.copy()

#Mask
thresh = .7
msk_time = np.zeros(east.shape)
msk_space= np.zeros(east.shape)
msk_veg =np.zeros(east.shape)
msk_time[temporalCoherence_a>thresh] += 1
msk_space[avgSpatialCoh_a>thresh] += 1
msk_water=waterMask_a.copy()
msk_concom=maskConnComp_a.copy()
msk_veg[lc==2] = 1
msk_veg[lc==3] = 1 
msk_veg[lc==4] = 1 
msk_veg[lc==5] = 1 
msk_veg[lc==6] = 1 
msk_veg[lc==12] = 1 
msk_veg[lc==19] = 1 
msk_veg[lc==20] = 1 

msk_sum = msk_time + msk_space + msk_water + msk_concom
fig,ax = plt.subplots(2,2,figsize=(7,7))
ax[0,0].imshow(msk_time,cmap='magma');ax[0,0].set_title('msk_time')
ax[0,1].imshow(msk_space,cmap='magma');ax[0,1].set_title('msk_space')
ax[1,0].imshow(msk_veg,cmap='magma');ax[1,0].set_title('msk_vegetation')
ax[1,1].imshow(msk_concom,cmap='magma');ax[1,1].set_title('msk_concom')

plt.figure()
plt.imshow(msk_sum,cmap='magma');plt.title('mask sum')

msk = np.ones(vert.shape)
msk[msk_sum<3] = 0

vmin,vmax = -20,20

from skimage.morphology import remove_small_objects
minimumPixelsInRegion=6000
mask = msk==1
msk_trim = remove_small_objects(mask, minimumPixelsInRegion, connectivity=1)
plt.figure();plt.imshow(msk_concom);plt.title('Conn comp')
plt.figure();plt.imshow(msk_trim);plt.title('Final mask')
plt.show()
# eastm[msk_sum<3] = np.nan
# vertm[msk_sum<3] = np.nan
# velocity_a[msk_sum<3] = np.nan
# vertm[msk_veg==0] = np.nan

vertm[msk_trim==0] = np.nan



# velocity_d_resamp[msk_sum<3] = np.nan
vmin,vmax = -10,10
kernel = Gaussian2DKernel(x_stddev=.2)
vertm_filt = convolve(vertm, kernel)
vertm_filt[waterMask_a==0] = np.nan

plt.figure();plt.imshow(vertm,vmin=vmin,vmax=vmax)
plt.figure();plt.imshow(vertm_filt,vmin=vmin,vmax=vmax)

# #remove mean from each disconnected region
# minPix = 1000
# labels = util.getConCom(msk,minPix)
# fig,ax = plt.subplots(2,1,figsize=(5,6))
# ax[0].imshow(msk);ax[0].set_title('mask')
# ax[1].imshow(labels);ax[1].set_title('connected regions')

# vertm_filt_bridged = vertm_filt.copy()

# for ii in range(int(labels.max())):
#     if len(vertm_filt_bridged[labels==ii+1]) < minPix:
#         vertm_filt_bridged[labels==ii+1] = np.nan # mask out small islands of data
#         msk[labels==ii+1] = 0
#     else:
#         vertm_filt_bridged[labels==ii+1]-=np.nanmean(vertm_filt_bridged[labels==ii+1])


ps_a.minlat = lats_a.min()
ps_a.minlon = lons_a.min()
ps_a.maxlat = lats_a.max()
ps_a.maxlon = lons_a.max()

vertGeo = util.geocodeKM(vertm_filt,10,lons_a,lats_a,ps_a)
plt.figure();plt.imshow(vertGeo,vmin=-20,vmax=0)
util.writeGeotiff(vertGeo, (ps_a.minlat,ps_a.maxlat), (ps_a.minlon,ps_a.maxlon), 'geotifs/OahuVertical_1017.tif')
util.geotiff_to_kmz('geotifs/OahuVertical_1017.tif', 'geotifs/OahuVertical_1017.kmz',vmin=-10,vmax=10,average_elevation=3,colormap='RdBu')


fig,ax = plt.subplots(2,2)
ax[0,0].imshow(velocity_a,vmin=vmin,vmax=vmax);ax[0,0].set_title('asc')
# ax[0,1].imshow(velocity_d_resamp,vmin=vmin,vmax=vmax);ax[0,1].set_title('des')
ax[1,0].imshow(vertm,vmin=vmin,vmax=vmax);ax[1,0].set_title('vertical')
ax[1,1].imshow(eastm,vmin=vmin,vmax=vmax);ax[1,1].set_title('east')

# ymin,ymax,xmin,xmax = 500,1500,1000,4000
ymin,ymax,xmin,xmax = 0,ps_d.nyl,0,ps_d.nxl

bg = 'World_Imagery' #'World_Shaded_Relief'
zoom=15
pad=0.0
vmin,vmax=-6,6

plt.figure();plt.imshow(vertm_filt,vmin=vmin,vmax=vmax)
plt.figure();plt.imshow(vertm_filt_bridged,vmin=vmin,vmax=vmax)

# makeMap.mapBackground(bg, ps_a.minlon,  ps_a.maxlon,  ps_a.minlat,  ps_a.maxlat, 0, zoom, title, scalebar=1, borders=True,plotFaults=True)
# img_handle = plt.imshow(np.flipud(east),vmin=vmin,vmax=vmax,cmap='RdBu_r', origin='upper', zorder=12, transform=ccrs.PlateCarree(), extent=[ps_a.minlon-pad, ps_a.maxlon+pad, ps_a.minlat-pad, ps_a.maxlat+pad])
# plt.colorbar(img_handle,fraction=0.03, pad=0.05,orientation='horizontal')
# plt.scatter( -117.31, 35.752, s=6,color='black',transform=ccrs.PlateCarree())
# plt.scatter( -117.31, 35.752, s=1,color='white',transform=ccrs.PlateCarree())
from PyPS2 import makeMap
title='East'
makeMap.mapImg(eastm[ymin:ymax,xmin:xmax], lons_a[ymin:ymax,xmin:xmax], lats_a[ymin:ymax,xmin:xmax], vmin, vmax, pad,zoom, title, bg=bg, cm='RdBu_r', plotFaults= False,alpha=1,contour=False)
title='Vertical'
makeMap.mapImg(vertm_filt[ymin:ymax,xmin:xmax], lons_a[ymin:ymax,xmin:xmax], lats_a[ymin:ymax,xmin:xmax], vmin, vmax, pad,zoom, title, bg=bg, cm='RdBu_r', plotFaults= False,alpha=1,contour=False)
makeMap.mapImg(vertm_filt_bridged[ymin:ymax,xmin:xmax], lons_a[ymin:ymax,xmin:xmax], lats_a[ymin:ymax,xmin:xmax], vmin, vmax, pad,zoom, title, bg=bg, cm='RdBu_r', plotFaults= False,alpha=1,contour=False)


ds = gdal.Open('geotifs/OahuVertical_new.tif')
vel = ds.GetVirtualMemArray()
makeMap.mapImg(vel, lons_a, lats_a, vmin, vmax, pad,zoom, title, bg=bg, cm='RdBu_r', plotFaults= False,alpha=1,contour=False)


# plt.figure()
# plt.imshow(velocity_d_resamp,vmin=vmin,vmax=vmax)

# # Test a point for the vert-horz
# py,px = 968,3202
# va = velocity_a[py,px]
# vd = velocity_d_resamp[py,px]
# az_a_p = az_a[py,px]
# az_d_p = az_d_resamp[py,px]
# inc_a_p = inc_a[py,px]
# inc_d_p = inc_d_resamp[py,px]
# losA = util.getLOSvec(az_a_p,inc_a_p)
# losD = util.getLOSvec(az_d_p,inc_d_p)
# o = np.concatenate((np.array([va,vd]).T,np.array([0,0])),axis=0)
# mx = np.array([1,0,0])
# # my = np.array([0,1,0])
# mz = np.array([0,0,1])
# smooth=.4
# A = np.array([[np.dot(losA,mx), np.dot(losA,mz)],
#               [np.dot(losD,mx), np.dot(losD,mz)],
#               [smooth,          0              ],
#               [0,               smooth         ]])
# Aa = np.dot( np.linalg.inv(np.dot(A.T,A)), A.T)
# vertHor = np.dot(Aa,o)