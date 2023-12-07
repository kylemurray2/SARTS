#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 12:58:37 2022

Load and plot MintPy outputs

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


ds = gdal.Open('./merged/geom_reference/waterMask.rdr.full.crop.vrt')
waterMask = ds.GetVirtualMemArray()

# maskConnComp
filename = mintDir + 'maskConnComp.h5'
ds = h5py.File(filename,'r+')   
maskConnComp = np.asarray(ds['mask'])
ds.close()

# # demErr.h5  (this comes in meters?)
# filename = mintDir + 'demErr.h5' 
# ds = h5py.File(filename,'r+')   
# demErr = np.asarray(ds['dem'])  #convert to cm
# ds.close()

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
avgSpatialCoh[waterMask<1] = np.nan


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

phs_model = util.phaseElev(velocity, hgt,msk, 1000, 4000, 500, 2000,makePlot=True)
velocity_c = velocity-phs_model
velocity_c[msk==0] = np.nan
velocity[msk==0] = np.nan

vmin,vmax = -20,20
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

title = 'DEM Error'
demErr[waterMask==0] = np.nan
makeMap.mapImg(demErr[::decimate,::decimate], lons[::decimate,::decimate], lats[::decimate,::decimate], vmin,vmax, pad,zoomLevel, title, bg, cm='magma',plotFaults= False,alpha=1)

#Spatial Filter
kernel = Gaussian2DKernel(x_stddev=1)
velFilt = convolve(velocity, kernel)
velFilt[msk==0] = np.nan
title='velocity'
makeMap.mapImg(velFilt[::decimate,::decimate], lons[::decimate,::decimate], lats[::decimate,::decimate],vmin,vmax, pad,zoomLevel, title, bg, cm='RdBu_r',plotFaults= False,alpha=1)

# Geocode and ouput geotiff
velGeo = util.geocodeKM(velocity,20,lons,lats,ps)
plt.figure();plt.imshow(velGeo,vmin=-30,vmax=0)
bounds = ps.bounds.split(',')
util.writeGeotiff(velGeo, (ps.minlat,ps.maxlat),(ps.minlon,ps.maxlon), 'HIlos.tif')

    
# timeseries.h5
filename = mintDir + 'timeseries.h5'
dsts = h5py.File(filename,'r+')  

vmin,vmax = -20,20
plt.figure();plt.imshow(velocity,vmin=vmin,vmax=vmax)

py,px =  775,2600
plt.figure();plt.plot(ps.dec_year,dsts['timeseries'][:,py,px],'.');plt.title('South Flank')

py,px =  1454,1000
plt.figure();plt.plot(ps.dec_year,dsts['timeseries'][:,py,px],'.');plt.title('Mauna Loa')

py,px =  1166,2243
plt.figure();plt.plot(ps.dec_year,dsts['timeseries'][:,py,px],'.');plt.title('Kilauea Summit')


if os.path.isfile('./lola.py'):
    import lola # file with lon/lat points list. lola = [[lon,lat],[lon,lat]]
    lonlat = lola.lola

ii=1
py1,px1 = util.ll2pixel(lons, lats, lonlat[ii][0], lonlat[ii][1])
velocity[py1[0],px1[0]]

# Subsections
nw = [2000,3400, 100, 2400]
sw = [890,2000, 710, 2500]
se = [90, 1780, 2500,5500]
ne = [1850, 4010, 2400,4160]

ratem = velocity.copy()
ratem[avgSpatialCoh<.7] = np.nan
ratem[ratem==0] = np.nan
# bg = 'World_Shaded_Relief'
bg = 'World_Imagery'

zoom=16
# Southwest map
offset = -1
vmin,vmax=-15,15
makeMap.mapImg(offset+ratem[sw[0]:sw[1], sw[2]:sw[3]], lons[sw[0]:sw[1], sw[2]:sw[3]], lats[sw[0]:sw[1], sw[2]:sw[3]], vmin, vmax, 0, zoom,'Rates (mm/yr)',bg=bg, cm='RdBu_r')
plt.savefig('southwestRates.svg',dpi=400)

plt.figure();plt.plot(dec_year,100*dsts['timeseries'][:,1244,1403],'.');plt.title('Wet-n-wild');plt.ylabel('LOS displacement (cm)')
plt.figure();plt.plot(dec_year,100*dsts['timeseries'][:,1399,1125],'.');plt.title('Landfill');plt.ylabel('LOS displacement (cm)')
plt.figure();plt.plot(dec_year,100*dsts['timeseries'][:,1750,976],'.');plt.title('Landfill');plt.ylabel('LOS displacement (cm)')
plt.figure();plt.plot(dec_year,100*dsts['timeseries'][:,1400,2240],'.');plt.title('Field');plt.ylabel('LOS displacement (cm)')

# Northwest map
offset = -2
vmin,vmax=-8,8
makeMap.mapImg(offset+ratem[nw[0]:nw[1], nw[2]:nw[3]], lons[nw[0]:nw[1], nw[2]:nw[3]], lats[nw[0]:nw[1], nw[2]:nw[3]], vmin, vmax, 0, zoom,'Rates (mm/yr)',bg=bg, cm='RdBu_r')
# plt.savefig('northwestRates.svg',dpi=400)

plt.figure();plt.plot(dec_year,100*dsts['timeseries'][:,2331,713],'.');plt.title('Landslide?');plt.ylabel('LOS displacement (cm)')

# Northeast map
offset = -3
vmin,vmax=-8,8
makeMap.mapImg(offset+ratem[ne[0]:ne[1], ne[2]:ne[3]], lons[ne[0]:ne[1], ne[2]:ne[3]], lats[ne[0]:ne[1], ne[2]:ne[3]], vmin, vmax, 0, zoom,'Rates (mm/yr)',bg=bg, cm='RdBu_r')
# plt.savefig('northeastRates.svg',dpi=400)

# Southeast map
offset = 0
zoom=15
vmin,vmax=-10,10
makeMap.mapImg(offset+velFilt[se[0]:se[1], se[2]:se[3]], lons[se[0]:se[1], se[2]:se[3]], lats[se[0]:se[1], se[2]:se[3]], vmin, vmax, 0, zoom,'Rates (mm/yr)',bg=bg, cm='RdBu_r')
# plt.savefig('southeastRates.svg',dpi=400)
# makeMap.mapBackground('World_Imagery', minlon, maxlon, minlat, maxlat, 14, 'Imagery')

# Sand Island map
offset = -2.5
vmin,vmax=-8,8
zoom=15
bg='World_Imagery'
makeMap.mapImg(offset+ratem[630:780, 3190:3380], lons[630:780, 3190:3380], lats[630:780, 3190:3380], vmin, vmax, 0, zoom,'Rates (mm/yr)',bg=bg, cm='RdBu_r')




fig,ax=plt.subplots();ax.set_title('avgPhaseVelocity')
ax.imshow(avgPhaseVelocity);plt.show()

avgPhaseVelocity_geo = geocode(avgPhaseVelocity,lons,lats,3,method='linear')
avgPhaseVelocity_geo[avgSpatialCoh_geo<.6] = np.nan
avgPhaseVelocity_geo[np.isnan(avgSpatialCoh_geo)] = np.nan

makeMap.mapBackground(bg, minlon, maxlon, minlat, maxlat, zoomLevel, title, pad=0,scalebar=False)
img_handle = plt.imshow(avgPhaseVelocity_geo,vmin=-4,vmax=4,cmap='RdBu_r', alpha=1, origin='upper',zorder=12, transform=ccrs.PlateCarree(), extent=[minlon, maxlon, minlat, maxlat])



# temporalCoherence
filename = mintDir + 'temporalCoherence.h5'
ds = h5py.File(filename,'r+')   
temporalCoherence = np.asarray(ds['temporalCoherence'])
ds.close()
plotImg(temporalCoherence,'temporalCoherence')

# numInvIfgram
filename = mintDir + 'numInvIfgram.h5'
ds = h5py.File(filename,'r+')   
numInvIfgram = np.asarray(ds['mask'])
ds.close()
plotImg(numInvIfgram,'numInvIfgram')

# numTriNonzeroIntAmbiguity
filename = mintDir + 'numTriNonzeroIntAmbiguity.h5'
ds = h5py.File(filename,'r+')   
numTriNonzeroIntAmbiguity = np.asarray(ds['mask'])
ds.close()
plotImg(numTriNonzeroIntAmbiguity,'numTriNonzeroIntAmbiguity')

fns_ifgs_wrapped = glob.glob(workDir + '/Fringe/PS_DS/2*.int')
ds = gdal.Open(fns_ifgs_wrapped[1])
ifg = ds.GetVirtualMemArray()
util.show(np.angle(ifg),'ifg','jet')

for ii in range(concom.shape[0]):
    concomzero.append(concom[ii,3067,3333])

unwStack = np.asarray(ds['unwrapPhase'])
unwBridgeStack = np.asarray(ds['unwrapPhase_bridging'])
ds.close()

fig,ax=plt.subplots(1,2);ax[0].set_title('unwrap Phase');ax[1].set_title('Bridged Phase')
ax[0].imshow(unwStack[0,:,:]);ax[1].imshow(unwBridgeStack[0,:,:]);plt.show()
diff = unwStack[0,:,:]-unwBridgeStack[0,:,:]
diff[avgSpatialCoh<.7] = np.nan
fig,ax=plt.subplots();ax.set_title('unwrapPhase-bridge')
ax.imshow(diff);plt.show()

# demErr.h5
filename = mintDir + 'demErr.h5'
ds = h5py.File(filename,'r+')   
demErr = np.asarray(ds['dem'])
ds.close()
plotImg(demErr,'demErr')

# timeseries_ERA5.h5
filename = mintDir + 'timeseries_ERA5.h5'
ds = h5py.File(filename,'r+')   
timeseries_ERA5 = np.asarray(ds['timeseries'])
ds.close()
plotImg(timeseries_ERA5[5,:,:],'timeseries_ERA5')

# ERA5.h5
filename = mintDir + 'inputs/ERA5.h5'
ds = h5py.File(filename,'r+')   
ERA5 = np.asarray(ds['timeseries'])
ds.close()
plotImg(ERA5[5,:,:],'ERA5')

# timeseriesResidual.h5
filename = mintDir + 'timeseriesResidual.h5'
ds = h5py.File(filename,'r+')   
timeseriesResidual = np.asarray(ds['timeseries'])
ds.close()
plotImg(timeseriesResidual[10,:,:],'timeseriesResidual')
