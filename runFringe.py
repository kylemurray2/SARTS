#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 15:06:09 2023

Reads input from file ./Fringe/params.cfg
for example:
    # --bbox/-b  S N W E
    bounds = 2600 6300 17400 38100
    maxMem = 100000
    nx = 20700
    ny = 3700
    filterFlag      = True
    unwrap          = True 
    filterStrength  = '.4'
    fixImage        = False  
    azlooks = 1
    rglooks = 4

    waterMask = True
    # NLCD raster for landcover and watermask
    nlcd_in = /d/HI/Landcover/data/hi_oahu_2011_ccap_hr_land_cover20140619.img 
    nlcd_tif = /d/HI/Landcover/data/oahu.tif


    Run this from Fringe directory
    Run this before starting spyder: export LD_PRELOAD=/home/km/Software/miniconda3/envs/isce/lib/libmkl_core.so.2:/home/km/Software/miniconda3/envs/isce/lib/libmkl_sequential.so.2:/home/km/Software/miniconda3/envs/isce/lib/libmkl_avx512.so.2:/home/km/Software/miniconda3/envs/isce/lib/libmkl_def.so.2


@author: km
"""

import numpy as np
import os
import glob
from matplotlib import pyplot as plt
from osgeo import gdal
from PyPS2 import util
from mintpy.utils import readfile, writefile
import argparse
import sequential_PL 
import tops2vrt
import nmap
import adjustMiniStacks
import ampdispersion

ps = np.load('./ps.npy',allow_pickle=True).all()
inps = argparse.Namespace()

# For tops2vrt
inps.indir      = ps.mergeddir
inps.stackdir   = './Fringe/coreg_stack'
inps.geomdir    = './Fringe/geometry'
inps.outdir     = './Fringe/slcs'
inps.bbox       = np.array([ps.cropymin,ps.cropymax,ps.cropxmin,ps.cropxmax])
inps.geobbox    = None

# For nmap
inps.inputDS        = './Fringe/coreg_stack/slcs_base.vrt'
inps.outputDS       = './Fringe/KS2/nmap'
inps.countDS        = './Fringe/KS2/count'

if os.path.isfile('./merged/geom_reference/waterMask.rdr.full.crop.vrt'):
    inps.maskDS     ='./merged/geom_reference/waterMask.rdr.full.crop.vrt'
else:
    inps.maskDS     = None

inps.linesPerBlock  = 256
inps.memorySize  = ps.maxMem
inps.halfWindowX = 11
inps.halfWindowY = 5
inps.pValue      = 0.05
inps.method      = 'KS2'
inps.noGPU       = 'True'

# For sequential_PL
inps.inputDir       = ps.slcdir
inps.weightDS       = inps.outputDS
inps.outputDir      = './Fringe/Sequential'
inps.memorySize     = inps.memorySize
inps.minNeighbors   = 5
inps.miniStackSize  = 30
inps.forceprocessing= False

# For adjustMiniStacks
inps.slcDir = './Fringe/slcs/'
inps.miniStackDir = './Fringe/Sequential/miniStacks/'
inps.datumDir = './Fringe/Sequential/Datum_connection/'
inps.outDir = './Fringe/adjusted_wrapped_DS'
inps.unwrapped = False

# For ampdispersion
inps.meanampDS = './Fringe/ampDispersion/mean'
inps.refBand = 1
inps.outputAD = './Fringe/ampDispersion/ampdispersion' # outputDS is different for nmap, so define it here.

inps.nx = str(ps.nx)
inps.ny = str(ps.ny)

#______________________________________________________________________________
tops2vrt.main(inps)

nmap.main(inps)

sequential_PL.main(inps)

adjustMiniStacks.main(inps)

ampdispersion.main(inps)
#______________________________________________________________________________

# if ps.sensor=='ALOS':
#     os.system('phase_link.py -i ./Fringe/coreg_stack/slcs_base.vrt -w ./Fringe/KS2/nmap -o ./Fringe/PhaseLink -r 70000 -x 11 -y 5')
#     ds_SLCS = glob.glob('./Fringe/PhaseLink/*slc')
#     for fn_slc in ds_SLCS:
#         util.write_xml(fn_slc,nx,ny,1,dataType='CFLOAT',scheme='BIP')


# ds = gdal.Open('./Fringe/KS2/nmap')
# nmap = ds.GetVirtualMemArray()
# nmap = np.asarray(nmap,dtype=np.float32)
# plt.figure();plt.imshow(nmap[7,:,:])

#Output the ps pixels by using a threshold on ampdispersion
ampDispersionThreshold = '0.45'
os.system('imageMath.py -e="a<' + ampDispersionThreshold + '" --a=./Fringe/ampDispersion/ampdispersion  -o ./Fringe/ampDispersion/ps_pixels -t byte')

#make tcorrMean.bin_____________________________________________________________
miniStacks_tcorr_files = glob.glob('./Fringe/Sequential/miniStacks/*/EVD/tcorr.bin')
chunk_size = 400  # Number of time slices to load at once
num_chunks = ps.ny // chunk_size
tcorrsMean = np.zeros((ps.ny, ps.nx),dtype='float32')

# ii=0;start = ii * chunk_size ;end = start+chunk_size + 1

for ii in range(ps.ny):
    tc = np.zeros((len(miniStacks_tcorr_files),ps.nx))
    for jj,mtf in enumerate(miniStacks_tcorr_files):
        ds = gdal.Open(mtf)
        tc[jj,:] = ds.GetVirtualMemArray()[ii,:]
    tcorrsMean[ii,:] = np.nanmean(tc,axis=0)
    
plt.figure();plt.imshow(tcorrsMean,vmin=0,vmax=1,cmap='magma');plt.show()

# if ps.sensor=='ALOS':
#     ds = gdal.Open('Fringe/PhaseLink/tcorr.bin')
#     tcorrsMean = ds.GetVirtualMemArray()
#     plt.figure();plt.imshow(tcorrsMean,vmin=0,vmax=1,cmap='magma');plt.show()

# Write the average tcorr file 
fmt = "GTiff"
driver = gdal.GetDriverByName(fmt)
[cols, rows] = tcorrsMean.shape
outDataRaster = driver.Create("./Fringe/tcorrMean.bin", rows, cols, 1, gdal.GDT_Float32)
bnd = outDataRaster.GetRasterBand(1)
bnd.WriteArray(tcorrsMean)
bnd.FlushCache()


# Now use the script ifgs.py to make unwrapped ifgs
