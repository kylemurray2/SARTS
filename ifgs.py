#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 17:05:26 2023

Takes output from Fringe
creates:
    Delaunay pairs
    fine.int
    fine_lk.int
    filt_lk.int
    filt_lk.cor
    filt_lk.unw

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
    
@author: km
"""

import numpy as np
import os
from osgeo import gdal
from isce.applications import looks
import FilterAndCoherence
import integratePS
import multiprocessing
from SARTS import unwrap
from isce.components import isceobj

# from Scripts.Fringe import writeStackVRT
# from matplotlib import pyplot as plt
# from datetime import date
# import scipy.spatial
# from isce.components import isceobj
# import glob
# from SARTS import util
from mintpy.utils import readfile, isce_utils
import localParams

ps = localParams.getLocalParams()

dlunw = True
makeIfgs = True
makeVrts = True            
num_processes = 5 # For parallel unwrapping (memory intensive)

nx              = ps.nx
ny              = ps.ny
bounds          = ps.bounds
filterFlag      = True
dounwrap          = ps.unwrap # Usually better to leave False and use runSnaphu.py for more options and outputs
filterStrength  = float(ps.filtStrength)
ps.azlooks      = int(ps.alks)
ps.rglooks      = int(ps.rlks)

ymin = int(np.floor(float(bounds.split(',')[0])))
ymax = int(np.ceil(float(bounds.split(',')[1])))
xmin = int(np.floor(float(bounds.split(',')[2])))
xmax = int(np.ceil(float(bounds.split(',')[3])))



ps.nyl     = ny//ps.azlooks
ps.nxl     = nx//ps.rglooks


# Make ifgs from normal SLCs
def makeIfg(slc1_fn,slc2_fn,ifg_fn):
    ds1 = gdal.Open(slc1_fn)
    ds2 = gdal.Open(slc2_fn)
    slc1 = ds1.GetVirtualMemArray()[ps.cropymin:ps.cropymax,ps.cropxmin:ps.cropxmax]
    slc2 = ds2.GetVirtualMemArray()[ps.cropymin:ps.cropymax,ps.cropxmin:ps.cropxmax]
    ifg = np.multiply(slc1,np.conj(slc2))
    
    out = isceobj.createIntImage() # Copy the interferogram image from before
    out.dataType = 'CFLOAT'
    out.filename = ifg_fn
    out.width = ifg.shape[1]
    out.length = ifg.shape[0]
    out.dump(out.filename + '.xml') # Write out xml
    fid=open(out.filename,"wb+")
    fid.write(ifg)
    out.renderHdr()
    out.renderVRT()  
    fid.close()
    
    return ifg

if not os.path.isdir(ps.intdir):
    print('Making merged/interferograms directory')
    os.mkdir(ps.intdir)

for ii in range(len(ps.dates)-1):
    d1 = ps.dates[ii]
    d2 = ps.dates[ii+1]
    slc1_fn = os.path.join(ps.slcdir,d1,d1+'.slc.full')
    slc2_fn = os.path.join(ps.slcdir,d2,d2+'.slc.full')
    pair = d1 + '_' + d2
    print('creating ' + pair + '.int')
    ifg_fn = os.path.join(ps.intdir,pair,pair+'.int')
    if not os.path.isfile(ifg_fn):
        if not os.path.isdir(os.path.join(ps.intdir,pair)):
            os.mkdir(os.path.join(ps.intdir,pair))
        makeIfg(slc1_fn,slc2_fn,ifg_fn)
    else:
        print(ifg_fn + ' already exists')

def downlook(pair):
    # Downlook ifgs
    pairDir    = os.path.join(ps.intdir,pair )
    ps.infile  = os.path.join(pairDir, pair + '.int')
    ps.outfile = pairDir + '/fine_lk.int'
    cor_file_out = pairDir + '/filt_lk.cor'
    filt_file_out =  pairDir + '/filt_lk.int'

    if not os.path.isfile(ps.outfile):
        print('downlooking ' + pair)
        looks.main(ps)
        
        # Filter and coherence
        if not os.path.isfile(filt_file_out):
            print('\n making ' + pair)
            FilterAndCoherence.runFilter(ps.outfile,filt_file_out,filterStrength)
            FilterAndCoherence.estCoherence(filt_file_out, cor_file_out)
    else:
        print(pair + '/' + filt_file_out + ' is already file.')
            

# for pair in pairsList:    
wavelength = 0.056
defo_max=0
max_comp=32
init_only=True
init_method='MCF'
cost_mode='SMOOTH'
atr = {}
atr['WIDTH']=ps.nxl
atr['LENGTH']=ps.nyl
atr['ALOOKS']=1
atr['RLOOKS']=1

def unwrapsnaphu(pair):  
    pairDir =  ps.intdir + '/' + pair 
    if not os.path.isfile(pairDir + '/filt_lk.unw'):
        print(pair + ' unwrapping')
        cor_file = pairDir + '/filt_lk.cor'
        int_file =  pairDir + '/filt_lk.int'
        unw_file =  pairDir + '/filt_lk.unw'        
        unwrap.unwrap_snaphu(int_file, cor_file, unw_file,atr, wavelength=wavelength, defo_max=defo_max, max_comp=max_comp, init_only=init_only, init_method=init_method, cost_mode=cost_mode)
    else:
        print(pair + ' already unwrapped.')

# for pair in ps.pairs:
#     pairdash = pair.replace('_','-')
#     if os.path.isdir(fringeDir + 'interferograms/' + pairdash):
#         if not os.path.isdir(fringeDir + 'PS_DS/sequential1/' + pair):
#             os.mkdir(fringeDir + 'PS_DS/sequential1/' + pair)
#         os.system('cp ./Fringe/interferograms/' + pairdash + '/* ./Fringe/PS_DS/sequential1/' + pair + '/')
#         os.system('fixImageXml.py -i ./Fringe/PS_DS/sequential1/' + pair + '/*.xml -f' )


# pair=ps.pairs[50]
# unwrapsnaphu(pair)


if dlunw:
    # for pair in ps.pairs:
    #     downlook(pair)
    # for pair in ps.pairs:
    #     unwrapsnaphu(pair)

    pool = multiprocessing.Pool(processes=num_processes)
    pool.map(downlook, ps.pairs)
    pool.close()
    pool.join()
    
    pool = multiprocessing.Pool(processes=num_processes)
    pool.map(unwrapsnaphu, ps.pairs)
    pool.close()
    pool.join()


# for pair in ps.pairs2:
#     pairDir = ps.intdir + '/' + pair
#     if not os.path.isdir(pairDir):
#         os.system('mkdir -p ' + pairDir)
#     if not os.path.isfile(pairDir + '/filt_lk.int'):
#         ifgOutName = pairDir + '/fine.int'
#         ifgLkName = pairDir + '/fine_lk.int'
#         print('making ' + pair)
#         d1 = pair.split('_')[0]
#         d2 = pair.split('_')[1]
#         fn_slc1 = ps.slcdir +'/' + d1 + '/' + d1 +  '.slc.full.vrt'
#         fn_slc2 = ps.slcdir +'/' + d2 + '/' + d2 +  '.slc.full.vrt'
#         ds1 = gdal.Open(fn_slc1)
#         ds2 = gdal.Open(fn_slc2)
#         slc1 = ds1.GetVirtualMemArray()
#         slc2 = ds2.GetVirtualMemArray()
#         ifg = np.multiply(slc1,np.conj(slc2))
#         out = isceobj.createImage() # Copy the interferogram image from before
#         out.dataType = 'CFLOAT'
#         out.filename = ifgOutName
#         out.width = ifg.shape[1]
#         out.length = ifg.shape[0]
#         out.dump(out.filename + '.xml') # Write out xml
#         fid=open(out.filename,"wb+")
#         fid.write(ifg)
#         out.renderHdr()
#         out.renderVRT()  
#         fid.close()
#         # Downlook
#         ps.infile = ifgOutName
#         ps.outfile = ifgLkName
#         looks.main(ps)
#         # Filter and coherence
#         cor_file = pairDir + '/filt_lk.cor'
#         int_file =  pairDir + '/filt_lk.int'
#         if not os.path.isfile(cor_file):
#             print('\n making ' + pair)
#             FilterAndCoherence.runFilter(ifgLkName,int_file,float(filterStrength))
#             FilterAndCoherence.estCoherence(int_file, cor_file)
            
# if makeVrts:            
#     # Load unw stack as vrt file
#     # You can write a stack vrt file with the script ~/Software/Scripts/Fringe/writeStackVRT.py
#     inps = argparse.Namespace()
#     inps.outdir = fringeDir + 'PS_DS/sequential1/unw'
#     inps.stackdir = fringeDir
#     inps.stacklist = glob.glob(fringeDir + 'PS_DS/sequential1/*/filt_lk.unw')
#     inps.indir = './merged/'
#     inps.outFn = 'unwStack.vrt'
#     # run the code
#     writeStackVRT.main(inps)
    
#     # Same for cor files
#     inps = argparse.Namespace()
#     inps.outdir = fringeDir + 'PS_DS/sequential1/cor'
#     inps.stackdir = fringeDir
#     inps.stacklist = glob.glob(fringeDir + 'PS_DS/sequential1/*/*.cor')
#     inps.indir = './merged/'
#     inps.outFn = 'corStack.vrt'
#     # Run the code
#     writeStackVRT.main(inps)
#     # Load cor stack
    
    
# # Load that in
# ds = gdal.Open(fringeDir + 'unwStack.vrt')
# unwStack = ds.GetVirtualMemArray()
# plt.figure();plt.imshow(unwStack[51,:,:])

# ds = gdal.Open(fringeDir + 'corStack.vrt')
# corStack = ds.GetVirtualMemArray()
# plt.figure();plt.imshow(corStack[10,:,:])

#idx = 51
# unw = np.asarray(unwStack[idx,:,:],dtype=np.float32).copy()
# cor = np.asarray(corStack[idx,:,:]).copy()


# if ps.waterMask:
#     wmds = gdal.Open('./merged/geom_reference/waterMask_lk.rdr.vrt')
#     wm = wmds.GetVirtualMemArray()


# ds = gdal.Open('./Fringe2/PS_DS/sequential1/' + ps.pairs[51] + '/filt_lk.unw.vrt')
# unw = ds.GetVirtualMemArray().copy()

# ds = gdal.Open('./Fringe2/PS_DS/sequential1/' + ps.pairs[51] + '/filt_lk.cor.vrt')
# cor = ds.GetVirtualMemArray().copy()
# unw[wm==0]= np.nan
# cor[wm==0]= np.nan

# unwm = unw.copy()
# unwm[cor<.7] = np.nan

# util.show(unwm);util.show(cor)
# from astropy.convolution import Gaussian2DKernel
# from scipy.signal import convolve as scipy_convolve
# from astropy.convolution import convolve
# kernel = Gaussian2DKernel(x_stddev=5)

# # Mask unw
# unwFilt = convolve(unwm, kernel)
# unwFilt[wm==0]= np.nan
# util.show(unwFilt)

# #Now add back in the good data
# unwFilt[cor>.7] = unw[cor>.7]