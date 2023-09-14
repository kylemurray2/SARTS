#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert land cover from NLCD to rdr coordinates.

Created on Fri Sep 23 16:41:59 2022

@author: km
"""

import numpy as np
import os
from osgeo import gdal
from matplotlib import pyplot as plt
from isce.components import isceobj
# from mintpy.utils import writefile
import argparse


def cmdLineParser():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser(
        description='Converts nlcd landcover img file to geotiff and regrids to the radar coordinates of your stack.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input_file', type=str, dest='fileName',
        required=True, help='Input nlcd img path/file name')
    return parser.parse_args()


def convert_land_cover(fileName, plot_flag=False):
    if not os.path.isdir('./Fringe'):
        os.mkdir('./Fringe')
    
    latFile = './merged/geom_reference/lat.rdr.full'
    lonFile = './merged/geom_reference/lon.rdr.full' 
    rdr_outfile = './merged/geom_reference/landCover.rdr.full'
    
    nlcd_tif = fileName + '.tif'
    os.system('gdalwarp ' +  fileName + ' ' + nlcd_tif + ' -t_srs "+proj=longlat +ellps=WGS84"')
    
    decimate = 1
    ds = gdal.Open(nlcd_tif, gdal.GA_ReadOnly) 
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    LCimage = arr[0:arr.shape[0]:decimate, 0:arr.shape[1]:decimate].copy(order='C')
    
    LCimage[LCimage == 0] = 21
    
    transform = ds.GetGeoTransform()
    startLat = transform[3]
    deltaLat = transform[5]
    startLon = transform[0]
    deltaLon = transform[1]
    
    latim = isceobj.createImage()
    latim.load(latFile + '.xml')
    lat = latim.memMap()[:,:,0]
    lonim = isceobj.createImage()
    lonim.load(lonFile + '.xml')
    lon = lonim.memMap()[:,:,0]
    
    lati = np.clip(((lat - startLat) / deltaLat).astype(int), 0, LCimage.shape[0] - 1)
    loni = np.clip(((lon - startLon) / deltaLon).astype(int), 0, LCimage.shape[1] - 1)
    cropped = (LCimage[lati, loni] + 1)
    cropped = np.reshape(cropped, (latim.coord2.coordSize, latim.coord1.coordSize))
    cropped.tofile(rdr_outfile)
    
    croppedim = isceobj.createImage()
    croppedim.initImage(rdr_outfile, 'read', cropped.shape[1], 'BYTE')
    croppedim.renderHdr()
    

    lc_fn = './Fringe/LCimage.vrt'
    os.system('rm ' + lc_fn)
    atr = {
        'WIDTH': LCimage.shape[1],
        'LENGTH': LCimage.shape[0],
        'FILE_TYPE': 'landCover'
    }
    # writefile.write(LCimage, out_file=lc_fn, metadata=atr)
    
    waterMask = np.zeros(cropped.shape, dtype=np.uint8)
    waterMask[cropped != 22] = 1
    
    if plot_flag:
        plt.figure()
        plt.imshow(cropped, cmap='magma')
        plt.title('Land Cover')
        
        plt.figure()
        plt.imshow(waterMask, cmap='magma')
        plt.title('Water Mask')
    else:
        print('Plotting turned off.')
        
        
    outName = './merged/geom_reference/waterMask.rdr.full'
    im2 = croppedim.clone()
    im2.filename = outName
    im2.dump(outName + '.xml')
    im2.renderHdr()
    im2.renderVRT()
    waterMask.tofile(outName)
    im2.finalizeImage()

    import imageio
    imageio.imwrite('watermask.tif',waterMask)

if __name__ == '__main__':
    inps = cmdLineParser()
    convert_land_cover(inps.fileName)
