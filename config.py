#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 15:46:44 2023

@author: km
"""

import numpy as np
import os,re,yaml,argparse,sys
import importlib.util
from osgeo import gdal
from pathlib import Path

def load_yaml_to_namespace(yaml_file):
    # Load the YAML file into a dictionary
    with open(yaml_file, 'r') as yaml_in:
        yaml_dict = yaml.safe_load(yaml_in)

    # Create a namespace from the dictionary
    namespace = argparse.Namespace(**yaml_dict)
    
    return namespace


def getPS(directory='.'):
    
    # Load the params from the yaml file
    yaml_file = os.path.join(directory,'params.yaml')
    lp_file = os.path.join(directory,'localParams.py')
    
    if os.path.isfile(yaml_file):
        print('Parsing yaml file and updating ps namespace...')
        params = load_yaml_to_namespace(yaml_file)
        # Load the ps namespace
        if os.path.isfile(os.path.join(directory,'ps.npy')):
            ps = np.load(os.path.join(directory,'ps.npy'),allow_pickle=True).all()
        else:
            ps = params
        
        # Update ps with any changes to params
        for attr in dir(params):
            if not attr.startswith('_'):
                # print(attr)
                setattr(ps, attr, getattr(params, attr))
    
# Set up some additional variables 
        ps.workdir = os.getcwd()
        ps.sensor = ps.sat
        ps.startDate               = ps.start[0:10]
        ps.stopDate                = ps.end[0:10]
        miny,maxy,minx,maxx        = ps.bounds.split(sep=',')
        ps.bbox                    = miny +' '+ maxy  +' '+ minx +' '+  maxx #demBounds[0] + ' ' + demBounds[1] + ' ' + demBounds[2] + ' ' + demBounds[3] #SNWE
        

        if ps.sat=='ALOS':
            ps.slcDir = ps.slc_dirname
            ps.inputDir = ps.slc_dirname



        if ps.sat=='SENTINEL-1':
            ps.swath_num = str(ps.swath_num)

            if ps.numProcess=='auto':
                ps.numProcess = os.cpu_count()
                ps.numProcess4topo = int(ps.numProcess/3) 


        ps.mergeddir= ps.workdir + '/merged'
        ps.intdir   = ps.mergeddir + '/interferograms'
        ps.reference_date = str(ps.reference_date)
        ps.slcdir   = ps.mergeddir + '/SLC'

        ps.ps_output = Path(os.path.join(ps.dolphin_work_dir, 'PS', 'ps_pixels.tif'))
        ps.amp_mean_file = Path(os.path.join(ps.dolphin_work_dir, 'PS', 'amp_mean.tif'))
        ps.amp_dispersion_file = Path(os.path.join(ps.dolphin_work_dir, 'PS', 'amp_dispersion.tif'))
        ps.block_shape=(ps.block_shape_x,ps.block_shape_y)
        
        if not hasattr(ps,'reference_date'):
            ps.reference_date=None
        elif ps.reference_date=='None':
            ps.reference_date=None

        ps.tsdir    = ps.workdir + '/TS'
        if 'nx' in ps.__dict__.keys():
            ps.nxl = ps.nx//ps.rlks
            ps.nyl = ps.ny//ps.alks

    elif os.path.isfile(lp_file):
        print('Using localParams.py...  This will be depricated in future versions. Use a yaml file instead.')       
        # Load the module
        spec = importlib.util.spec_from_file_location('localParams', lp_file)
        localParams = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(localParams)
        
        ps = localParams.getLocalParams(directory)
    
    else:
        print('No params file was found.')
        sys.exit(1)
    
    # Can remove this later..
    ps.geom=None
    
    # Save the updated ps namespace 
    np.save(os.path.join(directory, 'ps.npy'),ps)
    
    return ps







# def radarGeometryTransformer(latfile, lonfile, epsg=4326):
#     '''
#     Create a coordinate transformer to convert map coordinates to radar image line/pixels.
#     '''

#     driver = gdal.GetDriverByName('VRT')
#     inds = gdal.OpenShared(latfile, gdal.GA_ReadOnly)
#     tempds = driver.Create('', inds.RasterXSize, inds.RasterYSize, 0)
#     inds = None

#     tempds.SetMetadata({'SRS' : 'EPSG:{0}'.format(epsg),
#                         'X_DATASET': lonfile,
#                         'X_BAND' : '1',
#                         'Y_DATASET': latfile,
#                         'Y_BAND' : '1',
#                         'PIXEL_OFFSET' : '0',
#                         'LINE_OFFSET' : '0',
#                         'PIXEL_STEP' : '1',
#                         'LINE_STEP' : '1'}, 'GEOLOCATION')

#     trans = gdal.Transformer( tempds, None, ['METHOD=GEOLOC_ARRAY'])

#     return trans

# def lonlat2pixeline(lonFile, latFile, lon, lat):

#     trans = radarGeometryTransformer(latFile, lonFile)

#     ###Checkour our location of interest
#     success, location = trans.TransformPoint(1, lon, lat, 0.)
#     if not success:
#         print('Location outside the geolocation array range')

#     return location


# def getLinePixelBbox(geobbox, latFile, lonFile):

#     south,north, west, east = geobbox

#     se = lonlat2pixeline(lonFile, latFile, east, south)
#     nw = lonlat2pixeline(lonFile, latFile, west, north)

#     ymin = np.int16(np.round(np.min([se[1], nw[1]])))
#     ymax = np.int16(np.round(np.max([se[1], nw[1]])))

#     xmin = np.int16(np.round(np.min([se[0], nw[0]])))
#     xmax = np.int16(np.round(np.max([se[0], nw[0]])))

#     print("x min-max: ", xmin, xmax)
#     print("y min-max: ", ymin, ymax)

#     return ymin, ymax, xmin, xmax


# def update_yaml_key(file_path, key, new_value):
#     with open(file_path, "r") as f:
#         lines = f.readlines()

#     with open(file_path, "w") as f:
#         for line in lines:
#             # Try to match a YAML key-value pair line
#             match = re.match(rf"(\s*)({key}:\s*)(.+)", line)
#             if match:
#                 # Replace the value while preserving leading whitespaces and the key
#                 line = f"{match.group(1)}{match.group(2)}{new_value}\n"
#             f.write(line)
            

#     if 'geobbox' in ps.__dict__.keys():
#         if ps.geobbox is not None:
#             #keep the radar coordinate crop values if they are non zero
#             if np.sum([ps.cropymin,ps.cropymax,ps.cropxmin,ps.cropxmax]) == 0:
#                 # get crop bounds from geobbox
#                 latFile = os.path.join(ps.mergeddir,'geom_reference','lat.rdr.full')
#                 lonFile = os.path.join(ps.mergeddir,'geom_reference','lon.rdr.full')
#                 # if the bounding box in geo-coordinate is given, this has priority
#                 print("finding bbox based on geo coordinates of {} ...".format(ps.geobbox))
#                 ps.cropymin,ps.cropymax,ps.cropxmin,ps.cropxmax = getLinePixelBbox(ps.geobbox, latFile, lonFile)    
    
#                 update_yaml_key("params.yaml", "cropymin", ps.cropymin)
#                 update_yaml_key("params.yaml", "cropymax", ps.cropymax)
#                 update_yaml_key("params.yaml", "cropxmin", ps.cropxmin)
#                 update_yaml_key("params.yaml", "cropxmax", ps.cropxmax)