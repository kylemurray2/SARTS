#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 11:02:56 2022

Get the copernicus DEM 30m

@author: km
"""
import os, requests

def getDEM(bounds,demtype='COP30',apiKey='b4615a3e690e63f77c49d353b6797095'):
    '''
    bounds 'S,N,W,E'
    '''

    baseurl = "https://portal.opentopography.org/API/globaldem?"
    miny,maxy,minx,maxx = miny,maxy,minx,maxx = bounds.split(sep=',')
    
    data = dict(demtype='COP30',
      south=miny,
      north=maxy,
      west=minx,
      east=maxx,
      outputFormat='GTiff',
      API_Key=apiKey)
    
    if not os.path.isdir('./DEM'):
        os.mkdir('DEM')
    
    r = requests.get(baseurl, params=data, timeout=100, allow_redirects=True)
    open('DEM/dem.tif', 'wb').write(r.content)
    os.system('gdal_translate DEM/dem.tif -Of ISCE DEM/cop_dem_glo30_wgs84.dem')
    
# bounds = '35.63,35.81,-117.44,-117.23'#'SNWE
# getDEM(bounds)

