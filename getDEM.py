#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 11:02:56 2022

Get the copernicus DEM 30m

Either include your api key in the function, or make a file in your home 
directory called '~/.otkey' containing just the key on the first line.

@author: km
"""
import os, requests

def getDEM(bounds,demtype='COP30',apiKey=None):
    '''
    bounds 'S,N,W,E'
    '''
    home_dir = os.environ['HOME']

    if apiKey == None:
        print("Didn't inlcude API key. Reading the file ~/.otkey for open topography API key...")
        with open(home_dir + '/.otkey') as f:
            apiKey = f.readline().rstrip()
        
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

# Example usage
# bounds = '35.63,35.81,-117.44,-117.23'#'SNWE
# getDEM(bounds)

