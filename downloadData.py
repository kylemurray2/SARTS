#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 15:10:53 2023

Download SAR SLCs and Orbit files and DEM

@author: km
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import glob
import time
import requests
import stackSentinel
import zipfile
import asfQuery
import localParams
import getDEM
import pandas as pd
ps = localParams.getLocalParams()

dlSlc = True
dlOrbs = True
searchData = True
# def setupStack(makeStack=False,dlSlc=False):

if searchData:
    slcUrls,gran,dates = asfQuery.getGran(ps.path,ps.frame,ps.start,ps.end,ps.sat,ps.bounds,ps.point,ps.poly)

else:
    slcUrls = pd.read_csv('out.csv')['URL']
    gran = pd.read_csv('out.csv')['Granule Name']
    dates = pd.read_csv('out.csv')["Acquisition Date"]

dates.sort()
zips = glob.glob(ps.slc_dirname + '*zip')

# Make directories and download the slcs and orbits and move them to directories
if dlOrbs:
    orbUrls = []
    for g in gran:
        orbUrls.append(asfQuery.get_orbit_url(g))

    if not os.path.isdir(ps.slc_dirname):
        os.mkdir(ps.slc_dirname)
    if not os.path.isdir('orbits'):
        os.mkdir('orbits')

    for url in orbUrls:
        if not os.path.isfile('orbits/' + url[39:]):
            print('downloading orbit ' + url[39:])
            os.system('wget -P  ./orbits -nc -c ' + url + ' >> log') # Downloads the orbit files

if dlSlc:
    for ii, url in enumerate(slcUrls):
        if not os.path.isfile(ps.slc_dirname + gran[ii] + '.zip'):
            print('downloading orbit ' + url)
            os.system('wget -P ' + ps.slc_dirname + ' -nc -c ' + url + ' >> log') # Downloads the zip SLC files


def getRectBounds(minlat,maxlat,minlon,maxlon):
    lon_bounds = np.array([minlon,maxlon,maxlon,maxlon,maxlon,minlon,minlon,minlon])
    lat_bounds = np.array([maxlat,maxlat,maxlat,minlat,minlat,minlat,minlat,maxlat])
    return lon_bounds,lat_bounds

# Figure out what bounds to use for the DEM
minlats,maxlats,minlons,maxlons = [],[],[],[]
minlatssub,maxlatssub,minlonssub,maxlonssub = [],[],[],[]

for z in zips:
    safe = stackSentinel.sentinelSLC(z)
    safe.get_lat_lon_v2()
    # S,N,W,E = safe.SNWE[0],safe.SNWE[1],safe.SNWE[2],safe.SNWE[3]
    minlats.append(safe.SNWE[0])
    maxlats.append(safe.SNWE[1])
    minlons.append(safe.SNWE[2])
    maxlons.append(safe.SNWE[3])

    boundsAll = safe.getkmlQUAD(z)
    blons,blats = [],[]
    for b in boundsAll:
        blons.append(float(b.split(sep=",")[0]))
        blats.append(float(b.split(sep=",")[1]))

    blons,blats = np.asarray(blons), np.asarray(blats)

    blons.sort()
    blats.sort()

    minlatssub.append(blats[1])
    maxlatssub.append(blats[2])
    minlonssub.append(blons[1])
    maxlonssub.append(blons[2])

plt.figure()
for g in zips:
    safe = stackSentinel.sentinelSLC(g)
    safe.get_lat_lon_v2()
    # S,N,W,E = safe.SNWE[0],safe.SNWE[1],safe.SNWE[2],safe.SNWE[3]

    boundsAll = safe.getkmlQUAD(g)
    blons,blats = [],[]
    for b in boundsAll:
        blons.append(float(b.split(sep=",")[0]))
        blats.append(float(b.split(sep=",")[1]))

    blons,blats = np.asarray(blons), np.asarray(blats)

    blons.sort()
    blats.sort()

    lon_bounds,lat_bounds = getRectBounds(safe.SNWE[0],safe.SNWE[1],safe.SNWE[2],safe.SNWE[3])
    lon_bounds_sub,lat_bounds_sub = getRectBounds(blats[1],blats[2],blons[1],blons[2])

    plt.plot(lon_bounds,lat_bounds,linewidth=2,color='red',zorder=10)
    plt.plot(lon_bounds_sub,lat_bounds_sub,linewidth=2,color='blue',zorder=10)

if not os.path.isdir('Figs'):
    os.mkdir('Figs')
if not os.path.isdir('Npy'):
    os.mkdir('Npy')
plt.savefig('Figs/bounds.png')

minlat = min(minlats)
maxlat = max(maxlats)
minlon = min(minlons)
maxlon = max(maxlons)

#demBounds = [str(int(np.floor(minlat))), str(int(np.ceil(maxlat))), str(int(np.floor(minlon))), str(int(np.ceil(maxlon)))]
demBounds = str(int(np.floor(minlat))) +','+ str(int(np.ceil(maxlat)))+ ','+ str(int(np.floor(minlon))) +','+ str(int(np.ceil(maxlon)))

#os.system('mv run_files run_files_o')
   # Download dem if it doesn't exist
if not os.path.isdir('./DEM'):
    getDEM.getDEM(demBounds)
    DEM = glob.glob(ps.workdir + '/DEM/*wgs84.dem')[0]
    # Updating DEM’s wgs84 xml to include full path to the DEM
else:
    DEM = glob.glob(ps.workdir + '/DEM/*wgs84.dem')[0]

# os.system('fixImageXml.py -f -i ' + DEM + ' >> log')

minx,maxx,miny,maxy = ps.bounds.split(sep=',')

#setupParams = argparse.Namespace()
ps.dem_bounds  = demBounds # Get the DEM and define dem location
ps.dem  = DEM
np.save('ps.npy',ps)