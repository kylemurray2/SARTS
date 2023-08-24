#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 16:32:18 2022

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
from SARTS import asfQuery,getDEM
import localParams
import pandas as pd
ps = localParams.getLocalParams()

makeStack = False
dlSlc = False
dlOrbs = True
searchData = True
setupStack = True
# def setupStack(makeStack=False,dlSlc=False):

# Delete stuff. Careful! This can only be done manually now
if False:
    os.system('rm -r out.csv baselines coarse* configs coreg* log* ESD geom* inter* master reference seconda* merged misreg slaves run_files stack isce.log SAFE_files.txt')

if searchData:
    slcUrls,gran,dates = asfQuery.getGran(ps.path,ps.frame,ps.start,ps.end,ps.sat,ps.bounds,ps.point,ps.poly)

else:
    slcUrls = pd.read_csv('out.csv')['URL']
    gran = pd.read_csv('out.csv')['Granule Name']

dates=[]
for l in gran:
    dates.append(l[17:25])

dates.sort()


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

#check to make sure all the files are big enough and zipfile is valid
zips = glob.glob(ps.slc_dirname + '*zip')
zips.sort()
for z in zips:
    zipSize = os.stat(z).st_size
    if zipSize < 1e9:
        print('removing ' + z + ' because it is too small.')
        # os.remove(z)

for z in zips:
    try:
        x = zipfile.ZipFile(z)
        print("%s opened ok" % z)
        x.close()
    except:
        print("File %s is corrupt. Deleting." % z)
       # os.remove(z)
        continue

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

os.system('fixImageXml.py -f -i ' + DEM + ' >> log')

minx,maxx,miny,maxy = ps.bounds.split(sep=',')

#setupParams = argparse.Namespace()
ps.dem_bounds  = demBounds # Get the DEM and define dem location
ps.dem  = DEM
np.save('ps.npy',ps)

if setupStack:
    stackSentinel.main(ps)


runScripts = glob.glob(ps.workdir + '/run_files/run*')
runScripts.sort()

if makeStack:
    if not os.path.isdir('./logFiles'):
        os.mkdir('./logFiles')
    startT = time.time()
    for ii in np.arange(0,len(runScripts)):
        print('running ' + runScripts[ii].split('/')[-1])
        os.system('bash ' + runScripts[ii] + ' > logFiles/runlog_' + str(ii+1) )
#        if ii==5:
#            os.system('rm coreg_secondarys/*/overlap/IW?/*off*')
#        if ii==6:
#            os.system('rm -r ESD coarse_interferograms')
#        if ii==9:
#            os.system('rm coreg_secondarys/*/IW*/*off*')
#            os.system('mv SLCS SLCS_o')
#            os.system('mkdir SLCS')
#            os.system('mv SLCS_o/*' + ps.reference_date + '* SLCS/')
#            os.system('rm SLCS_o')

    endT = time.time()
    elapsed = endT-startT
    print('This run took ' + str(np.round(elapsed/60)) + ' minutes.')
