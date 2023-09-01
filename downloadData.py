#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 15:10:53 2023

Download SAR SLCs, Orbit files, and DEM

@author: km
"""

import os
import argparse
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from stackSentinel import sentinelSLC
from SARTS import asfQuery, getDEM
import localParams

ps = localParams.getLocalParams()

def cmdLineParser():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser(
        description='Download SLCs, orbits, and DEM for stack processing',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--search-data', type=bool, dest='searchData_flag', default=True)
    parser.add_argument('-d', '--download-slc', type=bool, dest='dlSlc_flag', default=True)
    parser.add_argument('-o', '--download-orbits', type=bool, dest='dlOrbs_flag', default=True)

    return parser.parse_args()


def searchData():
    slcUrls, gran, _,response = asfQuery.getGran(ps.path, ps.frame, ps.start, ps.end, ps.sat, ps.bounds, ps.point, ps.poly)
    return slcUrls, gran


def dlOrbs(gran):
    # Make directories and download the slcs and orbits and move them to directories
    orbUrls = [asfQuery.get_orbit_url(g) for g in gran]

    if not os.path.isdir(ps.slc_dirname):
        os.mkdir(ps.slc_dirname)
    if not os.path.isdir('orbits'):
        os.mkdir('orbits')

    for url in orbUrls:
        orbit_filename = url[39:]
        if not os.path.isfile(os.path.join('orbits', orbit_filename)):
            print('Downloading orbit ' + orbit_filename)
            os.system(f'wget -P ./orbits -nc -c {url} >> log')  # Downloads the orbit files



def dlSlc(slcUrls, gran):
    for ii, url in enumerate(slcUrls):
        fname = os.path.join(ps.slc_dirname, gran[ii] + '.zip')
        if not os.path.isfile(fname):
            print('Downloading orbit ' + url)
            #os.system(f'wget -P {ps.slc_dirname} -nc -c {url} >> log')  # Downloads the zip SLC files with wget
            os.system(f'curl -o {fname} -C - {url}')  # Downloads the zip SLC files




def dlSlc(url, gr):
        fname = os.path.join('test', gr + '.zip')
        dir = 'test'
        if not os.path.isfile(fname):
            print('Downloading orbit ' + url)
            os.system(f'wget -P {dir} -nc -c {url} >> log')  # Downloads the zip SLC files with wget
            # os.system(f'curl -o {fname} -C - {url}')  # Downloads the zip SLC files

with concurrent.futures.ThreadPoolExecutor() as exector : 
   exector.map(download_zip, urls)

print('Using an existing out.csv file.')
df = pd.read_csv('out.csv')
slcUrls, gran = df['URL'], df['Granule Name']



import requests
from concurrent.futures import ThreadPoolExecutor



response = requests.get(url,stream=True,allow_redirects=True)
# Get the content of the response (the file data)
file_data = response.content


def dl(url,outname):
    import requests
    response = requests.get(url,stream=True,allow_redirects=True)
    
    # Open the local file for writing
    with open(outname, 'wb') as file:
        # Iterate through the content and write to the file
        for data in response.iter_content(chunk_size=int(2**14)):
            file.write(data)
      
url = slcUrls[0]
gr = gran[0]

outNames = []
for ii in range(len(gran)):
    outNames.append(os.path.join('test', gran[ii] + '.zip'))

# dl(url,outname)   
import concurrent.futures
     
with concurrent.futures.ThreadPoolExecutor(max_workers=os.cpu_count) as executor:  # Adjust max_workers as needed
    futures = [executor.submit(dl, url, outName) for url, outName in zip(slcUrls, outNames)]
    concurrent.futures.wait(futures)


# Open the local file for writing
with open('./test/file.zip', 'wb') as file:
    file.write(response.content)


# dlSlc(url, gr)


with ThreadPoolExecutor(max_workers=15) as exector : 
   exector.map(dlSlc, slcUrls,gran)

   
asf_username='kylemurray2'
asf_password='Kujira8900'

















def dlDEM():
    zips = glob.glob(os.path.join(ps.slc_dirname, '*zip'))
    # Figure out what bounds to use for the DEM
    minlats, maxlats, minlons, maxlons = [], [], [], []

    for z in zips:
        safe = sentinelSLC(z)
        safe.get_lat_lon_v2()
        minlats.append(safe.SNWE[0])
        maxlats.append(safe.SNWE[1])
        minlons.append(safe.SNWE[2])
        maxlons.append(safe.SNWE[3])

    if not os.path.isdir('Figs'):
        os.mkdir('Figs')
    if not os.path.isdir('Npy'):
        os.mkdir('Npy')

    plt.savefig('Figs/bounds.png')

    minlat = min(minlats)
    maxlat = max(maxlats)
    minlon = min(minlons)
    maxlon = max(maxlons)

    demBounds = f"{int(np.floor(minlat))},{int(np.ceil(maxlat))},{int(np.floor(minlon))},{int(np.ceil(maxlon))}"

    # Download dem if it doesn't exist
    if not os.path.isdir('./DEM'):
        getDEM.getDEM(demBounds)
        DEM = glob.glob(os.path.join(ps.workdir, 'DEM', '*wgs84.dem'))[0]
        # Updating DEM’s wgs84 xml to include the full path to the DEM
    else:
        DEM = glob.glob(os.path.join(ps.workdir, 'DEM', '*wgs84.dem'))[0]
    
    os.system(f'fixImageXml.py -f -i {DEM} >> log')
    return demBounds, DEM


def main(inps):
    if inps.searchData_flag:
        slcUrls, gran  = searchData()
    else:
        print('Using an existing out.csv file.')
        df = pd.read_csv('out.csv')
        slcUrls, gran = df['URL'], df['Granule Name']

    dates = [l[17:25] for l in gran]
    dates.sort()
    print(dates)

    if inps.dlOrbs_flag:
        dlOrbs(gran)

    if inps.dlSlc_flag:
        dlSlc(slcUrls, gran)

    demBounds, DEM = dlDEM()
    ps.dem_bounds = demBounds  # Get the DEM and define dem location
    ps.dem = DEM
    np.save('ps.npy', ps)


if __name__ == '__main__':
    '''
    Main driver.
    '''
    inps = cmdLineParser()
    main(inps)
