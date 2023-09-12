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
from SARTS import asfQuery, getDEM, setupStack, config
import concurrent.futures
import requests
import zipfile,re

nproc = int(os.cpu_count())


def cmdLineParser():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser(
        description='Download SLCs, orbits, and DEM for stack processing',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--search-data', action='store_true', dest='searchData_flag', help='Search ASF for data and output to out.csv')
    parser.add_argument('-d', '--download-slc', action='store_true', dest='dlSlc_flag', help='download SLCs from ASF')
    parser.add_argument('-o', '--download-orbits', action='store_true', dest='dlOrbs_flag', help='download orbit files')
    parser.add_argument('-srtm', '--get-srtm', action='store_true', dest='get_srtm', help='Use SRTM dem instead of copernicus')

    return parser.parse_args()


def searchData(ps):
    slcUrls, gran, _,_ = asfQuery.getGran(ps.path, ps.start, ps.end, ps.sat, ps.bounds, ps.poly)
    return slcUrls, gran

# Example url
# url = 'https://datapool.asf.alaska.edu/SLC/SA/S1A_IW_SLC__1SDV_20221226T163241_20221226T163300_046505_059280_1BDE.zip'

def checkSizes(slcUrls,ps):
    bad= []
    for url in slcUrls:
        with requests.get(url, stream=True) as response:
            file_size_remote = int(response.headers['Content-Length'])
            
        fn = ps.slc_dirname + url.split('/')[-1]
        local_file_size = os.path.getsize(fn)
        if file_size_remote != local_file_size:
            print('bad file ' + url)
            bad.append(url)


def dl(url,outname):
    
    # Parse .netrc...This shouldn't be necessary, should change in future
    # homeDir = os.path.expanduser("~")
    # with open(os.path.join(homeDir,'.netrc'), 'r') as f:
    #     lines = f.readlines()
    # uname = lines[1].split()[1]
    # pword = lines[2].split()[1]
    #,auth=(uname, pword)
    #______________________________________
    
    response = requests.get(url,stream=True,allow_redirects=True)
    
    # Open the local file for writing
    with open(outname, 'wb') as file:
        # Iterate through the content and write to the file
        for data in response.iter_content(chunk_size=int(2**14)):
            file.write(data)
 
    
def dlOrbs(gran,outdir):

    # Make directories and download the slcs and orbits and move them to directories
    # orbUrls = [asfQuery.get_orbit_url(g) for g in gran]

    # Create an empty list to store the returned URLs
    orbUrls = []
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=nproc) as executor:  # Adjust max_workers as needed
        futures = [executor.submit(asfQuery.get_orbit_url, g) for g in gran]
    
        for future in concurrent.futures.as_completed(futures):
            try:
                url = future.result()
                orbUrls.append(url)
            except Exception as e:
                print(f"An exception occurred: {e}")
    
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    with open('orbits/orbUrls.txt', 'w') as file:
        for item in orbUrls:
            file.write(f"{item}\n")
    

    outNames = []
    dlorbs = []
    for url in orbUrls:
        fname = os.path.join(outdir,url.split('/')[-1])
        if not os.path.isfile(fname):
            outNames.append(fname)
            dlorbs.append(url)
        else:
            print('already exists ' + fname)

    # Download urls in parallel and in chunks
    with concurrent.futures.ThreadPoolExecutor(max_workers=nproc) as executor:  # Adjust max_workers as needed
        futures = [executor.submit(dl, url, outName) for url, outName in zip(dlorbs, outNames)]
        concurrent.futures.wait(futures)

    # Or do it in series with wget
    # for url in orbUrls:
    #     orbit_filename = url[39:]
    #     if not os.path.isfile(os.path.join('orbits', orbit_filename)):
    #         print('Downloading orbit ' + orbit_filename)
    #         os.system(f'wget -P ./orbits -nc -c {url} >> log')  # Downloads the orbit files


# The old way to download slcs (in series with wget or curl)
# def dlSlc(slcUrls, gran):
#     for ii, url in enumerate(slcUrls):
#         fname = os.path.join(ps.slc_dirname, gran[ii] + '.zip')
#         if not os.path.isfile(fname):
#             print('Downloading SLC ' + url)
#             #os.system(f'wget -P {ps.slc_dirname} -nc -c {url} >> log')  # Downloads the zip SLC files with wget
#             os.system(f'curl -o {fname} -C - {url}')  # Downloads the zip SLC files




def dlSlc(slcUrls, gran,outdir):

    # Create a list of file path/names
    outNames = []
    dlSLCs = []
    for ii in range(len(gran)):
        fname = os.path.join(outdir, gran[ii] + '.zip')
        if not os.path.isfile(fname):
            outNames.append(os.path.join(outdir, gran[ii] + '.zip'))
            dlSLCs.append(slcUrls[ii])
            
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
        
    # Download urls in parallel and in chunks
    
    print('Downloading the following files:')
    print(dlSLCs)
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=nproc) as executor:  # Adjust max_workers as needed
        futures = [executor.submit(dl, url, outName) for url, outName in zip(dlSLCs, outNames)]
        concurrent.futures.wait(futures)


def check_aux_cal(dir_path):
    # Check if directory exists
    if not os.path.exists(dir_path):
        print(f"The directory {dir_path} does not exist.")
        return False

    # Count the number of files that end with 'SAFE'
    safe_files_count = sum(1 for filename in os.listdir(dir_path) if filename.endswith('SAFE'))

    # Check if the directory has more than 80 such files
    if safe_files_count > 80:
        print(f"The directory {dir_path} exists and has more than 80 files ending in 'SAFE'.")
        return True
    else:
        print(f"The directory {dir_path} exists but does not have more than 80 files ending in 'SAFE'.")
        return False
    

def get_download_links(base_url, num_pages):
    all_links = []
    for page in range(1, num_pages + 1):
        url = f"{base_url}?page={page}"
        response = requests.get(url)
        if response.status_code == 200:
            # Using regular expression to find all download links
            download_links = re.findall(r'href="(/download/[^"]+)', response.text)
            # Construct the full URLs
            full_urls = [f"https://sar-mpc.eu{link}" for link in download_links]
            all_links.extend(full_urls)
        else:
            print(f"Failed to fetch page {page}")
    return all_links
    

def dlAuxCal(aux_cal_out_dir):
    '''
    Download aux_cal files from ESA website. 
    '''
    if not os.path.isdir(aux_cal_out_dir):
        print('did not find aux_cal directory. Creating new one, and downloading files.')
        os.mkdir(aux_cal_out_dir)
        
    # Base URL and number of pages to scrape
    base_url = "https://sar-mpc.eu/adf/aux_cal/"
    num_pages = 5
    
    download_links = get_download_links(base_url, num_pages)
    
    # aux_cal_out_dir = '/d/S1/aux'
    auxUrls = []
    outNames = []
    for url in download_links:
        response = requests.head(url)
        location = response.headers['Location']
        outname = location.split('/')[-1].split('?')[0]
        outFile = os.path.join(aux_cal_out_dir, outname)
        outNames.append(outFile)
        auxUrls.append(location)
    

    with concurrent.futures.ThreadPoolExecutor(max_workers=nproc) as executor:  # Adjust max_workers as needed
        futures = [executor.submit(dl, url, outName) for url, outName in zip(auxUrls, outNames)]
        concurrent.futures.wait(futures)
        
    # Loop through each file in the directory
    for filename in os.listdir(aux_cal_out_dir):
        if filename.endswith('.zip'):
            # Construct the full path of the zip file and the folder to unzip to
            zip_path = os.path.join(aux_cal_out_dir, filename)
            unzip_dir = os.path.join(aux_cal_out_dir, filename[:-4])  # Removes '.zip'
    
            # Create a new directory to unzip files into
            if not os.path.exists(unzip_dir):
                os.makedirs(unzip_dir)
    
            # Unzip the zip file
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(unzip_dir)
    
            # Remove the zip file
            os.remove(zip_path)


def dlDEM(ps):
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
        if inps.get_srtm:
            getDEM.getDEM(demBounds,srtm=True)
            DEM = glob.glob(os.path.join(ps.workdir, 'DEM', '*wgs84'))[0]
        else:
            getDEM.getDEM(demBounds)
            DEM = glob.glob(os.path.join(ps.workdir, 'DEM', '*wgs84.dem'))[0]
    else:
        if inps.get_srtm:
            DEM = glob.glob(os.path.join(ps.workdir, 'DEM', '*wgs84'))[0]
        else:
            DEM = glob.glob(os.path.join(ps.workdir, 'DEM', '*wgs84.dem'))[0]

    # Updating DEM’s wgs84 xml to include the full path to the DEM
    os.system(f'fixImageXml.py -f -i {DEM} >> log')

    if len(DEM) == 0:
        print('Error: DEM does not exists.')

    return demBounds, DEM


def main(inps):
    
    ps = config.getPS()

    if inps.searchData_flag==True:
        slcUrls, gran  = searchData(ps)
    else:
        print('Using an existing out.csv file.')
        df = pd.read_csv('out.csv')
        slcUrls, gran = df['URL'], df['Granule Name']

    dates = [l[17:25] for l in gran]
    dates.sort()
    print(dates)

    if inps.dlOrbs_flag:
        print('Downloading orbits')
        # dlOrbs(gran,ps.orbit_dirname)
    
        # Check if aux_cal files exist:
        result = check_aux_cal(ps.aux_dirname)
        if not result:
            dlAuxCal(ps.aux_dirname)

    if inps.dlSlc_flag:
        # Check for current SLCs and remove any bad ones
        zips = glob.glob(os.path.join(ps.slc_dirname,'*.zip'))
        if len(zips)>0:
            flag = setupStack.checkSLC(ps)

        dlSlc(slcUrls, gran, ps.slc_dirname)

    demBounds, DEM = dlDEM(ps)
    ps.dem_bounds = demBounds  # Get the DEM and define dem location
    ps.dem = DEM
    np.save('ps.npy', ps)


if __name__ == '__main__':
    '''
    Main driver.
    '''
    # For debugging
    inps = argparse.Namespace()
    inps.searchData_flag = True
    inps.dlSlc_flag = True
    inps.dlOrbs_flag = True
    inps.get_srtm = False

    inps = cmdLineParser()
    main(inps)
