#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 15:10:53 2023

Download SAR SLCs, Orbit files, and DEM

@author: Kyle Murray
"""

import os, argparse, glob, zipfile, re, requests,sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from stackSentinel import sentinelSLC
from SARTS import asfQuery, getDEM, setupStack, config
import concurrent.futures
from datetime import datetime,timedelta


nproc = int(os.cpu_count())

print('downloading with '  + str(nproc) + ' cpus')

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
    response = requests.get(url,stream=True,allow_redirects=True)
    
    # Open the local file for writing
    with open(outname, 'wb') as file:
        # Iterate through the content and write to the file
        for data in response.iter_content(chunk_size=int(2**14)):
            file.write(data)


def dlOrbs(gran,outdir):
    # Create an empty list to store the returned URLs
    sat_dates=[]
    for g in gran:
        string = g.split('_')[0] + g.split('_')[5][0:8]
        sat_dates.append(string)
    
    sat_dates = np.unique(sat_dates)
    sat_dates.sort()
    
    orbUrls = asfQuery.get_orbit_url(sat_dates)
 
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    if not os.path.isdir('./orbits'):
        os.system('ln -s ' + outdir + ' ./')
        
    orbUrls = np.unique(orbUrls)
    orbUrls.sort()
    
    print('Searching for orbits...')

    outNames = []
    dlorbs = []
    for url in orbUrls:
        fname = os.path.join(outdir,url.split('/')[-1])
        if not os.path.isfile(fname):
            print('No file ' + fname)
            outNames.append(fname)
            dlorbs.append(url)
        else:
            if os.path.getsize(fname) < 1024: 
                print('Overwriting ' + fname + ' because it was too small...')
                outNames.append(fname)
                dlorbs.append(url)
    print('Dowloading orbit files...')

    print(dlorbs)

    # Download urls in parallel and in chunks
    with concurrent.futures.ThreadPoolExecutor(max_workers=nproc) as executor:  # Adjust max_workers as needed
        futures = [executor.submit(dl, url, outName) for url, outName in zip(dlorbs, outNames)]
        concurrent.futures.wait(futures)

    outNames = []
    dlorbs = []
    redflag = False

    for url in orbUrls:
        fname = os.path.join(outdir,url.split('/')[-1])
        if not os.path.isfile(fname):
            print('No file ' + fname)
            outNames.append(fname)
            dlorbs.append(url)
            redflag=True
        else:
            if os.path.getsize(fname) < 1024: 
                print('Deleting ' + fname + ' because it was too small...')
                outNames.append(fname)
                os.remove(fname)
                dlorbs.append(url)
                redflag=True
            # else:
                # print('Downloaded OK ' + fname)

    
    if redflag:
        print('Some orbit files may have not been properly downloaded. Please try again.')
        #sys.exit(1)

def dlSlc(slcUrls, gran,outdir,sat='S1'):

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

    # check the files
    for ii in range(len(gran)):
        fname = os.path.join(outdir, gran[ii] + '.zip')
        if not os.path.isfile(fname):
            print('Warning: File does not exist ' + fname)
            sys.exit(1)
        else:
            if sat!='ALOS':
                if os.path.getsize(fname) < 2**30: # If it's smaller than 1 Gb
                    print('Warning: ' + fname + ' is too small. Try again.')
                    sys.exit(1)

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
        os.makedirs(aux_cal_out_dir)
        
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
    if ps.sat=='SENTINEL-1':
        for z in zips:
            safe = sentinelSLC(z)
            safe.get_lat_lon_v2()
            minlats.append(safe.SNWE[0])
            maxlats.append(safe.SNWE[1])
            minlons.append(safe.SNWE[2])
            maxlons.append(safe.SNWE[3])

        minlat = min(minlats)
        maxlat = max(maxlats)
        minlon = min(minlons)
        maxlon = max(maxlons)
    else:
        minlat = np.floor( float(ps.bounds.split(',')[0]) )
        maxlat = np.ceil(float(ps.bounds.split(',')[1]))
        minlon = np.floor(float(ps.bounds.split(',')[2]))
        maxlon = np.ceil(float(ps.bounds.split(',')[3]))

    if not os.path.isdir('Figs'):
        os.mkdir('Figs')
    if not os.path.isdir('Npy'):
        os.mkdir('Npy')

    plt.savefig('Figs/bounds.png')



    demBounds = f"{int(np.floor(minlat))},{int(np.ceil(maxlat))},{int(np.floor(minlon))},{int(np.ceil(maxlon))}"

    # Download dem if it doesn't exist
    file_list = glob.glob('./DEM/*wgs84*')
    if not file_list:
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


def update_yaml_key(file_path, key, new_value):
    new_value = new_value + ' '
    with open(file_path, "r") as f:
        lines = f.readlines()
    
    with open(file_path, "w") as f:
        for line in lines:
            # The regular expression is updated to accommodate varying whitespaces and possible comments
            match = re.match(rf"^(\s*{key}\s*:\s*)([^#]*)(.*)$", line)
            if match:
                # Replace the value while preserving leading whitespaces, the key, and comments
                line = f"{match.group(1)}{new_value}{match.group(3)}\n"
            f.write(line)
        


def main(inps):
    
    ps = config.getPS()

    if inps.searchData_flag==True:
        slcUrls, gran, _,_ = asfQuery.getGran(ps.path, ps.start, ps.end, ps.sat, ps.bounds, ps.poly)
    else:
        print('Using an existing out.csv file.')
        df = pd.read_csv('out.csv')
        slcUrls, gran = df['URL'], df['Granule Name']

    dates = [l[17:25] for l in gran]
    dates = np.unique(dates)
    dates.sort()
    print(dates)

    if inps.dlOrbs_flag:
        
        dlOrbs(gran,ps.orbit_dirname)
        
        # Check if aux_cal files exist:
        result = check_aux_cal(ps.aux_dirname)
        if not result:
            dlAuxCal(ps.aux_dirname)

    if ps.reference_date:
        print(f'Reference date is ({ps.reference_date})')
    else:
        print(f'setting reference date to the first date ({dates[0]})')
        update_yaml_key('params.yaml', 'reference_date', str(dates[0]))
        ps.reference_date = dates[0]


        
        
        
    if inps.dlSlc_flag:
        # Check for current SLCs and remove any bad ones
        zips = glob.glob(os.path.join(ps.slc_dirname,'*.zip'))
        if len(zips)>0:
            if ps.sat=='SENTINEL-1':
                flag = setupStack.checkSLC(ps)

        dlSlc(slcUrls, gran, ps.slc_dirname,ps.sat)
        
        # make sure we have the reference date
        matching_files = [filename for filename in zips if ps.reference_date in filename]
        if len(matching_files) ==0:
            print('attempting to download slc for reference date')
            date_obj1 = datetime.strptime(ps.reference_date, '%Y%m%d')
            date_obj2 =  date_obj1 + timedelta(days=1)
            date_obj1 -= timedelta(days=1)
            formatted_date_str1 = date_obj1.strftime('%Y-%m-%dT%H:%M:%SZ')
            formatted_date_str2 = date_obj2.strftime('%Y-%m-%dT%H:%M:%SZ')
            try:
                slcUrls, gran, _,_ = asfQuery.getGran(ps.path, formatted_date_str1, formatted_date_str2, ps.sat, ps.bounds, ps.poly)
                dlSlc(slcUrls, gran, ps.slc_dirname,ps.sat)
            except:
                print('failed to find SLC for reference date')
                sys.exit(1)


    demBounds, DEM = dlDEM(ps)
    ps.dem_bounds = demBounds  # Get the DEM and define dem location
    ps.dem = DEM
    np.save('ps.npy', ps)

    print('Next step is setupStack.py')
    
if __name__ == '__main__':
    '''
    Main driver.
    '''
    inps = cmdLineParser()

    # For debugging
    # inps = argparse.Namespace()
    # inps.searchData_flag = True
    # inps.dlSlc_flag = True
    # inps.dlOrbs_flag = True
    # inps.get_srtm = False

    main(inps)
