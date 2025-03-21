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
import shapely


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


def dlOrbs(gran, outdir):
    """
    Download orbit files for the given granules.
    
    Parameters
    ----------
    gran : list
        List of granule names
    outdir : str
        Output directory for orbit files
    """
    # Create the output directory if it doesn't exist
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    if not os.path.isdir('./orbits'):
        os.system('ln -s ' + outdir + ' ./')
    
    # Extract unique dates from the granule names
    unique_dates = set()
    for g in gran:
        # Check if this is a burst format granule (S1_322422_IW3_20200229T010154_VH_C8E8-BURST)
        if "-BURST" in g:
            try:
                parts = g.split('_')
                # Find the part that looks like a date (starts with year)
                date_part = next((part for part in parts if len(part) >= 8 and part[:4].isdigit()), None)
                
                if date_part and 'T' in date_part:
                    # Format is like 20200229T010154
                    date_str = date_part[:8]  # Extract YYYYMMDD part
                    
                    # Determine the satellite (S1A or S1B) - for bursts we may need to check acquisition date
                    # Since burst names don't always contain A/B designator, default to S1A if unsure
                    sat = "S1A"
                    
                    # Combine satellite identifier and date
                    sat_date = sat + date_str
                    unique_dates.add(sat_date)
                    print(f"Extracted date {date_str} from burst {g}")
            except Exception as e:
                print(f"Could not parse date from burst granule {g}: {e}")
        else:
            # Standard SLC granule format (S1A_IW_SLC__1SDV_20171117T015310_...)
            try:
                sat = g.split('_')[0][:3]  # Extract S1A or S1B
                date_str = g.split('_')[5][0:8]  # Extract YYYYMMDD part
                sat_date = sat + date_str
                unique_dates.add(sat_date)
                print(f"Extracted date {date_str} from standard granule {g}")
            except Exception as e:
                print(f"Could not parse date from standard granule {g}: {e}")
    
    # Convert to list and sort
    sat_dates = list(unique_dates)
    sat_dates.sort()
    
    if len(sat_dates) == 0:
        print("No valid dates could be extracted from granule names. Cannot download orbit files.")
        return
    
    print(f"Found {len(sat_dates)} unique acquisition dates for orbit files")
    print(sat_dates)
    
    # Get orbit URLs for each satellite/date combination
    orbUrls = asfQuery.get_orbit_url(sat_dates)
    
    # Filter out None values from orbit URLs
    orbUrls = [url for url in orbUrls if url is not None]
    
    if len(orbUrls) == 0:
        print("No orbit files found for the given granules.")
        return
    
    # Remove duplicates and sort
    orbUrls = list(set(orbUrls))
    orbUrls.sort()
    
    print('Searching for orbits...')

    outNames = []
    dlorbs = []
    for url in orbUrls:
        fname = os.path.join(outdir, url.split('/')[-1])
        if not os.path.isfile(fname):
            print('No file ' + fname)
            outNames.append(fname)
            dlorbs.append(url)
        else:
            if os.path.getsize(fname) < 1024: 
                print('Overwriting ' + fname + ' because it was too small...')
                outNames.append(fname)
                dlorbs.append(url)
    
    if len(dlorbs) == 0:
        print("All orbit files already exist and are valid.")
        return
    
    print('Downloading orbit files...')
    print(dlorbs)

    # Download urls in parallel and in chunks
    with concurrent.futures.ThreadPoolExecutor(max_workers=nproc) as executor:
        futures = [executor.submit(dl, url, outName) for url, outName in zip(dlorbs, outNames)]
        concurrent.futures.wait(futures)

    # Verify files were downloaded correctly
    outNames = []
    dlorbs = []
    redflag = False

    for url in orbUrls:
        fname = os.path.join(outdir, url.split('/')[-1])
        if not os.path.isfile(fname):
            print('No file ' + fname)
            outNames.append(fname)
            dlorbs.append(url)
            redflag = True
        else:
            if os.path.getsize(fname) < 1024: 
                print('Deleting ' + fname + ' because it was too small...')
                outNames.append(fname)
                os.remove(fname)
                dlorbs.append(url)
                redflag = True
    
    if redflag:
        print('Some orbit files may have not been properly downloaded. Please try again.')

def dlSlc(slcUrls, gran, outdir, sat='S1'):
    """
    Download SLC or burst files from the provided URLs.
    
    Parameters
    ----------
    slcUrls : list
        List of URLs to download
    gran : list
        List of granule names corresponding to URLs
    outdir : str
        Output directory for downloads
    sat : str
        Satellite identifier (default: 'S1')
    """
    # Create output directory if it doesn't exist
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    
    # Determine if we're dealing with burst data
    is_burst = any("-BURST" in g for g in gran)
    
    if is_burst:
        print("\n=== BURST DATA DOWNLOAD INFORMATION ===")
        print("The URLs provided are for individual Sentinel-1 burst TIFF files.")
        print("ASF requires authentication and a different download mechanism for burst data.")
        print("To download burst data, you have three options:")
        
        print("\nOption 1: Use ASF's API for authenticated programmatic access:")
        print("1. Create an account at https://urs.earthdata.nasa.gov/")
        print("2. Create a .netrc file in your home directory with:")
        print("   machine urs.earthdata.nasa.gov")
        print("   login YOUR_USERNAME")
        print("   password YOUR_PASSWORD")
        
        print("\nOption 2: Use ASF Vertex to download burst data:")
        print("1. Visit https://search.asf.alaska.edu/")
        print("2. Search for your area and time range")
        print("3. Filter for Sentinel-1 Burst Products")
        print("4. Select and download the data manually")
        
        print("\nOption 3: Use the ASF Python API with authenticated credentials:")
        print("1. Install the asf_search package: pip install asf_search")
        print("2. Use your Earthdata credentials with the package")
        print("3. Example script available at: https://github.com/asfadmin/Discovery-asf_search/")
        
        print("\n=== CREATING PLACEHOLDER FILES ===")
        print("Creating placeholder files for the burst data to allow workflow to continue.")
        print("You'll need to replace these with actual data later.")
        
        # Create placeholder files for now
        for i, (url, name) in enumerate(zip(slcUrls, gran)):
            placeholder_path = os.path.join(outdir, f"{name}.zip")
            
            # Get info about which burst this is
            burst_parts = url.split('/')
            scene_id = burst_parts[-4]  # S1B_IW_SLC__1SDV_...
            swath = burst_parts[-3]     # IW3
            pol = burst_parts[-2]       # VH or VV
            burst_num = burst_parts[-1].replace('.tiff', '')  # 6
            
            # Create a placeholder with information for manual download
            with open(placeholder_path, 'w') as f:
                f.write(f"PLACEHOLDER FOR SENTINEL-1 BURST DATA\n")
                f.write(f"URL: {url}\n")
                f.write(f"Granule: {name}\n")
                f.write(f"Scene ID: {scene_id}\n")
                f.write(f"Swath: {swath}\n")
                f.write(f"Polarization: {pol}\n")
                f.write(f"Burst Number: {burst_num}\n")
                f.write("\nThis is a placeholder file. You need to download the actual burst data.\n")
                f.write("See instructions printed in the console output.\n")
            
            print(f"Created placeholder for {name}")
        
        print("\n=== CONTINUING WITH WORKFLOW ===")
        print("The workflow will continue with orbit files and DEM, but burst data is missing.")
        print("You'll need to manually download the burst data before proceeding to processing steps.")
        return
    
    # For regular SLC files, continue with normal download process
    # Check for existing files and determine what we need to download
    print("Checking for existing files and preparing download list...")
    outNames = []
    dlSLCs = []
    existing_files = []
    
    for ii in range(len(gran)):
        fname = os.path.join(outdir, gran[ii] + '.zip')
        
        # If the file doesn't exist, add it to download list
        if not os.path.isfile(fname):
            outNames.append(fname)
            dlSLCs.append(slcUrls[ii])
        else:
            # Check if the existing file is a valid size
            if os.path.getsize(fname) < 1024*1024:  # At least 1 MB for all files
                print(f"File exists but is too small, re-downloading: {fname}")
                outNames.append(fname)
                dlSLCs.append(slcUrls[ii])
            else:
                print(f"Skipping existing valid file: {fname}")
                existing_files.append(fname)
    
    if not dlSLCs:
        print("All files already exist and appear to be valid.")
        return
    
    # Print download summary
    print(f"Downloading {len(dlSLCs)} files out of {len(gran)} total files")
    print("URLs to download:")
    for url in dlSLCs:
        print(f"  {url}")
    
    min_expected_size = 1024*1024*1024  # 1 GB for SLCs
    
    print(f"Detected file type: SLC")
    print(f"Minimum expected file size: {min_expected_size/(1024*1024):.1f} MB")
    
    # Download files with progress tracking
    download_failed = []
    print(f"Starting downloads with {nproc} parallel threads...")
    
    # First check if we can access the first URL to determine if authentication is needed
    try:
        test_url = dlSLCs[0]
        print(f"Testing download access with URL: {test_url}")
        response = requests.head(test_url, allow_redirects=True, timeout=30)
        
        if response.status_code == 401 or response.status_code == 403:
            print("Authentication required for ASF downloads.")
            print("You may need to add credentials to your .netrc file or use ASF's API for authenticated downloads.")
            print("See: https://asf.alaska.edu/how-to/data-tools/data-tool-use/")
            return
        
        expected_size = int(response.headers.get('Content-Length', 0))
        print(f"Expected file size from headers: {expected_size/(1024*1024):.1f} MB")
        
        if expected_size < 1024*1024:  # Less than 1 MB
            print("Warning: The server reports a very small file size. This may indicate:")
            print("1. The download requires authentication")
            print("2. The URL points to a redirect or error page")
            print("3. The actual data is not available at this URL")
    except Exception as e:
        print(f"Warning: Could not check URL access: {e}")
        
    # Now proceed with the downloads
    with concurrent.futures.ThreadPoolExecutor(max_workers=nproc) as executor:
        futures = [executor.submit(dl, url, outName) for url, outName in zip(dlSLCs, outNames)]
        concurrent.futures.wait(futures)
    
    # Verify the downloads
    print("Verifying downloaded files...")
    for ii, (url, outName) in enumerate(zip(dlSLCs, outNames)):
        # Check if file exists
        if not os.path.isfile(outName):
            print(f"Error: File not downloaded: {outName}")
            download_failed.append(url)
            continue
        
        # Check file size
        file_size = os.path.getsize(outName)
        if file_size < 1024:  # Less than 1 KB
            print(f"Error: Downloaded file is too small (only {file_size} bytes): {outName}")
            print("This could be an error page or empty file.")
            download_failed.append(url)
            continue
            
        if file_size < min_expected_size:
            print(f"Warning: File is smaller than expected ({file_size/(1024*1024):.1f} MB): {outName}")
            try:
                # Try to open as zip to verify
                with open(outName, 'rb') as f:
                    header = f.read(4)
                if header != b'PK\x03\x04':  # Zip file signature
                    print(f"Error: File is not a valid zip file: {outName}")
                    download_failed.append(url)
            except Exception as e:
                print(f"Error checking file: {e}")
                download_failed.append(url)
    
    # Report download status
    if download_failed:
        print(f"\nWARNING: {len(download_failed)} files failed to download properly:")
        for url in download_failed:
            print(f"  {url}")
        print("\nPossible reasons for download failure:")
        print("1. Authentication is required (check if you need ASF credentials)")
        print("2. The URLs may be incorrect or data is not available")
        print("3. Network issues or server-side problems")
        print("\nSuggested solutions:")
        print("1. Check the ASF website to verify the data is downloadable")
        print("2. Try manually downloading one file from the ASF website")
        print("3. Set up a .netrc file with your ASF credentials if needed")
        
        # Don't exit with error, so the script can continue if possible
        print("\nContinuing with available data...")
    else:
        print(f"\nSuccessfully downloaded {len(dlSLCs)} files")
        
    return

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
    
    # If we don't have SLC zip files, use the bounds from the config
    if not zips:
        print("No SLC zip files found. Using bounds from configuration.")
        if ps.bounds:
            try:
                minlat = np.floor(float(ps.bounds.split(',')[0]))
                maxlat = np.ceil(float(ps.bounds.split(',')[1]))
                minlon = np.floor(float(ps.bounds.split(',')[2]))
                maxlon = np.ceil(float(ps.bounds.split(',')[3]))
                print(f"Using bounds from config: {minlat}, {maxlat}, {minlon}, {maxlon}")
            except Exception as e:
                print(f"Error parsing bounds from config: {e}")
                print("Using default bounds for Colorado")
                # Default bounds for Colorado if all else fails
                minlat, maxlat = 36.0, 42.0
                minlon, maxlon = -110.0, -102.0
                print(f"Using default bounds: {minlat}, {maxlat}, {minlon}, {maxlon}")
        elif ps.poly:
            try:
                # Extract approximate bounds from polygon
                print("Extracting bounds from polygon")
                geom = shapely.wkt.loads(ps.poly)
                minlon, minlat, maxlon, maxlat = geom.bounds
                minlat = np.floor(minlat)
                maxlat = np.ceil(maxlat)
                minlon = np.floor(minlon)
                maxlon = np.ceil(maxlon)
                print(f"Bounds from polygon: {minlat}, {maxlat}, {minlon}, {maxlon}")
            except Exception as e:
                print(f"Error extracting bounds from polygon: {e}")
                print("Using default bounds for Colorado")
                # Default bounds for Colorado if all else fails
                minlat, maxlat = 36.0, 42.0
                minlon, maxlon = -110.0, -102.0
                print(f"Using default bounds: {minlat}, {maxlat}, {minlon}, {maxlon}")
        else:
            print("No bounds or polygon found in config. Using default bounds for Colorado.")
            # Default bounds for Colorado if all else fails
            minlat, maxlat = 36.0, 42.0
            minlon, maxlon = -110.0, -102.0
            print(f"Using default bounds: {minlat}, {maxlat}, {minlon}, {maxlon}")
    else:
        # Figure out what bounds to use for the DEM from the SLC zip files
        minlats, maxlats, minlons, maxlons = [], [], [], []
        if ps.sat=='SENTINEL-1':
            for z in zips:
                try:
                    safe = sentinelSLC(z)
                    safe.get_lat_lon_v2()
                    minlats.append(safe.SNWE[0])
                    maxlats.append(safe.SNWE[1])
                    minlons.append(safe.SNWE[2])
                    maxlons.append(safe.SNWE[3])
                    print(f"Extracted bounds from {z}: {safe.SNWE}")
                except Exception as e:
                    print(f"Error extracting bounds from {z}: {e}")
                    continue

            if minlats and maxlats and minlons and maxlons:
                minlat = min(minlats)
                maxlat = max(maxlats)
                minlon = min(minlons)
                maxlon = max(maxlons)
                print(f"Computed bounds from SLCs: {minlat}, {maxlat}, {minlon}, {maxlon}")
            else:
                # Fallback if we couldn't extract bounds from any SLCs
                print("Could not extract bounds from SLCs. Using bounds from config.")
                if ps.bounds:
                    minlat = np.floor(float(ps.bounds.split(',')[0]))
                    maxlat = np.ceil(float(ps.bounds.split(',')[1]))
                    minlon = np.floor(float(ps.bounds.split(',')[2]))
                    maxlon = np.ceil(float(ps.bounds.split(',')[3]))
                    print(f"Using bounds from config: {minlat}, {maxlat}, {minlon}, {maxlon}")
                else:
                    print("No bounds in config. Using default bounds for Colorado.")
                    minlat, maxlat = 36.0, 42.0 
                    minlon, maxlon = -110.0, -102.0
                    print(f"Using default bounds: {minlat}, {maxlat}, {minlon}, {maxlon}")
        else:
            # For non-Sentinel-1 satellites
            if ps.bounds:
                minlat = np.floor(float(ps.bounds.split(',')[0]))
                maxlat = np.ceil(float(ps.bounds.split(',')[1]))
                minlon = np.floor(float(ps.bounds.split(',')[2]))
                maxlon = np.ceil(float(ps.bounds.split(',')[3]))
                print(f"Using bounds from config for non-Sentinel-1: {minlat}, {maxlat}, {minlon}, {maxlon}")
            else:
                print("No bounds in config for non-Sentinel-1. Using default bounds.")
                minlat, maxlat = 36.0, 42.0
                minlon, maxlon = -110.0, -102.0
                print(f"Using default bounds: {minlat}, {maxlat}, {minlon}, {maxlon}")

    # Create directories if they don't exist
    if not os.path.isdir('Figs'):
        os.mkdir('Figs')
    if not os.path.isdir('Npy'):
        os.mkdir('Npy')

    # Create a simple figure showing the bounds
    try:
        plt.figure(figsize=(8, 6))
        plt.plot([minlon, maxlon, maxlon, minlon, minlon], [minlat, minlat, maxlat, maxlat, minlat], 'r-')
        plt.title(f"DEM Bounds: {minlat:.2f}, {maxlat:.2f}, {minlon:.2f}, {maxlon:.2f}")
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        plt.grid(True)
        plt.savefig('Figs/bounds.png')
        plt.close()
        print("Saved bounds figure to Figs/bounds.png")
    except Exception as e:
        print(f"Error creating bounds figure: {e}")
        # Continue even if figure creation fails

    # Format the DEM bounds string
    demBounds = f"{int(np.floor(minlat))},{int(np.ceil(maxlat))},{int(np.floor(minlon))},{int(np.ceil(maxlon))}"
    print(f"DEM bounds: {demBounds}")

    # Create DEM directory if it doesn't exist
    if not os.path.isdir('DEM'):
        os.mkdir('DEM')

    # Download dem if it doesn't exist
    file_list = glob.glob('./DEM/*wgs84*')
    if not file_list:
        print(f"No DEM found. Downloading DEM with bounds: {demBounds}")
        if inps.get_srtm:
            getDEM.getDEM(demBounds, srtm=True)
            DEM = glob.glob(os.path.join(ps.workdir, 'DEM', '*wgs84'))[0]
        else:
            getDEM.getDEM(demBounds)
            DEM = glob.glob(os.path.join(ps.workdir, 'DEM', '*wgs84.dem'))[0]
    else:
        print(f"DEM already exists: {file_list[0]}")
        if inps.get_srtm:
            DEM = glob.glob(os.path.join(ps.workdir, 'DEM', '*wgs84'))[0]
        else:
            DEM = glob.glob(os.path.join(ps.workdir, 'DEM', '*wgs84.dem'))[0]

    # Updating DEM's wgs84 xml to include the full path to the DEM
    print(f"Fixing XML for DEM: {DEM}")
    os.system(f'fixImageXml.py -f -i {DEM} >> log')

    if os.path.exists(DEM):
        print(f"DEM ready at: {DEM}")
    else:
        print(f'Error: DEM does not exist at expected path: {DEM}')
        sys.exit(1)

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
