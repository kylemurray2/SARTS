#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 15:24:40 2023

Check that the SLCs can be opened and set up the run_files and configs

@author: km
"""

import numpy as np
import os,glob,sys,zipfile
import stackSentinel
from stackSentinel import sentinelSLC
from SARTS import config

flag = False
doRemove = True

def checkSLC(ps,doRemove=True):
    flag = True

    # Check to make sure all the files are big enough and the zip files are valid
    zips = glob.glob(ps.slc_dirname + '*zip')
    
    # Check if we're dealing with burst data
    is_burst_data = any("-BURST" in os.path.basename(z) for z in zips)
    
    if is_burst_data:
        print("Detected Sentinel-1 burst data files")
        print("Using special validation for burst products")
        
        for z in zips:
            file_basename = os.path.basename(z)
            # For burst placeholder files, check if they have enough information
            try:
                with open(z, 'r') as f:
                    content = f.read()
                    
                if "PLACEHOLDER FOR SENTINEL-1 BURST DATA" in content:
                    print(f"{file_basename} is a placeholder file - you need to download the actual burst data")
                    flag = False
                elif len(content) < 100:  # Very small files are likely invalid
                    print(f"{file_basename} appears to be too small or corrupt")
                    if doRemove:
                        print(f"Removing invalid file: {z}")
                        os.remove(z)
                    flag = False
                else:
                    # This might be a valid burst TIFF wrapped in a zip or a proper burst product
                    print(f"{file_basename} appears to be a burst data file")
            except UnicodeDecodeError:
                # This might be a binary file - could be a real downloaded burst
                try:
                    # Try to validate if it's a zip file
                    with zipfile.ZipFile(z, 'r') as zf:
                        # It's a valid zip file
                        file_sizes = [info.file_size for info in zf.infolist()]
                        total_size = sum(file_sizes)
                        
                        if total_size < 1024 * 1024:  # Less than 1 MB
                            print(f"{file_basename} is a valid zip but seems too small for burst data")
                            if doRemove:
                                print(f"Removing potentially corrupt file: {z}")
                                os.remove(z)
                            flag = False
                        else:
                            print(f"{file_basename} validated as a zip file with {total_size/1024/1024:.2f} MB content")
                except zipfile.BadZipFile:
                    # Not a valid zip, might be a direct TIFF or other format
                    file_size = os.path.getsize(z)
                    if file_size > 5 * 1024 * 1024:  # Greater than 5 MB
                        print(f"{file_basename} is not a zip but has substantial size ({file_size/1024/1024:.2f} MB)")
                        # We'll consider this potentially valid
                    else:
                        print(f"{file_basename} is neither a valid zip nor of substantial size")
                        if doRemove:
                            print(f"Removing invalid file: {z}")
                            os.remove(z)
                        flag = False
        
        if not flag:
            print("\n===== BURST DATA VALIDATION FAILED =====")
            print("You need to download the actual burst data from ASF.")
            print("Options:")
            print("1. Use ASF Vertex: https://search.asf.alaska.edu/")
            print("2. Use asf_search Python API with your Earthdata Login credentials")
            print("\nExample code for asf_search:")
            print("""
import asf_search as asf
session = asf.ASFSession().auth_with_creds('your_username', 'your_password')
results = asf.search(
    platform='SENTINEL-1',
    processingLevel='BURST',
    relativeOrbit=151,
    start='2020-01-01T00:00:00UTC',
    end='2020-03-01T23:59:59UTC',
    intersectsWith='POLYGON((-105.2751 39.8385,-105.0987 39.8385,-105.0987 39.9177,-105.2751 39.9177,-105.2751 39.8385))'
)
results.download(path="./SLCS")
            """)
        else:
            print("\n===== BURST DATA VALIDATION PASSED =====")
            print("All burst files appear to be valid")
            
        return flag
    
    # Standard SLC checking for non-burst products
    for z in zips:
        zipSize = os.stat(z).st_size
        if zipSize < 1e9:
            print(f'May want to delete {z} because it is too small.')
            if doRemove:
                print(f"File {z} is corrupt. Deleting...")
                os.remove(z)
            else:
                print(f"File {z} is corrupt. Recommend deleting.")
            flag = False
            
    zips = glob.glob(ps.slc_dirname + '*zip')
    for z in zips:
        try:
            safe = sentinelSLC(z)
            safe.get_lat_lon_v2()
            # x = zipfile.ZipFile(z)
            # x.close()
            print(f"{z} opened ok")
        except:
            if doRemove:
                print(f"File {z} is corrupt. Deleting...")
                os.remove(z)
            else:
                print(f"File {z} is corrupt. Recommend deleting.")
            flag = False
            continue
    
    return flag

def main():
    ps = config.getPS()
    inps=ps
    
    # Check if we're working with burst data
    zips = glob.glob(ps.slc_dirname + '*zip')
    is_burst_data = any("-BURST" in os.path.basename(z) for z in zips)
    
    if is_burst_data:
        print("\n======= SENTINEL-1 BURST DATA DETECTED =======")
        print("The standard stackSentinel workflow is designed for full SLC products, not individual bursts.")
        print("You have two options for processing burst data:")
        
        print("\n1. DOWNLOAD FULL SLC PRODUCTS:")
        print("   Instead of using burst data, download the full Sentinel-1 SLC scenes that contain these bursts.")
        print("   This will work with the current stackSentinel workflow.")
        print("   You can find the corresponding SLCs by searching at https://search.asf.alaska.edu/")
        
        print("\n2. USE SPECIALIZED BURST PROCESSING:")
        print("   ISCE has dedicated tools for working with burst data, but they use a different workflow.")
        print("   For burst processing with ISCE, you'll need to use the following approach:")
        print("   a. Download the burst data using ASF search tools")
        print("   b. Use topsApp.py with appropriate XML configuration instead of stackSentinel")
        print("   c. Set up proper burst selection in your topsApp.xml file")
        
        print("\nFor burst processing documentation, see:")
        print("https://github.com/isce-framework/isce2/blob/main/examples/input_files/topsApp.xml.burst_processing")
        
        print("\nExample workflow for burst processing:")
        print("""
# Example workflow for burst processing with ISCE:
1. Install asf_search: pip install asf_search
2. Download burst data:
   ```python
   import asf_search as asf
   session = asf.ASFSession().auth_with_creds('your_username', 'your_password')
   results = asf.search(
       platform='SENTINEL-1',
       processingLevel='BURST',
       relativeOrbit=151,
       start='2020-01-01T00:00:00UTC',
       end='2020-03-01T23:59:59UTC',
       intersectsWith='POLYGON((-105.2751 39.8385,-105.0987 39.8385,-105.0987 39.9177,-105.2751 39.9177,-105.2751 39.8385))'
   )
   results.download(path="./bursts")
   ```
3. Create a topsApp.xml file with burst parameters
4. Run: topsApp.py topsApp.xml
        """)
        return
    
    if ps.sat=='SENTINEL-1':
        flag = checkSLC(ps)
    else:
        print('Skipping SLC check because not processing Sentinel-1')
        flag=True
        
    if flag:
        # stackSentinel.main(ps)
        if os.path.exists(os.path.join(ps.work_dir, 'run_files')):
            print('')
            print('**************************')
            print('run_files folder exists.')
            print(os.path.join(ps.work_dir, 'run_files'), ' already exists.')
            print('Please remove or rename this folder and try again.')
            print('')
            print('**************************')
            sys.exit(1)

        acquisitionDates, stackReferenceDate, secondaryDates, safe_dict, updateStack = stackSentinel.checkCurrentStatus(ps)
        # selecting pairs for interferograms / correlation / offset workflows
        pairs = stackSentinel.selectNeighborPairs(acquisitionDates, stackReferenceDate, secondaryDates, ps.num_connections, updateStack)

        print ('*****************************************')
        print ('Coregistration method: ', ps.coregistration )
        print ('Workflow: ', ps.workflow)
        print ('*****************************************')

        i = stackSentinel.slcStack(ps, acquisitionDates, stackReferenceDate, secondaryDates, safe_dict, updateStack, mergeSLC=True)

        #Checks presence of ion parameter file. If it exists, do ionosphere estimation.
        if ps.param_ion is None:
            print("Ion parameter file is not specified. Ionospheric estimation will not be done.")
        elif not os.path.isfile(ps.param_ion):
            print("Ion parameter file is missing. Ionospheric estimation will not be done.")
        else:
            dateListIon, pairs_same_starting_ranges_update, pairs_diff_starting_ranges_update, safe_dict = stackSentinel.checkCurrentStatusIonosphere(ps)
            i = stackSentinel.ionosphereStack(ps, dateListIon, stackReferenceDate, pairs_same_starting_ranges_update, pairs_diff_starting_ranges_update, safe_dict, i)

        print('Next step is to make the co-registered SLC stack with runISCE.py')

    else:
        print("Failed file check. Make sure all files are fully downloaded or delete corrupt files.")
        
if __name__ == '__main__':
    '''
    Main driver.
    '''
    # inps = cmdLineParser()
    main()