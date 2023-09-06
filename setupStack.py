#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 15:24:40 2023

Check that the SLCs can be opened and set up the run_files and configs

@author: km
"""

import numpy as np
import os
import glob
import zipfile
import stackSentinel
import localParams
from stackSentinel import sentinelSLC

ps = localParams.getLocalParams()
flag = False
doRemove = True

def checkSLC(doRemove=True):
    flag = True

    # Check to make sure all the files are big enough and the zip files are valid
    zips = glob.glob(ps.slc_dirname + '*zip')
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
    flag = checkSLC()

    if flag:
        stackSentinel.main(ps)
    else:
        print("Failed file check. Make sure all files are fully downloaded or delete corrupt files.")

if __name__ == '__main__':
    '''
    Main driver.
    '''
    # inps = cmdLineParser()
    main()