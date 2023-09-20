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
    flag = checkSLC(ps)

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