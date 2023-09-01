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

ps = localParams.getLocalParams()
flag = True

# Check to make sure all the files are big enough and the zip files are valid
zips = glob.glob(ps.slc_dirname + '*zip')
for z in zips:
    zipSize = os.stat(z).st_size
    if zipSize < 1e9:
        print(f'May want to delete {z} because it is too small.')
        # os.remove(z)
        flag = False

for z in zips:
    try:
        x = zipfile.ZipFile(z)
        print(f"{z} opened ok")
        x.close()
    except:
        print(f"File {z} is corrupt. Recommend deleting.")
        flag = False
        # os.remove(z)
        continue

if flag:
    stackSentinel.main(ps)
else:
    print("Failed file check. Make sure all files are fully downloaded or delete corrupt files.")
