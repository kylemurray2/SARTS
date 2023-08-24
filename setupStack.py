#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 15:24:40 2023

Check that the SLCs can be opened and setup the run_files and configs

@author: km
"""

import numpy as np
import os
import glob
import stackSentinel
import zipfile

ps = np.load('./ps.npy',allow_pickle=True).all()

#check to make sure all the files are big enough and zipfile is valid
zips = glob.glob(ps.slc_dirname + '*zip')
for z in zips:
    zipSize = os.stat(z).st_size
    if zipSize < 1e9:
        print('May want to delete ' + z + ' because it is too small.')
        # os.remove(z)

for z in zips:
    try:
        x = zipfile.ZipFile(z)
        print("%s opened ok" % z)
        x.close()
    except:
        print("File %s is corrupt. Recommend deleting." % z)
       # os.remove(z)
        continue
    
stackSentinel.main(ps)
