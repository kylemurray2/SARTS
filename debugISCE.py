#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 13:20:16 2022

@author: km
"""

import os
import numpy as np
from glob import glob

runMissing=False

ps = np.load('./ps.npy',allow_pickle=True).all()
coreg_secondarys_dirnames = glob('./coreg_secondarys/*')
coreg_secondarys_dirnames.sort()
dates = coreg_secondarys_dirnames[0].split('/')[2]


for ii in range(len(coreg_secondarys_dirnames)):
    '''
    before you can run run_10_fullBurst_resample these all need to exist
    you may need to run the last script (run_09_fullBurst_geo2rdr)
        e.g., SentinelWrapper.py -c ./configs/config_fullBurst_geo2rdr_20150611
    
    '''
    if not os.path.isfile(coreg_secondarys_dirnames[ii] + '/IW2/range_01.off.xml'):
        print(coreg_secondarys_dirnames[ii])
        if runMissing:
            date = coreg_secondarys_dirnames[ii].split('/')[2]
            os.system('SentinelWrapper.py -c ./configs/config_fullBurst_geo2rdr_' + date)



for d in coreg_secondarys_dirnames:
    '''
    before you can run run_11_extract_stack_valid_region these all need to exist
    you may need to run the last script (run_10_fullBurst_resample)
        e.g., SentinelWrapper.py -c ./configs/config_fullBurst_resample_20200215
    '''
    if not os.path.isfile(d + '/IW2.xml'):
        print(d)
        

