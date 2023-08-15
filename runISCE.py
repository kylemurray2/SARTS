#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 16:16:14 2023

@author: km
"""

import numpy as np
from matplotlib import pyplot as plt
import os
import glob
import time

ps = np.load('./ps.npy',allow_pickle=True).all()

delFiles = False # Incrementally delete files to save space. (can be dangerous)
start = 0
stop = 13

runScripts = glob.glob(ps.workdir + '/run_files/run*')
runScripts.sort()

if not os.path.isdir('./logFiles'):
    os.mkdir('./logFiles')
    
startT = time.time()

for ii in np.arange(start,stop):
    print('running ' + runScripts[ii].split('/')[-1])
    os.system('bash ' + runScripts[ii] + ' > logFiles/runlog_' + str(ii+1) )
    
    if delFiles:
        if ii==5:
            os.system('rm coreg_secondarys/*/overlap/IW?/*off*')
        if ii==6:
            os.system('rm -r ESD coarse_interferograms')
        if ii==9:
            os.system('rm coreg_secondarys/*/IW*/*off*')
            os.system('mv SLCS SLCS_o')
            os.system('mkdir SLCS')
            os.system('mv SLCS_o/*' + ps.reference_date + '* SLCS/')
            os.system('rm SLCS_o')

endT = time.time()
elapsed = endT-startT
print('This run took ' + str(np.round(elapsed/60)) + ' minutes.')
