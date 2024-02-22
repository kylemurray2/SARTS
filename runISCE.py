#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 16:16:14 2023

Author: km
"""

import os
import glob
import time
import numpy as np
from SARTS import config

ps = config.getPS()
# Incrementally delete files to save space (use with caution)
delFiles = False


# Find and sort run scripts
runScripts = sorted(glob.glob(os.path.join(ps.workdir, 'run_files', 'run*')))

# Start and stop indices for running scripts
start = 0
stop = len(runScripts)

# Create a directory for log files if it doesn't exist
if not os.path.isdir('./logFiles'):
    os.mkdir('./logFiles')

startT = time.time()




def append_wait_to_files(file_list):
    """
    Append the string "wait" to the end of each file in the provided list.
    This is needed so each parallel process finishes before the next file is run
    :param file_list: List of file paths to be processed.
    """
    for filename in file_list:
        with open(filename, 'a') as f:
            f.write("\nwait\n")

append_wait_to_files(runScripts)

for ii in range(start, stop):
    script_path = runScripts[ii]
    script_name = os.path.basename(script_path)
    print('Running ' + script_name)
    
    # Run the script and redirect the output to a log file
    log_path = os.path.join('./logFiles', 'runlog_' + str(ii + 1))
    os.system(f'bash {script_path} > {log_path}')
    
    if delFiles:
        if ii == 5:
            os.system('rm coreg_secondarys/*/overlap/IW?/*off*')
        elif ii == 6:
            os.system('rm -r ESD coarse_interferograms')
        elif ii == 9:
            os.system('rm coreg_secondarys/*/IW*/*off*')
            os.system('mv SLCS SLCS_o')
            os.system('mkdir SLCS')
            os.system(f'mv SLCS_o/*{ps.reference_date}* SLCS/')
            os.system('rm -r SLCS_o')

endT = time.time()
elapsed = endT - startT
print(f'This run took {np.round(elapsed/60)} minutes.')
