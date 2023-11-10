#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 12:14:22 2022

To add more dates later, we need to keep the last six dates (at least for the
    safe files) as well as the reference date (ps.reference_date)


To add more data to the stack you need the following:
    SLCS
    orbits
    DEM
    coreg_secondarys
    geom_reference
    secondarys
    reference

Delete:
    run_files


@author: km
"""

import os
import time
import glob
import numpy as np
from SARTS import config

ps = config.getPS()

# It will keep the last six dates in secondarys and coreg_secondarys
delList = ['SLCS',
           'coarse_interferograms',
           'coreg_secondarys',
           'secondarys',
           'ESD',
           'misreg',
           'orbits']

if False:
    coregSLCS = 'merged/SLC/2*/*.full'
    os.system('rm -r ' + coregSLCS)

print('WARNING: About to delete: ')
print(delList)
countdown = 5
for ii in range(countdown):
    print('\r' + str(countdown), end=' ')
    countdown-=1
    time.sleep(1)


dates = []
datesFN = glob.glob(os.path.join(ps.mergeddir,'SLC/*'))
for d in datesFN:
    dates.append(d.split('/')[-1])
dates.sort()


numDatesToKeep = 6 #len(ps.dates) - ps.dates.index(ps.reference_date) + 6
numDatesToKeep+=1
print('keeping dates: ')
print(dates[-numDatesToKeep:])


if 'coreg_secondarys' in delList:
    print('removing coreg_secondarys')
    for d in dates[:-numDatesToKeep]:
        if d != ps.reference_date: #Don't delete the reference date
            # print('removing ./coreg_secondarys/' + d)
            os.system('rm -r ./coreg_secondarys/' + d + '/*')

if 'secondarys' in delList:
    print('removing secondarys')
    for d in dates[:-numDatesToKeep]:
        if d != ps.reference_date: #Don't delete the reference date
            # print('removing ./coreg_secondarys/' + d)
            os.system('rm -r ./secondarys/' + d + '/*')

d_keep = dates[-numDatesToKeep:]

if 'SLCS' in delList:
    print('removing safe files')
    safeList = glob.glob(ps.slc_dirname + '*zip') # make a list of all the safe files
    for ii,d in enumerate(dates[:-numDatesToKeep]):  # Loop through the dates
        for safe in safeList:           # Delete all safe files with that date in the name
            if d in safe:
                if d != ps.reference_date: #Don't delete the reference date
                    # print('removing ' + safe)
                    os.remove(safe)

if 'orbits' in delList:
    print('removing orbits')
    orbList = glob.glob('./orbits/*EOF') # make a list of all the safe files
    for ii,d in enumerate(dates[:-numDatesToKeep]):  # Loop through the dates
        for orb in orbList:           # Delete all safe files with that date in the name
            if d in orb:
                if d != ps.reference_date: #Don't delete the reference date
                    # print('removing ' + orb)
                    os.remove(orb)
