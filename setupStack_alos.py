#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 15:24:40 2023

Check that the SLCs can be opened and set up the run_files and configs

if you get this error: 
TypeError: can only concatenate str (not "int") to str

stripmapStack/Stack line 69 and 70 change to:
        self.f.write('alks : ' + str(self.alks) +'\n')
        self.f.write('rlks : ' + str(self.rlks) +'\n')

@author: km
"""

from stripmapStack import stackStripMap
from SARTS import config

ps = config.getPS()
inps=ps
del(inps.crop)
stackStripMap.main(inps)

