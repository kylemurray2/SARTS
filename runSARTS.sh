#!/usr/bin/bash
# -*- coding: utf-8 -*-

downloadData.py --search-data True --download-slc True --download-orbits True
setupStack.py 
runISCE.py 