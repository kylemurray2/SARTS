#!/usr/bin/bash
# -*- coding: utf-8 -*-

# Search data, download slcs, download orbits
downloadData.py -sdo
# Check stack and write configs and run_files
setupStack.py 
# Run the files in run_files and output to log files
runISCE.py
# Downlook and Crop geom files, replace existing, plot on or off
adjustGeom.py -dcrp
# Run Fringe
runFringe.py
# Make ifgs, downlook, filter, unwrap
ifgs.py -dum -n 8



