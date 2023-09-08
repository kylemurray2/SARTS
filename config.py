#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 15:46:44 2023

@author: km
"""

import numpy as np
import os
import yaml
import argparse
import sys

def load_yaml_to_namespace(yaml_file):
    # Load the YAML file into a dictionary
    with open(yaml_file, 'r') as yaml_in:
        yaml_dict = yaml.safe_load(yaml_in)

    # Create a namespace from the dictionary
    namespace = argparse.Namespace(**yaml_dict)
    
    return namespace

def getPS():
    
    # Load the configs from the yaml file
    yaml_file = 'configs.yaml'
    
    
    if os.path.isfile(yaml_file):
        print('Parsing yaml file and updating ps namespace...')
        configs = load_yaml_to_namespace(yaml_file)
        # Load the ps namespace
        if os.path.isfile('./ps.npy'):
            ps = np.load('./ps.npy',allow_pickle=True).all()
        else:
            ps = configs
        
        # Update ps with any changes to configs
        for attr in dir(configs):
            if not attr.startswith('_'):
                print(attr)
                setattr(ps, attr, getattr(configs, attr))
    
# Set up some additional variables 
        ps.workdir = os.getcwd()
        ps.sensor = ps.sat
        ps.startDate               = ps.start[0:10]
        ps.stopDate                = ps.end[0:10]
        miny,maxy,minx,maxx        = ps.bounds.split(sep=',')
        ps.bbox                    = miny +' '+ maxy  +' '+ minx +' '+  maxx #demBounds[0] + ' ' + demBounds[1] + ' ' + demBounds[2] + ' ' + demBounds[3] #SNWE
        
        if ps.numProcess=='auto':
            ps.numProcess = os.cpu_count()
            ps.numProcess4topo = int(ps.numProcess/3) 

        ps.mergeddir= ps.workdir + '/merged'
        ps.intdir   = ps.mergeddir + '/interferograms'
        ps.tsdir    = ps.workdir + '/TS'
        ps.slcdir   = ps.mergeddir + '/SLC'
        if 'nx' in ps.__dict__.keys():
            ps.nxl = ps.nx//ps.rlks
            ps.nyl = ps.ny//ps.alks


    elif os.path.isfile('localParams.py'):
        print('Using localParams.py...  This will be depricated in future versions. Use a yaml file instead.')
        import localParams
        ps = localParams.getLocalParams()
    
    else:
        print('No config file was found.')
        sys.exit(1)
            
    # Save the updated ps namespace 
    np.save('ps.npy',ps)
    
    return ps