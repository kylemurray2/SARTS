#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: km

Use this if you want to change the number of looks and keep the originals

moves dolphin/interferograms directory to interferograms_alks_rlks
makes new interferograms dir
links all of the full res ifgs from the original to the new in respective date pair dirs
updates params.yaml with new rlks/alks if given

Then you should be able to just rerun ifgs.py -du to make new dl ifgs and unw

"""

import os,argparse
from SARTS import config,util


def cmdLineParser():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser(
        description='Crop and downlook geom files. Save parameters',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-a', '--new-alks',type=int, dest='alks_new', default=None, help='Value to replace alks')
    parser.add_argument('-r', '--new-rlks',type=int, dest='rlks_new', default=None, help='Value to replace rlks')

    return parser.parse_args()

def link_ifgs(ps, alks_orig, rlks_orig):
    source_base_dir = os.path.join(ps.workdir,ps.dolphin_work_dir, f'interferograms_{alks_orig}_{rlks_orig}')
    print(source_base_dir)
    target_base_dir = os.path.join(ps.workdir,ps.dolphin_work_dir, 'interferograms')
    
    for p in ps.pairs_seq:
        source_dir = os.path.join(source_base_dir, p)
        target_dir = os.path.join(target_base_dir, p)
        os.makedirs(target_dir, exist_ok=True)  # Ensure the target directory exists
    
        for file_ext in ['.int.xml', '.int.vrt', '.int']:
            source_file = os.path.join(source_dir, p + file_ext)
            target_file = os.path.join(target_dir, p + file_ext)
    
            if os.path.exists(source_file) and not os.path.exists(target_file):
                os.symlink(source_file, target_file)
                print(f"Link created: {target_file}")
            else:
                print(f"File not found or link already exists: {source_file} -> {target_file}")
    
    # now update the yaml params file with new alks and rlks values
    print('updating looks in params.yaml')


def main(inps):
    ps = config.getPS()
    alks_orig = ps.alks
    rlks_orig = ps.rlks
    link_ifgs(ps, alks_orig, rlks_orig)
    
    if inps.alks_new:
        util.update_yaml_key('params.yaml', 'alks', inps.alks_new)
    if inps.rlks_new:
        util.update_yaml_key('params.yaml', 'rlks', inps.rlks_new)
    
    
if __name__ == '__main__':
    '''
    Main driver.
    '''
    inps = cmdLineParser()
    main(inps)