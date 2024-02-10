#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: km

Use this if you want to change the number of looks and keep the originals

moves dolphin/interferograms directory to interferograms_alks_rlks
makes new interferograms dir
links all of the full res ifgs from the original to the new in respective date pair dirs
makes a new merged/geom_reference directory and links files.
updates params.yaml with new rlks/alks if given

Then you should be able to just rerun ifgs.py to make new dl ifgs and unw
But first rerun adjustGeom.py

adjustGeom.py -dr 
ifgs.py -du
prep_mintpy.py -a # -r #

Then you can make another MintPy directory 

"""

import os,argparse,glob
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
    target_base_dir = os.path.join(ps.workdir,ps.dolphin_work_dir, 'interferograms')
    
    if not os.path.isdir(source_base_dir):
        os.system('mv ' + target_base_dir + ' ' + source_base_dir)
    
    for p in ps.pairs:
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


def link_geom(ps, alks_orig, rlks_orig):
    source_base_dir = os.path.join(ps.mergeddir,f'geom_reference_{alks_orig}_{rlks_orig}')
    target_base_dir = os.path.join(ps.mergeddir,f'geom_reference')
    
    if not os.path.isdir(source_base_dir):
        os.system('mv ' + target_base_dir + ' ' + source_base_dir)
        os.mkdir(target_base_dir)

    crop_pns = glob.glob(source_base_dir+'/*')

    for pn in crop_pns: 
        fn = pn.split('/')[-1]
        source_file = os.path.join(source_base_dir, fn)
        target_file = os.path.join(target_base_dir, fn)

        if os.path.exists(source_file) and not os.path.exists(target_file):
            os.symlink(source_file, target_file)
            print(f"Link created: {target_file}")
        else:
            print(f"File not found or link already exists: {source_file} -> {target_file}")
    


def main(inps):
    ps = config.getPS()
    alks_orig = ps.alks
    rlks_orig = ps.rlks
    link_ifgs(ps, alks_orig, rlks_orig)
    link_geom(ps, alks_orig, rlks_orig)
    if inps.alks_new:
        util.update_yaml_key('params.yaml', 'alks', inps.alks_new)
    if inps.rlks_new:
        util.update_yaml_key('params.yaml', 'rlks', inps.rlks_new)
    util.update_yaml_key('params.yaml', 'doCropSlc', 'False')
    
if __name__ == '__main__':
    '''
    Main driver.
    '''
    inps = cmdLineParser()
    main(inps)