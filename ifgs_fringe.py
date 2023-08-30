#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 17:05:26 2023

Takes output from Fringe
creates:

    <date1>_<date2>.int
    fine_lk.int
    filt_lk.int
    filt_lk.cor
    filt_lk.unw

@author: km
"""

import numpy as np
import os
from isce.applications import looks
import FilterAndCoherence
import integratePS
import multiprocessing
from SARTS import unwrap
import argparse


def cmdLineParser():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser(
        description='Crop and downlook geom files. Save parameters',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--downlook', type=bool, dest='dlunw', default=True)
    parser.add_argument('-m', '--make-ifgs', type=bool, dest='makeIfgs', default=True)
    parser.add_argument('-n', '--nproc', type=int, dest='num_processes', default=4)
    parser.add_argument('-p', '--parallel', type=bool, dest='parallelProcess', default=True)

    return parser.parse_args()


def main(inps):
    
    ps = np.load('./ps.npy',allow_pickle=True).all()
    fringeDir = './Fringe/'
    
    ps.azlooks      = int(ps.alks)
    ps.rglooks      = int(ps.rlks)
    
    ps.intdir  = fringeDir + 'PS_DS/' + ps.networkType
    if ps.sensor=='ALOS':
        ps.slcdir  = fringeDir + 'PhaseLink'
        ps.dsStackDir = fringeDir + 'PhaseLink'
    
    else:
        ps.slcdir  = fringeDir + 'adjusted_wrapped_DS'
        ps.dsStackDir = fringeDir + 'adjusted_wrapped_DS'
    
    
    ps.slcStack       = fringeDir + 'coreg_stack/slcs_base.vrt'
    ps.tcorrFile      = fringeDir + 'tcorrMean.bin'
    ps.psPixelsFile   = fringeDir + 'ampDispersion/ps_pixels'
    ps.outDir         = fringeDir + 'PS_DS/' + ps.networkType
    ps.coregSlcDir    = './merged/SLC'
    ps.pairs          = ps.pairs
    ps.unwrapMethod   = None
    #______________________
    if inps.makeIfgs:
        integratePS.main(ps)
    #______________________
    
    
    
    def downlook(pair):
        # Downlook ifgs
        pairDir         = os.path.join(ps.outDir, pair )
        ps.infile       = os.path.join(pairDir, f"{pair}.int")
        ps.outfile      = os.path.join(pairDir, 'fine_lk.int')
        cor_file_out    = os.path.join(pairDir, 'filt_lk.cor')
        filt_file_out   = os.path.join(pairDir, 'filt_lk.int')
    
        if not os.path.isfile(ps.outfile):
            print(f"Downlooking {pair}")
            looks.main(ps)
            
            # Filter and coherence
            if not os.path.isfile(filt_file_out):
                print(f"Filtering and computing coherence for {pair}")
                FilterAndCoherence.runFilter(ps.outfile,filt_file_out,.4)
                FilterAndCoherence.estCoherence(filt_file_out, cor_file_out)
        else:
            print(pair + '/' + filt_file_out + ' is already file.')
                  
    
    def unwrapsnaphu(pair):  
        pairDir =  ps.outDir + '/' + pair 
        if not os.path.isfile(pairDir + '/filt_lk.unw'):
            print(f"Unwrapping {pair}")
            cor_file = os.path.join(pairDir, 'filt_lk.cor')
            int_file = os.path.join(pairDir, 'filt_lk.int')
            unw_file = os.path.join(pairDir, 'filt_lk.unw')        
            unwrap.unwrap_snaphu(int_file, cor_file, unw_file, ps)
        else:
            print(f"{pair} is already unwrapped.")
    
    
    if inps.dlunw:
        
        if inps.parallelProcess:
            pool = multiprocessing.Pool(processes=inps.num_processes)
            pool.map(downlook, ps.pairs)
            pool.close()
            pool.join()
            
            pool = multiprocessing.Pool(processes=inps.num_processes)
            pool.map(unwrapsnaphu, ps.pairs)
            pool.close()
            pool.join()
            
        else:
            for pair in ps.pairs:
                downlook(pair)
            for pair in ps.pairs:
                unwrapsnaphu(pair)
                
if __name__ == '__main__':
    '''
    Main driver.
    '''
    inps = cmdLineParser()
    main(inps)
