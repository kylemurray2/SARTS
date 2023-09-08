#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 17:05:26 2023
@author: km

Takes output from Fringe
creates:

    <date1>_<date2>.int
    fine_lk.int
    filt_lk.int
    filt_lk.cor
    filt_lk.unw

"""

import os
from isce.applications import looks
import FilterAndCoherence
import integratePS
import multiprocessing
from SARTS import unwrap, config
import argparse


def cmdLineParser():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser(
        description='Crop and downlook geom files. Save parameters',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--downlook', action='store_true', dest='downlook', default=True,help='Downlook interferograms')
    parser.add_argument('-d', '--unwrap', action='store_true', dest='unwrap', default=True,help='Unwrap interferograms')
    parser.add_argument('-m', '--make-ifgs', action='store_true', dest='makeIfgs', default=True,help='Make the interferograms')
    parser.add_argument('-n', '--nproc', type=int, dest='num_processes', default=5, help='Number of parallel processes. Use 1 for no parallelization')

    return parser.parse_args()


def downlook(pair,ps):
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
              

def unwrapsnaphu(pair,ps):  
    pairDir =  ps.outDir + '/' + pair 
    if not os.path.isfile(pairDir + '/filt_lk.unw'):
        print(f"Unwrapping {pair}")
        cor_file = os.path.join(pairDir, 'filt_lk.cor')
        int_file = os.path.join(pairDir, 'filt_lk.int')
        unw_file = os.path.join(pairDir, 'filt_lk.unw')        
        unwrap.unwrap_snaphu(int_file, cor_file, unw_file, ps)
    else:
        print(f"{pair} is already unwrapped.")


def main(inps):
    ps = config.getPS()
    
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
    
    if inps.num_processes>1:
        if inps.downlook:
            pool = multiprocessing.Pool(processes=inps.num_processes)
            pool.map(downlook, ps.pairs,ps)
            pool.close()
            pool.join()
        
        if inps.unwrap:
            pool = multiprocessing.Pool(processes=inps.num_processes)
            pool.map(unwrapsnaphu, ps.pairs,ps)
            pool.close()
            pool.join()
        
    else:
        if inps.downlook:
            for pair in ps.pairs:
                downlook(pair,ps)
        if inps.unwrap:
            for pair in ps.pairs:
                unwrapsnaphu(pair,ps)
            
if __name__ == '__main__':
    '''
    Main driver.
    '''
    inps = cmdLineParser()
    main(inps)