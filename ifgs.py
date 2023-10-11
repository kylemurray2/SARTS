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
from isce.components import isceobj
import FilterAndCoherence
import integratePS
import multiprocessing
from SARTS import unwrap, config
import argparse
from osgeo import gdal
import numpy as np

def cmdLineParser():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser(
        description='Crop and downlook geom files. Save parameters',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--downlook', action='store_true', dest='downlook',help='Downlook interferograms')
    parser.add_argument('-u', '--unwrap', action='store_true', dest='unwrap',help='Unwrap interferograms')
    parser.add_argument('-m', '--make-ifgs', action='store_true', dest='makeIfgs',help='Make the interferograms')
    parser.add_argument('-n', '--nproc', type=int, dest='num_processes', default=3, help='Number of parallel processes. Use 1 for no parallelization')
    parser.add_argument('-f', '--noFringe', action='store_true', dest='noFringe', help='Use this flag if you are not using Fringe psds.')

    return parser.parse_args()


# Make ifgs from normal SLCs
def makeIfg(slc1_fn,slc2_fn,ifg_fn,ps):
    ds1 = gdal.Open(slc1_fn)
    ds2 = gdal.Open(slc2_fn)
    slc1 = ds1.GetVirtualMemArray()[ps.cropymin:ps.cropymax,ps.cropxmin:ps.cropxmax]
    slc2 = ds2.GetVirtualMemArray()[ps.cropymin:ps.cropymax,ps.cropxmin:ps.cropxmax]
    ifg = np.multiply(slc1,np.conj(slc2))
    
    out = isceobj.createIntImage() # Copy the interferogram image from before
    out.dataType = 'CFLOAT'
    out.filename = ifg_fn
    out.width = ifg.shape[1]
    out.length = ifg.shape[0]
    out.dump(out.filename + '.xml') # Write out xml
    fid=open(out.filename,"wb+")
    fid.write(ifg)
    out.renderHdr()
    out.renderVRT()  
    fid.close()
    
    return ifg


def downlook(args):
    pair, ps = args
    # Downlook ifgs
    pairDir         = os.path.join(ps.outDir, pair )
    ps.infile       = os.path.join(pairDir, f"{pair}.int")
    ps.outfile      = os.path.join(pairDir, 'fine_lk.int')
    cor_file_out    = os.path.join(pairDir, 'filt_lk.cor')
    filt_file_out   = os.path.join(pairDir, 'filt_lk.int')

    # Downlook
    if not os.path.isfile(ps.outfile):
        print(f"Downlooking {pair}")
        looks.main(ps)
    # Filter
    if not os.path.isfile(filt_file_out):
        print(f"Filtering {pair}")
        FilterAndCoherence.runFilter(ps.outfile,filt_file_out,ps.filterStrength)
    # coherence
    if not os.path.isfile(cor_file_out):
        print(f"Computing coherence for {pair}")
        FilterAndCoherence.estCoherence(filt_file_out, cor_file_out)
    # else:
    #     print(pair + '/' + filt_file_out + ' is already file.')
              

def unwrapsnaphu(args):  
    pair, ps = args
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
    
    ps.azlooks      = int(ps.alks)
    ps.rglooks      = int(ps.rlks)
    ps.coregSlcDir    = './merged/SLC'
    ps.unwrapMethod   = None
    
    if not inps.noFringe:
        fringeDir = ps.stackdir.split('/')[1]
        ps.intdir  = os.path.join( fringeDir, 'PS_DS', ps.networkType)
        if ps.sensor=='ALOS':
            ps.slcdir  =  os.path.join(fringeDir,'PhaseLink')
            ps.dsStackDir =  os.path.join(fringeDir, 'PhaseLink')
        
        else:
            ps.slcdir  =  os.path.join(fringeDir, 'adjusted_wrapped_DS')
            ps.dsStackDir =  os.path.join(fringeDir, 'adjusted_wrapped_DS')
        
        ps.slcStack       =  os.path.join(fringeDir, 'coreg_stack/slcs_base.vrt')
        ps.tcorrFile      =  os.path.join(fringeDir,'tcorrMean.bin')
        ps.psPixelsFile   =  os.path.join(fringeDir, 'ampDispersion/ps_pixels')
        ps.outDir         =  os.path.join(fringeDir, 'PS_DS', ps.networkType)
    else:
        print("Not using Fringe for IFG formation")
        ps.outDir     = ps.intdir
    
    #__Make IFGS____________________
    if not inps.noFringe:
        # Make ifgs with fringe
        if inps.makeIfgs:
            integratePS.main(ps)
    else:
        # Make ifgs without fringe
        if inps.makeIfgs:

            if not os.path.isdir(ps.intdir):
                print('Making merged/interferograms directory')
                os.mkdir(ps.intdir)
    
            for ii in range(len(ps.dates)-1):
                d1 = ps.dates[ii]
                d2 = ps.dates[ii+1]
                slc1_fn = os.path.join(ps.slcdir,d1,d1+'.slc.full')
                slc2_fn = os.path.join(ps.slcdir,d2,d2+'.slc.full')
                pair = d1 + '_' + d2
                print('creating ' + pair + '.int')
                ifg_fn = os.path.join(ps.intdir,pair,pair+'.int')
                if not os.path.isfile(ifg_fn):
                    if not os.path.isdir(os.path.join(ps.intdir,pair)):
                        os.mkdir(os.path.join(ps.intdir,pair))
                    makeIfg(slc1_fn,slc2_fn,ifg_fn)
                else:
                    print(ifg_fn + ' already exists')
    #______________________
    

    args_list = [(pair, ps) for pair in ps.pairs]


    if inps.num_processes>1:
        if inps.downlook:
            pool = multiprocessing.Pool(processes=inps.num_processes)
            pool.map(downlook, args_list)
            pool.close()
            pool.join()
        
        if inps.unwrap:
            pool = multiprocessing.Pool(processes=inps.num_processes)
            pool.map(unwrapsnaphu,args_list)
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
    # inps.makeIfgs = True
    main(inps)