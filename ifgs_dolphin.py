#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 17:05:26 2023
@author: km

Takes output from dolphin
creates:

    <date1>_<date2>.int
    fine_lk.int
    filt_lk.int
    filt_lk.cor
    filt_lk.unw

"""

import os,glob
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
    parser.add_argument('-f', '--nodolphin', action='store_true', dest='nodolphin', help='Use this flag if you are not using dolphin psds.')

    return parser.parse_args()


# Make ifgs from normal SLCs
def makeIfg(slc1_fn,slc2_fn,ifg_fn,ps):
    ds1 = gdal.Open(slc1_fn)
    ds2 = gdal.Open(slc2_fn)
    slc1 = ds1.GetVirtualMemArray()
    slc2 = ds2.GetVirtualMemArray()
    ifg = np.multiply(slc2,np.conj(slc1))
    
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
        
    # Downlook the cor file
    ps.infile       = os.path.join(pairDir, f"{pair}.cor")
    ps.outfile      = os.path.join(pairDir, 'filt_lk.cor')
    if not os.path.isfile(ps.outfile):
        print(f"Downlooking {pair}")
        looks.main(ps)
        
    # coherence
    # if not os.path.isfile(cor_file_out):
    #     print(f"Computing coherence for {pair}")
    #     FilterAndCoherence.estCoherence(filt_file_out, cor_file_out)
    # else:
    #     print(pair + '/' + filt_file_out + ' is already file.')
              

def unwrapsnaphu(args):  
    pair, ps = args
    pairDir =  ps.outDir + '/' + pair 
    if not os.path.isfile( os.path.join( pairDir, pair + '.unw')):
        print(f"Unwrapping {pair}")
        cor_file =  os.path.join( pairDir, 'filt_lk.cor')
        int_file =  os.path.join( pairDir, 'fine_lk.int')
        unw_file =  os.path.join( pairDir, 'filt_lk.unw')        
        unwrap.unwrap_snaphu(int_file, cor_file, unw_file, ps)
    else:
        print(f"{pair} is already unwrapped.")


def main(inps):
    ps = config.getPS()
    
    ps.azlooks      = int(ps.alks)
    ps.rglooks      = int(ps.rlks)
    ps.coregSlcDir    = './merged/SLC'
    ps.unwrapMethod   = None
    

    dolphinDir = ps.slcStackVrt.split('/')[1]
    ps.intdir  = os.path.join( dolphinDir, 'interferograms')

    ps.slcdir  =  os.path.join(dolphinDir, 'linked_phase')
        
    ps.tcorrFile      =  glob.glob(ps.slcdir + '/temporal_coherence*tif')[0]
    ps.psPixelsFile   =  os.path.join(dolphinDir, 'PS/ps_pixels.tif')
    ps.outDir         =  os.path.join(dolphinDir, 'interferograms', ps.networkType)

    
    
    #organize the ifgs into pairs directories inside the ps.networkType dir
    if not os.path.isdir( os.path.join(ps.intdir,ps.networkType)):
        os.mkdir(os.path.join(ps.intdir,ps.networkType))
    for pair in ps.pairs:
        if not os.path.isdir( os.path.join(ps.intdir,ps.networkType,pair) ):
            os.mkdir(os.path.join(ps.intdir,ps.networkType,pair))
            
            # if not os.path.isfile( os.path.join(ps.intdir,ps.networkType,pair,pair + '.cor')):
            os.system('mv ' + ps.intdir + '/' + pair + '.* ' + ps.intdir + '/' + ps.networkType + '/' + pair + '/'  )
            
        ifg_fn = ps.intdir + '/' + ps.networkType + '/' + pair + '/' + pair + '.int'
        cor_fn = ps.intdir + '/' + ps.networkType + '/' + pair + '/' + pair + '.cor'

        os.system('gdal2isce_xml.py -i ' + ifg_fn )
        os.system('gdal2isce_xml.py -i ' + cor_fn )
    
    args_list = [(pair, ps) for pair in ps.pairs]
       
    if inps.num_processes>1:
        if inps.downlook:
            pool = multiprocessing.Pool(processes=4)#inps.num_processes)
            pool.map(downlook, args_list)
            pool.close()
            pool.join()
        
        if inps.unwrap:
            pool = multiprocessing.Pool(processes=5)#inps.num_processes)
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
            
            
    # link files to sequential1 if that is not the chosen network type
    if not inps.nodolphin:
        if ps.networkType != 'sequential1':
            seq1Dir = os.path.join(dolphinDir,'PS_DS','sequential1')
            if not os.path.isdir(seq1Dir):
                os.mkdir(seq1Dir)
            for p in ps.pairs_seq:
                pairdir = os.path.join(ps.outDir,p)
                if os.path.isdir(pairdir):
                    dest =os.path.join(seq1Dir,p)
                    # dest_abs = os.path.abspath(os.path.join(seq1Dir,p))
                    # source =os.path.join(ps.outDir,p)
                    source_rel = os.path.join('..',ps.networkType,p)

                    if not os.path.isdir(dest):
                        os.symlink(source_rel,dest)
                    else:
                        print(pairdir + ' Already exists in sequential1')
                else:
                    print(pairdir + ' Does not exist. sequential network is disconnected')


if __name__ == '__main__':
    '''
    Main driver.
    '''
    inps = argparse.Namespace()
    
    # inps = cmdLineParser()
    # # inps.makeIfgs = True
    # main(inps)
    
    
