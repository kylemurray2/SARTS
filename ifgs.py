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


# Mask filt_lk.int files with temp coherence
# import h5py
# # Temporal Coh
# filename = 'MintPy_dol_seq3_filt_1_3/temporalCoherence.h5'
# ds = h5py.File(filename,'r+')   
# temporalCoherence = np.asarray(ds['temporalCoherence'])
# ds.close()

# ps = config.getPS()
# dolphinDir = os.path.join(ps.workdir, ps.dolphin_work_dir)
# ps.intdir  = os.path.join( dolphinDir, 'interferograms')

# for pair in ps.pairs:
#     pairDir         = os.path.join(ps.intdir, pair )
#     filt_file_out   = os.path.join(pairDir, 'filt_lk.int')
#     intImage = isceobj.createIntImage()
#     intImage.load(filt_file_out + '.xml')
#     ifg = intImage.memMap()[:,:,0]
#     ifg = ifg.copy()
#     ifg[temporalCoherence<.3] = 0
    
#     fidc=open(filt_file_out,"wb")
#     fidc.write(ifg)

    
#     intImage.dump(filt_file_out + '.xml') # Write out xml
#     intImage.renderHdr()
#     intImage.renderVRT()
    

def downlook(args):
    '''
    downlook ifgs
    '''
    pair, ps = args
    pairDir         = os.path.join(ps.intdir, pair )
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
        
    # # Downlook the cor file
    # ps.infile       = os.path.join(pairDir, f"{pair}.cor")
    # ps.outfile      = os.path.join(pairDir, 'filt_lk.cor')
    # if not os.path.isfile(ps.outfile):
    #     print(f"Downlooking {pair}")
    #     looks.main(ps)
        
    # coherence
    if not os.path.isfile(cor_file_out):
        print(f"Computing coherence for {pair}")
        FilterAndCoherence.estCoherence(filt_file_out, cor_file_out)
    else:
        print(pair + '/' + filt_file_out + ' is already file.')
              

def unwrapsnaphu(args):  
    '''
    unwrap ifgs
    '''
    pair, ps = args
    pairDir = os.path.join(ps.intdir, pair )

    if not os.path.isfile( os.path.join( pairDir, 'filt_lk.unw')):
        print(f"Unwrapping {pair}")
        cor_file =  os.path.join( pairDir, 'filt_lk.cor')
        int_file =  os.path.join( pairDir, 'filt_lk.int')
        unw_file =  os.path.join( pairDir, 'filt_lk.unw')        
        unwrap.unwrap_snaphu(int_file, cor_file, unw_file, ps)
    else:
        print(f"{pair} is already unwrapped.")


def makePSDS(ds_slc1_fn, ds_slc2_fn, slc1_fn, slc2_fn, ps_mask_fn, out_fn):
    ds1 = gdal.Open(slc1_fn)
    ds2 = gdal.Open(slc2_fn)
    slc1 = ds1.GetVirtualMemArray()
    slc2 = ds2.GetVirtualMemArray()
    ifg_raw = np.multiply(slc1,np.conj(slc2))
    # del(slc1,slc2)
    
    ds1 = gdal.Open(ds_slc1_fn)
    ds2 = gdal.Open(ds_slc2_fn)
    ds_slc1 = ds1.GetVirtualMemArray()
    ds_slc2 = ds2.GetVirtualMemArray()
    ifg_ds_ps = np.multiply(ds_slc1,np.conj(ds_slc2))
    # del(ds_slc1,ds_slc2)
    
    ds = gdal.Open(ps_mask_fn)
    psPixels = ds.GetVirtualMemArray()

    # get the data for PS pixels
    ifg_ds_ps[psPixels == 1] = ifg_raw[psPixels == 1]

    out = isceobj.createIntImage() # Copy the interferogram image from before
    out.dataType = 'CFLOAT'
    out.filename = out_fn
    out.width = ifg_ds_ps.shape[1]
    out.length = ifg_ds_ps.shape[0]
    out.dump(out.filename + '.xml') # Write out xml
    fid=open(out.filename,"wb+")
    fid.write(ifg_ds_ps)
    out.renderHdr()
    out.renderVRT()  
    fid.close()

    return None


def main(inps):
    ps = config.getPS()
    ps.azlooks      = int(ps.alks)
    ps.rglooks      = int(ps.rlks)
    ps.coregSlcDir    = './merged/SLC'
    ps.unwrapMethod   = None
    
    dolphinDir = os.path.join(ps.workdir, ps.dolphin_work_dir)

    if inps.nodolphin:
        ps.intdir  = os.path.join( ps.mergeddir, 'interferograms')
    else:
        ps.intdir  = os.path.join( dolphinDir, 'interferograms')
        

    pix_spacing_range = round(2.3 * ps.rlks,2)
    pix_spacing_az = round(14.1 * ps.alks,2)
    print(f'Az looks is set to {ps.alks}')
    print(f'Az pixel spacing will be ~{pix_spacing_az} m')
    print(f'Rg looks is set to {ps.rlks}')
    print(f'Rg pixel spacing will be ~{pix_spacing_range} m')


    #__Make IFGS____________________
    if not inps.nodolphin:
        # Make ifgs with dolphin
        print('Using PS_DS slcs')
        if inps.makeIfgs:
            # integrate PS and DS
            for pair in ps.pairs:
                out_fn = os.path.join(dolphinDir,'interferograms',pair,pair + '.int')

                if not os.path.isfile(out_fn):
                    print('making ' + pair + ' PSDS ifg...')
                    d1,d2 = pair.split('_')
                    if ps.crop:
                        slc1_fn = os.path.join('merged','SLC',d1,d1+'.slc.full.crop.vrt')
                        slc2_fn = os.path.join('merged','SLC',d2,d2+'.slc.full.crop.vrt')
                    else:
                        slc1_fn = os.path.join('merged','SLC',d1,d1+'.slc.full.vrt')
                        slc2_fn = os.path.join('merged','SLC',d2,d2+'.slc.full.vrt')

                    ds_slc1_fn = os.path.join(dolphinDir,'linked_phase',d1+'.slc.tif')
                    ds_slc2_fn = os.path.join(dolphinDir,'linked_phase',d2+'.slc.tif')
                    ps_mask_fn = os.path.join(dolphinDir,'PS','ps_pixels.tif')
                    if not os.path.isdir(os.path.join(dolphinDir,'interferograms',pair)):
                        os.makedirs(os.path.join(dolphinDir,'interferograms',pair))
                    makePSDS(ds_slc1_fn, ds_slc2_fn, slc1_fn, slc2_fn, ps_mask_fn, out_fn)
                else:
                    print(pair + ' PSDS already exists.. skipping')
    else:
        # Make ifgs without dolphin PSDS
        print('Using original slcs because nodolphin flag is set')
        print(ps.intdir)
        if inps.makeIfgs:

            if not os.path.isdir(ps.intdir):
                print('Making merged/interferograms directory')
                os.makedirs(ps.intdir)
    
            for pair in ps.pairs:
                d1 = pair.split('_')[0]
                d2 = pair.split('_')[1]
                if ps.crop:
                    slc1_fn = os.path.join(ps.slcdir,d1,d1+'.slc.full.crop.vrt')
                    slc2_fn = os.path.join(ps.slcdir,d2,d2+'.slc.full.crop.vrt')
                else:
                    slc1_fn = os.path.join(ps.slcdir,d1,d1+'.slc.full.vrt')
                    slc2_fn = os.path.join(ps.slcdir,d2,d2+'.slc.full.vrt') 

                pair = d1 + '_' + d2
                print('creating ' + pair + '.int')
                ifg_fn = os.path.join(ps.intdir,pair,pair+'.int')
                if not os.path.isfile(ifg_fn):
                    if not os.path.isdir(os.path.join(ps.intdir,pair)):
                        os.makedirs(os.path.join(ps.intdir,pair))
                    makeIfg(slc1_fn,slc2_fn,ifg_fn,ps)
                else:
                    print(ifg_fn + ' already exists')
    #______________________
    # #organize the ifgs into pairs directories inside the ps.networkType dir
    # if not os.path.isdir( os.path.join(ps.intdir,ps.networkType)):
    #     os.mkdir(os.path.join(ps.intdir,ps.networkType))
    # for pair in ps.pairs:
    #     if not os.path.isdir( os.path.join(ps.intdir,ps.networkType,pair) ):
    #         os.mkdir(os.path.join(ps.intdir,ps.networkType,pair))
        
    #     suffixList = ['.int','.int.aux.xml','.int.hdr','.int.vrt','.cor','.cor.aux.xml','.cor.hdr']
    #     if not os.path.isfile( os.path.join(ps.intdir,ps.networkType,pair,pair + '.cor')):
    #         for suff in suffixList:
    #             print('ln -s ../../' + pair + suff +' ' + ps.intdir + '/' + ps.networkType + '/' + pair + '/'  )
    #             os.system('ln -s ../../' + pair + suff +' ' + ps.intdir + '/' + ps.networkType + '/' + pair + '/' + pair + suff   )

            
    #     ifg_fn = ps.intdir + '/' + ps.networkType + '/' + pair + '/' + pair + '.int'
    #     cor_fn = ps.intdir + '/' + ps.networkType + '/' + pair + '/' + pair + '.cor'

    #     os.system('gdal2isce_xml.py -i ' + ifg_fn )
    #     os.system('gdal2isce_xml.py -i ' + cor_fn )
    
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
            
            
    # # link files to sequential1 if that is not the chosen network type
    # if not inps.nodolphin:
    #     if ps.networkType != 'sequential1':
    #         seq1Dir = os.path.join(dolphinDir,'PS_DS','sequential1')
    #         if not os.path.isdir(seq1Dir):
    #             os.mkdir(seq1Dir)
    #         for p in ps.pairs_seq:
    #             pairdir = os.path.join(ps.outDir,p)
    #             if os.path.isdir(pairdir):
    #                 dest =os.path.join(seq1Dir,p)
    #                 # dest_abs = os.path.abspath(os.path.join(seq1Dir,p))
    #                 # source =os.path.join(ps.outDir,p)
    #                 source_rel = os.path.join('..',ps.networkType,p)

    #                 if not os.path.isdir(dest):
    #                     os.symlink(source_rel,dest)
    #                 else:
    #                     print(pairdir + ' Already exists in sequential1')
    #             else:
    #                 print(pairdir + ' Does not exist. sequential network is disconnected')


if __name__ == '__main__':
    '''
    Main driver.
    '''
    # inps = argparse.Namespace()
    
    inps = cmdLineParser()
    # inps.makeIfgs = True
    main(inps)