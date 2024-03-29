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
from SARTS import unwrap, config, filter_ifg
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
    parser.add_argument('-q', '--filt_full', action='store_true', dest='filt_full', help='Filter the full res ifg')

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
    if inps.filt_full:
        print('using filtered int')
        ps.infile       = os.path.join(pairDir, f"{pair}_filt.int")
    else:
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
  

def unwrapsnaphu_test(pair,ps):  
    '''
    unwrap ifgs
    '''
    pairDir = os.path.join(ps.intdir, pair )

    if not os.path.isfile( os.path.join( pairDir, 'filt_lk.unw')):
        print(f"Unwrapping {pair}")
        cor_file =  os.path.join( pairDir, 'filt_lk.cor')
        int_file =  os.path.join( pairDir, 'filt_lk.int')
        unw_file =  os.path.join( pairDir, 'filt_lk2.unw')        
        unwrap.unwrap_snaphu(int_file, cor_file, unw_file, ps)
    else:
        print(f"{pair} is already unwrapped.")
          
# unwrapsnaphu_test(ps.pairs[0],ps)
        

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
    # ds1 = gdal.Open(slc1_fn)
    # ds2 = gdal.Open(slc2_fn)
    # slc1 = ds1.GetVirtualMemArray()
    # slc2 = ds2.GetVirtualMemArray()
    # ifg_raw = np.multiply(slc1,np.conj(slc2))
    # del(slc1,slc2)
    
    ds1 = gdal.Open(ds_slc1_fn)
    ds2 = gdal.Open(ds_slc2_fn)
    ds_slc1 = ds1.GetVirtualMemArray()
    ds_slc2 = ds2.GetVirtualMemArray()
    ifg_ds_ps = np.multiply(ds_slc1,np.conj(ds_slc2))
    # del(ds_slc1,ds_slc2)
    
    # ds = gdal.Open(ps_mask_fn)
    # psPixels = ds.GetVirtualMemArray()

    # get the data for PS pixels
    # ifg_ds_ps[psPixels == 1] = ifg_raw[psPixels == 1]

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


def link_full_ifgs(dir_name_target,dir_name_source,base_dir,newpairs):
    '''
    This is for when you want to make a new network type, and link the available
    full IFGs from a previous network rather than making them all again.
    
    example:
        dolphinDir = os.path.join(ps.workdir, ps.dolphin_work_dir)     
        dir_name_target = 'interferograms'
        dir_name_source = 'interferograms_seq3'
        link_full_ifgs(dir_name_target,dir_name_source,dolphinDir,ps.pairs)
    
    '''
    
    dir_target = os.path.join(base_dir,dir_name_target)
    dir_source = os.path.join(base_dir,dir_name_source)
    
    exts = ['.int','.int.xml','.int.vrt']
    
    if not os.path.isdir(dir_target):
        os.mkdir(dir_target)
    
    for p in newpairs:
        pairdir_source = os.path.join(dir_source,p)
        pairdir_target = os.path.join(dir_target,p)
        if os.path.isdir(pairdir_source):
            if not os.path.isdir(pairdir_target):
                os.mkdir(pairdir_target)
            for ext in exts:
                dest =os.path.join(dir_target,p,p+ext)
                source_rel = os.path.join('../..',dir_name_source,p,p+ext)
        
                if not os.path.isdir(dest):
                    os.symlink(source_rel,dest)
                else:
                    print(f'{pairdir_source} Already exists in {dir_name_source}')
        else:
            print(f'{pairdir_source} Does not exist.')


def mask_ifgs(ps,msk,file_name='filt_lk.int'):
    # Loop through filt_lk.int and mask out water
    ds = gdal.Open('merged/geom_reference/waterMask_lk.rdr.vrt')
    msk = ds.GetVirtualMemArray()
    
    
    for pair in ps.pairs:
        pairDir         =  os.path.join(ps.intdir,pair)
        f_out   = os.path.join(pairDir,file_name)
        intImage = isceobj.createIntImage()
        intImage.load(f_out + '.xml')
        ifg = intImage.memMap()[:,:,0]
        ifg = ifg.copy()
        ifg[msk==0] = 0
        
        fidc=open(f_out,"wb")
        fidc.write(ifg)
    
        intImage.dump(f_out + '.xml') # Write out xml
        intImage.renderHdr()
        intImage.renderVRT()

def filt_full_ifgs(args):
    
    pair, ps = args
    pairDir         =  os.path.join(ps.intdir,pair)
    f_in   = os.path.join(pairDir,pair +'.int')
    f_out   = os.path.join(pairDir,pair +'_filt.int')
    if not os.path.isfile(f_out):
        intImage = isceobj.createIntImage()
        intImage.load(f_in + '.xml')
    
        ifg = intImage.memMap()[:,:,0]
        ifg = ifg.copy()
        real = np.real(ifg)
        imag = np.imag(ifg)
        realf = filter_ifg.gamma_map_filter(real, kernel_size=7, cu=0.4)
        imagf = filter_ifg.gamma_map_filter(imag, kernel_size=7, cu=0.4)
        ifgf    = realf + 1j * imagf
        
        intOut = intImage.clone()
        intOut.filename = f_out
        fidc=open(f_out,"wb")
        fidc.write(ifgf)
    
        intOut.dump(f_out + '.xml') # Write out xml
        intOut.renderHdr()
        intOut.renderVRT()
    
def make_ifgs_dolphin(args):

    pair, ps = args
    # integrate PS and DS
    out_fn = os.path.join(ps.dolphinDir,'interferograms',pair,pair + '.int')
    if not os.path.isfile(out_fn):
        print('making ' + pair + ' PSDS ifg...')
        d1,d2 = pair.split('_')
        if ps.crop:
            slc1_fn = os.path.join('merged','SLC',d1,d1+'.slc.full.crop.vrt')
            slc2_fn = os.path.join('merged','SLC',d2,d2+'.slc.full.crop.vrt')
        else:
            slc1_fn = os.path.join('merged','SLC',d1,d1+'.slc.full.vrt')
            slc2_fn = os.path.join('merged','SLC',d2,d2+'.slc.full.vrt')

        ds_slc1_fn = os.path.join(ps.dolphinDir,'linked_phase',d1+'.slc.tif')
        ds_slc2_fn = os.path.join(ps.dolphinDir,'linked_phase',d2+'.slc.tif')
        ps_mask_fn = os.path.join(ps.dolphinDir,'PS','ps_pixels.tif')
        if not os.path.isdir(os.path.join(ps.dolphinDir,'interferograms',pair)):
            os.makedirs(os.path.join(ps.dolphinDir,'interferograms',pair))
        makePSDS(ds_slc1_fn, ds_slc2_fn, slc1_fn, slc2_fn, ps_mask_fn, out_fn)
    else:
        print(pair + ' PSDS already exists.. skipping')
    # Filter full res if specified

            
def make_ifgs(args):
    pair, ps = args
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

def main(inps):
    ps = config.getPS()
    ps.azlooks      = int(ps.alks)
    ps.rglooks      = int(ps.rlks)
    ps.coregSlcDir    = './merged/SLC'
    ps.unwrapMethod   = None
    
    ps.dolphinDir = os.path.join(ps.workdir, ps.dolphin_work_dir)

    if inps.nodolphin:
        ps.intdir  = os.path.join( ps.mergeddir, 'interferograms')
    else:
        ps.intdir  = os.path.join( ps.dolphinDir, 'interferograms')
        

    if ps.sat=='SENTINEL-1':
        print('ALOS data')
        pix_spacing_range = round(2.3 * ps.rlks,2)
        pix_spacing_az = round(14.1 * ps.alks,2)
    else:
        print('ALOS data')
        pix_spacing_range = round(9.7 * ps.rlks,2)
        pix_spacing_az = round(3.1 * ps.alks,2)
        
    print(f'Az looks is set to {ps.alks}')
    print(f'Az pixel spacing will be ~{pix_spacing_az} m')
    print(f'Rg looks is set to {ps.rlks}')
    print(f'Rg pixel spacing will be ~{pix_spacing_range} m')

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
        #__Make IFGS____________________

        if inps.makeIfgs:
            if not inps.nodolphin:
                print('Using PS_DS slcs')
                pool = multiprocessing.Pool(processes=inps.num_processes)
                pool.map(make_ifgs_dolphin, args_list)
                pool.close()
                pool.join()
                
            else:
                print('Using original slcs because nodolphin flag is set')
                print(ps.intdir)
                if not os.path.isdir(ps.intdir):
                    print('Making merged/interferograms directory')
                    os.makedirs(ps.intdir)
                pool = multiprocessing.Pool(processes=inps.num_processes)
                pool.map(make_ifgs, args_list)
                pool.close()
                pool.join() 
                
        if inps.filt_full:
            print('filtering full res ifgs')
            pool = multiprocessing.Pool(processes=inps.num_processes)
            pool.map(filt_full_ifgs, args_list)
            pool.close()
            pool.join()  
                
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
                downlook((pair,ps))
        if inps.unwrap:
            for pair in ps.pairs:
                unwrapsnaphu((pair,ps))


if __name__ == '__main__':
    '''
    Main driver.
    '''
    # inps = argparse.Namespace()
    
    inps = cmdLineParser()
    # inps.makeIfgs = True
    main(inps)