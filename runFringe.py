#!/usr/bin/env python3
'''
Wrapper for running Fringe in SARTS workflow

KM

'''
import os,glob,argparse
import numpy as np
from matplotlib import pyplot as plt
from osgeo import gdal
import tops2vrt,nmap,sequential_PL,adjustMiniStacks,ampdispersion
from SARTS import config,util


ds = gdal.Open('Fringe/coreg_stack/slcs_base.vrt')
stack = ds.GetVirtualMemArray()

plt.figure();plt.plot( np.angle(stack[:,10,10]))
np.where(np.angle(stack[:,10,10])==0)

def cmdLineParser():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser(
        description='Run Fringe sequential estimator with phase linking',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('-t','--tops2vrt', action='store_true', help='Run tops2vrt function.')
    parser.add_argument('-n','--nmap', action='store_true', help='Run nmap function.')
    parser.add_argument('-s','--sequential_PL', action='store_true', help='Run sequential_PL function.')
    parser.add_argument('-a','--adjustMiniStacks', action='store_true', help='Run adjustMiniStacks function.')
    parser.add_argument('-d','--ampdispersion', action='store_true', help='Run ampdispersion function.')
    parser.add_argument('-c','--makeTcorr', action='store_true', help='Make Tcorr average')
    
    return parser.parse_args()


def run_sequential_PL(ps):
    sequential_PL.main(ps)
    ds_SLCS = glob.glob(ps.outDir + '/*slc')
    for fn_slc in ds_SLCS:
        util.write_xml(fn_slc,ps.nx,ps.ny,1,dataType='CFLOAT',scheme='BIP')
    ds_tcorrs= glob.glob(ps.outDir + '/*bin')
    for fn_slc in ds_tcorrs:
        util.write_xml(fn_slc,ps.nx,ps.ny,1,dataType='CFLOAT',scheme='BIP')
    ds_SLCS = glob.glob('./Fringe/Sequential/Datum_connection/EVD/*slc')
    for fn_slc in ds_SLCS:
        util.write_xml(fn_slc,ps.nx,ps.ny,1,dataType='CFLOAT',scheme='BIP')

def run_adjustMiniStacks(ps):
    adjustMiniStacks.main(ps)
    slcFns = glob.glob( ps.outDir + '/*slc')
    for  f in slcFns:
        util.write_xml(f,ps.nx,ps.ny,1,'CFLOAT','BSQ')
        
def run_ampdispersion(ps):
    ampdispersion.main(ps)
    util.write_xml(ps.outputAD,ps.nx,ps.ny,1,'FLOAT','BSQ')
    util.write_xml(ps.meanampDS,ps.nx,ps.ny,1,'FLOAT','BSQ')
    
    #Output the ps pixels by using a threshold on ampdispersion
    os.system('imageMath.py -e="a<' + str(ps.ampDispersionThreshold) + '" --a=./Fringe/ampDispersion/ampdispersion  -o ./Fringe/ampDispersion/ps_pixels -t byte')


def makeTcorrMean(ps):
    #make tcorrMean.bin_____________________________________________________________
    miniStacks_tcorr_files = glob.glob('./Fringe/Sequential/miniStacks/*/EVD/tcorr.bin')
    tcorrsMean = np.zeros((int(ps.ny), int(ps.nx)),dtype='float32')
    # ii=0;start = ii * chunk_size ;end = start+chunk_size + 1
    
    for ii in range(int(ps.ny)):
        tc = np.zeros((len(miniStacks_tcorr_files),int(ps.nx)))
        for jj,mtf in enumerate(miniStacks_tcorr_files):
            ds = gdal.Open(mtf)
            tc[jj,:] = ds.GetVirtualMemArray()[ii,:]
        tcorrsMean[ii,:] = np.nanmean(tc,axis=0)
        
    plt.figure();plt.imshow(tcorrsMean,vmin=0,vmax=1,cmap='magma');plt.show()

    # if ps.sensor=='ALOS':
    #     ds = gdal.Open('Fringe/PhaseLink/tcorr.bin')
    #     tcorrsMean = ds.GetVirtualMemArray()
    #     plt.figure();plt.imshow(tcorrsMean,vmin=0,vmax=1,cmap='magma');plt.show()

    # Write the average tcorr file 
    fmt = "GTiff"
    driver = gdal.GetDriverByName(fmt)
    [cols, rows] = tcorrsMean.shape
    outDataRaster = driver.Create("./Fringe/tcorrMean.bin", rows, cols, 1, gdal.GDT_Float32)
    bnd = outDataRaster.GetRasterBand(1)
    bnd.WriteArray(tcorrsMean)
    bnd.FlushCache()

def main(flags):
    
    ps = config.getPS()   
    ps.indir      = ps.mergeddir
    ps.bbox       = np.array([ps.cropymin,ps.cropymax,ps.cropxmin,ps.cropxmax])
    ps.memorySize = ps.maxMem
    ps.inputDir   = ps.slcdir
    ps.weightDS   = ps.outputDS
    ps.nx = str(ps.nx)
    ps.ny = str(ps.ny)
    ps.geobbox = None
    
    if flags.tops2vrt:
        tops2vrt.main(ps)

    if flags.nmap:
        nmap.main(ps)

    if flags.sequential_PL:
        run_sequential_PL(ps)
        
    if flags.adjustMiniStacks:
        run_adjustMiniStacks(ps)

    if flags.ampdispersion:
        run_ampdispersion(ps)
    
    if flags.makeTcorr:
        makeTcorrMean(ps)

    # Now use the script ifgs.py to make unwrapped ifgs
if __name__ == '__main__':
    '''
    Main driver.
    '''
    
    flags = cmdLineParser()

    # flags = argparse.Namespace()
    # flags.tops2vrt = True
    # flags.nmap = True
    # flags.sequential_PL = True
    # flags.adjustministacks = True
    # flags.ampdispersion = True
    # flags.makeTcorr = True

    main(flags)
