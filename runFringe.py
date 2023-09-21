import os,glob,yaml,argparse
import numpy as np
from matplotlib import pyplot as plt
from osgeo import gdal
import tops2vrt,nmap,sequential_PL,adjustMiniStacks,ampdispersion
from SARTS import config,util

def cmdLineParser():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser(description='Run fringe sequential estimator.')
    parser.add_argument('-t','--tops2vrt', action='store_true', help='Run tops2vrt function.')
    parser.add_argument('-n','--nmap', action='store_true', help='Run nmap function.')
    parser.add_argument('-s','--sequential_PL', action='store_true', help='Run sequential_PL function.')
    parser.add_argument('-a','--adjustMiniStacks', action='store_true', help='Run adjustMiniStacks function.')
    parser.add_argument('-d','--ampdispersion', action='store_true', help='Run ampdispersion function.')
    
    return parser.parse_args()


def run_sequential_PL(inps):
    sequential_PL.main(inps)
    ds_SLCS = glob.glob(inps.outDir + '/*slc')
    for fn_slc in ds_SLCS:
        util.write_xml(fn_slc,inps.nx,inps.ny,1,dataType='CFLOAT',scheme='BIP')
    ds_tcorrs= glob.glob(inps.outDir + '/*bin')
    for fn_slc in ds_tcorrs:
        util.write_xml(fn_slc,inps.nx,inps.ny,1,dataType='CFLOAT',scheme='BIP')
    ds_SLCS = glob.glob('./Fringe/Sequential/Datum_connection/EVD/*slc')
    for fn_slc in ds_SLCS:
        util.write_xml(fn_slc,inps.nx,inps.ny,1,dataType='CFLOAT',scheme='BIP')

def run_adjustMiniStacks(inps):
    adjustMiniStacks.main(inps)
    slcFns = glob.glob( inps.outDir + '/*slc')
    for  f in slcFns:
        util.write_xml(f,inps.nx,inps.ny,1,'CFLOAT','BSQ')
        
def run_ampdispersion(inps):
    ampdispersion.main(inps)
    util.write_xml(inps.outputAD,inps.nx,inps.ny,1,'FLOAT','BSQ')
    util.write_xml(inps.meanampDS,inps.nx,inps.ny,1,'FLOAT','BSQ')
    
def main(inps):
    
    ps = config.getPS()   
    inps.indir      = ps.mergeddir
    inps.bbox       = np.array([ps.cropymin,ps.cropymax,ps.cropxmin,ps.cropxmax])
    inps.memorySize = ps.maxMem
    inps.inputDir   = ps.slcdir
    inps.weightDS   = inps.outputDS
    inps.nx = str(ps.nx)
    inps.ny = str(ps.ny)
    
    if inps.tops2vrt:
        tops2vrt.main(inps)

    if inps.nmap:
        nmap.main(inps)

    if inps.sequential_PL:
        run_sequential_PL(inps)
        
    if inps.adjustMiniStacks:
        run_adjustMiniStacks(inps)

    if inps.ampdispersion:
        run_ampdispersion(inps)
        
        
    # if ps.sensor=='ALOS':
    #     os.system('phase_link.py -i ./Fringe/coreg_stack/slcs_base.vrt -w ./Fringe/KS2/nmap -o ./Fringe/PhaseLink -r 70000 -x 11 -y 5')
    #     ds_SLCS = glob.glob('./Fringe/PhaseLink/*slc')
    #     for fn_slc in ds_SLCS:
    #         util.write_xml(fn_slc,nx,ny,1,dataType='CFLOAT',scheme='BIP')
    
    #Output the ps pixels by using a threshold on ampdispersion
    os.system('imageMath.py -e="a<' + int(ps.ampDispersionThreshold) + '" --a=./Fringe/ampDispersion/ampdispersion  -o ./Fringe/ampDispersion/ps_pixels -t byte')

    #make tcorrMean.bin_____________________________________________________________
    miniStacks_tcorr_files = glob.glob('./Fringe/Sequential/miniStacks/*/EVD/tcorr.bin')
    tcorrsMean = np.zeros((ps.ny, ps.nx),dtype='float32')
    # ii=0;start = ii * chunk_size ;end = start+chunk_size + 1
    
    for ii in range(ps.ny):
        tc = np.zeros((len(miniStacks_tcorr_files),ps.nx))
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


    # Now use the script ifgs.py to make unwrapped ifgs
if __name__ == '__main__':
    '''
    Main driver.
    '''
    # inps = argparse.Namespace()
    # inps.tops2vrt = True
    # inps.nmap = True
    # inps.sequential_PL = True
    # inps.adjustministacks = True
    # inps.ampdispersion = True

    inps = cmdLineParser()
    main(inps)
