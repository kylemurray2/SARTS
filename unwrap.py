''' modified from MintPy'''

from contrib.Snaphu.Snaphu import Snaphu
import time
import numpy as np
from SARTS import util


def unwrap_snaphu(int_file, cor_file, unw_file, ps):
    '''Unwrap interferograms using SNAPHU via isce2.

    Modified from ISCE-2/topsStack/unwrap.py
    Notes from Piyush:
        SNAPHU is an iterative solver, starting from the initial solution. It can get
            stuck in an infinite loop.
        The initial solution is created using MCF or MST method. The MST initial solution
            typically require lots of iterations and may not be a good starting point.
        DEFO cost mode requires geometry info for the program to interpret the coherence
            correctly and setup costs based on that. DEFO always sounds more theoretical
            to me, but I haven not fully explored it. TOPO cost mode requires spatial baseline.
            SMOOTH cost mode is purely data driven.
        Amplitude of zero is a mask in all cost modes. For TOPO mode, amplitude is used to find
            layover; for SMOOTH mode, only non-zero amplitude matters.

    Default configurations in ISCE-2/topsStack:
        init_only = True
        init_method = 'MCF'
        cost_mode = 'SMOOTH'
    Default configurations in FRInGE:
        init_only = False
        init_method = 'MST'
        cost_mode = 'DEFO'

    Parameters: int_file    - str, path to the wrapped interferogram file
                cor_file    - str, path to the correlation file
                unw_file    - str, path to the output unwrapped interferogram file
                defo_max    - float, maximum number of cycles for the deformation phase
                max_comp    - int, maximum number of connected components
                init_only   - bool, initlize-only mode
                init_method - str, algo used for initialization: MCF, MST
                cost_mode   - str, statistical-cost mode: TOPO, DEFO, SMOOTH, NOSTATCOSTS
    Returns:    unw_file    - str, path to the output unwrapped interferogram file
    '''

    init_only=True            
    start_time = time.time()
    width = ps.nxl
    length = ps.nyl
    
    # altitude = float(atr['HEIGHT'])
    # earth_radius = float(atr['EARTH_RADIUS'])
    # wavelength = float(atr['WAVELENGTH'])
    wavelength = ps.lam
    
    rg_looks = 1
    az_looks = 1
    # corr_looks = float(atr.get('NCORRLOOKS', rg_looks * az_looks / 1.94))

    altitude = 800000.0
    earth_radius = 6371000.0
    rg_looks =1
    az_looks =1
    corr_looks =1

    ## setup SNAPHU
    # Use these defaults if they weren't defined
    if 'init_method' not in vars(ps):
        print('params not found in ps namespace. Using defaults')
        ps.init_method = 'MCF'
        ps.cost_mode = 'SMOOTH'
        ps.max_comp = 32
        ps.defo_max = 0
    # https://web.stanford.edu/group/radar/softwareandlinks/sw/snaphu/snaphu.conf.full
    # https://github.com/isce-framework/isce2/blob/main/contrib/Snaphu/Snaphu.py
    print('phase unwrapping with SNAPHU ...')
    print(f'SNAPHU cost mode: {ps.cost_mode}')
    print(f'SNAPHU init only: {init_only}')
    print(f'SNAPHU init method: {ps.init_method}')
    print(f'SNAPHU max number of connected components: {ps.max_comp}')
    


    snp = Snaphu()

    # file IO
    snp.setInput(int_file)
    snp.setOutput(unw_file)
    snp.setCorrfile(cor_file)
    snp.setCorFileFormat('FLOAT_DATA')
    snp.setWidth(width)

    
    # runtime options
    snp.setCostMode(ps.cost_mode)
    snp.setInitOnly(init_only)
    snp.setInitMethod(ps.init_method)

    # geometry parameters
    # baseline info is not used in deformation mode, but is very important in topography mode
    snp.setAltitude(altitude)
    snp.setEarthRadius(earth_radius)
    snp.setWavelength(wavelength)
    snp.setRangeLooks(rg_looks)
    snp.setAzimuthLooks(az_looks)
    snp.setCorrLooks(corr_looks)

    # deformation mode parameters
    snp.setDefoMaxCycles(ps.defo_max)

    # connected component control
    # grow connectedc components if init_only is True
    # https://github.com/isce-framework/isce2/blob/main/contrib/Snaphu/Snaphu.py#L413
    snp.setMaxComponents(ps.max_comp)

    ## run SNAPHU
    snp.prepare()
    snp.unwrap()
    print('finished SNAPHU running')

    # mask out wired values from SNAPHU
    # based on https://github.com/isce-framework/isce2/pull/326
    flag = np.fromfile(int_file, dtype=np.complex64).reshape(length, width)
    data = np.memmap(unw_file, dtype='float32', mode='r+', shape=(length*2, width))
    data[0:length*2:2, :][np.nonzero(flag == 0)] = 0
    data[1:length*2:2, :][np.nonzero(flag == 0)] = 0

    ## render metadata
    print(f'write metadata file: {unw_file}.xml')
    util.write_xml(unw_file,width,length,2,'FLOAT','BIL')
    if snp.dumpConnectedComponents:
        print(f'write metadata file: {unw_file}.conncomp.xml')
        util.write_xml(f'{unw_file}.conncomp',width,length,1,'BYTE','BIP')

    # time usage
    m, s = divmod(time.time() - start_time, 60)
    print(f'time used: {m:02.0f} mins {s:02.1f} secs.')

    return unw_file
