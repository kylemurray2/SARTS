import os
import numpy as np


def getLocalParams():
    import argparse

    if os.path.isfile('./ps.npy'):
        ps = np.load('./ps.npy',allow_pickle=True).all()
    else:
        ps = argparse.Namespace()

    # FOR PyPS:________________________________________________________________
    ps.workdir = os.getcwd() # Use current directory as working directory
     # working directory (should be where merged is)
    ps.skip = 1
    ps.alks = int(1) # number of looks in azimuth
    ps.rlks = int(4) # number of looks in range
    ps.networkType = 'delaunay' # sequential1, sequential#,

    if 'nx' in ps.__dict__.keys():
        ps.nxl = ps.nx//ps.rlks
        ps.nyl = ps.ny//ps.alks

    ps.seaLevel = -200
    ps.ifgMode = False

    ps.crop=True
    ps.cropymin = 0
    ps.cropymax = 0
    ps.cropxmin = 0
    ps.cropxmax = 0

    ps.mergeddir    = ps.workdir + '/merged'
    ps.intdir       = ps.mergeddir + '/interferograms'
    ps.tsdir        = ps.workdir + '/TS'
    ps.slcdir       = ps.mergeddir + '/SLC'


    # The following are only necessary to run setupStack (not setupPyPS.py)
    #________________________change these______________________________________
    ps.path         = ''  # '###'
    ps.fl           = ''  #'Ascending' or 'Descending
    ps.bounds       = '' #'minlat,maxlat,minlon,maxlon'  # it won't use bursts outside of this area
    ps.poly         = ''
    ps.swath_num    = '' #default='1 2 3'
    ps.reference_date = '' # YYYYMMDD default=None
    #__________________________________________________________________________

    ps.start        = '2014-01-01T00:00:00Z'    #'2010-05-01T00:00:00Z'
    ps.end          = '2028-01-25T23:59:00Z'   #'2020-02-01T23:59:00Z'
    ps.sat          ='SENTINEL-1'
    ps.sensor       = ps.sat
    ps.point        = None # lon,lat
    ps.frame        = None #'472,473,474,475'   #'####'
    ps.slc_dirname  = './SLCS/'
    ps.lam          = 0.56

    # Fringe___________________________________
    ps.waterMask = True
    # NLCD raster for landcover and watermask
    ps.nlcd_in = '/d/HI/Landcover/data/hi_oahu_2011_ccap_hr_land_cover20140619.img'
    ps.nlcd_tif = '/d/HI/Landcover/data/oahu.tif'
    ps.maxMem = 100
    ps.filterFlag      = True
    ps.filterStrength  = 0.3
    ps.unwrap = True
    #__________________________________________________________________________

    # FOR ISCE STACK PROCESSOR_________________________________________________
    ps.orbit_dirname           = './orbits'
    ps.aux_dirname             = './aux_cal/'
    ps.work_dir                = './'
    ps.exclude_dates           = None #default=None
    ps.include_dates           = None #default=None
    ps.polarization            = 'vv' #default='vv'
    ps.coregistration          = 'NESD' #default='NESD'
    ps.virtualMerge            = True #default=None
    ps.useGPU                  = True # default=False
    ps.numProcess              = 1
    ps.numProcess4topo         = 1
    ps.azimuthLooks            = '1'
    ps.rangeLooks              = '1'
    ps.filtStrength            = '0.1'
    ps.esdCoherenceThreshold   = '0.85'
    ps.startDate               = ps.start[0:10]
    ps.stopDate                = ps.end[0:10]
    ps.num_connections         = '1'
    ps.num_overlap_connections = '3'
    miny,maxy,minx,maxx        = ps.bounds.split(sep=',')
    ps.bbox                    = miny +' '+ maxy  +' '+ minx +' '+  maxx #demBounds[0] + ' ' + demBounds[1] + ' ' + demBounds[2] + ' ' + demBounds[3] #SNWE
    ps.workflow                = 'slc' #default='interferogram'
    ps.text_cmd                =''
    ps.snrThreshold            ='10'
    ps.unwMethod               ='snaphu'
    ps.rmFilter                =False
    ps.param_ion               =None
    ps.num_connections_ion     ='3'

    np.save('ps.npy',ps)

    return ps
