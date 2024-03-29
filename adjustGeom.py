#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 15:46:09 2023
@author: km

First, change the crop values in localParams
Read from ps.npy. this should have all of the processing params
Crop and downlook geom files (used in mintpy)

"""
import numpy as np
import os,sys,glob,argparse,time,re
from datetime import date
import isce.components.isceobj as isceobj
from mroipac.looks.Looks import Looks
from SARTS import util,config,network
from osgeo import gdal
import rasterio

def cmdLineParser():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser(
        description='Crop and downlook geom files. Save parameters',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--downlook', action='store_true', dest='doDownlook',help='Downlook geometry files.')
    parser.add_argument('-c', '--crop', action='store_true', dest='doCrop',help='Crop geometry files.')
    parser.add_argument('-r', '--replace', action='store_true', dest='replace',help='Overwrite cropped/downlooked geometry files')
    parser.add_argument('-f', '--fix-images', action='store_true', dest='fixImages',help='Fix file path in xml files (use if files were moved to different directory')
    parser.add_argument('-p', '--plot-off', action='store_false', dest='plot',help='Turn plotting off')

    return parser.parse_args()


def getbl(d):
    bl_file = os.path.join('./merged/baselines', d, d+'.vrt')
    ds = gdal.Open(bl_file)
    bl = ds.GetVirtualMemArray()
    bl = np.nanmean(bl)
    return bl


def update_yaml_key(file_path, key, new_value):
    with open(file_path, "r") as f:
        lines = f.readlines()

    with open(file_path, "w") as f:
        for line in lines:
            # Try to match a YAML key-value pair line
            match = re.match(rf"({key}\s*:\s*)(\S+)", line)
            if match:
                # Replace the value while preserving the key and any surrounding whitespace
                line = f"{match.group(1)}{new_value}\n"
            f.write(line)

          
def geo2rdr(ps):
    print('Getting approximate crop bounds from geobbox')
    print('Loading lon.rdr.full and lat.rdr.full')
    if ps.sat=='ALOS':
        infile = os.path.join(ps.mergeddir,'geom_reference','lon.rdr')
    else:
        infile= os.path.join(ps.mergeddir,'geom_reference','lon.rdr.full')
    imgi = isceobj.createImage()
    imgi.load(infile+'.xml')
    lon_full = imgi.memMap().copy()
    lon_full = lon_full.astype(float)  # Convert to float
    lon_full = np.squeeze(lon_full)
    lon_full[lon_full==0] = np.nan
    
    if ps.sat=='ALOS':
        infile = os.path.join(ps.mergeddir,'geom_reference','lat.rdr')
    else:
        infile= os.path.join(ps.mergeddir,'geom_reference','lat.rdr.full')
    
    imgi = isceobj.createImage()
    imgi.load(infile+'.xml')
    lat_full = imgi.memMap().copy()
    lat_full = lat_full.astype(float)  # Convert to float
    lat_full = np.squeeze(lat_full)
    lat_full[lat_full==0] = np.nan
    print(ps.geobbox)
    latmin,latmax,lonmin,lonmax = ps.geobbox
    y_ll,x_ll = util.ll2pixel(lon_full,lat_full,lonmin,latmin)
    y_lr,x_lr = util.ll2pixel(lon_full,lat_full,lonmax,latmin)
    y_ul,x_ul = util.ll2pixel(lon_full,lat_full,lonmin,latmax)
    y_ur,x_ur = util.ll2pixel(lon_full,lat_full,lonmax,latmax)
    
    ps.cropymin = np.nanmin([y_ll,y_lr,y_ul,y_ur])
    ps.cropymax = np.nanmax([y_ll,y_lr,y_ul,y_ur])
    ps.cropxmin = np.nanmin([x_ll,x_lr,x_ul,x_ur])
    ps.cropxmax = np.nanmax([x_ll,x_lr,x_ul,x_ur])
    
    update_yaml_key("params.yaml", "cropymin", ps.cropymin)
    update_yaml_key("params.yaml", "cropymax", ps.cropymax)
    update_yaml_key("params.yaml", "cropxmin", ps.cropxmin)
    update_yaml_key("params.yaml", "cropxmax", ps.cropxmax)


def main(inps):
    ps = config.getPS()

    # Update crop bounds if geobbox exists
    if 'geobbox' in ps.__dict__.keys() and  ps.geobbox is not None:
        if np.sum([ps.cropymin,ps.cropymax,ps.cropxmin,ps.cropxmax]) == 0:
            geo2rdr(ps)
        else:
            print('crop coordinates found')
    else:
        print('geobbox not found')
    
    
    if ps.crop:
        if np.sum([int(ps.cropymin),int(ps.cropymax),int(ps.cropxmin),int(ps.cropxmax)]) ==0:
            print('Crop is set to True but values are all zero. Set crop bounds and rerun.')
        else:
            print('Crop bounds found.')
    else:
        print('Crop set to False... not cropping')
        time.sleep(5)
    if inps.plot:
        import matplotlib.pyplot as plt
        import logging
        # Set the logging level for matplotlib to suppress DEBUG messages
        logging.getLogger('matplotlib').setLevel(logging.INFO)  
        plt.close('all')
    
    if inps.replace:
        print('replacing files...')
        os.system('rm ./merged/geom_reference/*crop*')
        os.system('rm ./merged/geom_reference/*lk*')
        os.system('rm ./merged/SLC/*/*crop*')
    
    # Make directories
    if not os.path.isdir(ps.tsdir):
        os.mkdir(ps.tsdir)
    if not os.path.isdir('Npy'):
        os.mkdir('Npy')
    if not os.path.isdir(ps.workdir + '/Figs'):
        os.mkdir(ps.workdir + '/Figs')
    
    if inps.fixImages:
        if ps.sensor == 'ALOS':
            geomList = glob.glob(ps.mergeddir + '/geom_reference/*rdr')
            slcList = glob.glob(ps.slcdir + '/*/.slc')
            blList = glob.glob(ps.mergeddir + '/baselines/????????/????????')
        else:
            geomList = glob.glob(ps.mergeddir + '/geom_reference/*full')
            slcList = glob.glob(ps.slcdir + '/*/*full')
            blList = glob.glob(ps.mergeddir + '/baselines/????????/????????')
        for fname in slcList:
            os.system('fixImageXml.py -i ' + fname + ' -f')
        for fname in geomList:
            os.system('fixImageXml.py -i ' + fname + ' -f')
        for fname in blList:
            os.system('fixImageXml.py -i ' + fname + ' -f')
    
    # Get the watermask using the nlcd land cover map        
    if ps.sat == 'ALOS':
        waterMaskFn = ps.mergeddir + '/geom_reference/waterMask.rdr'
    else:
        waterMaskFn = ps.mergeddir + '/geom_reference/waterMask.rdr.full'
    if ps.waterMask:
        if not os.path.isfile(waterMaskFn):
            from SARTS.landCover2rdr import convert_land_cover
            convert_land_cover(ps.nlcd_in,ps)
    
    # Get the list of geometry files to work on
    if ps.sat == 'ALOS':       
        geomList = glob.glob(ps.mergeddir + '/geom_reference/*rdr')
        slcList = glob.glob(ps.slcdir + '/*/*slc')


    else:
        geomList = glob.glob(ps.mergeddir + '/geom_reference/*full')
        slcList = glob.glob(ps.slcdir + '/*/*full')


    geomList = [item for item in geomList if '_lk' not in item]
    slcList.sort()
    # Get the acquisition dates
    flist = glob.glob( os.path.join(ps.slcdir, '2*'))
    dates = []
    for f in flist:
        dates.append(f[-8:])
    dates.sort()
    
    ghosts = []
    print('searching for ghosts...')

    for ii in range(len(dates)):
        da = dates[ii]
        if ps.sat == 'ALOS':
            fn = os.path.join(ps.slcdir,da,da + '.slc')
        else:
            fn = os.path.join(ps.slcdir,da,da + '.slc.full')
        if not os.path.isfile(fn):
            ghosts.append(da)
            print('Warning: ' + fn + ' was not found.')
    
    if len(ghosts)>0:
        print('Missing SLCs.  Either processes them or delete the directory in merged/SLC/')
        sys.exit(1)
    else:
        print('No missing files found')


    networkObj = network.Network()
    networkObj.dateList = dates
    networkObj.baselineDict[ps.reference_date] = 0.0
    
    
    def getbl(secDir):
        print(secDir)
        bl_file =secDir  + '/' + secDir.split('/')[-1] + '.vrt'
        ds = gdal.Open(bl_file)
        bl = ds.GetVirtualMemArray()
        bl = np.nanmean(bl)
        return bl

    bls = []
    for d in networkObj.dateList:
        if d == ps.reference_date:
            bls.append(0)
        else:
            secDir = './merged/baselines/' + d
            baseline = getbl(secDir)
            networkObj.baselineDict[d] = baseline
            bls.append(float(baseline))

    if ps.networkType=='singleMaster':
        networkObj.single_master()
    elif ps.networkType=='delaunay':
        networkObj.delaunay()
    elif ps.networkType=='maxbandwidth':
        networkObj.maxbandwidth(ps.bandwidth)
    else:
        print('choose valid networkType in ps.networkType')
    
    #print(networkObj.pairsDates)
    # for pair in networkObj.pairsDates:
        # print(pair)
    print(str(len(networkObj.pairsDates)) + ' pairs')


    def is_leap_year(year):
        return (year % 4 == 0 and year % 100 != 0) or year % 400 == 0
    dec_year = []
    dn = []
    for d in networkObj.dateList:
        yr, mo, day = int(d[0:4]), int(d[4:6]), int(d[6:8])
        
        # Convert to date object and get ordinal
        current_date = date(yr, mo, day)
        dt = current_date.toordinal()
        
        # Ordinal of the first day of the year
        d0 = date(yr, 1, 1).toordinal()
        
        # Calculate day of the year
        doy = dt - d0 + 1
        
        # Determine if it's a leap year
        is_leap = is_leap_year(yr)
        
        # Adjust the divisor for leap years
        days_in_year = 366 if is_leap else 365
        
        # Calculate decimal year
        dec_year.append(yr + (doy - 1) / days_in_year)
        dn.append(dt)
    
    # If needed, convert dec_year to a NumPy array
    dec_year = np.array(dec_year)
    dn = np.asarray(dn)
    dn0 = dn-dn[0] # make relative to first date
    
    # if inps.plot or ps.networkType=='delaunay':
    #     bls = []
    #     for d in networkObj.dateList:
    #         if d == ps.reference_date:
    #             bls.append(0)
    #         else:
    #             baseline = getbl(d)
    #             networkObj.baselineDict[d] = baseline
    #             bls.append(float(baseline))

    if inps.plot:
        plt.figure()
        plt.scatter(dec_year,bls)
        plt.xlabel('Time (yrs)')
        plt.ylabel('Perpindicular baseline (m)')
    
        for p in networkObj.pairsDates:
            id1 = networkObj.dateList.index(p.split('_')[0])
            id2 = networkObj.dateList.index(p.split('_')[1])
            bl1 = bls[id1]
            bl2 = bls[id2]
            decy1 =dec_year[id1]
            decy2 = dec_year[id2]
            if inps.plot:
                plt.plot([decy1,decy2],[bl1,bl2])

    pairs_seq = []
    for ii,d in enumerate(dates[0:-1]):
        pairs_seq.append(dates[ii] + '_' + dates[ii+1])
    
    
    # Now make pairs2
    pairs = networkObj.pairsDates
    nd = len(pairs)

    # Get width and length
    if ps.sat == 'ALOS':
        f_lon = ps.mergeddir + '/geom_reference/lon.rdr'
        ref_slc = os.path.join(ps.slcdir, dates[0],dates[0]+'.slc')
        os.system('gdal2isce_xml.py -i ' + ref_slc)
    else:
        f_lon = ps.mergeddir + '/geom_reference/lon.rdr.full'
        
    gImage = isceobj.createIntImage()
    gImage.load(f_lon + '.xml')
    nyf = gImage.length
    nxf = gImage.width
    
    if ps.crop:
        ny = ps.cropymax-ps.cropymin
        nx = ps.cropxmax-ps.cropxmin
    else:
        ny = gImage.length
        nx = gImage.width
        ps.cropxmin=0
        ps.cropxmax=nx
        ps.cropymin=0
        ps.cropymax=ny
    
    
    nxl = nx//int(ps.rlks)
    nyl = ny//int(ps.alks)
    ps.ny =         ny
    ps.nx =         nx
    ps.nxl =        nxl
    ps.nyl =        nyl
    ps.nxf =        nxf
    ps.nyf =        nyf


    if inps.doCrop and ps.crop:
        for infile in geomList:
            if os.path.isfile(infile):
                if not os.path.isfile(infile+'.crop'):
                    imgi = isceobj.createImage()
                    imgi.load(infile+'.xml')
                    geomIm = imgi.memMap()
                    
                    # Adapt to variable dimension order from isce images
                    shape = geomIm.shape
                    y_dim = shape.index(ps.nyf)
                    x_dim = shape.index(ps.nxf)
                    # z_dim = 3 - y_dim - x_dim
                    # Perform cropping dynamically based on identified dimensions
                    slices = [slice(None)] * 3  # Default slice (full range) for all dimensions
                    slices[y_dim] = slice(ps.cropymin, ps.cropymax)
                    slices[x_dim] = slice(ps.cropxmin, ps.cropxmax)
                    geomIm = geomIm[tuple(slices)]

                    # Write out cropped file
                    imgo = imgi.clone()
                    imgo.filename = infile+'.crop'
                    imgo.width  = ps.nx #ps.cropxmax-ps.cropxmin
                    imgo.length = ps.ny #ps.cropymax-ps.cropymin
                    geomIm.tofile(imgo.filename)
                    imgo.dump(imgo.filename+'.xml')
                    imgo.finalizeImage()
                    del(geomIm)
            else:
                print('no file ' + infile)
    else:
        print('Skipping cropping geom')
    

    if inps.doCrop and ps.crop and ps.doCropSlc:
        for infile in slcList:
            if os.path.isfile(infile):
                if not os.path.isfile(infile+'.crop'):
                    imgi = isceobj.createImage()
                    imgi.load(infile+'.xml')
                    slcIm = imgi.memMap()
                    
                    # Adapt to variable dimension order from isce images
                    shape = slcIm.shape
                    y_dim = shape.index(ps.nyf)
                    x_dim = shape.index(ps.nxf)
                    # Perform cropping dynamically based on identified dimensions
                    slices = [slice(None)] * 3  # Default slice (full range) for all dimensions
                    slices[y_dim] = slice(ps.cropymin, ps.cropymax)
                    slices[x_dim] = slice(ps.cropxmin, ps.cropxmax)
                    slcIm = slcIm[tuple(slices)]

                    # Write out cropped file
                    imgo = isceobj.createImage()
                    imgo.filename = infile+'.crop'
                    imgo.width  = ps.nx #ps.cropxmax-ps.cropxmin
                    imgo.length = ps.ny #ps.cropymax-ps.cropymin
                    imgo.dataType='CFLOAT'
                    slcIm.tofile(imgo.filename)
                    imgo.dump(imgo.filename+'.xml')
                    imgo.finalizeImage()
                    del(slcIm)
            else:
                print('no file ' + infile)
    else:
        print('Skipping cropping slc')
    
    if inps.doDownlook:
        if ps.crop:
            fList = glob.glob(ps.mergeddir + '/geom_reference/*crop')
        else:
            fList = glob.glob(ps.mergeddir + '/geom_reference/*rdr.full')

        def downLook(infile, outfile,alks,rlks):
            inImage = isceobj.createImage()
            inImage.load(infile + '.xml')
            inImage.filename = infile
            lkObj = Looks()
            lkObj.setDownLooks(alks)
            lkObj.setAcrossLooks(rlks)
            lkObj.setInputImage(inImage)
            lkObj.setOutputFilename(outfile)
            lkObj.looks()
            
        for infile in fList:
            f = infile.split('/')[-1].split('.')[0]
            if os.path.isfile(infile):
                outfile = ps.mergeddir + '/geom_reference/' + f + '_lk.rdr'
    
                if not os.path.isfile(outfile):
                    print('downlooking ' + f)
                    downLook(infile, outfile, ps.alks, ps.rlks)
                    
                else:
                    print(outfile + ' already exists')
    
            else:
                print('no file ' + infile)
    else:
        print('skipping downlooking')
    
    
    fList = glob.glob(ps.mergeddir + '/geom_reference/*lk.rdr')
    geom_data = {}
    
    
    for infile in fList:
        
        imgi = isceobj.createImage()
        imgi.load(infile+'.xml')
        geomIm = imgi.memMap().copy()
        geomIm = geomIm.astype(float)  # Convert to float

        geomIm = np.squeeze(geomIm)
        geomIm[geomIm==0] = np.nan
        var_name = infile.split('/')[-1].split('.')[0]
        geom_data[var_name] = geomIm


    # Separate az and los from los_lk
    axis = geom_data['los_lk'].shape.index(2)
    los_ifg, az_ifg = np.split(geom_data['los_lk'], 2, axis=axis)
    los_ifg = np.ascontiguousarray(np.squeeze(los_ifg)).astype(np.float32)
    az_ifg =  np.ascontiguousarray(np.squeeze(az_ifg)).astype(np.float32)
    
    _, inc_ifg = np.split(geom_data['incLocal_lk'], 2, axis=axis)
    inc_ifg = np.ascontiguousarray(np.squeeze(inc_ifg)).astype(np.float32)

    # Write out a new az file
    # azOutname = ps.mergeddir + '/geom_reference/az_lk.rdr'
    # util.writeISCEimg(az_ifg.astype(np.float32), azOutname, 1, nxl, nyl, 'FLOAT')
    azOutname = ps.mergeddir + '/geom_reference/az_lk.rdr'
    fidc=open(azOutname,"wb")
    fidc.write(az_ifg)
    #write out an xml file for it
    out = isceobj.createIntImage() # Copy the interferogram image from before
    out.dataType = 'FLOAT'
    out.bands = 1
    out.filename = azOutname
    out.width = nxl
    out.length = nyl
    out.dump(azOutname + '.xml') # Write out xml
    out.renderHdr()
    out.renderVRT()
   
    # Write out a new inc file
    incOutname = ps.mergeddir + '/geom_reference/inc_lk.rdr'
    fidc=open(incOutname,"wb")
    fidc.write(inc_ifg)
    #write out an xml file for it
    out = isceobj.createIntImage() # Copy the interferogram image from before
    out.dataType = 'FLOAT'
    out.bands = 1
    out.filename = incOutname
    out.width = nxl
    out.length = nyl
    out.dump(incOutname + '.xml') # Write out xml
    out.renderHdr()
    out.renderVRT()

    if not os.path.isdir('dolphin'):
        os.mkdir('dolphin')

    if ps.waterMask:
        # Make a waterMask tif for dolhpin
        if not os.path.isfile('dolphin/nodata_mask.tif') or inps.replace:
            if ps.crop:
                if ps.sat=='ALOS':
                    ds = gdal.Open(ps.mergeddir + '/geom_reference/waterMask.rdr.crop.vrt')
                else:
                    ds = gdal.Open(ps.mergeddir + '/geom_reference/waterMask.rdr.full.crop.vrt') 
            else:
                ds = gdal.Open(ps.mergeddir + '/geom_reference/waterMask.rdr.full.vrt')

            wm = ds.GetVirtualMemArray()
            with rasterio.open('dolphin/nodata_mask.tif', 'w', driver='GTiff',
                    height=wm.shape[0], width=wm.shape[1],
                    count=1, dtype=np.uint8) as dst:
                dst.write(wm, 1)  # Writing data to the first band
        else:
            print('dolphin/nodata_mask.tif watermask already exists.. skipping')
    
    if inps.plot:
        cmap = 'Spectral_r'
        fig,ax = plt.subplots(3,2,figsize=(9,9))
        ax[0,0].imshow(geom_data['lon_lk'],cmap=cmap);ax[0,0].set_title('lon_ifg')
        ax[0,1].imshow(geom_data['lat_lk'],cmap=cmap);ax[0,1].set_title('lat_ifg')
        ax[1,0].imshow(geom_data['hgt_lk'],cmap=cmap);ax[1,0].set_title('hgt_ifg')
        ax[1,1].imshow(los_ifg,cmap=cmap);ax[1,1].set_title('los_ifg')
        ax[2,0].imshow(geom_data['shadowMask_lk'],cmap=cmap);ax[2,0].set_title('shm_ifg')
        ax[2,1].imshow(inc_ifg,cmap=cmap);ax[2,1].set_title('inc_ifg')
        plt.savefig(ps.workdir + '/Figs/geom.svg',transparent=True,dpi=100 )
        plt.show()
    
    ps.dates =      dates
    ps.pairs =      pairs
    ps.pairs_seq =  pairs_seq
    ps.dec_year =   dec_year
    ps.dn =         dn
    ps.dn0 =        dn0
    ps.nd =         nd
    ps.minlon =     np.nanmin(geom_data['lon_lk'])
    ps.maxlon =     np.nanmax(geom_data['lon_lk'])
    ps.minlat =     np.nanmin(geom_data['lat_lk'])
    ps.maxlat =     np.nanmax(geom_data['lat_lk'])
    
    # Save the namespace
    np.save('./ps.npy',ps)
    np.save('./geom.npy',geom_data)

if __name__ == '__main__':
    '''
    Main driver.
    '''
    # inps = argparse.Namespace()
    # inps.doDownlook = True
    # inps.doCrop = True
    # inps.replace = False
    # inps.fixImages = False
    # inps.plot = False
    inps = cmdLineParser()
    main(inps)
