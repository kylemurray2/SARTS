#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 15:46:09 2023

First, change the crop values in localParams
Read from ps.npy. this should have all of the processing params
Crop and downlook geom files (used in mintpy)


@author: km
"""



import numpy as np
import glob
import os
from datetime import date
import isce.components.isceobj as isceobj
from mroipac.looks.Looks import Looks
import localParams
from SARTS import util
from Network import Network
from osgeo import gdal
import argparse


def cmdLineParser():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser(
        description='Crop and downlook geom files. Save parameters',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', '--plot', type=bool, dest='plot', default=True)
    parser.add_argument('-d', '--downlook', type=bool, dest='doDownlook', default=True)
    parser.add_argument('-r', '--replace', type=bool, dest='replace', default=False)
    parser.add_argument('-f', '--fix-images', type=bool, dest='fixImages', default=False)

    return parser.parse_args()

def getbl(secDir):
    bl_file =secDir  + '/' + secDir.split('/')[-1] + '.vrt'
    ds = gdal.Open(bl_file)
    bl = ds.GetVirtualMemArray()
    bl = np.nanmean(bl)
    return bl

def main(inps):
    ps = localParams.getLocalParams()

    if inps.plot:
        import matplotlib.pyplot as plt
        plt.close('all')
    
    if inps.replace:
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
    if ps.waterMask:
        if not os.path.isfile(ps.mergeddir + '/geom_reference/waterMask.rdr.full'):
            from SARTS.landCover2rdr import convert_land_cover
            convert_land_cover(ps.nlcd_in)
    
    # Get the list of geometry files to work on
    geomList = glob.glob(ps.mergeddir + '/geom_reference/*full')
    geomList = [item for item in geomList if '_lk' not in item]
    
    # Get the acquisition dates
    flist = glob.glob(ps.slcdir + '/2*')
    dates = []
    for f in flist:
        dates.append(f[-8:])
    dates.sort()
    
    networkObj = Network()
    networkObj.dateList = dates
    networkObj.baselineDict[ps.reference_date] = 0.0
    
    
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
    elif ps.networkType=='sequential1':
        networkObj.sequential1()
    elif ps.networkType=='sequential2':
        networkObj.sequential2()
    elif ps.networkType=='sequential3':
        networkObj.sequential3()
    elif ps.networkType=='sequential4':
        networkObj.sequential4()
    elif ps.networkType=='sequential5':
        networkObj.sequential5()
    else:
        print('choose valid networkType in ps.networkType')
    
    print(networkObj.pairsDates)
    for pair in networkObj.pairsDates:
        print(pair)
    print(len(networkObj.pairsDates))

    dec_year = []
    dn = []
    for d in networkObj.dateList:
        yr = d[0:4]
        mo = d[4:6]
        day = d[6:8]
        dt = date.toordinal(date(int(yr), int(mo), int(day)))
        dn.append(dt)
        dt = date.toordinal(date(int(yr), int(mo), int(day)))
        d0 = date.toordinal(date(int(yr), 1, 1))
        doy = np.asarray(dt)-d0+1
        dec_year.append(float(yr) + (doy/365.25))
    dn = np.asarray(dn)
    dn0 = dn-dn[0] # make relative to first date
    
    if inps.plot:
        plt.figure()
        plt.scatter(dec_year,bls)
        plt.xlabel('Time (yrs)')
        plt.ylabel('Perpindicular baseline (m)')
    
    for p in networkObj.pairsDates:
        id1 = networkObj.dateList.index(p.split('_')[0])
        id2 = networkObj.dateList.index(p.split('_')[1])
        bl1 = bls[id1]
        bl2 =bls[id2]
        decy1 =dec_year[id1]
        decy2 = dec_year[id2]
        if inps.plot:
            plt.plot([decy1,decy2],[bl1,bl2])

    pairs = []
    for ii,d in enumerate(dates[0:-1]):
        pairs.append(dates[ii] + '_' + dates[ii+1])
    
    
    # Now make pairs2
    pairs2 = networkObj.pairsDates
    nd = len(pairs)

    # Get width and length
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

    
    if ps.crop:
        for infile in geomList:
            if os.path.isfile(infile):
                if not os.path.isfile(infile+'.crop'):
                    imgi = isceobj.createImage()
                    imgi.load(infile+'.xml')
                    if infile.split('/')[-1] in ['los.rdr']:
                        imgi.scheme = 'BSQ'
                        imgi.imageType = 'bsq'
                    # print(imgi.memMap().shape)
                    # Rearrange axes order from small to big
                    geomIm = util.orderAxes(imgi.memMap(),ps.nxf,ps.nyf)
                    geomIm = geomIm[:,ps.cropymin:ps.cropymax,ps.cropxmin:ps.cropxmax]
                    # geomIm = geomIm[:,ps.cropymin:ps.cropymax,ps.cropxmin:ps.cropxmax]
                    imgo = imgi.clone()
                    imgo.filename = infile+'.crop'
                    imgo.width  = ps.nx #ps.cropxmax-ps.cropxmin
                    imgo.length = ps.ny #ps.cropymax-ps.cropymin
                    if infile.split('/')[-1] in ['incLocal.rdr']:
                        imgo.scheme = 'BSQ'
                    geomIm.tofile(imgo.filename)
                    imgo.dump(imgo.filename+'.xml')
                    imgo.finalizeImage()
                    del(geomIm)
            else:
                print('no file ' + infile)
    
    
    
    if inps.doDownlook:
        if ps.crop:
            fList = glob.glob(ps.mergeddir + '/geom_reference/*crop')
        else:
            fList = glob.glob(ps.mergeddir + '/geom_reference/*rdr')
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
    
    
    # Get bounding coordinates (Frame)
    f_lon_lk = ps.mergeddir + '/geom_reference/lon_lk.rdr'
    f_lat_lk = ps.mergeddir + '/geom_reference/lat_lk.rdr'
    f_hgt_lk = ps.mergeddir + '/geom_reference/hgt_lk.rdr'
    f_los_lk = ps.mergeddir + '/geom_reference/los_lk.rdr'
    f_shm_lk = ps.mergeddir + '/geom_reference/shadowMask_lk.rdr'
    f_inc_lk = ps.mergeddir + '/geom_reference/incLocal_lk.rdr'
    
    # LON --------------
    Image = isceobj.createImage()
    Image.load(f_lon_lk + '.xml')
    lon_ifg = util.orderAxes(Image.memMap(),nxl,nyl)[0,:,:]
    lon_ifg = lon_ifg.copy().astype(np.float32)
    lon_ifg[lon_ifg==0]=np.nan
    Image.finalizeImage()
    
    # LAT --------------
    Image = isceobj.createImage()
    Image.load(f_lat_lk + '.xml')
    lat_ifg =util.orderAxes(Image.memMap(),nxl,nyl)[0,:,:]
    lat_ifg = lat_ifg.copy().astype(np.float32)
    lat_ifg[lat_ifg==0]=np.nan
    Image.finalizeImage()
    
    # HGT --------------
    Image = isceobj.createImage()
    Image.load(f_hgt_lk + '.xml')
    hgt_ifg = util.orderAxes(Image.memMap(),nxl,nyl)[0,:,:]
    hgt_ifg = hgt_ifg.copy().astype(np.float32)
    hgt_ifg[hgt_ifg==-500]=np.nan
    Image.finalizeImage()
    
    # LOS --------------
    Image = isceobj.createImage()
    Image.load(f_los_lk + '.xml')
    # Image.bands=2
    # Image.scheme='BIP'
    los_ifg = util.orderAxes(Image.memMap(),nxl,nyl)[0,:,:]
    los_ifg = los_ifg.copy()
    util.show(los_ifg)
    az_ifg = util.orderAxes(Image.memMap(),nxl,nyl)[1,:,:]
    az_ifg = az_ifg.copy()
    Image.finalizeImage()
    
    # Write out a new los file
    losOutname = ps.mergeddir + '/geom_reference/los2_lk.rdr'
    fidc=open(losOutname,"wb")
    fidc.write(los_ifg)
    #write out an xml file for it
    out = isceobj.createIntImage() # Copy the interferogram image from before
    out.dataType = 'FLOAT'
    out.bands = 1
    out.filename = losOutname
    out.width = nxl
    out.length = nyl
    out.dump(losOutname + '.xml') # Write out xml
    out.renderHdr()
    out.renderVRT()
    
    
    # Write out a new az file
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
    
    # if you want to save these to geom
    los_ifg = los_ifg.copy().astype(np.float32)
    los_ifg[los_ifg==0]=np.nan
    az_ifg = az_ifg.copy().astype(np.float32)
    az_ifg[az_ifg==0]=np.nan
    
    Image = isceobj.createImage()
    Image.load(f_shm_lk + '.xml')
    Image.bands=1
    shm_ifg = util.orderAxes(Image.memMap(),nxl,nyl)[0,:,:]
    shm_ifg = shm_ifg.copy().astype(np.float32)
    shm_ifg[np.isnan(hgt_ifg)]=np.nan
    Image.finalizeImage()
    
    Image = isceobj.createImage()
    Image.load(f_inc_lk + '.xml')
    Image.bands=2
    Image.scheme='BSQ'
    inc = Image.memMap()
    #plt.figure();plt.imshow(Image.memMap()[0,:,:]);plt.show()
    # inc_ifg1 = Image.memMap()[0,:,:] # relative to the local plane of the ground
    inc_ifg = util.orderAxes(Image.memMap(),nxl,nyl)[1,:,:]# relative to surface normal vector (this is the one we want I think)
    inc_ifg = inc_ifg.copy()
    
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
    
    inc_ifg = inc_ifg.copy().astype(np.float32)
    inc_ifg[inc_ifg==0]=np.nan
    Image.finalizeImage()
    
    if inps.plot:
        cmap = 'Spectral_r'
        fig,ax = plt.subplots(3,2,figsize=(9,9))
        ax[0,0].imshow(lon_ifg,cmap=cmap);ax[0,0].set_title('lon_ifg')
        ax[0,1].imshow(lat_ifg,cmap=cmap);ax[0,1].set_title('lat_ifg')
        ax[1,0].imshow(hgt_ifg,cmap=cmap);ax[1,0].set_title('hgt_ifg')
        ax[1,1].imshow(los_ifg,cmap=cmap);ax[1,1].set_title('los_ifg')
        ax[2,0].imshow(shm_ifg,cmap=cmap);ax[2,0].set_title('shm_ifg')
        ax[2,1].imshow(inc_ifg,cmap=cmap);ax[2,1].set_title('inc_ifg')
        plt.savefig(ps.workdir + '/Figs/geom.svg',transparent=True,dpi=100 )
        plt.show()
    
    ps.dates =      dates
    ps.pairs =      pairs
    ps.pairs2 =     pairs2
    ps.dec_year =   dec_year
    ps.dn =         dn
    ps.dn0 =        dn0
    ps.nd =         nd
    
    ps.minlon =     np.nanmin(lon_ifg)
    ps.maxlon =     np.nanmax(lon_ifg)
    ps.minlat =     np.nanmin(lat_ifg)
    ps.maxlat =     np.nanmax(lat_ifg)
    
    # Save the namespace
    np.save('./ps.npy',ps)


if __name__ == '__main__':
    '''
    Main driver.
    '''
    inps = cmdLineParser()
    main(inps)
