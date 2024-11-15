#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 21:54:12 2020

@author: kdm95
"""
import numpy as np

def getTime(path,frame):
    '''
     Figure out what time the aquisition was
    '''
    import os
    import requests
    import pandas as pd
    path = str(path)
    frame = str(frame)
    # path = os.getcwd().split('/')[-2]
    # frame= os.getcwd().split('/')[-1]

    start='2014-05-01T00:00:00Z'
    end='2099-06-01T00:00:00Z'
    asfUrl = 'https://api.daac.asf.alaska.edu/services/search/param?platform=SENTINEL-1&processinglevel=SLC&output=CSV'
    call = asfUrl + '&relativeOrbit=' + path + '&frame=' + frame + '&start=' + start + '&end=' + end
    # Here we'll make a request to ASF API and then save the output info to .CSV file
    if not os.path.isfile('out.csv'):
        r =requests.get(call,timeout=100)
        with open('out.csv','w') as j:
            j.write(r.text)
    # Open the CSV file and get the URL and File names
    hour = pd.read_csv('out.csv')["Start Time"][0][11:13]
    minute = pd.read_csv('out.csv')["Start Time"][0][14:16]
    
    Lon = pd.read_csv('out.csv')["Near Start Lon"][0]
    Lat = pd.read_csv('out.csv')["Near Start Lat"][0]
    return int(hour),int(minute), Lon, Lat

def read_isce_image(file_name):
    from isce.components import isceobj
    imgi = isceobj.createImage()
    imgi.load(file_name+'.xml')
    outIM = imgi.memMap()
    
    return outIM


def writeGeotiff(array, lat_bounds, lon_bounds, output_file,epsg=4326):
    '''
    no data is zero
    '''
    from osgeo import gdal,osr
    rows, cols = array.shape

    # Define geotransform parameters (top-left corner coordinates, pixel size)
    x_min, y_max = lon_bounds[0], lat_bounds[1]
    x_max, y_min = lon_bounds[1], lat_bounds[0]
    x_pixel_size = (x_max - x_min) / float(cols)
    y_pixel_size = (y_max - y_min) / float(rows)
    geotransform = (x_min, x_pixel_size, 0, y_max, 0, -y_pixel_size)

    # Create the GeoTIFF file
    driver = gdal.GetDriverByName('GTiff')
    dataset = driver.Create(output_file, cols, rows, 1, gdal.GDT_Float32)

    # Set the geotransform parameters and coordinate system (WGS84)
    dataset.SetGeoTransform(geotransform)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)  # WGS84
    dataset.SetProjection(srs.ExportToWkt())

    # Write the NumPy array to the GeoTIFF file and set NoData value
    raster_band = dataset.GetRasterBand(1)
    raster_band.WriteArray(array)
    raster_band.SetNoDataValue(np.nan)  # Set zero values as NoData (transparent)

    dataset.FlushCache()
    dataset = None  # Close the file
    
    
def reproject_to_wgs84(input_filename, output_filename, resolution_reduction_factor=1):
    '''
    Reprojects a single GeoTIFF file to EPSG:4326.
    
    Parameters:
        input_filename (str): Path to the input GeoTIFF file.
        output_filename (str): Path where the output GeoTIFF will be saved.
        resolution_reduction_factor (int, optional): Factor by which the resolution is reduced. Default is 1 (no reduction).
    '''
    import rasterio
    from rasterio.warp import calculate_default_transform, reproject, Resampling

    dst_crs = 'EPSG:4326'  # Destination CRS

    # Open the source dataset
    with rasterio.open(input_filename) as src:
        src_crs = src.crs  # Source CRS
        src_transform = src.transform  # Source affine transformation matrix

        # Calculate the transformation and dimensions for the output image
        dst_transform, width, height = calculate_default_transform(
            src_crs, dst_crs, src.width, src.height, *src.bounds,
            dst_width=src.width // resolution_reduction_factor,
            dst_height=src.height // resolution_reduction_factor
        )

        # Update the metadata for the destination file
        dst_kwargs = src.meta.copy()
        dst_kwargs.update({
            'crs': dst_crs,
            'transform': dst_transform,
            'width': width,
            'height': height,
            'driver': 'GTiff'
        })

        # Create the destination file and reproject
        with rasterio.open(output_filename, 'w', **dst_kwargs) as dst:
            reproject(
                source=rasterio.band(src, 1),
                destination=rasterio.band(dst, 1),
                src_transform=src_transform,
                src_crs=src_crs,
                dst_transform=dst_transform,
                dst_crs=dst_crs,
                resampling=Resampling.nearest
            )

    return output_filename


def geotiff_to_kmz(geotiff_path, kmz_path,vmin,vmax, average_elevation, colormap='RdBu'):
    import rasterio,simplekml,zipfile
    import matplotlib.pyplot as plt
    # Step 1: Read the bounds and data from the GeoTIFF file
    with rasterio.open(geotiff_path) as src:
        bounds = src.bounds
        data = src.read(1)  # Assuming single-band image

    # Create an alpha (transparency) channel, and set it to fully opaque (255)
    alpha_channel = np.ones(data.shape, dtype=np.uint8) * 255
    
    # Set alpha to 0 (fully transparent) where data is NaN
    alpha_channel[np.isnan(data)] = 0
    
    # Rescale data to Byte (values between 0 and 255)
    data_min, data_max = vmin,vmax
    norm_data = (data - data_min) / (data_max - data_min)  # Normalize data to 0-1 range
    
    # Apply colormap to normalized data
    cm = plt.get_cmap(colormap)
    colored_data = (cm(norm_data)[:, :, :3] * 255).astype(np.uint8)  # Convert colormap float values to uint8 for RGB image
    
    # Combine RGB and alpha channel to get an RGBA image
    rgba_data = np.dstack((colored_data, alpha_channel))
    
    # Step 2: Save the RGBA data to PNG
    png_path = "temp_image.png"
    with rasterio.open(png_path, 'w', driver='PNG', width=data.shape[1], height=data.shape[0], count=4, dtype=np.uint8) as dst:
        for k, channel in enumerate([rgba_data[:, :, i] for i in range(4)]):  # Split RGBA
            dst.write(channel, k + 1)

    # Step 3: Create a KML representation of the raster image
    kml = simplekml.Kml()
    # Assuming average_elevation contains the average elevation in meters
    ground = kml.newgroundoverlay(name="GroundOverlay")
    ground.icon.href = png_path
    ground.altitudemode = simplekml.AltitudeMode.relativetoground
    ground.altitude = 10000#average_elevation  # Set the altitude
    ground.latlonbox.north = bounds.top
    ground.latlonbox.south = bounds.bottom
    ground.latlonbox.east = bounds.right
    ground.latlonbox.west = bounds.left
    kml_path = "temp_kml.kml"
    kml.save(kml_path)

    # Step 4: Package the PNG and KML into a KMZ file
    with zipfile.ZipFile(kmz_path, 'w') as kmz:
        kmz.write(png_path)
        kmz.write(kml_path)

    # Optionally delete temporary files (PNG and KML)
    import os
    os.remove(png_path)
    os.remove(kml_path)


def readGeotiff(fileName):
    '''
    Reads a geotiff file
    Outputs the array and lat/lon bounds
    '''
    from osgeo import gdal
    dataset = gdal.Open(fileName, gdal.GA_ReadOnly)
    array = dataset.ReadAsArray()
    
    geotransform = dataset.GetGeoTransform()
    x_min = geotransform[0]
    y_max = geotransform[3]
    x_max = x_min + geotransform[1] * dataset.RasterXSize
    y_min = y_max + geotransform[5] * dataset.RasterYSize
    
    # Optional: If the image is rotated or sheared, you may need to calculate the bounds differently
    # Check the values of geotransform[2] and geotransform[4] to handle rotations/shearing if necessary
    dataset = None
    return array,x_min,x_max,y_min,y_max


def reGridStack(inputImageStack,outputX,outputY):
    '''
    inputImageStack: (m-rows,n-colums,k-stack)
    outputX: output coordinates (probably lon_ifg) (m-out,n-out)
    outputY: output coordiantes (probably lat_ifg) (m-out,n-out)
    refer to https://stackoverflow.com/questions/20915502/speedup-scipy-griddata-for-multiple-interpolations-between-two-irregular-grids
    '''
    import scipy.interpolate as spint
    import scipy.spatial.qhull as qhull
    import itertools
    def interp_weights(xy, uv,d=2):
        tri = qhull.Delaunay(xy)
        simplex = tri.find_simplex(uv)
        vertices = np.take(tri.simplices, simplex, axis=0)
        temp = np.take(tri.transform, simplex, axis=0)
        delta = uv - temp[:, d]
        bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
        return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))
    
    def interpolate(values, vtx, wts):
        return np.einsum('nj,nj->n', np.take(values, vtx), wts)
    
    m,n,k = inputImageStack.shape
    mi,ni = outputX.shape    
    [Y,X]   =np.meshgrid(np.linspace(0,1,n),np.linspace(0,2,m))
    [Yi,Xi] =np.meshgrid(np.linspace(0,1,ni),np.linspace(0,2,mi))
    
    xy=np.zeros([X.shape[0]*X.shape[1],2])
    xy[:,0]=Y.ravel()
    xy[:,1]=X.ravel()
    uv=np.zeros([Xi.shape[0]*Xi.shape[1],2])
    uv[:,0]=Yi.ravel()
    uv[:,1]=Xi.ravel()
    
    #Computed once and for all !
    vtx, wts = interp_weights(xy, uv)
    outputStack = np.zeros((mi,ni,k))

    for kk in np.arange(k):
        values=inputImageStack[:,:,kk]
        valuesi=interpolate(values.ravel(), vtx, wts)
        valuesi=valuesi.reshape(Xi.shape[0],Xi.shape[1])
        outputStack[:,:,kk] = valuesi
        
    return outputStack
    

def improfile(z, x0, y0, x1, y1):
    """
    Get a profile
    Captures 1d profile values from point A to B specified as indices. 
    Inputs: im, x0, y0, x1, y1
    Outputs: z (vector of values along the profile) 
    """
    length = int(np.hypot(x1-x0, y1-y0))
    x, y = np.linspace(x0, x1, length), np.linspace(y0, y1, length)
    
    # Extract the values along the line
    zi = z[y.astype(int), x.astype(int)]
    return zi

def ll2pixel(lon_ifg, lat_ifg, lon, lat):
    """
    Output the pixels (radar coords) given the lat/lon matrices and lat/lon points.
    input: lon, lat
    output: y, x
    
    """
    if np.isscalar(lon) and np.isscalar(lat):
        if np.nanmean(lon_ifg) * lon < 0:
            print('WARNING: you may need to subtract 360')
        
        a = abs(lat_ifg - lat)
        b = abs(lon_ifg - lon)
        c = a + b
        # Replace NaN values with a large number so they don't interfere with finding the minimum
        c_no_nan = np.where(np.isnan(c), np.inf, c)
        # Find the absolute difference from zero
        abs_diff = np.abs(c_no_nan)
        # Find the index of the minimum value in the flattened array
        min_index_flat = np.nanargmin(abs_diff)
        # Convert the flat index to 2D index
        y, x = np.unravel_index(min_index_flat, c.shape)
        return y,x
    

# phase elevation model
def phaseElev(img, hgt,msk, ymin, ymax, xmin, xmax,makePlot=True):
    '''
    Take the ifg or rate map and the dem and mask and outputs phs/elev dependence.
    Use ymin, xmin, etc. if you want to only use a subset of the image.
      otherwise, put the image len/width for those values.
    '''
#     img[np.isnan(img)] = 0
    
#     hgt[np.isnan(hgt)] = 0
    p = img[ymin:ymax, xmin:xmax].copy()
    z = hgt[ymin:ymax, xmin:xmax].copy()

        
    p = p[msk[ymin:ymax, xmin:xmax]!=0] 
    z = z[msk[ymin:ymax, xmin:xmax]!=0]
    G = np.vstack([z.ravel(), np.ones((len(z.ravel()),1)).flatten()]).T
    Gg = np.dot( np.linalg.inv(np.dot(G.T,G)), G.T)
    moda = np.dot(Gg,p.ravel())
    phs_model = moda[0] * hgt.ravel() + moda[1]
    phs_model = phs_model.reshape(img.shape)
    
    elevs = np.arange(z.min(),z.max(),10)
    bestLine =elevs*moda[0] + moda[1]
    
    if makePlot:
        from matplotlib import pyplot as plt
        plt.figure(figsize=(10,5));plt.scatter(z.ravel(),p.ravel(),.1,color='black')
        plt.plot(elevs,bestLine,color='red')
        plt.title('Phase-elevation dependence')
        plt.xlabel('Elevation (m)')
        plt.ylabel('Displacement (mm)')
        plt.savefig('Figs/phs_elev.png',dpi=300)
    return phs_model


def px2ll(x, y, lon_ifgm,lat_ifgm):
    lon = lon_ifgm[y,x]
    lat = lat_ifgm[y,x]
    return lon,lat


def fitLong(image,order,mask):
    from astropy.convolution import Gaussian2DKernel,convolve
    kernel = Gaussian2DKernel(x_stddev=1) # For smoothing and nan fill
    # image = convolve(image,kernel)
    # image[np.isnan(image)] = 0

    image[~mask] = np.nan
    
    ny,nx = image.shape
    X,Y = np.meshgrid(range(nx),range(ny))
    X1,Y1 = X.ravel(),Y.ravel()
    image_raveled = image.ravel()
    # Create a mask for valid (non-nan) data points
    valid_mask = ~np.isnan(image_raveled)
    
    # Select valid data points
    X1_valid, Y1_valid = X1[valid_mask], Y1[valid_mask]
    G_valid = np.array([np.ones(X1_valid.shape), X1_valid, Y1_valid]).T
    image_valid = image_raveled[valid_mask]
    
    if order == 1:  # Plane
        # Perform the least squares fit on the valid data points only
        Gg_valid = np.dot(np.linalg.pinv(G_valid.T @ G_valid), G_valid.T)
        mod = np.dot(Gg_valid, image_valid)
        
        # Create synthetic image using the model
        # Note that we calculate synth using all X1, Y1 points, not just the valid ones
        synth = mod[0] + mod[1] * X1 + mod[2] * Y1
        synth = synth.reshape(ny, nx)
    
        # Where the original image was nan, set the synthetic image to nan as well
        # synth[np.isnan(image)] = np.nan
            
    if order==2: # Quadratic
        G_valid = np.array([np.ones(X1_valid.shape), X1_valid, Y1_valid,X1_valid**2, Y1_valid**2]).T
        Gg_valid = np.dot(np.linalg.pinv(G_valid.T @ G_valid), G_valid.T)
        mod   = np.dot(Gg_valid,image_valid)
        synth = mod[0] + mod[1] * X1 + mod[2] * Y1 + mod[3] * X1**2 + mod[4] * Y1**2 
        synth = synth.reshape(ny,nx)

    if order==3:
        G  = np.array([np.ones((len(X1),)), X1, Y1, X1**2, Y1**2,X1**3, Y1**3]).T
        Gg = np.dot( np.linalg.inv(np.dot(G.T,G)), G.T)
        mod   = np.dot(Gg,image.ravel())
        synth = mod[0] + mod[1] * X1 + mod[2] * Y1 + mod[3] * X1**2 + mod[4] * Y1**2 + mod[5] * X1**3 + mod[6] * Y1**3
        synth = synth.reshape(ny,nx)

    if order==4:
        G  = np.array([np.ones((len(X1),)), X1, Y1, X1**2, Y1**2,X1**3, Y1**3,X1**4, Y1**4]).T
        Gg = np.dot( np.linalg.inv(np.dot(G.T,G)), G.T)
        mod   = np.dot(Gg,image.ravel())
        synth = mod[0] + mod[1] * X1 + mod[2] * Y1 + mod[3] * X1**2 + mod[4] * Y1**2 + mod[5] * X1**3 + mod[6] * Y1**3 + mod[7] * X1**4 + mod[8] * Y1**4
        synth = synth.reshape(ny,nx)


    return synth



def json2bbox(file):
    '''
    Takes a json file (geojson) and returns the bounds of the associated rectangle
    You can make a geojson file here: http://geojson.io/
    (minlon, minlat, maxlon, maxlat). This is the format for stac
    '''
    import json
    import numpy as np
    f = open(file,)
    coords = json.loads(f.read())['features'][0]['geometry']['coordinates'][0]
    lons=[]
    lats=[]
    
    for coord in coords:
        lons.append(coord[0])
        lats.append(coord[1])
        
    bbox = [np.min(lons), np.min(lats), np.max(lons), np.max(lats)]    
    return bbox,lons,lats


import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def struct_fun(data, ny, nx, tot=600, lengthscale=600, plot_flag=0, binwidth=20, fun=None):
    """
    Calculate structure function from an unwrapped interferogram matrix (data).
    nan values in data are okay.
    
    Parameters:
    - data: 2D NumPy array representing the unwrapped interferogram.
    - ny, nx: dimensions of data.
    - tot: total number of points to generate.
    - lengthscale: scale of the structure function.
    - plot_flag: if nonzero, plot the results.
    - binwidth: width of bins for averaging.
    - fun: function to fit ('exp' or 'spherical').
    
    Returns:
    - sqrt(S/2): square root of half the structure function values.
    - S2: unused, set to 0 for compatibility.
    - dists: distances for each pair of points.
    - allnxy: array of the same size as dists, for compatibility.
    - yd: square root of half the binned structure function values.
    - yd_std: unused, set to 0 for compatibility.
    - xd: distances corresponding to the binned values.
    - sf_fit: fitted function values.
    """
    
    xx, yy = np.arange(nx), np.arange(ny)
    X, Y = np.meshgrid(xx, yy, sparse=False, indexing='ij')
    
    xd, yd = np.meshgrid(
        [0, 1, 2, 5, 10, 15, 20, 25, 35, (lengthscale - binwidth), lengthscale],
        [-lengthscale, (-lengthscale + binwidth), -35, -25, -20, -15, -10, -5, -2, -1, 0, 1, 2, 5, 10, 15, 20, 25, 35, (lengthscale - binwidth), lengthscale],
        sparse=False, indexing='ij'
    )
    
    tx, ty = np.random.randint(1, lengthscale, size=(2, tot))
    ty[::2] = -ty[::2]
    
    q = np.vstack([tx, ty]).T
    _, ids = np.unique(q, axis=0, return_index=True)
    tx, ty = tx[ids].astype(int), ty[ids].astype(int)
    
    tx, ty = np.append(tx, xd.flatten()), np.append(ty, yd.flatten())
    
    aty = np.abs(ty)
    S = np.empty(len(tx))
    allnxy = np.empty(len(tx))
    
    for i in range(len(tx)):
        if ty[i] >= 0:
            A = data[1 : ny - ty[i], tx[i] : nx - 1]
            B = data[ty[i] : ny - 1, 1 : nx - tx[i]]
        else:
            A = data[aty[i] : ny - 1, tx[i] : nx - 1]
            B = data[1 : ny - aty[i], 1 : nx - tx[i]]
        C = A - B  # Difference, nan values will propagate
        
        # Calculating square of differences and then mean while ignoring nan
        S[i] = np.nanmean(C**2)
        allnxy[i] = np.count_nonzero(~np.isnan(C.flatten()))  # Count non-nan entries for weights
    
    
    dists = np.sqrt(tx**2 + ty**2)
    bins = np.arange(0, dists.max(), binwidth, dtype=int)
    S_bins, dist_bins = [], []
    
    for bin_min in bins:
        bin_ids = np.where((dists < (bin_min + binwidth)) & (dists >= bin_min))[0]
        w = allnxy[bin_ids]
        if len(w) == 0:
            S_bins.append(np.nan)
            dist_bins.append(np.nan)
        else:
            S_bins.append(np.average(S[bin_ids], weights=w))
            dist_bins.append(np.nanmean(dists[bin_ids]))
    

    
    xd = np.asarray(dist_bins)
    yd = np.sqrt(np.asarray(S_bins, dtype=np.float32) / 2)
    yd[np.isnan(yd)] = 0
    
    sf_fit = np.zeros_like(xd)  # Default to zero if no fitting function specified

    if fun == 'exp':
        mask = ~np.isnan(xd) & ~np.isnan(yd)  # Ensure both xd and yd are not nan
        sf_fit = fit_log(xd[mask], yd[mask]) if np.any(mask) else np.zeros_like(xd)
    elif fun == 'spherical':
        mask = ~np.isnan(xd) & ~np.isnan(yd)
        sf_fit = fit_spherical(xd[mask], yd[mask]) if np.any(mask) else np.zeros_like(xd)
    else:
        sf_fit = np.zeros_like(xd)
        
    
    if plot_flag:
        plt.figure(figsize=(14, 10))
        plt.subplot(221)
        plt.title("Image")
        plt.imshow(data,vmin=-10,vmax=10)
        plt.colorbar()
        
        plt.subplot(222)
        plt.title("sqrt(S) vs. position")
        plt.scatter(tx, ty, c=np.sqrt(S))
        plt.scatter(-tx, -ty, c=np.sqrt(S))
        plt.xlabel('east')
        plt.ylabel('north')
        plt.colorbar()
        
        plt.subplot(212)
        plt.title("S vs. distance, colored by num points")
        plt.scatter(dists[1:], np.sqrt(S[1:]/2), c=allnxy[1:])

        # Ensure that dist_bins and sf_fit are not empty or all nan
        if len(dist_bins) > 0 and not np.all(np.isnan(sf_fit)):
            valid_mask = ~np.isnan(dist_bins) & ~np.isnan(sf_fit)
            plt.scatter(np.array(dist_bins)[valid_mask], np.sqrt(np.array(S_bins)[valid_mask]/2), label='Binned Data', alpha=0.7)
            plt.plot(np.array(dist_bins)[valid_mask], sf_fit[valid_mask], 'r-', label='Best Fit')
        else:
            print("No valid data for plotting best fit function.")
        
        plt.xlabel('Distance')
        plt.ylabel('sqrt(S/2)')

        plt.colorbar()
        plt.show()
    
    return np.sqrt(S / 2), 0, dists, allnxy, yd, 0, xd, sf_fit

def fit_log(xd, yd):
    """Fit logarithmic function."""
    def fit_function(x, a, b, c):
        return a * np.log(b * x) + c
    popt, _ = curve_fit(fit_function, xd, yd)
    return fit_function(xd, *popt)

def fit_spherical(xd, yd):
    """Fit spherical model."""
    def spherical(x, a, b):
        return b * (1.5 * x / a - 0.5 * (x / a) ** 3.0)
    popt, _ = curve_fit(spherical, xd, yd)
    return spherical(xd, *popt)





def estimate_dem_error(ts0, G0, tbase, date_flag=None, phase_velocity=False):
    """Estimate DEM error with least square optimization.
    Parameters: ts0            - 2D np.array in size of (numDate, numPixel), original displacement time-series
                G0             - 2D np.array in size of (numDate, numParam), design matrix in [G_geom, G_defo]
                tbase          - 2D np.array in size of (numDate, 1), temporal baseline
                date_flag      - 1D np.array in bool data type, mark the date used in the estimation
                phase_velocity - bool, use phase history or phase velocity for minimization
    Returns:    delta_z        - 2D np.array in size of (1,       numPixel) estimated DEM residual
                ts_cor         - 2D np.array in size of (numDate, numPixel),
                                    corrected timeseries = tsOrig - delta_z_phase
                ts_res         - 2D np.array in size of (numDate, numPixel),
                                    residual timeseries = tsOrig - delta_z_phase - defModel
    Example:    delta_z, ts_cor, ts_res = estimate_dem_error(ts, G, tbase, date_flag)
    """
    import scipy
    if len(ts0.shape) == 1:
        ts0 = ts0.reshape(-1, 1)
    if date_flag is None:
        date_flag = np.ones(ts0.shape[0], np.bool_)

    # Prepare Design matrix G and observations ts for inversion
    G = G0[date_flag, :]
    ts = ts0[date_flag, :]
    if phase_velocity:
        tbase = tbase[date_flag, :]
        G = np.diff(G, axis=0) / np.diff(tbase, axis=0)
        ts = np.diff(ts, axis=0) / np.diff(tbase, axis=0)

    # Inverse using L-2 norm to get unknown parameters X
    # X = [delta_z, constC, vel, acc, deltaAcc, ..., step1, step2, ...]
    # equivalent to X = np.dot(np.dot(np.linalg.inv(np.dot(G.T, G)), G.T), ts)
    #               X = np.dot(np.linalg.pinv(G), ts)
    X = scipy.linalg.lstsq(G, ts, cond=1e-15)[0]

    # Prepare Outputs
    delta_z = X[0, :]
    ts_cor = ts0 - np.dot(G0[:, 0].reshape(-1, 1), delta_z.reshape(1, -1))
    ts_res = ts0 - np.dot(G0, X)

    # for debug
    debug_mode = False
    if debug_mode:
        import matplotlib.pyplot as plt
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, figsize=(8, 8))
        ts_all = np.hstack((ts0, ts_res, ts_cor))
        ymin = np.min(ts_all)
        ymax = np.max(ts_all)
        ax1.plot(ts0, '.');           ax1.set_ylim((ymin, ymax)); ax1.set_title('Original  Timeseries')
        ax2.plot(ts_cor, '.');        ax2.set_ylim((ymin, ymax)); ax2.set_title('Corrected Timeseries')
        ax3.plot(ts_res, '.');        ax3.set_ylim((ymin, ymax)); ax3.set_title('Fitting Residual')
        ax4.plot(ts_cor-ts_res, '.'); ax4.set_ylim((ymin, ymax)); ax4.set_title('Fitted Deformation Model')
        plt.show()

    return delta_z, ts_cor, ts_res


def viewIFGstack(flip=True,chain=True):
    ''' look at all the ifgs with napari'''
    
    import numpy as np
    import isce.components.isceobj as isceobj
    import napari
    
    ps = np.load('./ps.npy',allow_pickle=True).all()
    
    if chain:
        pairs = ps.pairs
    else:
        pairs = ps.pairs2
    
    stack = np.zeros((len(pairs),ps.nyl,ps.nxl))
    for ii in range(len(pairs)):
        p = pairs[ii]
        f = './merged/interferograms/' + p + '/fine_lk_filt.int'
        intImage = isceobj.createIntImage()
        intImage.load(f + '.xml')
        ifg = intImage.memMap()[:,:,0] 
        ifgc = np.angle(ifg)
        if flip:
            stack[ii,:,:] = np.flipud(ifgc)
        else:
            stack[ii,:,:] = ifgc
    viewer = napari.view_image(stack,colormap='RdYlBu')


def viewUNWstack(flip=True,chain=True):
    ''' look at all the ifgs with napari'''
    
    import numpy as np
    import isce.components.isceobj as isceobj
    import napari
    import isce
    
    ps = np.load('./ps.npy',allow_pickle=True).all()
    gam = np.load('./Npy/gam.npy')
    if chain:
        pairs = ps.pairs
    else:
        pairs = ps.pairs2
    
    
    
    stack = np.zeros((len(pairs),ps.nyl,ps.nxl))
    for ii in range(len(pairs)):
        p = pairs[ii]
        f = './merged/interferograms/' + p + '/filt.unw'
        intImage = isceobj.createImage()
        intImage.dataType='FLOAT'
        intImage.load(f + '.xml')
        unw = intImage.memMap()[:,:,0]
        unw = unw.copy()
        unw[gam==0] = 0
        if flip:
            stack[ii,:,:] = np.flipud(unw)
        else:
            stack[ii,:,:] = unw
    viewer = napari.view_image(stack,colormap='RdYlBu')

    # fig,ax = plt.subplots(4,8)
    # kk=0
    # for a in ax.ravel():
    #     a.imshow(stack[kk,:,:])
    #     a.axes.xaxis.set_visible(False)
    #     a.axes.yaxis.set_visible(False)
    # plt.tight_layout()


def viewCORstack(flip=True,chain=True):
    ''' look at all the ifgs with napari'''
    
    import numpy as np
    import isce.components.isceobj as isceobj
    import napari
    
    ps = np.load('./ps.npy',allow_pickle=True).all()
    gam = np.load('./Npy/gam.npy')
    
    if chain:
        pairs = ps.pairs
    else:
        pairs = ps.pairs2
        
    stack = np.zeros((len(pairs),ps.nyl,ps.nxl))
    for ii in range(len(pairs)):
        p = pairs[ii]
        f = './merged/interferograms/' + p + '/cor.r4'
        intImage = isceobj.createImage()
        intImage.dataType='FLOAT'
        intImage.load(f + '.xml')
        unw = intImage.memMap()[:,:,0]
        unw = unw.copy()
        unw[gam==0] = 0
        if flip:
            stack[ii,:,:] = np.flipud(unw)
        else:
            stack[ii,:,:] = unw
    viewer = napari.view_image(stack[0,:,:],colormap='jet')
    
    
def getUNW(pair):
    import isce.components.isceobj as isceobj
    gam = np.load('Npy/gam.npy')
    f = './merged/interferograms/' + pair + '/filt.unw'
    intImage = isceobj.createImage()
    intImage.dataType='FLOAT'
    intImage.load(f + '.xml')
    unw = intImage.memMap()[:,:,0]
    unw = unw.copy()
    unw[gam==0] = np.nan
    return unw

def getConCom(msk, minimumPixelsInRegion=1000):
    '''
    Takes a binary input (like a mask) as input and outputs labels for
    regions greater than the given minimum pixels.
    '''
    
    import cv2
    ratesu8 = (msk*255).astype(np.uint8)
    num_labels, labels = cv2.connectedComponents(ratesu8)
    
    npix = []
    for ii in range(num_labels):
        npix.append(len(np.where(labels==ii)[0]))
    npix = np.asarray(npix)
    # concom = np.where(npix>minimumPixelsInRegion)[0]
    
    newLabels = np.zeros(msk.shape)
    
    for ii in range(len(npix)):
        lab = npix[ii]
        newLabels[labels==lab] = ii
        
    return labels



def coregister(img1,img2):
    """
    Coregister two images
    
    inputs
        img1: reference image you want to align img2 to
        img2: image you want to be aligned 
    outputs:
        img2_coreg: the coregistered version of img2
        H: the homography matrix
        
    Based on this tutorial:
        https://www.sicara.fr/blog/2019-07-16-image-registration-deep-learning
    
    The KAZE algorithm is written up here (AKAZE is a faster version of that):
        http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.304.4980&rep=rep1&type=pdf
        
    Written: 4/20/2022 @ 4:20:69

    """
    
    from matplotlib import pyplot as plt
    import cv2 as cv
    
    # img1 = cv.imread('/home/km/Pictures/Murray_175px.jpg', cv.IMREAD_GRAYSCALE)  # referenceImage
    # img2 = cv.imread('/home/km/Pictures/Murray2px.jpg', cv.IMREAD_GRAYSCALE)  # sensedImage
    
    img1_8 = cv.normalize(src=img1, dst=None, alpha=0, beta=255, norm_type=cv.NORM_MINMAX, dtype=cv.CV_8U)
    img2_8 = cv.normalize(src=img2, dst=None, alpha=0, beta=255, norm_type=cv.NORM_MINMAX, dtype=cv.CV_8U)
    
    
    # Initiate AKAZE detector
    akaze = cv.AKAZE_create()
    # Find the keypoints and descriptors with SIFT
    kp1, des1 = akaze.detectAndCompute(img1_8, None)
    kp2, des2 = akaze.detectAndCompute(img2_8, None)
    
    # BFMatcher with default params
    bf = cv.BFMatcher()
    matches = bf.knnMatch(des1, des2, k=2)
    
    # Apply ratio test
    good_matches = []
    for m,n in matches:
        if m.distance < 0.75*n.distance:
            good_matches.append([m])
            
    # Draw matches
    img3 = cv.drawMatchesKnn(img1_8,kp1,img2_8,kp2,good_matches,None,flags=cv.DrawMatchesFlags_NOT_DRAW_SINGLE_POINTS)
    cv.imwrite('matches.jpg', img3)
    
    # Select good matched keypoints
    ref_matched_kpts = np.float32([kp1[m[0].queryIdx].pt for m in good_matches])
    sensed_matched_kpts = np.float32([kp2[m[0].trainIdx].pt for m in good_matches])
    
    # Compute homography
    H, status = cv.findHomography(sensed_matched_kpts, ref_matched_kpts, cv.RANSAC,5.0) 
    
    # Warp image
    # warped_image = cv.warpPerspective(img2_8, H, (img2_8.shape[1], img2_8.shape[0]))
    img2_coreg = cv.warpPerspective(img2, H, (img1.shape[1], img1.shape[0]))
    
    fig,ax = plt.subplots(2,2)
    ax[0,0].imshow(img1);ax[0,0].set_title('img1')
    ax[0,1].imshow(img2);ax[0,1].set_title('img2')
    ax[1,0].imshow(img2_coreg);ax[1,0].set_title('img2_coreg')
    ax[1,1].imshow(img1-img2_coreg,vmin=-100,vmax=100,cmap='RdBu_r');ax[1,1].set_title('img1 - img2_coreg')
    
    meanRes = np.nanmean(abs(img1-img2_coreg))
    print('Mean Residual: ' + str(meanRes))

    return img2_coreg, H
    
def show(img,title=None,cmap='magma'):
    """
    just plots the image so you don't have to type as much. For quickly viewing.
    """
    from matplotlib import pyplot as plt
    plt.figure()
    plt.imshow(img,cmap=cmap)
    plt.show()
    if title:
        plt.title(title)


def gaussian_kernel(Sx, Sy, sig_x, sig_y):
    if np.mod(Sx,2) == 0:
        Sx = Sx + 1

    if np.mod(Sy,2) ==0:
            Sy = Sy + 1

    x,y = np.meshgrid(np.arange(Sx),np.arange(Sy))
    x = x + 1
    y = y + 1
    x0 = (Sx+1)/2
    y0 = (Sy+1)/2
    fx = ((x-x0)**2.)/(2.*sig_x**2.)
    fy = ((y-y0)**2.)/(2.*sig_y**2.)
    k = np.exp(-1.0*(fx+fy))
    a = 1./np.sum(k)
    k = a*k
    return k

def convolve(data, kernel):
    import cv2
    R = cv2.filter2D(data.real,-1,kernel)
    Im = cv2.filter2D(data.imag,-1,kernel)

    return R + 1J*Im


def butter(img,wavelength,nyq_freq=0.5,order=2):
    from scipy import signal
    
    #wavelength = 150  #Get rid of stuff happening over these spatial scales or longer
    #nyq_freq = .5 # 1 sample/pixel /2
    #order = 2
    cutoff_frequency = 1/(wavelength*2)
    
    def butterLow(cutoff, critical, order):
        normal_cutoff = float(cutoff) / critical
        b, a = signal.butter(order, normal_cutoff, btype='lowpass')
        return b, a
    
    def butterFilter(data, cutoff_freq, nyq_freq, order):
        b, a = butterLow(cutoff_freq, nyq_freq, order)
        y = signal.filtfilt(b, a, data)
        return y
    
    filt = butterFilter(img, cutoff_frequency, nyq_freq, order)
    
    return filt

def write_xml(filename,width,length,bands,dataType,scheme):
    import isce.components.isceobj as isceobj    
    img=isceobj.createImage()
    img.setFilename(filename)
    img.setWidth(width)
    img.setLength(length)
    img.setAccessMode('Read')
    img.bands=bands
    img.dataType=dataType
    img.scheme = scheme
    img.renderHdr()
    img.renderVRT()
    return

def writeISCEimg(img,outName,nBands,width,length,dtype):
    '''
    quick way to write a file with an xml file. Automatically checks datatype
    img: image to write to file with an xml
    outname: name of output file not including the .xml extension
    dtype: 'Float' or 'Complex' or 'int'
    ''' 
    import isce.components.isceobj as isceobj
    fidc=open(outName,"wb")
    fidc.write(img)
    #write out an xml file for it
    out = isceobj.createIntImage() # Copy the interferogram image from before
    
    if dtype=='Float':
        img = np.asarray(img,dtype=np.float32)
        out.dataType = 'FLOAT'
    elif dtype=='int':
        img = np.asarray(img,dtype=int)
        out.dataType = 'BYTE'
    elif dtype=='Complex':
        img = np.asarray(img,dtype=np.complex64)
        out.dataType = 'CFLOAT'
    
    out.scheme = 'BIP'
    out.bands = 1
    out.filename = outName
    out.width = width
    out.length = length
    out.dump(outName + '.xml') # Write out xml
    out.renderHdr()
    out.renderVRT()


def filtAndCoherence(infileIFG,filtFileOut,corFileOut,filterStrength):
    '''
    Runs filtering and coherence estimation
    '''
    
    import FilterAndCoherence as fc
    if filterStrength <= 0:
        print('Skipping filtering because filterStrength is 0')
    else:
        fc.runFilter(infileIFG, filtFileOut, filterStrength)    
    
    fc.estCoherence(filtFileOut, corFileOut)

def unwrap_snaphu(int_file,cor_file,unw_file,length, width,rlks,alks):
    '''
    Inputs:
        intfile
        corfile
        length,width: length and width of intfile
        rlks,alks: give rlks and alks just to record it in the xml file
    Outputs:
        unwfile: writes unw image to this file

    '''
    
    from contrib.Snaphu.Snaphu import Snaphu

    altitude = 800000.0
    earthRadius = 6371000.0
    wavelength = 0.056
    defomax = 2
    maxComponents = 20
    
    snp = Snaphu()
    snp.setInitOnly(False)
    snp.setInput(int_file)
    snp.setOutput(unw_file)
    snp.setWidth(width)
    snp.setCostMode('SMOOTH')
    snp.setEarthRadius(earthRadius)
    snp.setWavelength(wavelength)
    snp.setAltitude(altitude)
    snp.setCorrfile(cor_file)
    snp.setInitMethod('MCF')
    # snp.dumpConnectedComponents(True)
    snp.setMaxComponents(maxComponents)
    snp.setDefoMaxCycles(defomax)
    snp.setRangeLooks(rlks)
    snp.setAzimuthLooks(alks)
    snp.setCorFileFormat('FLOAT_DATA')
    snp.prepare()
    snp.unwrap()
    write_xml(unw_file, width, length, 2 , "FLOAT",'BIL')
    return

def getWaterMask(DEMfilename, lon_filename, lat_filename, outputfilename):
    import createWaterMask as wm
    bbox = wm.dem2bbox(DEMfilename)
    geo_file = wm.download_waterMask(bbox, DEMfilename, fill_value=0)
    wm.geo2radar(geo_file, outputfilename, lat_filename, lon_filename)
    
def tsFilt(alld, dec_year, N=5, desiredPeriod = 1):
    '''
    Temporal filter
    Inputs:
        alld: len(time) X n X m matrix
        N: Filter order
        desiredPeriod: roughly the cutoff period in years. (Anything shorter 
            than this value will be filtered out).
    Output:
        alldFilt
   
    Wn is the Cutoff frequency between 0 and 1.  0 is infinitely smooth and 1 is the original. 
        this is the frequency multiplied by the nyquist rate. 
        if we have 25 samples per year, then the nyquist rate would be ~12hz. So if we make Wn=.5
        we will have filtered to 6hz (letting signals with wavelengths of 2 months or longer).
        If we make wn=1/12 then we will filter to 1hz (letting only signals with wavelengths of 1 year).
    '''
    import scipy.signal as signal

    dec_year = np.asarray(dec_year)
    samplesPerYear = len(dec_year) / (dec_year.max()-dec_year.min())
    nyquistRate = samplesPerYear/2 #this is the highest freq we can resolve with our sampling rate
    Wn = 1/(desiredPeriod * nyquistRate)
    B, A = signal.butter(N, Wn, output='ba')
    
    alldFilt = signal.filtfilt(B,A, alld,axis=0)
    alldFilt[alldFilt==0]=np.nan
    
    return alldFilt

def tsFilt1d(ts, dec_year, N=5, desiredPeriod = 1):
    '''
    Temporal filter
    Inputs:
        ts: 1d time series corresponding to dec_year dates
        N: Filter order
        desiredPeriod: roughly the cutoff period in years. (Anything shorter 
            than this value will be filtered out).
    Output:
        alldFilt
   
    Wn is the Cutoff frequency between 0 and 1.  0 is infinitely smooth and 1 is the original. 
        this is the frequency multiplied by the nyquist rate. 
        if we have 25 samples per year, then the nyquist rate would be ~12hz. So if we make Wn=.5
        we will have filtered to 6hz (letting signals with wavelengths of 2 months or longer).
        If we make wn=1/12 then we will filter to 1hz (letting only signals with wavelengths of 1 year).
    '''
    import scipy.signal as signal

    dec_year = np.asarray(dec_year)
    samplesPerYear = len(dec_year) / (dec_year.max()-dec_year.min())
    nyquistRate = samplesPerYear/2 #this is the highest freq we can resolve with our sampling rate
    Wn = 1/(desiredPeriod * nyquistRate)
    B, A = signal.butter(N, Wn, output='ba')
    
    tsF = signal.filtfilt(B,A, ts,axis=0)
    tsF[tsF==0]=np.nan
    
    return tsF


def getLOSvec(psi,theta):
    '''
    Given the flight azimuth and the los azimuth angle, return LOS vector
    theta should be the angle from ground normal
    psi should be the angle from EAST (X) of the flight direction angle. (not the look direction)
    
    In ISCE:
    psi is from the second band in the los.rdr file (and also in the az_lk.rdr file from PyPS)
    theta is from the second band of incLocal.rdr file (which is the angle from the surface normal vector)
    
    '''
    psi+=90 # Add 90 degrees to get the look direction from the flight direction
    losA=np.zeros((3,))
    losA[0] =  np.sin(np.deg2rad(theta)) * np.cos(np.deg2rad(psi))
    losA[1] =  np.sin(np.deg2rad(theta)) * np.sin(np.deg2rad(psi))
    losA[2] =  np.cos(np.deg2rad(theta))
    return losA
    
    
def invertVertHor(asc,des,psi_a,theta_a,psi_d,theta_d,smooth):
    '''
    damped least squares inversion of ascending/descending LOS data
       returns horizontal (East-west) and Vertical displacements
    IN:
        asc: value of ascending ifg at a pixel
        des: value of descending ifg at a pixel
        psi_a/d: azimuth direction of asc/des 
        theta_a/d: incidence angle of asc/des
        smooth: damping factor 
    OUT:
        vertHor: east-west and vertical deformation 
    '''
    
    losA = getLOSvec(psi_a,theta_a)
    losD = getLOSvec(psi_d,theta_d)

    o = np.concatenate((np.array([asc,des]).T,np.array([0,0])),axis=0)
    # Define unit basis vectors
    mx = np.array([1,0,0])
    # my = np.array([0,1,0])
    mz = np.array([0,0,1])
    
    # A_nozero = np.array([[np.dot(losA,mx), np.dot(losA,mz)],
    #               [np.dot(losD,mx), np.dot(losD,mz)]])
    
    A = np.array([[np.dot(losA,mx), np.dot(losA,mz)],
                  [np.dot(losD,mx), np.dot(losD,mz)],
                  [smooth,          0              ],
                  [0,               smooth         ]])
    
    # A = np.array([[np.dot(losA,mx), np.dot(losA,my), np.dot(losA,mz)],
    #               [np.dot(losD,mx),np.dot(losD,my),np.dot(losD,mz)],
    #               [smooth,              0,         0              ],
    #               [0,              smooth,          0              ],
    #               [0,              0,              smooth         ]])
    
    Aa = np.dot( np.linalg.inv(np.dot(A.T,A)), A.T)
    vertHor = np.dot(Aa,o)
    
    # For L-curve (doesn't really make an L though)
    # solution_norm = np.linalg.norm(vertHor)
    # res_norm = np.linalg.norm( np.dot(A_nozero,vertHor) - o_nozero)
    
    return vertHor #, solution_norm, res_norm

def geocode(filename):
    ''' Geocodes the filename and outputs geocoded image in that directory 
        with .geo. I think filename can be a list of multiple names.
    '''
    
    import geocodeIsce
    import glob
    # setupParams = np.load('setupParams.npy',allow_pickle=True).item()
    
    ps = np.load('./ps.npy',allow_pickle=True).all()
    
    minlon = ps.lon_ifg.min()
    maxlon = ps.lon_ifg.max()
    minlat = ps.lat_ifg.min()
    maxlat = ps.lat_ifg.max()
    dem = glob.glob('./DEM/*wgs84')[0]
    bbox1 = [minlat,maxlat, minlon,maxlon]
    
    class inpsArgs():
        prodlist = filename
        bbox = bbox1
        demfilename = dem
        reference = ps.workdir + '/reference'
        secondary = ps.workdir + '/reference'
        numberRangeLooks = ps.rlks
        numberAzimuthLooks = ps.alks
        
    geocodeIsce.runGeocode(inpsArgs, inpsArgs.prodlist, inpsArgs.bbox, inpsArgs.demfilename, is_offset_mode=False)

def geocodeKM(img,resolution,lon_ifg,lat_ifg, method='linear'):
    '''
    resolution in meters (pixel)
    
    '''    
    from scipy.interpolate import griddata 
    
    minlat = lat_ifg.min()
    minlon = lon_ifg.min()
    maxlat = lat_ifg.max()
    maxlon = lon_ifg.max()
    
    delta_latitude = maxlat-minlat
    delta_longitude = maxlon-minlon
    avgLat = (minlat+maxlat)/2
    
    lonDist = 1000*delta_longitude * (40000 * np.cos( np.deg2rad(avgLat) ) / 360)
    latDist = 111111 * delta_latitude  
    
    nx = int(lonDist//resolution)
    ny = int(latDist//resolution)
    xx = np.linspace(minlon,maxlon,nx)
    yy = np.linspace(minlat,maxlat,ny)
    XX,YY = np.meshgrid(xx,yy)
    valid_mask = ~np.isnan(lon_ifg) | ~np.isnan(lat_ifg)


    imgRegrid = griddata((lon_ifg[valid_mask].ravel(),lat_ifg[valid_mask].ravel()), img[valid_mask].ravel(), (XX,YY), method=method)
    imgRegrid = np.flipud(imgRegrid)
    
    return imgRegrid,XX,YY

def rad2cm(input_vec,wavelength=.056,output='cm'):
    '''
    input: number or array in radians
    wavelength: in meters
        Sentinel-1: 0.056
        Alos-1: 0.23
    output: mm, cm, m
    
    The factor of 4pi comes from the fact that the phase difference represents 
    a round-trip change in distance (to the ground and back), and there are 2pi 
    radians in a full cycle. So a fringe = wavelength/2, which is ~2.78 cm LOS 
    displacement.
    '''
    if output=='m':
        factor=1
    if output=='cm':
        factor=100
    if output=='mm':
        factor=1000
    output_cm = (input_vec*0.056/(4*np.pi))*factor
    
    return output_cm

def orderAxes(inputArray,nx,ny):
    '''  Rearrange axes order from small to big '''
    imShape = np.asarray(inputArray.shape)
    smaA = np.where(imShape==imShape.min())[0][0]
    inputArray = np.moveaxis(inputArray,smaA,0)
    imShape = np.asarray(inputArray.shape)
    bigA = np.where(imShape==imShape.max())[0][0]
    if nx>ny:
        inputArray = np.moveaxis(inputArray,bigA,2)
    else:
        inputArray = np.moveaxis(inputArray,bigA,1)
    return inputArray

def fitSine1d(t, signal):
    """
    Fits a signal with a model consisting of a cosine and a sine component using the least squares method. 
   
    Parameters:
    t (array-like): Time vector.
    signal (array-like): Signal to fit.
   
    Returns:
    A tuple containing the frequency, amplitude, and phase shift of the fitted cosine and sine components.
    """
    from scipy.optimize import curve_fit

    # Define the function to fit
    def signal_model(t, freq_cos, freq_sin, amplitude_cos, amplitude_sin, phase_shift_cos, phase_shift_sin):
        return amplitude_cos * np.cos(2 * np.pi * freq_cos * t + phase_shift_cos) + amplitude_sin * np.sin(2 * np.pi * freq_sin * t + phase_shift_sin)
   
    # Set some initial parameter values
    freq_cos = 1
    freq_sin = 2
    amplitude_cos = np.mean(abs(signal))
    amplitude_sin = np.mean(abs(signal))/2
    phase_shift_cos = 0.0
    phase_shift_sin = 0.0
   
    # Fit the data
    popt, pcov = curve_fit(signal_model, t, signal, p0=[freq_cos, freq_sin, amplitude_cos, amplitude_sin, phase_shift_cos, phase_shift_sin],maxfev=5000)
    freq_cos, freq_sin, amplitude_cos, amplitude_sin, phase_shift_cos, phase_shift_sin = popt
    modeled_signal = amplitude_cos * np.cos(2 * np.pi * freq_cos * t + phase_shift_cos) + amplitude_sin * np.sin(2 * np.pi * freq_sin * t + phase_shift_sin)

        # try:
        #     popt, _ = curve_fit(signal_model, valid_years, valid_interpolated_values, p0=initial_guess, maxfev=10000)
        #     modeled_signal2 = signal_model(decimal_years_vector, *popt)  # Use decimal_years_vector for full length
        #     freq_cos2, freq_sin2, amplitude_cos2, amplitude_sin2, phase_shift_cos2, phase_shift_sin2 = popt
        # except RuntimeError as e:
        #     print(f"Error fitting sine model for station {station_name}: {e}")
        #     modeled_signal = np.full_like(decimal_years_vector, np.nan)
        #     freq_cos, freq_sin, amplitude_cos, amplitude_sin, phase_shift_cos, phase_shift_sin = [np.nan] * 6
            


    return modeled_signal,freq_cos, freq_sin, amplitude_cos, amplitude_sin, phase_shift_cos, phase_shift_sin

def vector2raster(shapefile_fn, raster_sample_fn, raster_out_fn, feature):
    '''
    Converts vector data from a shape file to raster data in geotiff file.
    
    shapefile_fn: File name for the shapefile that has the vector data
    
    raster_sample_fn: An example geotiff raster that you want the output raster 
        to look like.
        
    raster_out_fn: name of the output raster geotiff 
    
    feature: The column in the geopandas dataframe that will be the values of 
        the raster. Needs to be numbers/floats
        
    Returns: raster
    
    '''
    # Convert polygons to raster
    import geopandas as gpd
    import rasterio
    from rasterio import features
    import cartopy.crs as ccrs
    from osgeo import gdal
    
    projection = ccrs.PlateCarree()
    dataframe = gpd.read_file(shapefile_fn).to_crs(projection)
    
    # Crop like this
    #dataframe = dataframe.cx[minlon:maxlon,minlat:maxlat]

    rst = rasterio.open(raster_sample_fn)
    meta = rst.meta.copy()
    meta.update(compress='lzw')

    with rasterio.open(raster_out_fn, 'w+', **meta) as out:
        out_arr = out.read(1)
        # this is where we create a generator of geom, value pairs to use in rasterizing
        shapes = ((geom,value) for geom, value in zip(dataframe.geometry, dataframe[feature].astype(np.float32)))
        burned = features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)
        out.write_band(1, burned)
        
    print('raster was saved to ' + raster_out_fn)
        
    ds = gdal.Open(raster_out_fn)
    raster = ds.GetVirtualMemArray()
    return raster 


def update_yaml_key(file_path, key, new_value):
    import re
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
            
            
def date_string_2_dec_year(date_strings):
    '''
    date_strings: array of 'YYYYMMDD'
    '''
    from datetime import date

    def is_leap_year(year):
        return (year % 4 == 0 and year % 100 != 0) or year % 400 == 0
    dec_year = []
    for d in date_strings:
        yr, mo, day = int(d[0:4]), int(d[4:6]), int(d[6:8])
        current_date = date(yr, mo, day)
        dt = current_date.toordinal()
        d0 = date(yr, 1, 1).toordinal()
        doy = dt - d0 + 1
        is_leap = is_leap_year(yr)
        days_in_year = 366 if is_leap else 365
        dec_year.append(yr + (doy - 1) / days_in_year)
    
    # If needed, convert dec_year to a NumPy array
    dec_year = np.array(dec_year)
    return dec_year


def invert_ifgs(unw_stack,pairs,dates):
    '''
    unw_stack: (nifgs,ny,nx)
    pairs: array of 'YYYYMMDD_YYYYMMDD'
    dates: array of 'YYYYMMDD'
    returns time series of displacements : array (ndates,ny,nx)
    '''
    
    # Make G matrix for dates inversion
    G = np.zeros((len(pairs)+1,len(dates)))# extra row of zeros to make first date zero for reference
    for ii,pair in enumerate(pairs):
        a = dates.index(pair[0:8])
        b = dates.index(pair[9:17])
        G[ii,a] = 1
        G[ii,b] = -1
    G[-1,0]=1
    Gg = np.dot( np.linalg.inv(np.dot(G.T,G)), G.T)
    
    nxl = unw_stack.shape[2]
    nyl = unw_stack.shape[1]
    
    # Do dates inversion
    ts_full = np.zeros((len(dates),nxl*nyl))
    for ii in np.arange(0,nyl-1): #iterate through rows
        tmp = np.zeros((len(pairs)+1,nxl))
        for jj,pair in enumerate(pairs): #loop through each ifg and append to alld 
            tmp[jj,:] = unw_stack[jj,ii,:]
        ts_full[:,ii*nxl:nxl*ii+nxl] = np.dot(Gg, tmp)
    
    ts_full = ts_full.reshape(len(dates),nyl,nxl)
    ts_full = rad2cm(ts_full,wavelength=.056,output='cm')
    print('converted to cm')
    return ts_full


def fill_bad_geom(geom_path):
    '''
    geom_path: path to geometryRadar.h5 ('./MintPy/inputs/geometryRadar.h5')
    
    This will just fill the bad values with max lat/lon so when you geocode, 
    all of the values will collapse to a single pixel in the corner.  
    
    '''
    import h5py


    ds = h5py.File(geom_path,'r+')
    lons = np.asarray(ds['longitude'])
    lats = np.asarray(ds['latitude'])
    
    # fig,ax = plt.subplots(2,1)
    # ax[0].imshow(lons);ax[0].set_title('lon')
    # ax[1].imshow(lats);ax[1].set_title('lat')
    mask_geom = lons!=0
        
    lon_ramp = fitLong(lons, 1,mask_geom)
    lat_ramp = fitLong(lats, 1,mask_geom)
    
    # fig,ax = plt.subplots(2,1)
    # ax[0].imshow(lon_ramp);ax[0].set_title('lon_ramp')
    # ax[1].imshow(lat_ramp);ax[1].set_title('lat_ramp')
    
    lon_diff = abs(lon_ramp-lons)
    lat_diff = abs(lat_ramp-lats)
    
    # fig,ax = plt.subplots(2,1)
    # ax[0].imshow(lon_diff);ax[0].set_title('lon_diff')
    # ax[1].imshow(lat_diff);ax[1].set_title('lat_diff')
    
    lons_corrected = lons.copy()
    lats_corrected = lats.copy()
    
    lons_corrected[lon_diff>.5] = lon_ramp[lon_diff>.5]
    lats_corrected[lat_diff>.5] = lat_ramp[lat_diff>.5]
    
    # fig,ax = plt.subplots(2,1)
    # ax[0].imshow(lons_corrected[3500:,2500:]);ax[0].set_title('lon corrected')
    # ax[1].imshow(lats_corrected[3500:,2500:]);ax[1].set_title('lat corrected')
    
    
    lons_single_fill = lons_corrected.copy()
    lats_single_fill = lats_corrected.copy()
    
    lons_single_fill[lon_diff>.5] = np.nanmax(lons_corrected)
    lats_single_fill[lat_diff>.5] = np.nanmax(lats_corrected)
    
    # fig,ax = plt.subplots(2,1)
    # ax[0].imshow(lons_single_fill);ax[0].set_title('lon lons_single_fill')
    # ax[1].imshow(lats_single_fill);ax[1].set_title('lat lats_single_fill')
    
    ds['longitude'][...] = lons_single_fill
    ds['latitude'][...] = lats_single_fill
    
    
    ds.close()
