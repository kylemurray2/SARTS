#!/usr/bin/env python3
"""
InSAR LOS velocity decomposition into vertical and horizontal components.
Resamples descending velocities to ascending grid before decomposition.
"""

import numpy as np
import h5py 
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from SARTS import util
from astropy.convolution import Gaussian2DKernel, convolve
from skimage.morphology import remove_small_objects
import os

def load_mintpy_file(filename, dataset_name):
    """Load data from MintPy HDF5 file"""
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")
        
    try:
        with h5py.File(filename, 'r') as ds:
            return np.asarray(ds[dataset_name])
    except KeyError:
        raise KeyError(f"Dataset '{dataset_name}' not found in {filename}")

def load_geometry(mintpy_dir):
    """Load geometry data from MintPy geometry file"""
    filename = os.path.join(mintpy_dir, 'inputs/geometryRadar.h5')
    with h5py.File(filename, 'r') as ds:
        return {
            'longitude': np.asarray(ds['longitude']),
            'latitude': np.asarray(ds['latitude']),
            'azimuth': np.asarray(ds['azimuthAngle']),
            'incidence': np.asarray(ds['incidenceAngle']),
            'waterMask': np.asarray(ds['waterMask'])
        }

def load_masks(mintpy_dir):
    """Load coherence masks"""
    masks = {}
    
    # Temporal coherence
    masks['temporal'] = load_mintpy_file(
        os.path.join(mintpy_dir, 'temporalCoherence.h5'),
        'temporalCoherence'
    )
    
    # Connected components
    masks['connected'] = load_mintpy_file(
        os.path.join(mintpy_dir, 'maskConnComp.h5'),
        'mask'
    )
    
    return masks

def decompose_velocities(vel_asc, vel_desc, geom_asc, geom_desc, region=None, smooth=0.1):
    """
    Decompose LOS velocities into vertical and horizontal components
    
    Parameters
    ----------
    vel_asc : 2D array
        Ascending velocity map
    vel_desc : 2D array
        Descending velocity map (resampled to ascending grid)
    geom_asc : dict
        Ascending geometry (azimuth, incidence angles)
    geom_desc : dict
        Descending geometry (azimuth, incidence angles)
    region : tuple, optional
        (ymin, ymax, xmin, xmax) for subsetting
    smooth : float, optional
        Smoothing parameter for inversion
        
    Returns
    -------
    east : 2D array
        East-west velocity component
    vert : 2D array
        Vertical velocity component
    """
    # Validate inputs
    if vel_asc.shape != vel_desc.shape:
        raise ValueError("Ascending and descending velocities must have same shape")
    
    if not all(key in geom_asc for key in ['azimuth', 'incidence']):
        raise ValueError("Missing required geometry fields for ascending data")
    
    if not all(key in geom_desc for key in ['azimuth', 'incidence']):
        raise ValueError("Missing required geometry fields for descending data")
    
    if region is None:
        ymin, ymax = 0, vel_asc.shape[0]
        xmin, xmax = 0, vel_asc.shape[1]
    else:
        ymin, ymax, xmin, xmax = region
        
    shape = vel_asc[ymin:ymax, xmin:xmax].shape
    npix = shape[0] * shape[1]
    vert_hor = np.zeros((npix, 2))
    
    # Prepare flattened arrays
    asc_flat = vel_asc[ymin:ymax, xmin:xmax].ravel()
    desc_flat = vel_desc[ymin:ymax, xmin:xmax].ravel()
    az_a = geom_asc['azimuth'][ymin:ymax, xmin:xmax].ravel()
    az_d = geom_desc['azimuth'][ymin:ymax, xmin:xmax].ravel()
    inc_a = geom_asc['incidence'][ymin:ymax, xmin:xmax].ravel()
    inc_d = geom_desc['incidence'][ymin:ymax, xmin:xmax].ravel()
    
    # Process each pixel
    for ii in range(npix):
        if ii % 200000 == 0:
            print(f"\rProcessing: {100 * ii/npix:.1f}%", end='', flush=True)
            
        if np.isnan(asc_flat[ii]) or np.isnan(desc_flat[ii]):
            vert_hor[ii,:] = np.nan
            continue
            
        vert_hor[ii,:] = util.invertVertHor(
            asc_flat[ii], desc_flat[ii],
            az_a[ii], inc_a[ii],
            az_d[ii], inc_d[ii],
            smooth
        )
    
    print("\rProcessing: 100%")
    
    east = vert_hor[:,0].reshape(shape)
    vert = vert_hor[:,1].reshape(shape)
    
    return east, vert

def apply_masks_and_filter(data, masks, water_mask, temporal_thresh=0.7, 
                          kernel_std=0.6, min_region_size=6000):
    """
    Apply masks and spatial filtering
    
    Parameters
    ----------
    data : 2D array
        Input velocity map
    masks : dict
        Dictionary containing temporal coherence and connected component masks
    water_mask : 2D array
        Water mask
    temporal_thresh : float, optional
        Threshold for temporal coherence mask
    kernel_std : float, optional
        Standard deviation for Gaussian kernel. Set to 0 for no filtering.
    min_region_size : int, optional
        Minimum size of connected regions to keep
        
    Returns
    -------
    filtered_data : 2D array
        Masked and optionally filtered data
    mask : 2D array
        Binary mask showing valid pixels
    """
    # Validate input shapes
    if not (data.shape == masks['temporal'].shape == masks['connected'].shape == water_mask.shape):
        raise ValueError("All input arrays must have the same shape")
    
    # Create combined mask
    mask = (masks['temporal'] > temporal_thresh) & (masks['connected'] > 0)
    
    # Remove small regions
    mask = remove_small_objects(mask, min_region_size)
    
    # Apply mask
    masked_data = data.copy()
    masked_data[~mask] = np.nan
    masked_data[water_mask == 0] = np.nan
    
    # Apply spatial filter only if kernel_std > 0
    if kernel_std > 0:
        kernel = Gaussian2DKernel(x_stddev=kernel_std)
        filtered_data = convolve(masked_data, kernel)
        return filtered_data, mask
    
    return masked_data, mask

def main():
    # Change these paths and parameters_______
    base_dir = '/d/HI'
    proc_dir = 'Oahu'
    MintPy_dir = 'MintPy_1_4'

    # Optionally define a small region to process
    region = None  # Subset region
    # region = (760, 1260, 2120, 2530)  # Subset region example

    kernel_std = 0  # Set to 0 for no filtering
    # _________________________________________

    asc_dir = os.path.join(base_dir, 'Asc', proc_dir)
    desc_dir = os.path.join(base_dir, 'Des', proc_dir)
    output_dir = os.path.join('/d/HI/verthorz', proc_dir)
    
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, 'Npy'), exist_ok=True)
    os.chdir(output_dir)
    
    # Load geometry and velocities
    asc_geom = load_geometry(os.path.join(asc_dir, MintPy_dir))
    desc_geom = load_geometry(os.path.join(desc_dir, MintPy_dir))
    
    # Load velocities (convert to mm/yr)
    vel_asc = load_mintpy_file(
        os.path.join(asc_dir, MintPy_dir, 'velocity.h5'), 
        'velocity'
    ) * 1000
    
    vel_desc = load_mintpy_file(
        os.path.join(desc_dir, MintPy_dir, 'velocity.h5'),
        'velocity'
    ) * 1000
    
    # Resample descending to ascending grid
    vel_desc_resamp = griddata(
        (desc_geom['longitude'].ravel(), desc_geom['latitude'].ravel()),
        vel_desc.ravel(),
        (asc_geom['longitude'], asc_geom['latitude']),
        method='linear'
    )
    
    # Also resample geometry
    desc_geom_resamp = {}
    for key in ['azimuth', 'incidence']:
        desc_geom_resamp[key] = griddata(
            (desc_geom['longitude'].ravel(), desc_geom['latitude'].ravel()),
            desc_geom[key].ravel(),
            (asc_geom['longitude'], asc_geom['latitude']),
            method='linear'
        )
    
    # Load masks
    asc_masks = load_masks(os.path.join(asc_dir, MintPy_dir))
    
    # Decompose velocities
    east, vert = decompose_velocities(
        vel_asc, vel_desc_resamp, 
        asc_geom, desc_geom_resamp,
        region=region
    )
    
    # Apply masks and filtering (or no filtering with kernel_std=0)
    vert_masked, mask = apply_masks_and_filter(
        vert, 
        asc_masks,
        asc_geom['waterMask'],
        kernel_std=kernel_std  # Set to 0 for no filtering
    )
    
    # Save results
    os.makedirs('Npy', exist_ok=True)
    np.save('Npy/vert.npy', vert_masked)
    np.save('Npy/east.npy', east)
    
    # Basic visualization
    plt.figure(figsize=(12,5))
    plt.subplot(121)
    plt.imshow(vert_masked, cmap='RdBu', vmin=-10, vmax=10)
    plt.colorbar(label='mm/yr')
    plt.title('Vertical Velocity')
    
    plt.subplot(122)
    plt.imshow(east, cmap='RdBu', vmin=-10, vmax=10)
    plt.colorbar(label='mm/yr')
    plt.title('East-West Velocity')
    plt.show()

if __name__ == "__main__":
    main()