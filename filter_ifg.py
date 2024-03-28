#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 01:12:44 2024

@author: murray8
"""
import pywt
import numpy as np
import rasterio
import matplotlib.pyplot as plt
from skimage import io, img_as_float
from skimage.restoration import denoise_wavelet
from scipy.ndimage import uniform_filter, gaussian_filter

def normalize_data(data):
    min_val = np.min(data)
    max_val = np.max(data)
    return (data - min_val) / (max_val - min_val)


def wavelet_denoise_phase(image, wavelet='db4', level=None, threshold_scale=3):
    """
    Apply wavelet-domain regularization to denoise the phase of an image.

    Parameters:
    - image: 2D numpy array, input phase image.
    - wavelet: str, type of wavelet to use.
    - level: int, the level of wavelet decomposition. If None, max level is used.
    - threshold_scale: float, scale factor to apply to the threshold level.

    Returns:
    - denoised_image: 2D numpy array, the denoised phase image.
    
    Example usage:
    phs_f_wavelet = wavelet_denoise_phase(phs, wavelet='db4', level=None, threshold_scale=1)


    """
    # Determine the number of decomposition levels
    if level is None:
        level = pywt.dwt_max_level(data_len=min(image.shape), filter_len=pywt.Wavelet(wavelet).dec_len)
        
    print(f'number of levels: {level}')
    # Compute wavelet coefficients
    coeffs = pywt.wavedec2(image, wavelet, mode='symmetric', level=level)
    
    # Estimate the noise sigma from the median of the absolute deviation from the median
    # of the detail coefficients at the finest scale
    sigma = np.median(np.abs(coeffs[-1][-1])) / 0.6745
    
    # Calculate the threshold
    threshold = sigma * threshold_scale * np.sqrt(2 * np.log(image.size))
    print(f'threshold: {threshold}')

    # Apply soft thresholding to the detail coefficients
    new_coeffs = list(coeffs)
    for i in range(1, len(coeffs)):
        new_coeffs[i] = tuple(pywt.threshold(c, value=threshold, mode='soft') for c in coeffs[i])

    # Reconstruct the denoised image
    denoised_image = pywt.waverec2(new_coeffs, wavelet)
    
    return denoised_image

def wavelet_denoise_phase_with_levels(image, wavelet='db4', max_level=None, threshold_scale=3):
    """
    Apply wavelet-domain regularization to denoise the phase of an image and output images for each level.
    
    Parameters:
    - image: 2D numpy array, input phase image.
    - wavelet: str, type of wavelet to use.
    - max_level: int, the maximum level of wavelet decomposition. If None, max level is used.
    - threshold_scale: float, scale factor to apply to the threshold level.
    
    Returns:
    - denoised_images: List of 2D numpy arrays, the denoised phase image at each level.
    """
    # Compute wavelet coefficients
    coeffs = pywt.wavedec2(image, wavelet, mode='symmetric', level=max_level)
    max_level = pywt.dwt_max_level(data_len=min(image.shape), filter_len=pywt.Wavelet(wavelet).dec_len) if max_level is None else max_level
    
    # Initialize list to store denoised images at each level
    denoised_images = []
    
    for level in range(1, max_level+1):
        # Copy the coefficients
        temp_coeffs = list(coeffs)
        
        # Estimate noise sigma from the detail coefficients at the current level
        sigma = np.median(np.abs(temp_coeffs[level][-1])) / 0.6745
        threshold = sigma * threshold_scale * np.sqrt(2 * np.log(image.size))
        
        # Apply soft thresholding to the detail coefficients at and above the current level
        for i in range(level, len(temp_coeffs)):
            temp_coeffs[i] = tuple(pywt.threshold(c, value=threshold, mode='soft') for c in temp_coeffs[i])
        
        # Reconstruct the image from the modified coefficients
        level_image = pywt.waverec2(temp_coeffs, wavelet)
        denoised_images.append(level_image)
    
    return denoised_images


def lee_filter(img, kernel_size=7, damping_factor=1):
    """
    Apply Lee filter to an image.
    """
    img_mean = uniform_filter(img, size=kernel_size)
    img_sqr_mean = uniform_filter(img**2, size=kernel_size)
    img_variance = img_sqr_mean - img_mean**2

    overall_variance = np.var(img)

    img_weights = img_variance / (img_variance + overall_variance / damping_factor)
    filtered_img = img_mean + img_weights * (img - img_mean)
    return filtered_img

def frost_filter(img, kernel_size=7, damping_factor=2):
    """
    Apply Frost filter to an image.
    """
    img_mean = uniform_filter(img, size=kernel_size)
    local_variance = uniform_filter(img**2, size=kernel_size) - img_mean**2
    
    overall_variance = np.mean(local_variance)
    alpha = damping_factor * overall_variance
    
    def frost_kernel(distance):
        return np.exp(-distance / alpha)
    
    distances = np.arange(kernel_size) - kernel_size // 2
    kernel = frost_kernel(distances[:, None]**2 + distances[None, :]**2)
    
    filtered_img = gaussian_filter(img, sigma=1, mode='constant', cval=0, truncate=kernel_size / 6.0)
    return filtered_img

def gamma_map_filter(img, kernel_size=7, cu=0.25):
    """
    Apply Gamma MAP filter to an image.
    """
    local_mean = uniform_filter(img, size=kernel_size)
    local_var = uniform_filter(img**2, size=kernel_size) - local_mean**2
    overall_variance = np.mean(local_var)
    
    tao = cu * overall_variance
    filtered_img = local_mean + (np.maximum(local_var - tao, 0) / (local_var + tao)) * (img - local_mean)
    return filtered_img

def kuan_filter(img, kernel_size=7, cu=0.25):
    """
    Apply Kuan filter to an image.
    """
    local_mean = uniform_filter(img, size=kernel_size)
    local_var = uniform_filter(img**2, size=kernel_size) - local_mean**2
    overall_variance = np.mean(local_var)
    
    tao = cu * overall_variance
    filtered_img = local_mean + (local_var / (local_var + tao)) * (img - local_mean)
    return filtered_img


# def main():
#     ds = gdal.Open('/mnt/lustre/koa/scratch/murray8/HI/Asc/Oahu/merged/geom_reference/waterMask_lk.rdr.vrt')
#     wm = ds.GetVirtualMemArray()
#     mask = wm==1
#     plt.figure();plt.imshow(wm)
#     ifg_fn = 'dolphin/interferograms/20150111_20150429/fine_lk.int'
    
#     # Load the raster image
#     with rasterio.open(ifg_fn) as src:
#         image = src.read(1)  # Read the first band
    
#     real = np.real(image)
#     imag = np.imag(image)
#     phs = np.angle(image)
    
#     real/=np.max(real)
#     imag/=np.max(imag)
    
#     real[np.isnan(real)] = 0
#     real[np.isinf(real)] = 0
    
#     max_level = pywt.dwt_max_level(data_len=min(real.shape), filter_len=pywt.Wavelet('db4').dec_len)
#     chosen_level = max(1, min(max_level, 4))  # Ensures at least 1 and no more than 4 or the max level
    
    
#     # Apply wavelet denoising
#     wavelet='db4'
#     threshold_multiplier=1.5
#     level=4
#     real_f = wavelet_denoise(real,wavelet=wavelet,level=level,threshold_multiplier=threshold_multiplier)
#     imag_f = wavelet_denoise(imag,wavelet=wavelet,level=level,threshold_multiplier=threshold_multiplier)
    
#     phs_wave = denoise_wavelet(phs, wavelet='db4', mode='hard')
#     # Assuming 'phs' is your raster array to be filtered
#     # Apply each filter to the phs array
#     from scipy.ndimage import binary_dilation, binary_erosion
#     from scipy.ndimage import gaussian_filter
#     blurred_mask = gaussian_filter(wm.astype(float), sigma=10)
#     feathered_mask = np.clip(blurred_mask / blurred_mask.max(), 0, 1)
    
#     plt.figure;plt.imshow(feathered_mask);plt.show()
#     kernel = 3
#     cu=1
#     df = 4
#     phs_lee = lee_filter(phs,kernel_size=kernel,damping_factor=df)
#     phs_frost = frost_filter(phs,kernel_size=kernel,damping_factor=df*2)
#     phs_gamma_map = gamma_map_filter(phs, kernel_size=kernel, cu=cu)
#     phs_gamma_map2 = gamma_map_filter_masked(phs, smoothed_mask,kernel_size=kernel, cu=cu)
    
#     phs_kuan = kuan_filter(phs, kernel_size=kernel, cu=cu)
    
#     # Plot all four
#     fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    
#     axs[0, 0].imshow(phs_lee, cmap='jet')
#     axs[0, 0].set_title('Lee Filter')
#     axs[0, 0].axis('off')
    
#     axs[0, 1].imshow(phs_frost, cmap='jet')
#     axs[0, 1].set_title('Frost Filter')
#     axs[0, 1].axis('off')
    
#     axs[1, 0].imshow(phs_gamma_map2, cmap='jet')
#     axs[1, 0].set_title('Gamma Map mask Filter')
#     axs[1, 0].axis('off')
    
#     axs[1, 1].imshow(phs_gamma_map, cmap='jet')
#     axs[1, 1].set_title('gamma orig Filter')
#     axs[1, 1].axis('off')
#     plt.tight_layout()
#     plt.show()

# from scipy.ndimage import binary_dilation, binary_fill_holes, generic_filter

def fill_nan_with_nearest(image):
    """
    Custom filter function to replace NaN values with the nearest non-NaN value.
    """
    # Flatten the array and remove NaNs
    non_nan_values = image[~np.isnan(image)]
    if non_nan_values.size > 0:
        return np.nanmin(non_nan_values)  # Return the nearest non-NaN value (minimum in this case)
    else:
        return np.nan  # Return NaN if no non-NaN values are found

def grow_unmasked_regions(image, iterations=1, footprint_size=3):
    """
    Grow non-masked parts into NaN parts by a specified number of pixels (iterations).
    
    Parameters:
    - image: 2D numpy array, input image with NaNs as masked values.
    - iterations: int, the number of pixels to grow the non-masked regions.
    - footprint_size: int, the size of the footprint for the generic filter.
    
    Returns:
    - grown_image: 2D numpy array, the image with grown non-masked regions.
    """
    # Create a binary mask where True indicates a non-NaN value
    non_nan_mask = ~np.isnan(image)
    
    # Perform binary dilation on the mask
    dilated_mask = binary_dilation(non_nan_mask, iterations=iterations)
    
    # Fill holes inside the dilated mask if needed
    filled_mask = binary_fill_holes(dilated_mask)
    
    # Use generic_filter to replace NaNs in the dilated areas
    footprint = np.ones((footprint_size, footprint_size))
    grown_image = generic_filter(image, function=fill_nan_with_nearest, footprint=footprint, mode='nearest')
    
    # Ensure only the originally NaN areas that are now in the dilated region are updated
    final_image = np.where(~non_nan_mask & filled_mask, grown_image, image)
    
    return final_image,grown_image

# a,phs_dilate = grow_unmasked_regions(phs, iterations=3)
# plt.figure();plt.imshow(phs_dilate,cmap='jet')

# phs_gamma_map = gamma_map_filter(phs, kernel_size=7, cu=.5)
# fig, ax = plt.subplots(1,2, figsize=(10, 10))
# vmin,vmax = -3.4,3.4
# ax[0].imshow(phs[775:950,2400:2500],cmap='jet',vmin=vmin,vmax=vmax);ax[0].set_title('orig phase')

# ax[1].imshow(phs_gamma_map[775:950,2400:2500],cmap='jet',vmin=vmin,vmax=vmax);ax[1].set_title('filt phase')
# plt.show()


# cpxf    = real_f + 1j * imag_f

# plt.figure();plt.imshow(np.angle(cpxf)[775:950,2400:2500],cmap='jet');plt.title('filtered phase')
# plt.figure();plt.imshow(phs_f[775:950,2400:2500],cmap='jet');plt.title('filtered phase')
# plt.figure();plt.imshow(phs[775:950,2400:2500],cmap='jet');plt.title('orig phase')
# plt.figure();plt.imshow(phs[775:950,2400:2500]-phs_f[775:950,2400:2500],cmap='jet');plt.title('diff phase')

