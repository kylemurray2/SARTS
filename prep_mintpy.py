#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 14:13:32 2024

@author: km
"""
import os,shutil
from SARTS import  config
from osgeo import gdal
import numpy as np
from scipy.ndimage import uniform_filter
import matplotlib.pyplot as plt

ps = config.getPS()


sartPath = config.__file__
sartPath = os.path.dirname(sartPath)
mpPath = os.path.join(sartPath, 'docs', 'mintpy_template.inps')

if not os.path.isdir('./MintPy'):
    print("Setting up MintPy directory")
    os.mkdir('MintPy')
    shutil.copy(mpPath, './MintPy/inps.cfg')
else:
    print('MintPy directory already exists.')


# Find a box that has the highest coherence
# Read in a coherence file
pair = ps.pairs[-1]
#load the last 16 images
cor_stack = []
for pair in ps.pairs[-16:]:
    ds = gdal.Open(os.path.join(ps.dolphin_work_dir,'interferograms',pair,'filt_lk.cor'))
    image = ds.GetVirtualMemArray()
    image = image.copy()
    image[np.isnan(image)] = 0
    cor_stack.append(image)
    
cor_stack = np.asarray(cor_stack,dtype=np.float32)

avg_cor = np.nanmean(cor_stack,axis=0)

def find_highest_avg_box(original_image, box_size, downsample_factor):
    # Downsample the image
    downsampled_image = original_image[::downsample_factor, ::downsample_factor]
    # Compute the local mean on the downsampled image
    local_mean = uniform_filter(downsampled_image, size=box_size//downsample_factor)
    # Find the position of the highest average value in the downsampled image
    max_pos = np.unravel_index(np.argmax(local_mean), local_mean.shape)
    # Scale the coordinates back to the original image size
    y1, x1 = max_pos[0] * downsample_factor, max_pos[1] * downsample_factor
    y2, x2 = y1 + box_size, x1 + box_size
    return y1, y2, x1, x2


box_size = 300
downsample_factor =4
y1, y2, x1, x2 = find_highest_avg_box(avg_cor, box_size, downsample_factor)

print(f'The best {box_size}X{box_size} box is: ')
print(f'y1:y2, x1:x2 --> {y1}:{y2}, {x1}:{x2}')

x = [x1, x1, x2, x2, x1]
y = [y1, y2, y2, y1, y1]

plt.figure(figsize=(10, 10))
plt.imshow(avg_cor, cmap='magma')
plt.plot(x, y, 'black', linewidth=2)
plt.title('coherence')
plt.legend(['Best coherence box'])
plt.show()

