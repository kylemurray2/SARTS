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

mpDir = 'MintPy'


ps = config.getPS()


sartPath = config.__file__
sartPath = os.path.dirname(sartPath)
mpPath = os.path.join(sartPath, 'docs', 'mintpy_template.inps')

inps_fn = os.path.join(mpDir,'inps.cfg')

if not os.path.isdir(mpDir):
    print("Setting up MintPy directory")
    os.mkdir(mpDir)
    shutil.copy(mpPath, inps_fn)
else:
    print('MintPy directory already exists.')

def update_cfg_variable(cfg_file_path, variable_name, new_value):
    """
    Updates the value of a specific variable in a .cfg file.

    Parameters:
    - cfg_file_path: Path to the .cfg file.
    - variable_name: Name of the variable to update.
    - new_value: New value to assign to the variable.
    """
    # Read the file and store lines in memory
    with open(cfg_file_path, 'r') as file:
        lines = file.readlines()

    # Modify the line containing the variable
    with open(cfg_file_path, 'w') as file:
        for line in lines:
            if line.startswith(variable_name):
                # Split the line at the '=' sign, replace the value, and reassemble the line
                parts = line.split('=')
                if len(parts) == 2:
                    line = f"{parts[0]}= {new_value}\n"
            file.write(line)



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
original_image = avg_cor

def find_highest_avg_box(original_image, box_size, downsample_factor):
    # Downsample the image
    downsampled_image = original_image[::downsample_factor, ::downsample_factor]
    ds_box = box_size//downsample_factor
    # Compute the local mean on the downsampled image
    local_mean = uniform_filter(downsampled_image, size=ds_box)
    # Find the position of the highest average value in the downsampled image
    max_pos = np.unravel_index(np.argmax(local_mean), local_mean.shape)
    # Scale the coordinates back to the original image size
    y1, x1 = max_pos[0] * downsample_factor, max_pos[1] * downsample_factor
    
    y1 -=(box_size//2)
    x1 -=(box_size//2)
    
    
    y2, x2 = y1 + box_size, x1 + box_size
    return y1, y2, x1, x2


# plt.figure();plt.imshow(downsampled_image)
# plt.figure();plt.imshow(local_mean)
# plt.figure();plt.imshow(downsampled_image)


box_size = 300
downsample_factor =4
y1, y2, x1, x2 = find_highest_avg_box(avg_cor, box_size, downsample_factor)

print(f'The best {box_size}X{box_size} box is: ')
print(f'y1:y2, x1:x2 --> {y1}:{y2}, {x1}:{x2}')
x = [x1, x1, x2, x2, x1]
y = [y1, y2, y2, y1, y1]

plt.figure(figsize=(10, 10))
plt.imshow(avg_cor, cmap='magma')
plt.plot(x, y, 'red', linewidth=2)
plt.title('coherence')
plt.legend(['Best coherence box'])
plt.show()

# replace the value in inps.cfg
variable_name = 'mintpy.network.aoiYX'
new_value = f'{y1}:{y2}, {x1}:{x2}' 
update_cfg_variable(inps_fn, variable_name, new_value)

