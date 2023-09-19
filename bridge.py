#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 15:55:20 2023

Experimental
Takes a mask and an image (rate map or similar) and labels connected components
with a minimum number of pixels. Then an average is removed from each connected 
region. 

Note: It is often not ideal to remove an average if a significant number of pixels
have deformation.  

@author: km


"""

import numpy as np
from matplotlib import pyplot as plt
import cv2
from skimage.morphology import remove_small_objects

def getConCom(msk, minimumPixelsInRegion=3000):
    '''
    Takes a binary input (like a mask) as input and outputs labels for
    regions greater than the given minimum pixels.
    '''
    
    print('Removing small connected regions...')
    msk2 = remove_small_objects(msk, minimumPixelsInRegion, connectivity=1)

    ratesu8 = (msk2*255).astype(np.uint8)
    num_labels, labels = cv2.connectedComponents(ratesu8)
    
    return msk2,labels


def main(imageIn, mask, minPix=1000, plot=False):
    #remove mean from each disconnected region
    mask2,labels = getConCom(mask,minPix)
    
    fig,ax = plt.subplots(2,1,figsize=(5,6))
    ax[0].imshow(imageIn);ax[0].set_title('mask')
    ax[1].imshow(labels);ax[1].set_title('connected regions')
    
    imageOut = imageIn.copy()
    
    print('Re-referencing each connected region...')

    for ii in range(int(labels.max())):
        if len(imageOut[labels==ii+1]) < minPix:
            imageOut[labels==ii+1] = np.nan # mask out small islands of data
            imageIn[labels==ii+1] = 0
        else:
            imageOut[labels==ii+1]-=np.nanmean(imageOut[labels==ii+1])
    
    return imageOut, mask2



