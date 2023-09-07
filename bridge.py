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


def getConCom(msk, minimumPixelsInRegion=1000):
    '''
    Takes a binary input (like a mask) as input and outputs labels for
    regions greater than the given minimum pixels.
    '''
    
    ratesu8 = (msk*255).astype(np.uint8)
    num_labels, labels = cv2.connectedComponents(ratesu8)
    
    ## If the count of pixels less than a threshold, then set pixels to `0`.
    print('Removing small connected regions...')
    for i in range(1, num_labels+1):
        pts =  np.where(labels == i)
        if len(pts[0]) < minimumPixelsInRegion:
            labels[pts] = 0
            
    return labels


def main(imageIn, mask, minPix=1000, plot=False):
    #remove mean from each disconnected region
    labels = getConCom(mask,minPix)
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
    
    return imageOut


