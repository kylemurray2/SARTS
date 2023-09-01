#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 15:55:20 2023

@author: km


"""

import numpy as np
from matplotlib import pyplot as plt
from SARTS import util
import localParams

ps = localParams.getLocalParams()

def bridge(imageIn, mask, minPix=1000):
    #remove mean from each disconnected region
    minPix = 1000
    labels = util.getConCom(mask,minPix)
    fig,ax = plt.subplots(2,1,figsize=(5,6))
    ax[0].imshow(imageIn);ax[0].set_title('mask')
    ax[1].imshow(labels);ax[1].set_title('connected regions')
    
    imageOut = imageIn.copy()
    
    for ii in range(int(labels.max())):
        if len(imageOut[labels==ii+1]) < minPix:
            imageOut[labels==ii+1] = np.nan # mask out small islands of data
            imageIn[labels==ii+1] = 0
        else:
            imageOut[labels==ii+1]-=np.nanmean(imageOut[labels==ii+1])
    
    return imageOut


