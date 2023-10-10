#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 13:47:23 2023

@author: km
"""

import numpy as np
from matplotlib import pyplot as plt
import os,sys,glob,argparse,time,h5py
from datetime import date
import isce.components.isceobj as isceobj
from mroipac.looks.Looks import Looks
from SARTS import util,config
from Network import Network
from osgeo import gdal

from PyPS2 import noiseModel as nm


ps = config.getPS()
mintDir = './MintPy/'


filename = mintDir + 'timeseries.h5'
ds = h5py.File(filename, 'r+')
timeseries = ds['timeseries']
dates = np.asarray(ds['date'])
ts = np.asarray(ds['timeseries'][:, 1400, 1125])
ds.close()


fs=6 # Sampling interval in days (6 days probably for s1)

# Interpolate the time series to have even sampling
ts_interp,time_interp = nm.interpolate_timeseries(ps.dn0, ts, fs)
# Get decimal years for plotting time series
dec_year = []
yr0 = ps.dates[0][0:4]
dec_year_interp=[]
for dn in time_interp:
    yr = np.floor(dn/365) + int(yr0)
    doy = dn%365
    dec_year_interp.append(float(yr) + (doy/365.25))
dec_year_interp = np.asarray(dec_year_interp,dtype=np.float32)


# plt.figure();plt.imshow(timeseries[5,:,:],'magma')


# Compute uncertainties using 4 different methods

# lsq

mod, A, synth, residuals, C_model_white = nm.linear(ts_interp,time_interp)
plt.figure();plt.plot(dec_year_interp,ts_interp,'.')
plt.plot(dec_year_interp,synth)
m_uncertainty_white = (1.96*np.sqrt(np.diag(C_model_white)))

# bootstrap
mean_slope, std_slope, ci_slope, mean_intercept, std_intercept, ci_intercept = nm.bootstrap_linreg(time_interp, ts_interp,1000)

# spectral
m_uncertainty_white, m_uncertainty_color, rmse, spectral_index,sig = nm.get_uncertainties(residuals,A,fs,plot=True)

mod2 = [0.0021011/365, 1.35800236e-02]
mle_upper = np.dot(A,mod+ mod2)
mle_lower = np.dot(A,mod- mod2)

synth_upper_wh = np.dot(A, mod+m_uncertainty_white)
synth_lower_wh = np.dot(A, mod-m_uncertainty_white)
synth_upper_pl = np.dot(A, mod+m_uncertainty_color)
synth_lower_pl = np.dot(A, mod-m_uncertainty_color)



plt.figure()
plt.plot(dec_year_interp, ts_interp, '.',color='gray')
plt.plot(dec_year_interp, synth,'black')
plt.plot(dec_year_interp, synth_lower_wh,'g')
plt.plot(dec_year_interp, synth_upper_wh,'g')
plt.plot(dec_year_interp, synth_lower_pl,'--',color='purple')
plt.plot(dec_year_interp, synth_upper_pl,'--',color='purple')
# plt.plot(dec_year_interp, mle_upper,'--',color='red')
# plt.plot(dec_year_interp, mle_lower,'--',color='red')

plt.legend(['Data','mean rate','white','white','Powerlaw','powerlaw','mle'])
plt.show()

print('Boot strapped Rate uncertainty: ', str(np.round(1000*1.96*std_slope,3)))
print('White noise Rate uncertainty: ', str(np.round(1000*m_uncertainty_white[0],3)))
print('Colored noise Rate uncertainty: ', str(np.round(1000*m_uncertainty_color[0],3)))
print('spectral index: ', str(np.round(spectral_index,3)))
print('rate: ', str(np.round(mod[0],7)))


# mle
from scipy.optimize import minimize

def mle(C_diag_elements, r):
    from scipy.optimize import minimize

    # Construct diagonal matrix C from its diagonal elements
    C = np.diag(C_diag_elements)
    
    # Calculate log determinant of C
    log_det = np.log(np.linalg.det(C))
    
    # Calculate the inverse of C
    C_inv = np.linalg.inv(C)
    
    # Calculate the MLE
    term1 = np.dot(r, np.dot(C_inv, r))
    mle_value = log_det + term1 + len(r) * np.log(2 * np.pi)
    
    return mle_value
# Example data
r = residuals

# Initial guess for C's diagonal elements
C_init = np.ones_like(r)
# bounds = (np.diag(C_model_white)[0],np.diag(C_model_white)[0]*100)

result = minimize(mle, C_init, args=(r), method='L-BFGS-B', bounds=[(1e-12, 1e-10) for _ in r])

optimal_C_diag = result.x
optimal_C = np.diag(optimal_C_diag)