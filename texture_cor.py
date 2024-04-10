import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.interpolate import interp2d
from scipy.ndimage import gaussian_filter
import os
import h5py
from SARTS import config

dirName='Hawaii'
p2 = '291'

mp_dir_2 = 'MintPy'

baseDir =  os.path.join('/d/HI/Alos/Asc', dirName) 
os.chdir(baseDir)


dir2 = os.path.join(baseDir, p2)

ps_2 = config.getPS(dir2)

if not os.path.isdir('geotifs'):
    os.mkdir('geotifs')
if not os.path.isdir('Npy'):
    os.mkdir('Npy')
    
# cropxmin2,cropxmax2,cropymin2,cropymax2 = 139,ps_2.nxl,0,ps_2.nyl

mpdir_2 = os.path.join(dir2, mp_dir_2)



filename    = os.path.join(mpdir_2, 'inputs/geometryRadar.h5')
ds          = h5py.File(filename,'r+')   
hgt_2 = np.asarray(ds['height'])
lons_2      = np.asarray(ds['longitude'])
lats_2      = np.asarray(ds['latitude'])
az_2        = np.asarray(ds['azimuthAngle'])
inc_2       = np.asarray(ds['incidenceAngle'])
waterMask_2 = np.asarray(ds['waterMask'])
ds.close()


#Velocity
filename    = os.path.join(mpdir_2, 'velocity.h5')
ds          = h5py.File(filename,'r+')   
velocity_2  = np.asarray(ds['velocity']) *1000 #convert to mm
ds.close()



velocity_2[waterMask_2==0] = np.nan
velocity_2[velocity_2==0] = np.nan
fig,ax = plt.subplots(1,1)
ax[0].imshow(velocity_2);ax[0].set_title('velocity 2')
plt.show()


# Mask
spatial_coh_thresh  = 0.7
temporal_coh_thresh = 0.8
numtriplets_thresh_ratio  = 0.05
vel_std_thresh = 1

flipud=False
if not os.path.isfile(os.path.join(ps_2.workdir, p2, 'Npy/mask_full.npy')):
    mask_full_2 = mask_mp.get_mask(mintPy_dir=mpdir_2,spatial_coh_thresh  = spatial_coh_thresh,temporal_coh_thresh =temporal_coh_thresh,vel_std_thresh=vel_std_thresh,numtriplets_thresh_ratio=numtriplets_thresh_ratio,flipud=flipud)
    np.save(os.path.join(ps_2.workdir, p2, 'Npy/mask_full.npy'),mask_full_2)
else:
    mask_full_2 = np.load(os.path.join(ps_2.workdir, p2, 'Npy/mask_full.npy'))



# Load Data
N = len(dates)
N_ref = 2
Na_ref = 282
Nr_ref = 205
lamda = ps_2.lamda

# Acquisition time
Date = np.loadtxt('291/MintPy/timeseries.h5', dtype=str, delimiter=',', usecols=[0])
Date = [datetime.strptime(date, '%Y%m%d') for date in Date[:N]]
t = np.arange(N)



# Size
Na, Nr = hgt_2.shape

# Rate from mintpy
rate = np.loadtxt('../velocity.h5')  # unit: m/y
rate = rate.T * 100  # unit: cm/y

# Deformation mask
maskdef = 0
if maskdef == 1:
    mask_def = rate.copy()
    mask_def[np.abs(rate) > 1.5] = np.nan
    mask_def[np.logical_not(np.isnan(mask_def))] = 1
else:
    mask_def = np.ones_like(hgt_2)

# Coherence mask
mask_coh = np.loadtxt('../temporalCoherence.h5')
mask_coh[mask_coh < 0.87] = np.nan
mask_coh[np.logical_not(np.isnan(mask_coh))] = 1

# MASK
MASK = mask_def * mask_coh

# Phase time series
phase_ts_atm = np.transpose(np.loadtxt('../timeseries_ERA5_ramp_demErr.h5'))[:N]  # unit: m
phase_ts_atm = (4 * np.pi / lamda) * phase_ts_atm
phase_ts_atm -= phase_ts_atm[Na_ref, Nr_ref]

phase_ts_deramp = np.transpose(np.loadtxt('../timeseries_ramp_demErr.h5'))[:N]  # unit: m
phase_ts_deramp = (4 * np.pi / lamda) * phase_ts_deramp
phase_ts_deramp -= phase_ts_deramp[Na_ref, Nr_ref]
phase_ts_deramp[np.isnan(phase_ts_deramp)] = 0

# Local Slope Estimation Based on Texture Correlation
# Parameter setting
W = 101  # window size (square window, must be odd number)
r = 0.5  # overlap ratio

# Window segmentation
w = (W - 1) // 2  # half of window size (square window: 2*w+1)
overlap = round(r * (2 * w + 1))  # overlap pixels
Na_C = np.arange(w, Na, 2 * w - overlap)
Nr_C = np.arange(w, Nr, 2 * w - overlap)

# Slope estimation
k_HTC = slope_estimation(phase_ts_deramp, hgt_2, MASK, W, r, N_ref)

# Slope interpolation
Y, X = np.meshgrid(Na_C, Nr_C)
Xq, Yq = np.meshgrid(np.arange(Nr), np.arange(Na))
k_KTC_interp = np.zeros_like(phase_ts_deramp)
for n in range(N):
    k_HTC_filt = gaussian_filter(k_HTC[:, :, n], 7)
    k_KTC_interp[:, :, n] = interp2d(X, Y, k_HTC_filt, kind='cubic')(Xq, Yq)

# Intercept filtering
phase_ts_HTC_low = phase_ts_deramp.copy()
intercept = np.zeros_like(phase_ts_HTC_low)
for n in range(N):
    tmp = phase_ts_deramp[:, :, n] - k_KTC_interp[:, :, n] * hgt_2
    tmp_filt = gaussian_filter(tmp, 251)
    tmp -= tmp_filt
    phase_ts_HTC_low[:, :, n] = tmp
    intercept[:, :, n] = tmp_filt
phase_ts_HTC_low -= phase_ts_HTC_low[Na_ref, Nr_ref]

# High resolution correction
phase_ts_HTC_high, flag = high_resolution_correction(phase_ts_HTC_low, k_KTC_interp, MASK, t, 251, 251, Na_ref, Nr_ref)
phase_ts_HTC = phase_ts_HTC_high

# Velocity
tmp = phase_ts_HTC.copy()
tmp -= tmp[282, 205]
velocity = np.zeros_like(hgt_2)
for i in range(Na):
    for j in range(Nr):
        ts = tmp[i, j]
        coe = np.polyfit(t, ts, 1)
        velocity[i, j] = coe[0]
velocity = velocity * 365.25 * 100 / (4 * np.pi / lamda)  # unit: cm/y

plt.figure()
plt.imshow(velocity, cmap='jet', vmin=-2, vmax=2)
plt.colorbar(label='Velocity (cm/y)')
plt.scatter(Nr_ref, Na_ref, color='k', marker='o', s=20)
plt.axis('off')

# Time series-1D
Na_c = 372
Nr_c = 372

# Original
ts = phase_ts_deramp[Na_c, Nr_c] / (4 * np.pi / lamda) * 100  # unit: cm

# ERA5
ts0 = phase_ts_atm[Na_c, Nr_c] / (4 * np.pi / lamda) * 100  # unit: cm
coe0 = np.polyfit(t, ts0, 1)

# THC
ts2 = phase_ts_HTC[Na_c, Nr_c] / (4 * np.pi / lamda) * 100
coe2 = np.polyfit(t, ts2, 1)

# Slope & Intercept
k = k_KTC_interp[Na_c, Nr_c]
d = intercept[Na_c, Nr_c]

plt.figure()
plt.subplot(2, 1, 1)
plt.plot(Date, ts, ':o', markersize=3, color='black', markerfacecolor='black', label='Original')
plt.plot(Date, ts0, ':o', markersize=3, color='green', markerfacecolor='green', label='ERA5')
plt.plot(Date, ts2, ':>', markersize=3, color='red', markerfacecolor='red', label='Local linear fitting')
plt.plot(Date, coe0[0] * t + coe0[1], color='green')
plt.plot(Date, coe2[0] * t + coe2[1], color='red')
plt.legend()
plt.xlabel('Time Series')
plt.ylabel('Deformation (cm)')
plt.ylim(-2, 8)
plt.xticks(rotation=45)
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(Date, k)
plt.xlabel('Time Series')
plt.ylabel('Linear Coefficient')
plt.xticks(rotation=45)
plt.grid(True)

plt.show()
