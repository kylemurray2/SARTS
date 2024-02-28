from sklearn.decomposition import FastICA, PCA
import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
from sklearn.preprocessing import StandardScaler
from SARTS import config
from sklearn.linear_model import LinearRegression

plt.close('all')

ps = config.getPS()
ps.dn
ncomp = 4
# Load TS
mpdir = 'MintPy_1_4'
filename = os.path.join(mpdir, 'timeseries.h5')
ds = h5py.File(filename,'r+')   
data = np.asarray(ds['timeseries'][:,900:1150,2200:2600]) *1000 #convert to mm
ds.close()
z,y,x = data.shape
X = data.reshape(z, -1)  # Reshape to (z, y*x)

def getRates(ts,dn):
    '''
    ts is shape (time,npix)
    '''
    
    rates = np.zeros(ts.shape[1])
    dn_reshaped = ps.dn.reshape(-1, 1)
    
    # Perform linear regression for each time series
    for i in range(ts.shape[1]):
        # Create a linear regression object
        model = LinearRegression()
        model.fit(dn_reshaped, ts[:, i])
        rates[i] = model.coef_[0]
    
    rates_reshaped = rates.reshape(y,x)
    
    return rates_reshaped
    

# Compute ICA
ica = FastICA(n_components=ncomp, whiten="unit-variance")
S_ = ica.fit_transform(X)  # Reconstruct signals
A_ = ica.mixing_  # Get estimated mixing matrix
A_reshaped = A_.reshape(y,x,ncomp)
X_2 = np.dot(S_,A_.T)
X_2_reshaped =  X_2.reshape(z,y,x)

pca = PCA(n_components=ncomp)
H = pca.fit_transform(X)  # Reconstruct signals based on orthogonal components

plt.figure()
fig, ax = plt.subplots(3, 1)
ax[0].plot(X[:, 20000])
ax[0].set_title('Time Series')
# You should plot different components if S_ is a matrix with different independent components
ax[1].plot(S_[:, 0], label='IC 1')  # Plot the first independent component
ax[1].plot(S_[:, 1], label='IC 2')  # Plot the second independent component
ax[1].plot(S_[:, 2], label='IC 3')  # Plot the third independent component
ax[1].set_title('Independent Components')
ax[1].legend()  # Add the legend to the first subplot
ax[2].plot(X_2[:, 20000])
ax[2].set_title('Reconstructed Time Series')
plt.tight_layout()
plt.show()
vmin,vmax = -100,100


first_component = S_[:, 0].reshape(-1, 1)
first_mixing_vector = A_[:, 0].reshape(1, -1)  # This gets the entire third row, transpose it to column
second_component = S_[:, 1].reshape(-1, 1)
second_mixing_vector = A_[:, 1].reshape(1, -1)  # This gets the entire third row, transpose it to column
third_component = S_[:, 2].reshape(-1, 1)
third_mixing_vector = A_[:, 2].reshape(1, -1)  # This gets the entire third row, transpose it to column

X_reconstructed_first_component = np.dot(first_component, first_mixing_vector)#.reshape(z, y, x)
X_reconstructed_second_component = np.dot(second_component, second_mixing_vector)#.reshape(z, y, x)
X_reconstructed_third_component = np.dot(third_component, third_mixing_vector)#.reshape(z, y, x)

# Plot the three components form the mixing matrix
fig,ax = plt.subplots(3,1,figsize=(5,10))
vmin,vmax = None,None#-60,60
ax[0].imshow(A_reshaped[:,:,0],vmin=vmin,vmax=vmax);ax[0].set_title('ICA-1')
ax[1].imshow(A_reshaped[:,:,1],vmin=vmin,vmax=vmax);ax[1].set_title('ICA-2')
ax[2].imshow(A_reshaped[:,:,2],vmin=vmin,vmax=vmax);ax[2].set_title('ICA-3')
plt.show()

# Plot example time series
pty,ptx=22,215
# pty,ptx=103,298
plt.figure();
plt.plot(data[:,pty,ptx],'.')
plt.plot(X_2_reshaped[:,pty,ptx],'.')
plt.plot(X_reconstructed_second_component.reshape(z,y,x)[:,pty,ptx],'.')
plt.legend(['original ts','ica ts','ica ts 3rd'])

# Get the velocity maps for the original and three components
X_vel = getRates(X,ps.dn)
X_1_vel = getRates(X_reconstructed_first_component,ps.dn)
X_2_vel = getRates(X_reconstructed_second_component,ps.dn)
X_3_vel = getRates(X_reconstructed_third_component,ps.dn)

# Plot velocity maps
vmin,vmax = -.01,.01
fig,ax = plt.subplots(2,2)
ax[0,0].imshow(X_vel,vmin=vmin,vmax=vmax);ax[0,0].set_title('Original rates')
ax[0,1].imshow(X_1_vel,vmin=vmin,vmax=vmax);ax[0,1].set_title('1st ICA rates')
ax[1,0].imshow(X_2_vel,vmin=vmin,vmax=vmax);ax[1,0].set_title('2nd ICA rates')
ax[1,1].imshow(X_3_vel,vmin=vmin,vmax=vmax);ax[1,1].set_title('3rd ICA rates')
plt.show()

plt.figure()
plt.imshow(X_vel-X_3_vel)