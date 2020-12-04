import numpy as np
import linecache
from decimal import Decimal
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter

from scipy.interpolate import griddata
import scipy.ndimage.filters
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

import glob # matches filenames and directories knowing some part of the name

#voltage = '0.00'
voltage = '1.40'

# arrays to be filled in the PLDOS and their position along the trasnport direction
X, Y, Z, = np.array([]), np.array([]), np.array([])


for i in range(31):
    print(glob.glob('../../TS_{}/PLDOS_along_Z/dos_{}L_*.dat'.format(voltage, i+4))[0])
    f = np.loadtxt(glob.glob('../../TS_{}/PLDOS_along_Z/dos_{}L_*.dat'.format(voltage, i+4))[0], unpack=True)
    
    X_dat = f[0,:]
    Y_dat = np.linspace(9.52795032+(i*3.176), 9.52795032+((i+1)*3.176), num=400)
    #Y_dat = np.linspace(i, i+1, num=400)
    Z_dat = f[1,:]
    #print(X_dat[0],X_dat.shape)
    #print(Y_dat[0], Y_dat.shape)
    
    # Convert from pandas dataframes to numpy arrays
    #print(X.shape)
    #print(Y.shape)
    #print(Z.shape)
    for i in range(len(X_dat)):
            X = np.append(X, X_dat[i])
            Y = np.append(Y, Y_dat[i])
            Z = np.append(Z, Z_dat[i])
    #print(X[0],X.shape)
    #print(X.min(), X.max())
    #print(Y.min(), Y.max())
    
    
xi = np.linspace(-2, 2, num=21)
print(int(min(Y)), int(max(Y)))
yi = np.linspace(int(min(Y)), int(max(Y)), num=100)
#print(X[0].shape)
#print(Y.shape)
print(Z[0],Z.shape)
Z = scipy.ndimage.filters.gaussian_filter(Z, sigma=12)#, mode='constant')

z = griddata((X, Y), Z, (xi[None,:], yi[:,None]), method='nearest', 
fill_value=0.0, rescale=True)

# parameters for the plot
fig = plt.figure(figsize=[10, 6])
#plt.xlim(-2.0, 2.0)
#plt.ylim(0.0, 3.5)
plt.xlabel('Z (transport direction) (Ang)', fontsize=20) #($10^{-10}$)
plt.ylabel('Energy (eV)', fontsize=20)
plt.title('Pristine - Bias = 1.40 eV', fontsize=24, pad=25.0) 
#plt.xticks(np.range(0.0,1.51, step=0.1)

plt.axes().xaxis.set_major_locator(ticker.MultipleLocator(10))
plt.tick_params(axis='x', which='major', width=2.00, length=5.0, labelsize=20)
plt.axes().xaxis.set_minor_locator(ticker.MultipleLocator(2))
plt.tick_params(axis='x', which='minor', width=1.00, length=3.5, labelsize=14)
plt.axes().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
plt.axes().yaxis.offsetText.set_fontsize(20)
plt.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)
plt.axes().yaxis.set_major_locator(ticker.MultipleLocator(0.5))
plt.axes().yaxis.set_minor_locator(ticker.MultipleLocator(0.1   ))
plt.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)

#print(yi)
levels = MaxNLocator(nbins=110).tick_values(z.min(), z.max())

#cfff = plt.imshow(z.T, cmap=plt.cm.jet, interpolation="gaussian", aspect="auto")
#plt.colorbar(cfff)
#cf = plt.contourf(xi, yi, z, levels=levels, cmap=plt.cm.jet)
cf = plt.contourf(yi, xi, z.T, levels=levels, cmap=plt.cm.jet)#, interpolation="gaussian", aspect="auto")
plt.colorbar(cf)#, ax=axs[index])
#plt.set_title('{}'.format(lst_labels[index]), fontsize=28)

plt.savefig('PLDOS_along_trans_pris_{}.png'.format(voltage))

plt.show()
