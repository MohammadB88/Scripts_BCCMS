import numpy as np
import linecache
from decimal import Decimal
import math
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import scipy.ndimage.filters
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator


#layers = np.range(4,38)

# read in the PDOS at each layer
#f = np.loadtxt('dos_4L_73-96.dat', unpack=True)
#plt.plot(f[0,:], f[1,:], 'k.', ms=1)

X, Y, Z, = np.array([]), np.array([]), np.array([])

for i in range(31):
    #print(index)
    f = np.loadtxt('dos_{}L_{}-{}.dat'.format(i+4,(i+3)*24+1,(i+4)*24), unpack=True)
    
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
    
    
xi = np.linspace(-2, 2, 21)
print(int(min(Y)), int(max(Y)))
yi = np.linspace(int(min(Y)), int(max(Y)), 41)
#print(X[0].shape)
#print(Y.shape)
#print(Z[0],Z.shape)
Z = scipy.ndimage.filters.gaussian_filter(Z, 12)

z = griddata((X, Y), Z, (xi[None,:], yi[:,None]), method='nearest', 
fill_value=0.0) #rescale=True, 

#print(yi)
levels = MaxNLocator(nbins=110).tick_values(z.min(), z.max())

cf = plt.imshow(z, cmap=plt.cm.jet, interpolation="gaussian", aspect="auto")
#cf = plt.contourf(xi, yi, z, levels=levels, cmap=plt.cm.jet)#, interpolation="gaussian", aspect="auto")
plt.colorbar(cf)#, ax=axs[index])
#plt.set_title('{}'.format(lst_labels[index]), fontsize=28)

#plt.savefig('dis_{}.png'.format(geom))
#plt.show()
    
    
#plt.plot(f[0], f[1], color='black', label='pristine')

#plt.contourf(f[0], 4, f[1], 20, cmap='RdGy')
#plt.colorbar();

#handles, labels = ax1.get_legend_handles_labels()

#fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=4, fancybox=True, shadow=True, fontsize=23)

#plt.savefig('trans_DOS_Mo_vacancies_{}.png'.format(voltage))

plt.show()
