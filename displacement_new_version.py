import numpy as np
import linecache
from decimal import Decimal
import math

import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter

# storing the displacement for 1sv
#tmpdis = open("dis_LRE_s.txt", "wt")
## storing the displacement for 2sv-top
#tmpdis = open("displacement2svtop.txt", "wt")
#a1=np.zeros((147,4))
## storing the displacement for 2sv-par
#tmpdis = open("displacement2svpar.txt", "wt")
#a2=np.zeros((147,4))

i = 0
l = 0 

#Atoms_to_skip = [187, 703]
geoms = ['LRE187_703_s', 'LRE182_187_699_703_2stop','LRE174_187_703_715_2spar', 'LRE202_695_mo', 'LRE187_199_202_210_679_691_695_703_mo3s', 'LRE182_187_195_199_202_206_210_675_679_687_691_695_699_703_mo6s', 'LRE187_703_s_with_mo', 'LRE182_187_699_703_2stop_with_mo', 'LRE202_695_mo_with_s', 'LRE202_695_mo_with_2stop']


# Read, calculate, and write the displacement 
#for geom in geoms:
    #num_atoms = int(linecache.getline("not_optim/dev_w8_{}.fdf".format(geom), 8).rstrip().split()[1])
    #print(geom, 'number of atoms {}'.format(num_atoms))
    
    #a0=np.zeros((num_atoms,4))
    #for idx in range(0,num_atoms):
        #l_not_optim = linecache.getline("not_optim/dev_w8_{}.fdf".format(geom), idx+11)
        #linecache.clearcache()
        #coord_not_optim = l_not_optim.split()
        ##print(l_not_optim)
        ## calculating displacement for 1sv vacancy.
        #l_optim = linecache.getline("optim/dev_w8_{}.fdf".format(geom), idx+11)
        #linecache.clearcache()
        #coord_optim = l_optim.split()
        
        #displacement = math.sqrt(((float(coord_optim[0])-float(coord_not_optim[0]))**2)+((float(coord_optim[1])-float(coord_not_optim[1]))**2)+((float(coord_optim[2])-float(coord_not_optim[2]))**2))
        
        ##except IndexError:
            ##print("not_optim/dev_w8_{}.fdf".format(geom))
            ##print(linecache.getline("not_optim/comp_inp_dev_w8_{}.fdf".format(geom), idx+11))
            ##print(idx)
        #a0[idx,0] = float(coord_optim[0])
        #a0[idx,1] = float(coord_optim[1])
        #a0[idx,2] = float(coord_optim[2])
        #a0[idx,3] = displacement
        ##print(idx)
    #print(idx, '\n')
    #np.savetxt('dis_{}.txt'.format(geom), a0)
    
fig, axs = plt.subplots(1, 10, figsize=(20,10), sharex=True, squeeze=True)
lst_labels = ['$V_{S}$', '$V_{2S-top}$', '$V_{2S-par}$', '$V_{Mo}$', '$V_{Mo+3S}$', '$V_{Mo+6S}$', '$Mo_{S}$', '$Mo_{2S-top}$', '$S_{Mo}$', '$2S-top_{Mo}$']

# With a calculation, the min and max displacement has been calculated.
cb_min = 0.0
cb_max = 2.1874486506750372 
levels = MaxNLocator(nbins=600).tick_values(cb_min, cb_max)

# Set the parameters for the plots
for idx, geom in enumerate(geoms):
    axs[idx].tick_params(axis='x', which='major', width=2.00, length=5.0, labelsize=20, direction='in', bottom=True, top=True)
    
    axs[idx].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    axs[idx].yaxis.offsetText.set_fontsize(20)
    axs[idx].tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22, direction='in', left=True, right=True)
    axs[idx].yaxis.set_major_locator(ticker.MultipleLocator(20))
    
for idx, geom in enumerate(geoms):
    print(idx)
    dis = np.loadtxt('dis_{}.txt'.format(geom))

    X_dat = dis[:,0]
    Y_dat = dis[:,2]
    Z_dat = dis[:,3]
    print(X_dat.shape)
    print(Y_dat.shape)
    print(Z_dat.shape)
    
    # Convert from pandas dataframes to numpy arrays
    X, Y, Z, = np.array([]), np.array([]), np.array([])
    for i in range(len(X_dat)):
            X = np.append(X, X_dat[i])
            Y = np.append(Y, Y_dat[i])
            Z = np.append(Z, Z_dat[i])
            
    xi = np.linspace(int(X.min()), int(X.max()), 21)
    yi = np.linspace(int(Y.min()), int(Y.max()), 116)
    print(X.shape)
    print(Y.shape)
    print(Z.shape)
    
    z = griddata((X, Y), Z, (xi[None,:], yi[:,None]), method='linear', 
    rescale=True, fill_value=0.0)
    
    #print(geom)
    #print(z.min(), z.max(), '\n')
    #cb_min = 0.0
    #cb_max = 2.1874486506750372     
    #levels = MaxNLocator(nbins=600).tick_values(z.min(), z.max())

    cf = axs[idx].contourf(xi, yi, z, 50, levels=levels) #, cmap=plt.cm.jet
    axs[idx].set_title('{}'.format(lst_labels[idx]), fontsize=20)
    #cbar = fig.colorbar(cf, ax=axs[idx])#, aspect=20.0, pad=0.05)
    #cbar.ax.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=20)
    
    #plt.savefig('dis_{}.png'.format(geom))
    #plt.show()

# Hide x labels and tick labels for top plots and y ticks for right plots.
for idx, geom in enumerate(geoms):
    axs[idx].label_outer()
    
cbar = fig.colorbar(cf, ax=axs[-1], aspect=35.0, pad=0.10)
cbar.ax.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=20) #, direction='in', left=True, right=True)
#plt.subplots_adjust(right=1.00, wspace=0.00, hspace=0.0, bottom=0.0, top=0.005)
plt.tight_layout(rect=(0.0,0.00,1.00,1.0))
fig.savefig('dis_compar.png')
#plt.savefig('./compar_strains_1wv.png')
plt.show()
