import sisl
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

dr = 'vec_cur'
bias = 1.40
#energies = [-1.5, -1.25, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5]
energies = [-1.5, -0.5, 0.0, 0.5, 1.5] 

# read in the coordinates
x, y, z = np.loadtxt('{}/coords.txt'.format(dr))

# Create plot
fig , ax = plt.subplots(1, len(energies), figsize=(16,8), sharey=True)

#ax = fig.gca(projection='3d')
idx = 0
for eng in energies:
    print(idx, eng)
    # read in the vector currents from Left to Right electrode
    J_l_pos = np.loadtxt('{}/Left_{}_pos.txt'.format(dr, eng))
    J_l_neg = np.loadtxt('{}/Left_{}_neg.txt'.format(dr, eng))
    J_l_all = np.loadtxt('{}/Left_{}_all.txt'.format(dr, eng))
    # read in the vector currents from Right to Left electrode
    J_r_pos = np.loadtxt('{}/Right_{}_pos.txt'.format(dr, eng))
    J_r_neg = np.loadtxt('{}/Right_{}_neg.txt'.format(dr, eng))
    J_r_all = np.loadtxt('{}/Right_{}_all.txt'.format(dr, eng))
    
    print('Inputs have been read!')
    #ax.quiver(x, y, z, J[:, 0], J[:, 1], J[:,2]);
    
    # plot the vector currents from Left to Right electrode
    #ax[idx].quiver(x, z, J_l_pos[:, 0], J_l_pos[:,2], color='blue', width=0.009)#, lw=25, pivot='mid')
    #ax[idx].quiver(x, z, J_l_neg[:, 0], J_l_neg[:,2], color='red', width=0.009)#, lw=25, pivot='mid');
    ax[idx].quiver(x, z, J_l_all[:, 0], J_l_all[:,2], color='cyan', width=0.009)#, lw=25, pivot='mid');
    
    # plot the vector currents from Right to Left electrode
    #ax[idx].quiver(x, z, J_r_pos[:, 0], J_r_pos[:,2], color='blue', width=0.009)#, lw=25, pivot='mid')
    #ax[idx].quiver(x, z, J_r_neg[:, 0], J_r_neg[:,2], color='red', width=0.009)#, lw=25, pivot='mid');
    ax[idx].quiver(x, z, J_r_all[:, 0], J_r_all[:,2], color='red', width=0.009)#, lw=25, pivot='mid');
        
    # Add atoms as points to the plots
    ax[idx].scatter(x, z, color = 'black', s=5)
    
    #ax[idx].set_xlim()
    
    print('plots have been produced')
    ax[idx].set_title('At Energy {} eV'.format(eng))
    print('title is set')
    idx += 1
    
    #plt.pause(10)
    #plt.pause(0.0001)

plt.suptitle('pristine device at Bias {} eV'.format(bias))

plt.subplots_adjust(right=1.00, wspace=0.05, hspace=0.05, bottom=0.05, top=0.4)
plt.tight_layout()
plt.savefig('vector_current_B1.40.png')
plt.show()

# *********************************************
# plot bond currents on each atom in the device
# *********************************************

# read in the bond_current from sparse matrix
#J = sparse.load_npz('JJJ_pos.csr.npz')

#ia = 400
#print(J[ia, :].data[:])
#print(J[ia, :])

# read in the coordinates
#fh_geom = sisl.get_sile('geom.nc')
#geom = fh_geom.read_geometry()
#print(type(geom))
#print(geom)

# Get X and Y coordinates
#x, y, z = geom.xyz[:, :3].T
#print(z[ia])
#print(geom.axyz(J[ia, :].indices)[:, 2], geom.axyz(J[ia, :].indices)[:, 2] - z[ia])

# Add geometry
#ax.scatter(x, z)

#def r(N, x):
    #return np.repeat(x, N)

#scale = 1000000
#Jia = J[ia, :] #J[ia, :].data[:] > 0.2000E-10]
##print(ia, Jia)
##if Jia.nnz == 0:
    ### There are no bond-currents here
    ##continue
## Calculate vectors
#Rx = geom.axyz(Jia.indices)[:, 0] - x[ia]
##print(Rx)
##Ry = geom.axyz(Jia.indices)[:, 1] - y[ia]
#Rz = geom.axyz(Jia.indices)[:, 2] - z[ia]
##print(Rz)

##d = (Rx ** 2 + Ry ** 2) ** 0.5
#d = (Rx ** 2 + Rz ** 2) ** 0.5
#Rx *= Jia.data[:] / d * scale
##print(Rx)
##Ry *= Jia.data[:] / d * scale
#Rz *= Jia.data[:] / d * scale
##print(Rz)

#N = len(Rx)
##ax.quiver(r(N, x[ia]), r(N, y[ia]), Rx, Ry, Jia.data[:] / Jia.data.max())
#ax.quiver(r(N, x[ia]), r(N, z[ia]), Rx, Rz, Jia.data[:] / Jia.data.max())

## Add quiver plots
#scale = 0.001
#for ia in geom:
    #Jia = J[ia, :]
    ##print(ia, Jia)
    #if Jia.nnz == 0:
        ## There are no bond-currents here
        #continue
    ## Calculate vectors
    #Rx = geom.axyz(Jia.indices)[:, 0] - x[ia]
    ##Ry = geom.axyz(Jia.indices)[:, 1] - y[ia]
    #Rz = geom.axyz(Jia.indices)[:, 2] - z[ia]
    ##d = (Rx ** 2 + Ry ** 2) ** 0.5
    #d = (Rx ** 2 + Rz ** 2) ** 0.5
    #Rx *= Jia.data[:] / d * scale
    ##Ry *= Jia.data[:] / d * scale
    #Rz *= Jia.data[:] / d * scale

    #N = len(Rx)
    ##ax.quiver(r(N, x[ia]), r(N, y[ia]), Rx, Ry, Jia.data[:] / Jia.data.max())
    #ax.quiver(r(N, x[ia]), r(N, z[ia]), Rx, Rz, Jia.data[:] / Jia.data.max())

