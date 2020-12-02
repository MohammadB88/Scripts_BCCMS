import sisl
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# read in the coordinates
x, y, z = np.loadtxt('coords.txt')

# read in the vector currents
J = np.loadtxt('vect_cur_pos.txt')

# Create plot
fig = plt.figure()
ax = fig.add_subplot(111)
#ax = fig.gca(projection='3d')


#ax.quiver(x, y, z, J[:, 0], J[:, 1], J[:,2]);
#ax.scatter(x,y,z, color = 'black')
ax.quiver(x, z, J[:, 0], J[:,2], color='blue');
ax.scatter(x, z, color = 'black')

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

