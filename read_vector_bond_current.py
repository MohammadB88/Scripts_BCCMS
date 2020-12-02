import sisl 
import numpy as np
from scipy import sparse
#from sisl.io import get_sile, BaseSile

# Get TBT file
tbt = sisl.get_sile('MoS2_dev.TBT.nc')

dr = './vec_cur'
energies = [-1.5, -1.25, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5] 

for eng in energies:
    # read in and write the vector currents from Left to the Right electrode
    J = tbt.vector_current('Left', eng, only='+')
    np.savetxt("{}/Left_{}_pos.txt".format(dr, eng), J)

    J = tbt.vector_current('Left', eng, only='-')
    np.savetxt("{}/Left_{}_neg.txt".format(dr, eng), J)

    J = tbt.vector_current('Left', eng, only='all')
    np.savetxt("{}/Left_{}_all.txt".format(dr, eng), J)
    print(J[400, :]) 
    # read in and write the vector currents from Right to the Left electrode
    J = tbt.vector_current('Right', eng, only='+')
    np.savetxt("{}/Right_{}_pos.txt".format(dr, eng), J)

    J = tbt.vector_current('Right', eng, only='-')
    np.savetxt("{}/Right_{}_neg.txt".format(dr, eng), J)

    J = tbt.vector_current('Right', eng, only='all')
    np.savetxt("{}/Right_{}_all.txt".format(dr, eng), J)
    print(J[400, :]) 


# Get X and Y coordinates
geom = tbt.read_geometry()
coordinates = geom.xyz[:, :3].T
np.savetxt('{}/coords.txt'.format(dr), coordinates)


#geom = tbt.geometry
#output = sisl.get_sile('geom.nc', 'w')
#print(type(geom))
#w_geom = sisl.get_sile('geometry.nc')
#geom = tbt.read_geometry()
#output.write_geometry(geom)

# Retrieve bond-currents
# only = '+' only local currents are analyzed
# only = '-' only local currents are analyzed
# only = all transmission between a subset of orbitals is analyzed
#J = tbt.bond_current('Left', -1.00, only='all')
#J = tbt.bond_current('Left', -1.00)
#print(type(J))
#print(J[400, :].data[J[400, :].data[:] > 0.2000E-10]) 
#print(J[400, :])

#sparse.save_npz("JJJ_all.csr", J)

#J = tbt.bond_current('Left', -1.00, only='-')
#print(type(J))
#print(J[400, :].data[:]) 
#sparse.save_npz("JJJ_neg.csr", J)

#J = tbt.bond_current('Left', -1.00, only='all')
#print(type(J))
#print(J[400, :].data[:]) 

#J = tbt.bond_current('Left', -1.00, only='+')#, uc=True)
#print(type(J))
#print('Here is the next \n', J[400, :]) 
#J = tbt.bond_current('Left', -1.00, only='+', uc=True, kavg=)
#print(type(J))
#print(J[400, :].data[:]) 

#sparse.save_npz("JJJ_pos.csr", J)
