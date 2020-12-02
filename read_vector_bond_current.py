import sisl 
import numpy as np
from scipy import sparse
#from sisl.io import get_sile, BaseSile

# Get TBT file
tbt = sisl.get_sile('MoS2_dev.TBT.nc')

dr = './vec_cur'
biass = [-1.5, -1.0, -0.5, 0.0, 0.5, 1.4] 

for bias in biass:
    J = tbt.vector_current('Left', bias, only='+')
    print(J[400, :]) 
    np.savetxt("{}/EV_{}_pos.txt".format(dr, bias), J)

    J = tbt.vector_current('Left', bias, only='-')
    print(J[400, :]) 
    np.savetxt("{}/EV_{}_neg.txt".format(dr, bias), J)

    J = tbt.vector_current('Left', bias, only='all')
    print(J[400, :]) 
    np.savetxt("{}/EV_{}_all.txt".format(dr, bias), J)


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
