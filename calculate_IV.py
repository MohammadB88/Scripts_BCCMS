import sisl 
import matplotlib.pyplot as plt
import numpy as np

# ***********************************************************************************************
# This would plot the current obtained by calculating all the non-equilibrium Hamiltonians (NEHs)

V_NEHs = [-1.0, -0.5, -0.2, 0.0, 0.2, 0.5, 1.0]
I_NEHs = np.empty([len(V_NEHs)])
for i, v in enumerate(V_NEHs):
    I_NEHs[i] = sisl.get_sile('TS_{:.1f}/MoS2.TBT.nc'.format(v)).current()
    print(I_NEHs[i])

plt.plot(V_NEHs, I_NEHs * 1e9, label='non-equilibrium Hamiltonians')
plt.xlabel('Bias [V]'); plt.ylabel(r'Current [$\mathrm{nA}$]'); #plt.ylabel(r'Current [$\mu\mathrm{A}$]');
plt.legend(loc='upper right')

# ***********************************************************************************************
# This would plot the photocurrent obtained by calculating all the non-equilibrium Hamiltonians (NEHs)

I_ph = np.empty([len(V_NEHs)])
for i, v in enumerate(V_NEHs):
    I_ph[i] = sisl.get_sile('phcurrent_TBT_{:.1f}/MoS2.TBT.nc'.format(v)).current()
    print(I_ph[i])

plt.plot(V_NEHs, I_ph * 1e9, label='photocurrent')
plt.xlabel('Bias [V]'); plt.ylabel(r'Current [$\mathrm{nA}$]'); #plt.ylabel(r'Current [$\mu\mathrm{A}$]');
plt.legend(loc='upper right')

# ***********************************************************************************************
# This would plot the current obtained by interpolating the Hamiltonians for all the biases 
# using the calculated non-equilibrium Hamiltonians (NEHs) (InterH)

V_InterH = np.arange(-1.0,1.05,0.1)
I_InterH = np.empty([len(V_NEHs)])
for i, v in enumerate(V_InterH):
    I_InterH[i] = sisl.get_sile('current_interpolation/TBT_{:.1f}/MoS2.TBT.nc'.format(v)).current()
    print(I_InterH[i])

plt.plot(V_InterH, I_InterH * 1e9, label='interpolated current')
plt.xlabel('Bias [V]'); plt.ylabel(r'Current [$\mathrm{nA}$]'); #plt.ylabel(r'Current [$\mu\mathrm{A}$]');
plt.legend(loc='upper right')

plt.show()
