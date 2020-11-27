import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter

import numpy as np

# read in the IV inputs
f_pristine = np.loadtxt('pristine/MoS2_dev.TBT.AVTRANS_Left-Right', unpack=True)
f_s = np.loadtxt('int_LR_v_s/MoS2_dev.TBT.AVTRANS_Left-Right', unpack=True)
f_2stop = np.loadtxt('int_LR_v_2stop/MoS2_dev.TBT.AVTRANS_Left-Right', unpack=True)
f_2spar = np.loadtxt('int_LR_v_2spar/MoS2_dev.TBT.AVTRANS_Left-Right', unpack=True)
#f_mo = np.loadtxt('int_LR_v_mo/MoS2_dev.TBT.AVTRANS_Left-Right', unpack=True)
#f_mo3s = np.loadtxt('int_LR_v_mo3s/MoS2_dev.TBT.AVTRANS_Left-Right', unpack=True)
#f_mo6s = np.loadtxt('int_LR_v_mo6s/MoS2_dev.TBT.AVTRANS_Left-Right', unpack=True)

# parameters for the plot
fig = plt.figure(figsize=[20, 12])
plt.xlim(-2.0, 2.0)
plt.xlabel('Energy (eV)', fontsize=20)
plt.ylim(0.0, 3.5)
plt.ylabel('Transmission', fontsize=20) #($10^{-10}$)

#plt.xticks(np.range(0.0,1.51, step=0.1)

plt.axes().xaxis.set_major_locator(ticker.MultipleLocator(0.4))
plt.tick_params(axis='x', which='major', width=2.00, length=5.0, labelsize=20)
plt.axes().xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
plt.tick_params(axis='x', which='minor', width=1.00, length=3.5, labelsize=14)
plt.axes().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
plt.axes().yaxis.offsetText.set_fontsize(20)
plt.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)
#plt.axes().yaxis.set_major_locator(ticker.MultipleLocator(0.2))
plt.axes().yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
plt.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)


# plot the IV curves
#scale_factor = 1e10
#print(scale_factor)
plt.plot(f_pristine[0], f_pristine[1], color='black', label='pristine')
plt.plot(f_s[0], f_s[1], color='blue', label='$V_{S}$')
plt.plot(f_2stop[0], f_2stop[1], color='red', label='$V_{2Stop}$')
plt.plot(f_2spar[0], f_2spar[1], color='green', label='$V_{2Spar}$')
#plt.plot(f_mo[0], f_mo[1], color='magenta', label='$V_{Mo}$')
#plt.plot(f_mo3s[0], f_mo3s[1], color='cyan', label='$V_{Mo+3S}$')
#plt.plot(f_mo6s[0], f_mo6s[1], color='orange', label='$V_{Mo+6s}$')

plt.legend(fontsize=23, loc='upper right')
plt.savefig('comparision_transmission_S_vacancies.png')

plt.show()
