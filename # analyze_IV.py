# analyze_IV.py

from __future__ import print_function
import sisl
import numpy as np
import matplotlib.pyplot as plt



fh_pristine = open('cur_lr.txt', 'r')
fh_v_2svtop = open('cur_lr_v_2svtop.txt', 'r')
fh_v_mo = open('cur_lr_v_mo.txt', 'r')


fh_pr_lines = fh_pristine.readlines()
fh_2svtop_lines = fh_v_2svtop.readlines()
fh_mo_lines = fh_v_mo.readlines()
#print(fh_pr_lines)
#print(fh_2svtop_lines)
#print(fh_mo_lines)

V_I = []
for l in fh_lines[::2]:
    V_I.append((float(l.rstrip().split()[8]), float(l.rstrip().split()[11])))
    #print(v_i)
    #V_I.append(v_i)

V_I = sorted(V_I)
print(V_I)


np_vi = np.empty((9,2))
for idx in range(len(V_I)):
    np_vi[idx, 0] = V_I[idx][0]
    np_vi[idx, 1] = V_I[idx][1]
print(np_vi)


plt.plot(np_vi[:,0], np_vi[:,1]); #* 1e9
plt.xlabel('Bias [V]'); plt.ylabel(r'Current [$\mathrm{A}$]');


