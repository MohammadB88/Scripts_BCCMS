# analyze_IV.py
from __future__ import print_function
import sisl
import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')


fh_pristine = open('cur_lr.txt', 'r')
fh_v_2stop = open('cur_lr_v_2stop.txt', 'r')
fh_v_mo = open('cur_lr_v_mo.txt', 'r')
fh_v_s = open('cur_lr_v_s.txt', 'r')
fh_v_mo3s = open('cur_lr_v_mo3s.txt', 'r')


fh_pr_lines = fh_pristine.readlines()
fh_2stop_lines = fh_v_2stop.readlines()
fh_mo_lines = fh_v_mo.readlines()
fh_s_lines = fh_v_s.readlines()
fh_mo3s_lines = fh_v_mo3s.readlines()
#print(fh_pr_lines)
#print(fh_2svtop_lines)
#print(fh_mo_lines)
#print(fh_s_lines)
#print(fh_mo3s_lines)


V_I_pr = []
for l in fh_pr_lines[::2]:
    V_I_pr.append((float(l.rstrip().split()[8]), float(l.rstrip().split()[11])))
V_I_pr = sorted(V_I_pr)
#print(V_I_pr)
np_vi_pr = np.empty((len(V_I_pr),2))
for idx in range(len(V_I_pr)):
    np_vi_pr[idx, 0], np_vi_pr[idx, 1] = V_I_pr[idx][0], V_I_pr[idx][1]
print(np_vi_pr)


V_I_s= []
for l in fh_2stop_lines[::2]:
    V_I_s.append((float(l.rstrip().split()[8]), float(l.rstrip().split()[11])))
V_I_s = sorted(V_I_s)
#print(V_I_s)
np_vi_s = np.empty((len(V_I_s),2))
for idx in range(len(V_I_s)):
    np_vi_s[idx, 0], np_vi_s[idx, 1] = V_I_s[idx][0], V_I_s[idx][1]
print(np_vi_s)


V_I_2stop = []
for l in fh_2stop_lines[::2]:
    V_I_2stop.append((float(l.rstrip().split()[8]), float(l.rstrip().split()[11])))
V_I_2stop = sorted(V_I_2stop)
#print(V_I_2stop)
np_vi_2stop = np.empty((len(V_I_2stop),2))
for idx in range(len(V_I_2stop)):
    np_vi_2stop[idx, 0], np_vi_2stop[idx, 1] = V_I_2stop[idx][0], V_I_2stop[idx][1]
print(np_vi_2stop)


V_I_mo = []
for l in fh_mo_lines[::2]:
    V_I_mo.append((float(l.rstrip().split()[8]), float(l.rstrip().split()[11])))
V_I_mo = sorted(V_I_mo)
#print(V_I_mo)
np_vi_mo = np.empty((len(V_I_mo),2))
for idx in range(len(V_I_mo)):
    np_vi_mo[idx, 0], np_vi_mo[idx, 1] = V_I_mo[idx][0], V_I_mo[idx][1]
print(np_vi_mo)


V_I_mo3s = []
for l in fh_mo_lines[::2]:
    V_I_mo3s.append((float(l.rstrip().split()[8]), float(l.rstrip().split()[11])))
V_I_mo3s = sorted(V_I_mo3s)
#print(V_I_mo3s)
np_vi_mo3s = np.empty((len(V_I_mo3s),2))
for idx in range(len(V_I_mo3s)):
    np_vi_mo3s[idx, 0], np_vi_mo3s[idx, 1] = V_I_mo3s[idx][0], V_I_mo3s[idx][1]
print(np_vi_mo3s)


#plt.plot(np_vi_pr[:,0], np_vi_pr[:,1]* 1e15); #* 1e9
#plt.plot(np_vi_s[:,0], np_vi_s[:,1]); #* 1e9
#plt.plot(np_vi_2stop[:,0], np_vi_2stop[:,1]); #* 1e9
plt.plot(np_vi_mo[:,0], np_vi_mo[:,1]); #* 1e9
#plt.plot(np_vi_mo3s[:,0], np_vi_mo3s[:,1]); #* 1e9
#plt.xlabel('Bias [V]'); plt.ylabel(r'Current [$\mathrm{A}$]');
#plt.ylim(-0.2,0.2,0.002)

