import numpy as np
import fimport
import matplotlib.pyplot as plt
import sys

prof = sys.argv[1]
mass = sys.argv[2]

data = fimport.load_array('data{0}/prof{1}/summary_l1_prof{1}.txt'.format(mass, 
  prof), 6)
colors = ['m', 'c', 'green', 'orange', 'red', 'k']

mode_freqs = []

if mass == '1.5' and prof == '140':
  mode_freqs = [0.1369714754136950e3, # 415
      0.1385008651627E+003, # 416
      0.1393135703489933E+003, # 417
      0.1407849468287591E+003] # 418
  i_mode_freqs = np.searchsorted(data['Refreq'], mode_freqs)
  numax = 144
elif prof == '150' and mass == '1.5':
  numax = 77
elif prof == '138' and mass == '1.5':
  numax = 170
elif prof == '160' and mass == '1.5':
  numax = 52
elif prof == '122' and mass == '1.24':
  numax = 212
else:
  numax = 0

T = 1 / data['Refreq']

dT = (-T[1:] + T[:-1]) * 1e6

plt.plot(data['Refreq'][:-1], dT, '-', linewidth=2, zorder=1)
for i in range(len(mode_freqs)):
  plt.scatter([mode_freqs[i]], [dT[i_mode_freqs[i]]], c=colors[i], s=30, zorder=2)
plt.axvline(numax, linestyle=':', c='k', linewidth=2)

plt.xlim([numax - 70, numax + 70])
plt.xlabel(r'$\nu$ [$\mu$Hz]')
plt.ylabel('$\Delta$ P [s]')

plt.savefig('data{0}/plots/delta_P_{1}.pdf'.format(mass, prof))
