import numpy as np
import fimport
import sys
import matplotlib.pyplot as plt

prof = sys.argv[1]
mass = sys.argv[2]

data = fimport.load_array('data{}/measured_modes.txt'.format(mass), 1)

max_mode_i = np.max(data['mode_i'])

split_m = []
split_p = []
for i in range(int(max_mode_i + 1)):
  mode_data = data[np.where(data['mode_i'] == i)]
  split_m.append((mode_data['freq'][np.where(mode_data['m'] == 0)] - \
      mode_data['freq'][np.where(mode_data['m'] == -1)])[0])
  split_p.append((mode_data['freq'][np.where(mode_data['m'] == 1)] - \
      mode_data['freq'][np.where(mode_data['m'] == 0)])[0])

plt.plot(data['freq'][np.where(data['m'] == 0)], split_p, '.')
plt.plot(data['freq'][np.where(data['m'] == 0)], split_m, '.')
plt.show()

with open('data{}/measured_splits.txt'.format(mass),'w') as f:
  f.write('freq split_-1 split_1 split_ave\n')
  for i in range(len(split_m)):
    f.write('{} {} {} {}'.format(data['freq'][np.where(data['m'] == 0)][i], 
      split_m[i], split_p[i], (split_m[i] + split_p[i]) / 2) + '\n')

