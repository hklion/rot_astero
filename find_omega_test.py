import numpy as np
import sys
import fimport
import consts as c
from scipy.interpolate import InterpolatedUnivariateSpline
import scipy.optimize
from split_lib import *
import matplotlib.pyplot as plt

prof = sys.argv[1]

measured = {'g1':1.846938e-7, 'p2':1.3163265e-7}

data_fnames = {'mesa':'data/profile' + prof + '.data',
    'step':'data/prof' + prof + '_splittings_hestep.txt',
    'poly':'data/prof' + prof + '_splittings_envdiff_poly_1.0.txt'}

profile_names = ['step', 'poly']

# load in data
data = {}
head_data = {}

for f in data_fnames:
  if f == 'model':
    data[f] = fimport.load_array(data_fnames[f], 1)
  elif f in profile_names:
    data[f] = fimport.load_array(data_fnames[f], 1)
  else:
    file_data = fimport.load_array(data_fnames[f], 6)
    if f == 'mesa':
      file_data = np.flipud(file_data)
      header_row = 2
    else:
      header_row = 3
    data[f] = file_data
    file_head_data = fimport.load_header(data_fnames[f], header_row)
    head_data[f] = file_head_data

omega_e = data['mesa']['omega'][-1]

f_g = InterpolatedUnivariateSpline(data['poly']['ratio'],
        data['poly']['g1'] / omega_e, k=1)
f_p = InterpolatedUnivariateSpline(data['poly']['ratio'],
        data['poly']['p2'] / omega_e)
plt.axhline(measured['g1'])
plt.axhline(measured['p2'])

array = np.linspace(-5, 10)

plt.plot(array, f_g(array) * -1.19632e-6)
plt.plot(array, f_p(array) * -1.19632e-6)

plt.show()
