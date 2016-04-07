import numpy as np
import fimport
import sys
from split_lib import *

def find_r_last(data, field, val):
  for i in range(len(data[field])):
    if data[field][i] >= val:
      return data['radius'][i]

def find_r_first_eq(data, field, val):
  for i in range(len(data[field])):
    if data[field][len(data[field]) - i - 1] == val:
      return data['radius'][len(data[field]) - i - 1]

prof = sys.argv[1]
mass = sys.argv[2]
mass_prof = '{}_{}'.format(mass, prof)

data = fimport.load_array('data{}/profile{}.data'.format(mass, prof), 6)
max_g_data = fimport.load_array('g_max.txt', 1)

# he4 - find where drops below 99% of central fraction
r_he_core = find_r_last(data, 'he4', 0.999* data['he4'][-1])

r_h_burn = find_r_last(data, 'eps_nuc', 10.)

r_conv = find_r_first_eq(data, 'mixing_type', 1)

max_g = max_g_data[np.where(max_g_data['mass'] == strip_txt(mass))]['g_max']
print str(max_g_data['mass'])
print max_g

with open('data{}/model_data{}.txt'.format(mass, prof),'w') as f:
  f.write('he_core h_burn conv g_max\n')
  f.write('{} {} {} {}'.format(r_he_core, r_h_burn, r_conv, max_g[0]))
