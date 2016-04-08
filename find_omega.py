import numpy as np
import sys
import fimport
import consts as c
from scipy.interpolate import InterpolatedUnivariateSpline
import scipy.optimize
from split_lib import *
import matplotlib.pyplot as plt

prof = sys.argv[1]
mass = sys.argv[2]
p = sys.argv[3]
g = sys.argv[4]

if mass == '1.24':
  if p == 'p3' and g == 'g5':
    initial_guess = (20, 1e-7)
  if p == 'p1' and g == 'g2':
    initial_guess = (35, 2.3e-7)
elif mass == '1.5':
  if prof == '140':
    if p == 'p2' and g == 'g1':
      initial_guess = (30, 2e-7)
  elif prof == '138':
    if p == 'p2' and g == 'g1':
      initial_guess = (60, 8e-8)
else:
  initial_guess = (10, 1e-7)

measured_pairs = []
#g_vals = np.linspace(0.17e-6, 0.25e-6)
#p_vals = np.linspace(0.125e-6, 0.15e-6)
p_vals = np.linspace(0.1e-6, 0.17e-6)
g_vals = np.linspace(0.15e-6, 0.3e-6)
for g_val in g_vals:
  for p_val in p_vals:
    #measured_pairs.append({g:g_val, p:p_val})
    measured_pairs.append({g:g_val, p:p_val})

data_fnames = {'mesa':'data{}/profiles/profile{}.data'.format(mass, prof),
    'step':'data{}/prof{}_splittings_hestep.txt'.format(mass, prof),
    #'0.02step': 'data{}/prof{}_splittings_step_0.02.txt'.format(mass, prof),
    #'0.001step': 'data{}/prof{}_splittings_step_0.001.txt'.format(mass, prof),
    'poly':'data{}/prof{}_splittings_envdiff_poly_1.0.txt'.format(mass, prof)}

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
print omega_e

for profile in profile_names:
  with open('data{}/find_omega_{}_{}_{}_{}.txt'.format(mass, prof, profile, p, g), 'w') as fi:
    fi.write('{} {} ratio omega_e\n'.format(g, p))
    f_g = InterpolatedUnivariateSpline(data[profile]['ratio'],
        data[profile][g] / omega_e, k=1)
    f_p = InterpolatedUnivariateSpline(data[profile]['ratio'],
        data[profile][p] / omega_e)
    for measured in measured_pairs:
      def f(pair):
        ratio, surf = pair
        return (f_g(ratio) * surf - measured[g], f_p(ratio) * surf - measured[p])
      soln = scipy.optimize.fsolve(f, initial_guess)
      outstring = '{0} {1} {2} {3}\n'.format(measured[g], measured[p],
          soln[0], soln[1])
      fi.write(outstring)
