import numpy as np
import matplotlib.pyplot as plt
import sys
import fimport
import consts as c
import scipy.optimize
import matplotlib.transforms as mtransforms
from split_lib import *

usage = 'python angmomt.py prof mass prof_type {param}'

if len(sys.argv) < 4:
  print usage
  exit(1)

prof = sys.argv[1]
mass = sys.argv[2]
mass_prof = '{}_{}'.format(mass, prof)
prof_type = sys.argv[3]
if len(sys.argv) == 5:
  param = float(sys.argv[4])
else:
  param = 1.

outfile_names = {'hestep':'hestep', 'lin_conv':'envdiff_lin_conv',
    'rdiff':'envdiff_{}'.format(param), 'poly':'envdiff_poly_{}'.format(param),
    'step':'step_{}'.format(param)}
if prof_type not in outfile_names.keys():
  print 'prof_type must be one of {}'.format(outfile_names.keys())
  exit(1)

line_styles = ['-.', '--', '-',':']
mode_colors = ['m', 'c', 'green', 'orange', 'red', 'k']

mode_nums = {'1.5_138':{'g1':'330', 'g2':'331', 'p1':'332', 'p2':'333', 
      'g4':'334'},
    '1.5_140':{'g1':'415', 'p1':'416', 'g2':'417', 'p2':'418'},
    '1.24_122':{'g1':'306', 'g2':'307', 'p1':'308', 'p2':'309', 'g3':'310', 
      'g4':'311', 'g5':'312', 'p3':'313', 'g6':'314'}}

mode_keys_dict = {'1.5_138':['g1', 'g2', 'p1', 'p2', 'g3'],
    '1.5_140':['g1', 'p1', 'p2', 'g2'],
    '1.24_122':['g1', 'g2', 'p1', 'p2', 'g3', 'g4', 'g5', 'p3', 'g6']}

data_fnames = {'mesa':'data{}/profile{}.data'.format(mass, prof),
    'model':'data{}/model_data{}.txt'.format(mass, prof)}

for mode in mode_nums[mass_prof]:
  data_fnames[mode] = 'data{}/summary_l1_prof{}_00{}.txt'.format(mass, prof, mode_nums[mass_prof][mode])

'''data_fnames = {'mesa':'data/profile' + prof + '.data',
    'g1':'data/summary_l1_prof140_00415.txt',
    'p1':'data/summary_l1_prof140_00416.txt',
    'p2':'data/summary_l1_prof140_00417.txt',
    'g2':'data/summary_l1_prof140_00418.txt',
    'model':'data/model_data140.txt'}'''

'''data_fnames = {'mesa':"data/profile" + prof + ".data",
    'g1':"data/summary_l1_prof150_00500.txt",
    'g2':"data/summary_l1_prof150_00900.txt",
    'p1':"data/summary_l1_prof150_01100.txt",
    'p2':"data/summary_l1_prof150_01130.txt"}'''

#mode_keys = ['g1', 'g2', 'p1', 'p2']
#mode_keys = ['g1','p1', 'p2','g2']
#mode_keys = mode_nums[prof].keys()

mode_keys = mode_keys_dict[mass_prof]

# load in data
data = {}
head_data = {}
for f in data_fnames:
  if f == 'model':
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

he_core = data['model']['he_core'] * c.rsun
h_burn = data['model']['h_burn'] * c.rsun
conv_zone = data['model']['conv'] * c.rsun
eps = 1e-6 * c.rsun

r_mesa = 10. ** data['mesa']['logR'] * c.rsun
r_star = r_mesa[-1]
m_mesa = data['mesa']['q'] * c.msun * head_data['mesa']['star_mass']
x = data['g1']['x'] # Radial grid from GYRE, normalized to [0,1]
r = x * r_star
m = np.interp(r, r_mesa, m_mesa)
omega = data['g1']['Omega_rot'] / np.sqrt(r_star ** 3. / (c.ggrav * m[-1]))

min_l = angmomt(m, r, np.ones(len(m)) * omega[-1])
calc_l = angmomt(m, r, omega)
print 1 - min_l / calc_l

# Calculate fractional portion of K in each region
print integral(x, data['p2']['ReK'], False, [0., 1.])
print 'k core p', integral(x, data['p2']['ReK'], False, [0., he_core / r_star])
print 'k core g', integral(x, data['g1']['ReK'], False, [0., he_core / r_star])
print 'k env p', integral(x, data['p2']['ReK'], False, [he_core / r_star, 1])
print 'k env g', integral(x, data['g1']['ReK'], False, [he_core / r_star, 1])
print 'k shell p', integral(x, data['p2']['ReK'], False, [he_core / r_star, h_burn / r_star])
print 'k shell g', integral(x, data['g1']['ReK'], False, [he_core / r_star, h_burn / r_star])
print 'k through rcb p', integral(x, data['p2']['ReK'], False, [0., conv_zone / r_star])
print 'k through rcb g', integral(x, data['g1']['ReK'], False, [0., conv_zone / r_star])

k_Hshell = integral(x, data['p2']['ReK'], False, [he_core / r_star, h_burn / r_star])
k_Hcore = integral(x, data['p2']['ReK'], False, [h_burn / r_star, conv_zone / r_star])
k_env = integral(x, data['p2']['ReK'], False, [conv_zone / r_star, 1.])

#print k_core, k_Hshell, k_Hcore, k_env

Deltanu = 1. / (2. * integral(r_mesa, 1. / data['mesa']['csound'], False))
print Deltanu

'''omegas = {'mesa':omega,
    'hesharp':(omega[0], omega[-1], he_core, he_core + eps),
    'heh':(omega[0], omega[-1], he_core, h_burn),
    'heconv':(omega[0], omega[-1], he_core, conv_zone),
    'hsharp':(omega[0], omega[-1], h_burn, h_burn + eps),
    'hconv':(omega[0], omega[-1], h_burn, conv_zone),
    'csharp':(omega[0], omega[-1], conv_zone, conv_zone + eps),
    'outerup':(omega[0], omega[-1] * 10., he_core, h_burn),
    'surface':(omega[0], omega[-1], conv_zone, r_star),
    'lowin':(omega[0] * 0.1, omega[-1], he_core, he_core + eps)
    }'''

omegas = {}
omega_keys = []
if prof_type == ('hestep') or prof_type == 'step':
  if prof_type == 'hestep':
    step_at = he_core
  else:
    step_at = param * r_star
  for ratio in np.logspace(0., np.log10(omega[0] / omega[-1]), 20):
    key = str(np.around(ratio, 5))
    omegas[key] = (ratio * omega[-1], omega[-1], step_at, step_at + eps)
    omega_keys.append(key)
elif prof_type == 'lin_conv':
  for ratio in np.logspace(0., np.log10(omega[0] / omega[-1]), 20):
    key = str(np.around(ratio, 5))
    omegas[key] = (ratio * omega[-1], omega[-1], conv_zone, r_star - eps)
    omega_keys.append(key)
elif prof_type == 'lin_env_rdiff':
  angmomt_ratio = param
  for rdiff in np.logspace(np.log10(conv_zone), np.log10(r_star - eps), 20):
    def angmomt_rdiff(omc):
      return angmomt(m, r, omega_prof_lin(r, omc, omega[-1], conv_zone, rdiff)) \
          - angmomt(m, r, omega) * angmomt_ratio
    best_omc = scipy.optimize.broyden1(angmomt_rdiff, omega[-1] * 1.00001)
    print rdiff / r_star, best_omc / omega[-1], angmomt_rdiff(best_omc), \
        angmomt(m, r, omega_prof_lin(r, best_omc, omega[-1], he_core, rdiff))
    key = str(np.around(best_omc / omega[-1], 5))
    omegas[key] = (best_omc, omega[-1], conv_zone, rdiff)
    omega_keys.append(key)
elif prof_type == 'poly':
  for ratio in np.logspace(0., np.log10(omega[0] / omega[-1]), 20):
    key = str(np.around(ratio, 5))
    omegas[key] = (ratio * omega[-1], omega[-1], conv_zone * param)
    omega_keys.append(key)
else:
  print "Profile type unknown. Must be {hestep, lin_conv, rdiff, poly}"
  exit(1)

outstring = ''
for key in omega_keys:
  if key != 'mesa':
    outstring += outdata(key, x, m, r, [data[mode_key]['ReK'] for mode_key in mode_keys], [head_data[mode_key]['Rebeta'] for mode_key in mode_keys], arr=False, tp=omegas[key])
  else:
    outstring += outdata(key, x, m, r, [data[mode_key]['ReK'] for mode_key in mode_keys], [head_data[mode_key]['Rebeta'] for mode_key in mode_keys], arr=True, org_omega=omegas[key])

with open('data' + mass + '/prof' + prof + '_splittings_' + outfile_names[prof_type] + '.txt','w') as f:
  if prof_type == 'poly':
    f.write('ratio angmom ' + ' '.join(mode_keys) + ' r_in alpha\n')
  else:
    f.write('ratio angmom ' + ' '.join(mode_keys) + ' r_in r_out\n')
  f.write(outstring)
