import numpy as np
import matplotlib.pyplot as plt
import sys
import fimport
import consts as c
from split_lib import *
import astero_param as astero
#import seaborn as sns
#sns.set_style("ticks", {"xtick.direction": u'in', "ytick.direction": u"in"})

def label_gen(first, last):
  return {label:str(mode).zfill(4) for (label, mode) in (('g_{}'.format(i), i) for i in range(first, last))}

prof = sys.argv[1]
mass = sys.argv[2]
mass_prof = '{}_{}'.format(mass, prof)

mode_nums = {'1.24_122':{label:'0{}'.format(mode) for (label, mode) in (('g_{}'.format(i), i) for i in range(300, 325))},
    '1.5_138':{label:'0{}'.format(mode) for (label, mode) in (('g_{}'.format(i), i) for i in range(310, 340))},
    '1.5_150':label_gen(920, 950),
    '1.5_160':label_gen(1500, 1550)}

vel_e_all = {'1.24':np.linspace(2, 8, 5),
    '1.5':np.linspace(0.2, 1, 5)}

prof_obs = {'1.24':'122', '1.5':'138'}

data_fnames = {'mesa':'data{}/profiles/profile{}.data'.format(mass, prof),
    'mesa_obs':'data{}/profiles/profile{}.data'.format(mass, prof_obs[mass]),
    'model':'data{}/model_data/model_data{}.txt'.format(mass, prof)}

prof_types = ['hestep', 'poly']

for mode in mode_nums[mass_prof]:
  data_fnames[mode] = 'data{0}/prof{1}/summary_l1_prof{1}_0{2}.txt'.format(mass, prof, mode_nums[mass_prof][mode])

mode_keys = mode_nums[mass_prof].keys()
mode_keys.sort()

# load in data
data = {}
head_data = {}
for f in data_fnames:
  if f == 'model':
    data[f] = fimport.load_array(data_fnames[f], 1)
  else:
    file_data = fimport.load_array(data_fnames[f], 6)
    if f == 'mesa' or f == 'mesa_obs':
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
r_star_obs = (10. ** data['mesa_obs']['logR'] * c.rsun)[-1]
m_mesa = data['mesa']['q'] * c.msun * head_data['mesa']['star_mass']
x = data[mode_keys[0]]['x'] # Radial grid from GYRE, normalized to [0,1]
r = x * r_star
m = np.interp(r, r_mesa, m_mesa)

large_sep = astero.large_sep(data['mesa']) * 1e6
numax = astero.numax(data['mesa'], head_data['mesa'])

omega_c = omega_c_est(data['model']['g_max']) * 1e-6
vel_e = vel_e_all[mass] # km/s
omega_e = vel_e * 1e5 / r_star_obs
ratios = np.around(omega_c / omega_e, 5)

print 'ratios', omega_c / omega_e

# list of mode frequencies
freqs = [head_data[mode]['Refreq'] for mode in mode_keys]

omega_keys = [(ratio, prof_type) for ratio in ratios for prof_type in prof_types]
print omega_keys

splits = {}
for key in omega_keys:
  splits[key] = []
  (ratio, prof_type) = key
  if prof_type == 'hestep':
    omega_prof = omega_prof_lin(r,
        omega_c, omega_c / ratio, he_core, he_core + eps)
  elif prof_type == 'poly':
    omega_prof = omega_prof_pow(r,
        omega_c, omega_c / ratio, conv_zone)
  else:
    print prof_type
    print 'uhhhh'
  for mode in mode_keys:
    split_val = split(x, data[mode]['ReK'], head_data[mode]['Rebeta'], omega_prof)
    splits[key].append(split_val)

print freqs, large_sep

# Holding ratio constant...

colors = sns.color_palette()
tmp = colors[1]
colors[1] = colors[2]
colors[2] = tmp

for ratio in ratios:
  to_plot = [(ratio, 'hestep'), (ratio, 'poly')]


  for (i,key) in enumerate(to_plot):
    plt.plot(freqs - 1, splits[key] / max(splits[key]),'-o', c=colors[i], label=key[1])

  plt.ylim(0.3, 1.0)

  plt.legend(loc = 4)
  plt.title('mass {}, profile {}, ratio {}'.format(mass, prof, to_plot[0][0]))
  plt.savefig('data{}/plots/folded_splitting_{}_r{}.pdf'.format(mass, prof, np.around(ratio,2)))
  plt.clf()

  plt.axvline(numax, c='k',linestyle=":")
  for (i, key) in enumerate(to_plot):
    plt.plot(freqs, splits[key] / max(splits[key]),linestyle='-',marker='.', label=key[1],linewidth=1.5, c=colors[i])
  plt.ylim(0.3, 1.0)
  plt.legend(loc=4)
  plt.title('mass {}, profile {}, ratio {}'.format(mass, prof, to_plot[0][0]))
  plt.savefig('data{}/plots/splitting_{}_r{}.pdf'.format(mass, prof, np.around(ratio,2)))
  plt.clf()
