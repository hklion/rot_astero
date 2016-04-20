import numpy as np
import matplotlib.pyplot as plt
import sys
import fimport
import consts as c
from split_lib import *
import astero_param as astero
import seaborn.apionly as sns
#sns.set_style("ticks", {"xtick.direction": u'in', "ytick.direction": u"in"})



prof = sys.argv[1]
mass = sys.argv[2]
mass_prof = '{}_{}'.format(mass, prof)

mode_nums = {'1.24_122':{label:'0{}'.format(mode) for (label, mode) in (('g_{}'.format(i), i) for i in range(300, 325))},
    '1.24_133':label_gen(620,700),
    '1.24_148':label_gen(1360,1470),
    '1.5_138':{label:'0{}'.format(mode) for (label, mode) in (('g_{}'.format(i), i) for i in range(310, 348))},
    '1.5_150':label_gen(890, 970),
    '1.5_160':label_gen(1400, 1600),
    '1.5_129':label_gen(90, 110),
    '1.33_124':label_gen(260, 290),
    '1.33_123':label_gen(240, 260),
    '1.33_122':label_gen(210, 235),
    '1.33_1198':label_gen(210, 240),
    '1.33_1204':label_gen(220, 260),
    '1.33_1208':label_gen(230, 270),
    '1.33_1215':label_gen(250, 280),
    '1.33_norot_1208':label_gen(220,250),
    '1.33_norot_1214':label_gen(230,260),
    '1.33_norot_1220':label_gen(240,270),
    '1.33_norot_1227':label_gen(250,300)}

vel_e_all = {'1.24':np.linspace(2, 8, 5),
    '1.5':np.linspace(0.2, 1, 5),
    '1.33':np.linspace(2.5, 3.2, 5)}

prof_obs = {'1.24':'122', '1.5':'138', '1.33':'124'}

data_fnames = {'mesa':'data{}/profiles/profile{}.data'.format(mass, prof),
    'mesa_obs':'data{}/profiles/profile{}.data'.format(strip_txt(mass), 
      prof_obs[str(strip_txt(mass))]),
    'model':'data{}/model_data/model_data{}.txt'.format(mass, prof)}

prof_types = ['step', 'power law', '0.02R$_*$ step', 'conv step', 'core diff']

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
vel_e = vel_e_all[str(strip_txt(mass))] # km/s
omega_e = vel_e * 1e5 / r_star_obs
ratios = np.around(omega_c / omega_e, 5)

print 'ratios', omega_c / omega_e

# list of mode frequencies
freqs = [head_data[mode]['Refreq'] for mode in mode_keys]

omega_keys = [(ratio, prof_type) for ratio in ratios for prof_type in prof_types]

print he_core / r_star

splits = {}
for key in omega_keys:
  splits[key] = []
  (ratio, prof_type) = key
  if prof_type == 'step':
    omega_prof = omega_prof_lin(r,
        omega_c, omega_c / ratio, he_core, he_core + eps)
  elif prof_type == 'power law':
    omega_prof = omega_prof_pow(r,
        omega_c, omega_c / ratio, conv_zone)
  elif prof_type == '0.02R$_*$ step':
    omega_prof = omega_prof_lin(r,
        omega_c, omega_c / ratio, 0.02 * r_star, 0.02 * r_star + eps)
  elif prof_type == 'conv step':
    omega_prof = omega_prof_lin(r, omega_c, omega_c / ratio, conv_zone,
        conv_zone + eps)
  elif prof_type == 'core diff':
    omega_prof = omega_prof_pow(r, omega_c, omega_c / ratio, he_core / 2, he_core)
  elif prof_type == 'tanh':
    omega_prof = omega_prof_tanh(r, omega_c, omega_c / ratio, 0.02 * r_star, eps)
  else:
    print prof_type
    print 'uhhhh'
  for mode in mode_keys:
    split_val = split(x, data[mode]['ReK'], head_data[mode]['Rebeta'], omega_prof)
    splits[key].append(split_val)

# Holding ratio constant...

colors = sns.color_palette()
tmp = colors[1]
colors[1] = colors[2]
colors[2] = tmp

for ratio in ratios:
  conv_step_split_ratio = np.array(splits[(ratio, 'power law')]) / np.array(splits[(ratio, 'step')])
  print max(conv_step_split_ratio)

  to_plot = [(ratio, 'step'), (ratio, 'power law'), (ratio, '0.02R$_*$ step'),
      (ratio, 'conv step'), (ratio, 'core diff')]


  for (i,key) in enumerate(to_plot):
    plt.plot(np.mod(freqs - 1 * np.ones(len(freqs)), large_sep), splits[key] / max(splits[key]),'.', c=colors[i], label=key[1])

  plt.ylim(0.3, 1.0)

  plt.legend(loc = 4)
  plt.title('mass {}, profile {}, ratio {}'.format(mass, prof, to_plot[0][0]))
  plt.savefig('data{}/plots/folded_splitting_{}_r{}.pdf'.format(mass, prof, np.around(ratio,2)))
  plt.clf()

  plt.axvline(numax, c='k',linestyle=":", label=r'$\nu_\mathrm{max}$')
  for (i, key) in enumerate(to_plot):
    plt.plot(freqs, splits[key] / max(splits[key]),linestyle='-',
        marker=str((i % 4) + 1), markersize=10, label=key[1],linewidth=1, c=colors[i])
  if strip_txt(mass) == 1.33:
    measured_split = fimport.load_array('data{}/measured_splits.txt'.format(strip_txt(mass)), 1)
    plt.plot(measured_split['freq'],
        measured_split['split_ave'] / max(measured_split['split_ave']),
        label="measured", linewidth=1, c=colors[-2], linestyle='-', marker='.')
  plt.ylim(0.20, 1.0)
  plt.xlim(numax - 4 * large_sep, numax + 4 * large_sep)
  plt.legend(loc=3,ncol=2,fontsize=12)
  plt.title('mass {}, profile {}, ratio {}'.format(mass, prof, to_plot[0][0]))
  plt.savefig('data{}/plots/splitting_{}_r{}.pdf'.format(mass, prof, np.around(ratio,2)))
  plt.clf()
