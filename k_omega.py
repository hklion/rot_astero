import numpy as np
import matplotlib.pyplot as plt
import sys
import fimport
import consts as c
import scipy.optimize
import matplotlib.transforms as mtransforms
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from astero_param import *
from split_lib import *

prof = sys.argv[1]
mass = sys.argv[2]
mass_prof = '{}_{}'.format(mass, prof)

def mode_str(model, i):
  if model == '1.5_160':
    if i in [1534,1533]:
      m_type = 'p'
    else:
      m_type = 'g'
  elif model == '1.5_150':
    if i in [931,933]: #[939, 940]:
      m_type = 'p'
    else:
      m_type = 'g'
  elif model == '1.24_148':
    if i in [1403, 1405]:
      m_type = 'p'
    else:
      m_type = 'g'
  return '{}_{}'.format(m_type, i)

line_styles = ['-.', '--', '-',':']
mode_colors = {'p':'m', 'g':'c'}

mode_nums = {'1.5_138':{'g1':'0330', 'g2':'0331', 'p1':'0332', 'p2':'0333', 
      'g4':'0334'},
    '1.5_140':{'g1':'0415', 'p1':'0416', 'g2':'0417', 'p2':'0418'},
    '1.5_150':{label:'0{}'.format(mode) for (label, mode) in ((mode_str('1.5_150',i), i) for i in range(930, 940))},
    '1.5_160':{label:'{}'.format(mode) for (label, mode) in ((str(mode_str('1.5_160', i)), i) for i in range(1530, 1539))},
    '1.24_122':{'g1':'0306', 'g2':'0307', 'p1':'0308', 'p2':'0309', 'g3':'0310', 
      'g4':'0311', 'g5':'0312', 'p3':'0313', 'g6':'0314'},
    '1.5_129':label_gen(90, 110),
    '1.33_124':label_gen(260,300),
    '1.33_122':label_gen(210,240),
    '1.24_148':{label:'{}'.format(mode) for (label, mode) in ((str(mode_str('1.24_148', i)), i) for i in range(1400, 1406))}} # 1360, 1470

data_fnames = {'mesa':'data{}/profiles/profile{}.data'.format(mass, prof),
    'model':'data{}/model_data/model_data{}.txt'.format(mass, prof)}

for mode in mode_nums[mass_prof]:
  data_fnames[mode] = 'data{0}/prof{1}/summary_l1_prof{1}_0{2}.txt'.format(mass,
      prof, mode_nums[mass_prof][mode])

'''data_fnames = {'mesa':"data/profile" + prof + ".data",
    'g1':"data/summary_l1_prof150_00500.txt",
    'g2':"data/summary_l1_prof150_00900.txt",
    'p1':"data/summary_l1_prof150_01100.txt",
    'p2':"data/summary_l1_prof150_01130.txt"}'''

#mode_keys = ['g1', 'g2', 'p1', 'p2']
mode_keys = mode_nums[mass_prof].keys()

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

numax_val = numax(data['mesa'], head_data['mesa'])
print numax_val

r_mesa = 10. ** data['mesa']['logR'] * c.rsun
r_star = r_mesa[-1]
m_mesa = data['mesa']['q'] * c.msun * head_data['mesa']['star_mass']
x = data[mode_keys[0]]['x'] # Radial grid from GYRE, normalized to [0,1]
r = x * r_star
m = np.interp(r, r_mesa, m_mesa)
omega = data[mode_keys[0]]['Omega_rot'] / np.sqrt(r_star ** 3. / (c.ggrav * m[-1]))
brunt = np.sqrt(data['mesa']['brunt_N2']) * 1e6
brunt = np.nan_to_num(brunt)
lamb = np.sqrt(2.) * data['mesa']['csound'] / (r_mesa) * 1e6
print lamb

min_l = angmomt(m, r, np.ones(len(m)) * omega[-1])
calc_l = angmomt(m, r, omega)
print 1- min_l / calc_l

g_line = mlines.Line2D([],[], color='c', label="g")
p_line = mlines.Line2D([],[], color='m', label="p")
step_line = mlines.Line2D([],[],color='b', linestyle='-', label='step')
mesa_line = mlines.Line2D([],[],color='b', linestyle='--', label='MESA')
diff_line = mlines.Line2D([],[],color='b', linestyle='-.', label='envelope diff')
conv_patch = mpatches.Patch(color='gray', alpha=0.4, label='convection', linewidth=0)
eps_patch = mpatches.Patch(color='lightblue', alpha=0.6, label='$\epsilon_\mathrm{nuc} > 10$ erg/g/cm$^3$', linewidth=0)
g_patch = mpatches.Patch(color='orange', alpha=0.2, label='g modes',linewidth=0)
p_patch = mpatches.Patch(color='lightgreen', alpha=0.6, label='p modes', linewidth=0)

# Start plotting
fig, ax1 = plt.subplots()

omega_step = omega_prof_lin(r, omega[0], omega[-1], he_core, he_core + eps)
omega_pow = omega_prof_pow(r, omega[0], omega[-1], conv_zone)
omega_tanh = omega_prof_tanh(r, omega[0], omega[-1], he_core, 1e-3 * c.rsun)

plt.xlabel("R/R$_*$")
#ax1.loglog(x, omega, 'b', linewidth=1.5)
#ax1.loglog(x, omega_step, 'b-', linewidth=1.5)
#ax1.loglog(x, omega_pow, 'b-.', linewidth=1.5)
ax1.loglog(x, omega_tanh, 'r-.', linewidth=1.5)
ax1.set_ylabel("$\Omega$", color='b')

trans = mtransforms.blended_transform_factory(ax1.transData, ax1.transAxes)

# Shade p and g mode cavities
ax1.fill_between(r_mesa / r_star, 0, 1, where=np.all(np.vstack((lamb < numax_val, brunt < numax_val)), axis=0), facecolor='lightgreen', alpha=0.6, transform=trans, linewidth=0)

ax1.fill_between(r_mesa / r_star, 0, 1, where=np.all(np.vstack((lamb > numax_val, brunt > numax_val)), axis=0), facecolor='orange', alpha=0.2, transform=trans, linewidth=0)

# Shade eps_nuc > 10 erg/g/s
ax1.fill_between(r_mesa / r_star, 0, 1, where=data['mesa']['eps_nuc'] > 10, facecolor='lightblue', alpha=0.6, transform=trans, linewidth=0)

# Shade convective regions
ax1.fill_between(r_mesa / r_star, 0, 1, where=data['mesa']['mixing_type'] == 1, facecolor='gray', alpha=0.4, transform=trans, linewidth=0)

plt.xlim([5e-5,1.])
plt.ylim([1e-7,1e-4])

plt.legend(handles=[g_line, p_line, step_line, diff_line, conv_patch, eps_patch, g_patch, p_patch], fontsize=9.5, loc=(0.015, 0.63))

plt.savefig('data{}/plots/omega_prof_{}.pdf'.format(mass, prof))

ax2 = ax1.twinx()
for (i, key) in enumerate(mode_keys):
  # Choose color based on type of mode
  #if key[0] == 'g':
  #  color = 'k'
  #elif key[0] == 'p':
  #  color = 'r'
  #style = line_styles[int(key[1:]) % len(line_styles)]
  #color = mode_colors[i]
  color = mode_colors[key[0]]
  print color
  style = '-'
  ax2.loglog(x, integral(x, data[key]['ReK'], True), c=color, linestyle=style,
      linewidth=1.5)
  #ax2.loglog(x, data[key]['ReK'], c='b', linestyle=style, linewidth=1.5)
ax2.set_ylabel('$\int K_{n, \ell}$')
ax2.yaxis.set_label_coords(1.065, 0.5)

plt.xlim([5e-5,1.])
plt.ylim([5e-3,1.1])

plt.legend(handles=[g_line, p_line, step_line,  diff_line, conv_patch, eps_patch, g_patch, p_patch], fontsize=9.5, loc=(0.015, 0.63))

plt.title('{} M$_\odot$, profile {}'.format(mass, prof))

plt.savefig('k_omega{}_{}.pdf'.format(mass, prof))
