import numpy as np
import fimport
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import sys

prof = '138_140'

data_step_138 = fimport.load_array('data/prof138_splittings_hestep.txt', 1)
data_step_140 = fimport.load_array('data/prof140_splittings_step.txt', 1)
data_envdiff_poly_138 = fimport.load_array('data/prof138_splittings_envdiff_poly_1.0.txt', 1)
data_envdiff_poly_140 = fimport.load_array('data/prof140_splittings_envdiff_poly.txt', 1)

g_line = mlines.Line2D([],[], color='b', label="g")
p_line = mlines.Line2D([],[], color='r', label='p')
step_line_138 = mlines.Line2D([],[],ls='-',color='k', label='step 138', alpha=0.5)
step_line_140 = mlines.Line2D([],[],ls='-',color='k', label='step 140')
envdiff_poly_line_138 = mlines.Line2D([],[],ls='--',c='k',label='envelope diff 138', alpha=0.5)
envdiff_poly_line_140 = mlines.Line2D([],[],ls='--',c='k',label='envelope diff 140')
measured_patch = mpatches.Patch(color='black', alpha=0.4, label='Measured split', linewidth=0)

modes = ['g1', 'p2']
colors = ['b', 'r', 'orange', 'c']

for (i, key) in enumerate(modes):
  plt.plot(data_step_138['ratio'], data_step_138[key] * 1e6, ls='-', c=colors[i], alpha=0.5)
  plt.plot(data_step_140['ratio'], data_step_140[key] * 1e6, ls='-', c=colors[i])
  plt.plot(data_envdiff_poly_138['ratio'], data_envdiff_poly_138[key] * 1e6, '--', c=colors[i], alpha=0.5)
  plt.plot(data_envdiff_poly_140['ratio'], data_envdiff_poly_140[key] * 1e6, '--', c=colors[i])

plt.legend(handles=[g_line, p_line, step_line_138, step_line_140, envdiff_poly_line_138, envdiff_poly_line_140, measured_patch],loc=2)

plt.xlabel("$\Omega_c / \Omega_e$")
plt.ylabel("Split [$\mu$Hz]")

plt.fill((0, 0, 80, 80), (0.135+0.008, 0.135-0.008, 0.135-0.008, 0.135+0.008), facecolor='r', alpha=0.4, linewidth=0)
plt.axhline(0.135, c='r', linestyle='-', linewidth=0.5, alpha=0.6)
plt.fill((0, 0, 80, 80), (0.20, 0.25, 0.25, 0.20), facecolor='b', alpha=0.4, linewidth=0)
plt.axhline(0.225, c='b', linestyle='-', linewidth=0.5, alpha=0.6)


plt.savefig('split_v_ratio_comp_{}.pdf'.format(prof))
plt.xlim([0.6,10])
plt.ylim([0.03, 0.7])
plt.savefig('split_v_ratio_comp_zoom_{}.pdf'.format(prof))

plt.clf()

exit()

plt.plot(data_step['ratio'], data_step['g1'] / data_step['p2'], ls='-', c='k')
plt.plot(data_envdiff_poly['ratio'], data_envdiff_poly['g1'] / data_envdiff_poly['p2'], '--', c='k')
plt.fill((0,0,80,80), (0.25 / (0.135 - 0.008), 0.2 / (0.135 + 0.008), 0.2 / (0.135 + 0.008), 0.25 / (0.135 - 0.008)), facecolor='k', alpha=0.4, linewidth=0)
plt.axhline(0.225 / 0.135, linewidth = 0.5, alpha = 0.6, color='k')
plt.xlabel("$\Omega_c / \Omega_e$")
plt.ylabel(r"$\delta \nu_g / \delta \nu_p$")

plt.legend(handles=[step_line, envdiff_poly_line, measured_patch],loc=4)

plt.savefig('split_ratio_comp_{}.pdf'.format(prof))

plt.xlim([0.6,10])
plt.ylim([0.6, 1.6])
plt.savefig('split_ratio_comp_zoom_{}.pdf'.format(prof))

plt.clf()

plt.plot(data_step['ratio'], (data_step['g1'] - data_step['p2']) * 1e6, ls='-', c='k')
plt.plot(data_envdiff_poly['ratio'], (data_envdiff_poly['g1'] - data_envdiff_poly['p2']) * 1e6, '--', c='k')
plt.fill((0,0,80,80), (0.25 - (0.135 - 0.008), 0.2 -(0.135 + 0.008), 0.2 - (0.135 + 0.008), 0.25 - (0.135 - 0.008)), facecolor='k', alpha=0.4, linewidth=0)
plt.axhline(0.225 - 0.135, linewidth = 0.5, alpha = 0.6, color='k')
plt.xlabel("$\Omega_c / \Omega_e$")
plt.ylabel(r"$\delta \nu_g - \delta \nu_p$")

plt.legend(handles=[step_line, envdiff_poly_line, measured_patch],loc=2)

plt.savefig('split_diff_comp_{}.pdf'.format(prof))
plt.xlim([0.6,10])
plt.ylim([-0.1, 0.3])
plt.savefig('split_diff_comp_zoom_{}.pdf'.format(prof))
