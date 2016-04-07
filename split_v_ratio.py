import numpy as np
import fimport
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import sys

prof = sys.argv[1]
mass = sys.argv[2]
mass_prof = '{}_{}'.format(mass, prof)

file_tail = ''
if len(sys.argv) > 3:
  file_tail = sys.argv[3]

model = fimport.load_header('data{}/mass_data.txt'.format(mass), 1)
data_step = fimport.load_array('data{}/prof{}_splittings_hestep.txt'.format(mass, prof), 1)
#data_hburn_drop = fimport.load_array('data/prof' + prof + '_splittings_drop.txt', 1)
#data_hburn_step = fimport.load_array('data/prof' + prof + '_splittings_step_hburn.txt', 1)
data_envdiff_poly = fimport.load_array('data{}/prof{}_splittings_envdiff_poly_1.0.txt'.format(mass, prof), 1)
data_step_2em2 = fimport.load_array('data{}/prof{}_splittings_step_0.02.txt'.format(mass, prof), 1)
#data_envdiff_75_poly = fimport.load_array('data/prof' + prof + '_splittings_envdiff_0.75poly.txt',1)
#data_envdiff_50_poly = fimport.load_array('data/prof' + prof + '_splittings_envdiff_0.5poly.txt',1)
#data_envdiff_25_poly = fimport.load_array('data/prof' + prof + '_splittings_envdiff_g+dgpoly.txt',1)

g_line = mlines.Line2D([],[], color='b', label="g")
p_line = mlines.Line2D([],[], color='r', label='p')
p2_line = mlines.Line2D([],[], color='orange', label='less p')
step_line = mlines.Line2D([],[],ls='-',color='k', label='step')
envdiff_poly_line = mlines.Line2D([],[],ls='--',c='k',label='envelope diff')
step_2em2_line = mlines.Line2D([],[],ls='-.',c='k', label='step 0.02 R$_*$')
measured_patch = mpatches.Patch(color='black', alpha=0.4, label='Measured split', linewidth=0)

#modes = ['g1', 'p1', 'p2']
modes = ['g5', 'p3', 'g6', 'p1']
colors = ['b', 'r', 'orange', 'c','g']

p = model['p']
dp = model['dp']
g = model['g']
dg = model['dg']

for (i, key) in enumerate(modes):
  plt.plot(data_step['ratio'], data_step[key] * 1e6, ls='-', c=colors[i])
#  plt.plot(data_hburn_drop['ratio'], data_hburn_drop[key] * 1e6, ls=':', c=colors[i])
#  plt.plot(data_hburn_step['ratio'], data_hburn_step[key] * 1e6, ls='-.', c=colors[i])
  plt.plot(data_envdiff_poly['ratio'], data_envdiff_poly[key] * 1e6, '--', c=colors[i])
  plt.plot(data_step_2em2['ratio'], data_step_2em2[key] * 1e6, '-.', c=colors[i])
#  plt.plot(data_envdiff_75_poly['ratio'], data_envdiff_75_poly[key] * 1e6, '-.', c=colors[i])
#  plt.plot(data_envdiff_50_poly['ratio'], data_envdiff_50_poly[key] * 1e6, ':', c=colors[i])
#  plt.plot(data_envdiff_25_poly['ratio'], data_envdiff_25_poly[key] * 1e6, '-.', c=colors[i])

plt.legend(handles=[g_line, p_line, p2_line, step_line, envdiff_poly_line, step_2em2_line, measured_patch],loc=2)

plt.xlabel("$\Omega_c / \Omega_e$")
plt.ylabel("Split [$\mu$Hz]")

plt.fill((0, 0, 50, 50), (p+dp, p-dp, p-dp, p+dp), facecolor='r', alpha=0.4, linewidth=0)
plt.axhline(p, c='r', linestyle='-', linewidth=0.5, alpha=0.6)
plt.fill((0, 0, 50, 50), (g-dg, g+dg, g+dg, g-dg), facecolor='b', alpha=0.4, linewidth=0)
plt.axhline(g, c='b', linestyle='-', linewidth=0.5, alpha=0.6)


plt.savefig('split_v_ratio_step_{}_{}.pdf'.format(mass_prof, file_tail))
plt.xlim([0.6,10])
plt.ylim([0.03, 0.6])
plt.legend(handles=[g_line, p_line, p2_line, step_line, envdiff_poly_line, step_2em2_line, measured_patch],loc=4)
plt.savefig('split_v_ratio_step_zoom_{}_{}.pdf'.format(mass_prof, file_tail))

plt.clf()

plt.plot(data_step['ratio'], data_step['g1'] / data_step['p2'], ls='-', c='k')
plt.plot(data_envdiff_poly['ratio'], data_envdiff_poly['g1'] / data_envdiff_poly['p2'], '--', c='k')
plt.fill((0,0,80,80), (g+dg / (p - dp), g-dg / (p + dp), g-dg / (p + dp), g+dg / (p - dp)), facecolor='k', alpha=0.4, linewidth=0)
plt.axhline(g / p, linewidth = 0.5, alpha = 0.6, color='k')
plt.xlabel("$\Omega_c / \Omega_e$")
plt.ylabel(r"$\delta \nu_g / \delta \nu_p$")

plt.legend(handles=[step_line, envdiff_poly_line, measured_patch],loc=4)

plt.savefig('split_ratio_step_{}_{}.pdf'.format(mass_prof, file_tail))

plt.xlim([0.6,10])
plt.ylim([0.6, 1.6])
plt.savefig('split_ratio_step_zoom_{}_{}.pdf'.format(mass_prof, file_tail))

plt.clf()

plt.plot(data_step['ratio'], (data_step['g1'] - data_step['p2']) * 1e6, ls='-', c='k')
plt.plot(data_envdiff_poly['ratio'], (data_envdiff_poly['g1'] - data_envdiff_poly['p2']) * 1e6, '--', c='k')
plt.fill((0,0,80,80), (g+dg - (p - dp), g-dg -(p + dp), g-dg - (p + dp), g+dg - (p - dp)), facecolor='k', alpha=0.4, linewidth=0)
plt.axhline(g - p, linewidth = 0.5, alpha = 0.6, color='k')
plt.xlabel("$\Omega_c / \Omega_e$")
plt.ylabel(r"$\delta \nu_g - \delta \nu_p$")

plt.legend(handles=[step_line, envdiff_poly_line, measured_patch],loc=2)

plt.savefig('split_diff_step_{}_{}.pdf'.format(mass_prof, file_tail))
plt.xlim([0.6,10])
plt.ylim([-0.1, 0.3])
plt.savefig('split_diff_step_zoom_{}_{}.pdf'.format(mass_prof, file_tail))
