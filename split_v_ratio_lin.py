import numpy as np
import fimport
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import sys

prof = sys.argv[1]

data_step = fimport.load_array('data/prof' + prof + '_splittings_drop.txt', 1)
data_envdiff = fimport.load_array('data/prof' + prof + '_splittings_envdiff_1.0.txt', 1)
data_envdiff_101 = fimport.load_array('data/prof' + prof + '_splittings_envdiff_1.01.txt', 1)
data_envdiff_1001 = fimport.load_array('data/prof' + prof + '_splittings_envdiff_1.001.txt', 1)
data_envdiff_10001 = fimport.load_array('data/prof' + prof + '_splittings_envdiff_1.0001.txt', 1)
data_envdiff_09999 = fimport.load_array('data/prof' + prof + '_splittings_envdiff_0.9999.txt', 1)
data_envdiff_lin_conv = fimport.load_array('data/prof' + prof + '_splittings_envdiff_lin_conv.txt', 1)

#data_envdiff_0995 = fimport.load_array('data/prof' + prof + '_splittings_envdiff_0.999.txt', 1)

g_line = mlines.Line2D([],[], color='b', label="g")
p_line = mlines.Line2D([],[], color='r', label='p')
step_line = mlines.Line2D([],[],ls='-',color='k', label='step')
envdiff_lin_conv_line = mlines.Line2D([],[],ls=('-',(6,2,2,2,2,2)),c='k',label='fixed rdiff')
envdiff_line = mlines.Line2D([],[],ls=':', color='k', label='envelope diff')
envdiff_101_line = mlines.Line2D([],[],ls='--',c='k',label='$+10^{-2}$')
envdiff_1001_line = mlines.Line2D([],[],ls='-.',c='k',label='$+10^{-3}$')
envdiff_10001_line = mlines.Line2D([],[],ls=('-',(10, 2, 1, 2)),c='k',label='$+10^{-4}$')
envdiff_09999_line = mlines.Line2D([],[],ls=('-',(10, 2, 10, 2)),c='k',label='$-10^{-4}$')


modes = ['g1', 'p2']
colors = ['b', 'r', 'orange', 'c']

for (i, key) in enumerate(modes):
  plt.plot(data_step['ratio'], data_step[key] * 1e6, ls='-', c=colors[i])
  plt.plot(data_envdiff['ratio'], data_envdiff[key] * 1e6, ':', c=colors[i])
  plt.plot(data_envdiff_101['ratio'], data_envdiff_101[key] * 1e6, '--', c=colors[i])
  plt.plot(data_envdiff_1001['ratio'], data_envdiff_1001[key] * 1e6, '-.', c=colors[i])
  plt.plot(data_envdiff_10001['ratio'], data_envdiff_10001[key] * 1e6, ls=('-',(10, 2, 1, 2)), c=colors[i])
  plt.plot(data_envdiff_09999['ratio'], data_envdiff_09999[key] * 1e6, ls=('-',(10, 2, 10, 2)), c=colors[i])
  plt.plot(data_envdiff_lin_conv['ratio'], data_envdiff_lin_conv[key] * 1e6, ls=('-',(6,2,2,2,2,2)), c=colors[i])


plt.legend(handles=[g_line, p_line, step_line, envdiff_line, envdiff_101_line, envdiff_1001_line, envdiff_10001_line, envdiff_09999_line, envdiff_lin_conv_line],loc=2)


plt.xlabel("$\Omega_c / \Omega_e$")
plt.ylabel("Split [$\mu$Hz]")


plt.savefig('split_v_ratio.pdf')
plt.xlim([0.6,10])
plt.ylim([0.03, 0.6])
plt.savefig('split_v_ratio_zoom.pdf')
