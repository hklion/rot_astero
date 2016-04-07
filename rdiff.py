import numpy as np
import fimport
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

data_envdiff = fimport.load_array('data/prof140_splittings_envdiff_1.0.txt', 1)
data_envdiff_101 = fimport.load_array('data/prof140_splittings_envdiff_1.01.txt', 1)
data_envdiff_1001 = fimport.load_array('data/prof140_splittings_envdiff_1.001.txt', 1)
data_envdiff_10001 = fimport.load_array('data/prof140_splittings_envdiff_1.0001.txt', 1)
data_envdiff_09999 = fimport.load_array('data/prof140_splittings_envdiff_0.9999.txt', 1)

envdiff_line = mlines.Line2D([],[],ls=':', color='b', label='envelope diff')
envdiff_101_line = mlines.Line2D([],[],ls='--',c='b',label='$+10^{-2}$')
envdiff_1001_line = mlines.Line2D([],[],ls='-.',c='b',label='$+10^{-3}$')
envdiff_10001_line = mlines.Line2D([],[],ls=('-',(10, 2, 1, 2)),c='b',label='$+10^{-4}$')
envdiff_09999_line = mlines.Line2D([],[],ls=('-',(10, 2, 10, 2)),c='b',label='$-10^{-4}$')

plt.plot(data_envdiff['ratio'], data_envdiff['r_out'], ':',c='b')
plt.plot(data_envdiff_101['ratio'], data_envdiff_101['r_out'], '--',c='b')
plt.plot(data_envdiff_1001['ratio'], data_envdiff_1001['r_out'], '-.',c='b')
plt.plot(data_envdiff_10001['ratio'], data_envdiff_10001['r_out'], ls=('-',(10, 2, 1, 2)),c='b')
plt.plot(data_envdiff_09999['ratio'], data_envdiff_09999['r_out'], ls=('-',(10, 2, 10, 2)),c='b')

plt.legend(handles=[envdiff_line, envdiff_101_line, envdiff_1001_line, envdiff_10001_line, envdiff_09999_line])
plt.xlabel("$\Omega_c/ \Omega_e$")
plt.ylabel("Rdiff [R$_\odot$]")

plt.savefig('rdiff.pdf')
