import numpy as np
import fimport
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import sys

org_angmom = 6.13487288796e+49

data_step = fimport.load_array('data/prof140_splittings_drop.txt', 1)
data_envdiff_lin_conv = fimport.load_array('data/prof140_splittings_envdiff_lin_conv.txt', 1)

plt.semilogy(data_step['ratio'], abs(data_step['angmom'] / org_angmom - 1))
plt.semilogy(data_envdiff_lin_conv['ratio'], abs(data_envdiff_lin_conv['angmom'] / org_angmom - 1))

plt.xlabel("$\Omega_c / \Omega_e$")
plt.ylabel("Absolute val of fractional difference in angular momentum")

plt.legend(("Step", "Diff in envelope"), loc=5)

plt.savefig('angular_momentum_models.pdf')
