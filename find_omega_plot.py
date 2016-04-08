import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn.apionly as sns
import fimport
import consts as c

prof = sys.argv[1]
mass = sys.argv[2]
prof_type = sys.argv[3]
p = sys.argv[4]
g = sys.argv[5]

mass_prof = '{}_{}'.format(mass, prof)

lin = np.linspace(0, 1e-6, 100)

mode_ratio = {'1.24_122':1.7, '1.5_138':1.5, '1.5_140':1.35}

data = fimport.load_array('data{}/find_omega_{}_{}_{}_{}.txt'.format(mass, prof, prof_type, p, g), 1)
data_mesa = fimport.load_array('data{}/profiles/profile{}.data'.format(mass, prof), 6)
radius = data_mesa['radius'][0] * c.rsun

g_vals = np.unique(data[g])
p_vals = np.unique(data[p])

i = 0
while data[g][i] == data[g][0]:
  i += 1
j = len(data[g]) / i
assert int(j) == j

ratio = np.reshape(np.copy(data['ratio']), (j, i))
omega_e = np.reshape(np.copy(data['omega_e']), (j, i))

plt.pcolor(p_vals, g_vals, ratio, cmap=plt.get_cmap('cubehelix'), vmax=50, vmin=1)
plt.axis([p_vals.min(), p_vals.max(), g_vals.min(), g_vals.max()])
plt.fill_between(lin, 0, lin, color='white', alpha=0.4,linewidth=0)
plt.plot(lin, lin * mode_ratio[mass_prof], color='k')
plt.colorbar(label='$\Omega_c / \Omega_e$')
plt.xlabel(p)
plt.ylabel(g)
plt.title(prof + ' ' + prof_type)

plt.savefig('data{}/find_omega_ratio_{}_{}_{}_{}.pdf'.format(mass, prof, prof_type, p, g))

plt.clf()

plt.pcolor(p_vals, g_vals, omega_e * radius / 1e5, cmap=plt.get_cmap('cubehelix'), vmax=1, vmin=0)
plt.axis([p_vals.min(), p_vals.max(), g_vals.min(), g_vals.max()])
plt.fill_between(lin, 0, lin, color='white', alpha=0.4, linewidth=0)
plt.plot(lin, lin * mode_ratio[mass_prof], color='white')
#plt.plot(lin, lin * 1.4, color='white')
plt.colorbar(label='$v_e$ [km/s]')
plt.xlabel(p)
plt.ylabel(g)
plt.title(prof + ' ' + prof_type)

plt.savefig('data{}/plots/find_omega_ve_{}_{}_{}_{}.pdf'.format(mass, prof, prof_type, p, g))
