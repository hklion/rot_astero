import numpy as np
import matplotlib.pyplot as plt
import sys
import consts as c
import fimport
from astero_param import *

usage = "python propagation.py prof_number mass"

if len(sys.argv) != 3:
  print usage
  exit()

prof = int(sys.argv[1])
mass = sys.argv[2]

delta = 0.1
c.rsun = 6.96e10

mesa = fimport.load_array("data{}/profile{}.data".format(mass, prof), 6)
mesa = np.flipud(mesa)
#l0 = np.loadtxt("summary_l0_prof" + str(prof) + ".txt", skiprows=0)
#l1 = np.loadtxt("summary_l1_prof" + str(prof) + ".txt", skiprows=0)
#l2 = np.loadtxt("summary_l2_prof" + str(prof) + ".txt", skiprows=0)
#l3 = np.loadtxt("summary_l3_prof" + str(prof) + ".txt", skiprows=0)
head = fimport.load_header("data{}/profile{}.data".format(mass, prof), 2)

r = 10. ** mesa['logR']
brunt = np.sqrt(mesa['brunt_N2']) * 1e6
csound = mesa['csound']
lamb_l1 = np.sqrt(2.) * csound / (r * c.rsun) * 1e6
lamb_l2 = np.sqrt(6.) * csound / (r * c.rsun) * 1e6

plt.loglog(r / r[-1], brunt, c='b', linewidth=1.5)
plt.loglog(r / r[-1], lamb_l1, c='k', linewidth=1.5)
#plt.loglog(r, lamb_l2, c='g', linewidth=1.5)

numax = numax(mesa, head)

print 'numax', numax

plt.axhline(y=numax, c='r',linewidth=1.5)

plt.xlim([5e-5,1e0])
plt.legend(("Brunt-Vaisala", "Lamb $\ell = 1$", r'$\nu_\mathrm{max}$'))

plt.xlabel("Radius [R$_\odot$]")
plt.ylabel("Frequency [$\mu$ Hz]")

plt.title('{} M$_\odot$, profile {}'.format(mass, prof))

plt.savefig("data{}/propagation_{}.pdf".format(mass, prof))

exit()

def prop_diagram(r, brunt, lamb, freqs):
  plt.loglog(r, brunt, c='b', linewidth = 1.5)
  plt.loglog(r, lamb, c='k', linewidth = 1.5)
  for freq in freqs:
    p = []
    g = []
    for i in range(len(r)):
      if freq > brunt[i] and freq > lamb[i]:
        p.append(r[i])
      elif freq < brunt[i] and freq < lamb[i]:
        g.append(r[i])
    #print len(g), freq
    plt.loglog(g, freq * np.ones(len(g)), c='g')
    plt.loglog(p, freq * np.ones(len(p)), c='purple')
  plt.axhline(y=4000, c='r')
  plt.legend(("Brunt-Vaisala", "Lamb"), loc=2)
  plt.xlabel("Radius [R$_\odot$]")
  plt.ylabel("Frequency [$\mu$Hz]")

  plt.savefig('rsg_ell1.pdf')

def prop_diagram_ell_zero(r, brunt, freqs):
  plt.semilogy(r, brunt)
  for freq in freqs:
    p = [1.]
    for i in range(len(r)):
      if freq > brunt[i]:
        if r[i] > p[-1] - delta:
          p.append(r[i])
        else:
          plt.semilogy(p, freq * np.ones(len(p)), c='g')
          p = [r[i]]
    plt.semilogy(p, freq * np.ones(len(p)), c='g')
  plt.axhline(y=4000, c='r')
  plt.show()
  plt.xlabel("Radius [R$_\odot$]")
  plt.ylabel("Frequency [$\mu$Hz]")

  plt.savefig('rsg_ell0.pdf')

#prop_diagram_ell_zero(r, brunt, l0[:,3])
prop_diagram(r, brunt, lamb_l1, l1[:,3])
#prop_diagram(r, brunt, lamb_l2, l2[:,4])

#exit()
