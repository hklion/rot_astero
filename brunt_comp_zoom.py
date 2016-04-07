import numpy as np
import matplotlib.pyplot as plt
import sys

i = int(sys.argv[1])

data = np.loadtxt("data/profile" + str(i) + '.data', skiprows=6)

mass = data[:,19]
h1 = data[:,60]
brunt = data[:,91]
brunt_struct = data[:,92]
brunt_comp = data[:,93]

fig, ax1 = plt.subplots()

ax1.plot(mass, brunt * 1e2, 'b')
ax1.plot(mass, brunt_struct * 1e2, 'g')
ax1.plot(mass, brunt_comp * 1e2, 'purple')

ax1.set_xlabel("Mass coordinate")
ax1.set_ylabel("$N^2$ [$10^{-2}$ Hz$^2$]")

plt.legend(("Full", "Structural", "Composition"), loc=5)

ax1.set_ylim([-1e-1, 1.5])
ax1.set_xlim([0.14, 0.16])

ax2 = ax1.twinx()
ax2.plot(mass, h1, 'r')
ax2.set_ylabel("H mass fraction", color='r')
ax2.set_xlim([0.14, 0.16])


plt.savefig("brunt_comp" + str(i) + 'zoom.pdf')
