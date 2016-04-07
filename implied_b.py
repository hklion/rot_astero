import sys
import numpy as np
import matplotlib.pyplot as plt

i = int(sys.argv[1])

nu = 3e4

data = np.loadtxt("data/profile" + str(i) + ".data", skiprows=6)

q = data[:,19]
rho = 10. ** data[:,2]
omega = data[:,90]

b_r = np.sqrt(nu * 4. * np.pi * rho * omega)

plt.semilogy(q, omega)
plt.semilogy(q, rho)
plt.semilogy(q, b_r)

plt.legend(("Omega", "rho", "B_r"), loc=4)

plt.savefig("implied_b_" + str(i) + ".pdf")  
