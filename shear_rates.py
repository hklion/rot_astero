import numpy as np
import matplotlib.pyplot as plt
import sys

i = int(sys.argv[1])

data = np.loadtxt("data/profile" + str(i) + ".data", skiprows=6)

coul_log = 4.
sb = 5.67e-5
gamma = 5./3.
rsun = 6.96e10

mass = data[:,19]
brunt2 = data[:,91]
brunt_therm = np.sqrt(data[:,92])
brunt_comp = np.sqrt(data[:,93])
omega = data[:,90]
r = 10. ** data[:,4]
t = 10. ** data[:,1]
p = 10. ** data[:,3]
opacity = data[:,36]
rho = 10. ** data[:,2]

# Thermal diffusivity, Spruit 2002 sec 2.4
#kappa = 16. * sb * t ** 3 / (3. * opacity * rho ** 2. * cp)

# Resistivity, Menou 2004 eq 66 (also magnetic diffusivity)
eta = 5.2e11 * coul_log / t ** 1.5

# Radiative conductivity, Menou 2004 eq 67
chi_rad = 16. * t ** 3. * sb / (3. * opacity * rho) # Intermediary for rad diff.
# Radiative diffusivity, Menou 2004 eq 68
xi_rad = (gamma - 1.) / gamma * t / p * chi_rad

domega_dr = -(omega[1:] - omega[:-1]) / (r[1:] - r[:-1])
q = -r[:-1] / omega[:-1] * (omega[1:] - omega[:-1]) / (r[1:] - r[:-1])
q_other = -r[1:] / omega[1:] * (omega[1:] - omega[:-1]) / (r[1:] - r[:-1])

# Spruit 99 eq 29
q_cond = brunt2 * eta / (2 * omega ** 2. * xi_rad)

q0 = (brunt_comp / omega) ** (7./4.) * (eta / ((r * rsun) ** 2. * brunt_comp)) ** (1./4.)
q1 = (brunt_therm / omega) ** (7./4.) * \
    (eta / ((r* rsun) ** 2. * brunt_therm)) ** (1./4.) * (eta / xi_rad) ** (3./4.)

plt.plot(r, omega)
plt.plot(r[:-1], domega_dr)
#plt.plot(r[:-1], q)
plt.show()

plt.plot(mass[:-1], q)
#plt.plot(mass[1:], q_other)
plt.plot(mass, q_cond)
plt.plot(mass, q0)
plt.plot(mass, q1)
plt.plot(mass, q0 + q1)
#plt.plot(mass, brunt_therm)
#plt.plot(mass, brunt_comp)

plt.legend(("q", "q_cond", "q0", "q1", "q0 + q1"),loc=1)

plt.show()
