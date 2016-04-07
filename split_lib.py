import numpy as np
import sys
import fimport
import consts as c
import scipy.optimize

def pow_turning_point(r, center, surface, alpha):
  return (surface / center) ** (1. / alpha) * r[-1]

def pow_alpha(r, center, surface, turning_point1, turning_point2):
  return np.log(center / surface) / np.log(turning_point2 / turning_point1)

def integral(x, integrand, array, bounds=False):
  assert bounds==False or len(bounds) == 2
  dx = x[1:] - x[:-1]
  # Zero out unwanted parts of dx to bound integral
  if bounds:
    i_bounds = np.searchsorted(x, bounds)
    dx[:i_bounds[0]] = 0.
    dx[i_bounds[1]:] = 0.
  if not array:
    return sum(dx * integrand[1:])
  if array:
    arr = np.zeros(len(x))
    for i in range(len(dx)):
      arr[i + 1] = arr[i] + dx[i] * integrand[i]
    return arr

# Generate angular velocity profile that is uniform in the center, then falls off
# linearly to the envelope value, and is uniform until the surface
def omega_prof_lin(r, center, surface, turning_pt1, turning_pt2):
  # Find indices of turning points.
  i_turning_pt1 = np.searchsorted(r, turning_pt1)
  i_turning_pt2 = np.searchsorted(r, turning_pt2)
  # Calculate slope and intercept of linear portion
  slope = (center - surface) / (turning_pt1 - turning_pt2)
  inter = (center * turning_pt2 - surface * turning_pt1) / (turning_pt2 - turning_pt1)
  omega = np.ones(len(r))
  omega[:i_turning_pt1] *= center
  omega[i_turning_pt2:] *= surface
  omega[i_turning_pt1:i_turning_pt2] = slope * r[i_turning_pt1:i_turning_pt2] + inter
  return omega

def omega_prof_tanh(r, center, surface, turning_pt, sigma):
  return (center - surface) * 0.5 * (np.tanh((turning_pt - r) / sigma)) + (surface + center) / 2

def omega_prof_pow(r, center, surface, turning_pt1, turning_pt2 = 0.):
  if turning_pt2 == 0:
    turning_pt2 = r[-1]
  i_turning_pt1 = np.searchsorted(r, turning_pt1)
  i_turning_pt2 = np.searchsorted(r, turning_pt2)
  alpha = pow_alpha(r, center, surface, turning_pt1, turning_pt2)
  omega = np.ones(len(r))
  omega[:i_turning_pt1] *= center
  omega[i_turning_pt1:i_turning_pt2] = surface * turning_pt2 ** alpha / \
      r[i_turning_pt1:i_turning_pt2] ** alpha
  omega[i_turning_pt2:] *= surface
  return omega

# Calculate splitting for a single mode/rotational profile
def split(x, k, beta, omega):
  return integral(x, k * omega, False) * beta / (2. * np.pi)

def partial_split(x, k, beta, omega, lim):
  return integral(x, k * omega, lim) * beta / (2. * np.pi)

def angmomt(m, r, omega_profile):
  i = 0.4 * (m[1:] - m[:-1]) * (r[1:] ** 5 - r[:-1] ** 5) / (r[1:] ** 3 - r[:-1] ** 3)
  return np.sum(i * omega_profile[:-1])

# Generate data string for a single rotational configuration for a number of modes
def outdata(omega_name, x, m, r, k_all, beta, arr, org_omega=[], tp=()):
  pow_profile = False
  if arr:
    omega = org_omega
  else:
    if len(tp) == 4:
      omega = omega_prof_lin(r, tp[0], tp[1], tp[2], tp[3])
    else:
      omega = omega_prof_pow(r, tp[0], tp[1], tp[2])
      alpha = pow_alpha(r, tp[0], tp[1], tp[2])
      pow_profile = True
  l = angmomt(m, r, omega)
  splitting = []
  for (i, k_mode) in enumerate(k_all):
    splitting.append(split(x, k_mode, beta[i], omega))
  if pow_profile:
    splitting.extend([tp[2] / c.rsun, alpha])
  else:
    splitting.extend([tp[2] / c.rsun, tp[3] / c.rsun]) # internal and external turning points
  str_splitting = [str(val) for val in splitting]
  outstring = omega_name.ljust(8) + ' ' + str(l) + ' ' + ' '.join(str_splitting) + '\n'
  return outstring

# Estimate the central angular velocity assuming a g mode splitting with very
# little envelope contribution (-> beta = 1/2). g_max is in muHz
def omega_c_est(g_max):
  return g_max * (4 * np.pi)

def label_gen(first, last):
  return {label:str(mode).zfill(4) for (label, mode) in (('g_{}'.format(str(i).zfill(4)), i) for i in range(first, last))}

def strip_txt(mass):
  for (i, char) in enumerate(mass):
    if char == '_':
      return float(mass[:i])
  return float(mass)

'''def gen_rot_profile(model_data, prof_type, min_ratio, max_ratio, eps=1e-6, 
    param=1, n=20, **kwargs):
  omegas = {}
  omega_keys = []
  if prof_type == ('hestep') or prof_type == 'step':
    if prof_type == 'hestep':
      step_at = he_core
    else:
      step_at = param * r_star
    for ratio in np.logspace(0., np.log10(omega[0] / omega[-1]), 20):
      key = str(np.around(ratio, 5))
      omegas[key] = (ratio * omega[-1], omega[-1], step_at, step_at + eps)
      omega_keys.append(key)
  elif prof_type == 'lin_conv':
    for ratio in np.logspace(0., np.log10(omega[0] / omega[-1]), 20):
      key = str(np.around(ratio, 5))
      omegas[key] = (ratio * omega[-1], omega[-1], conv_zone, r_star - eps)
      omega_keys.append(key)
  elif prof_type == 'lin_env_rdiff':
    angmomt_ratio = param
    for rdiff in np.logspace(np.log10(conv_zone), np.log10(r_star - eps), 20):
      def angmomt_rdiff(omc):
        return angmomt(m, r, 
            omega_prof_lin(r, omc, omega[-1], conv_zone, rdiff)) \
          - angmomt(m, r, omega) * angmomt_ratio
      best_omc = scipy.optimize.broyden1(angmomt_rdiff, omega[-1] * 1.00001)
      print rdiff / r_star, best_omc / omega[-1], angmomt_rdiff(best_omc), \
          angmomt(m, r, omega_prof_lin(r, best_omc, omega[-1], he_core, rdiff))
      key = str(np.around(best_omc / omega[-1], 5))
      omegas[key] = (best_omc, omega[-1], conv_zone, rdiff)
      omega_keys.append(key)
  elif prof_type == 'poly':
    for ratio in np.logspace(0., np.log10(omega[0] / omega[-1]), 20):
      key = str(np.around(ratio, 5))
      omegas[key] = (ratio * omega[-1], omega[-1], conv_zone * param)
      omega_keys.append(key)
  else:
    print "Profile type unknown. Must be {hestep, step, lin_conv, rdiff, poly}"
    exit(1)
  return omega_keys, omegas'''
