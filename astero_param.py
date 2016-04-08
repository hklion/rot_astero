import numpy as np
import fimport
from split_lib import *
import consts as c

def large_sep(mesa_data):
  return 0.5 / integral(mesa_data['radius'] * c.rsun, 1. / mesa_data['csound'], False)

def numax(mesa_data, mesa_head_data):
  return mesa_head_data['star_mass'] * (mesa_head_data['Teff'] / 5777) ** 3.5 / (mesa_data['luminosity'][-1]) * 3090

if __name__ == "__main__":
  prof = sys.argv[1]
  mass = sys.argv[2]

  data = fimport.load_array('data{0}/profiles/profile{1}.data'.format(mass, prof),
      6)
  data = np.flipud(data)
  head_data = fimport.load_header('data{0}/profiles/profile{1}.data'.format(mass,
    prof), 2)

  print 'large sep',large_sep(data)
  print 'numax', numax(data, head_data)
