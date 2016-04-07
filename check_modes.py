import numpy as np
import fimport
import sys

usage = 'check_modes.py prof_number mass'

if len(sys.argv) != 3:
  print usage
  exit()

prof = sys.argv[1]
mass = sys.argv[2]

modes = fimport.load_array("data{}/summary_l1_prof{}.txt".format(mass, prof), 6)

winding = modes['n_p'] - modes['n_g']

err = False

for i in range(1, len(winding)):
  if winding[i] != winding[i-1] + 1:
    print modes[i]
    err = True

if err:
  print "Modes are missing."
else:
  print "All good."
