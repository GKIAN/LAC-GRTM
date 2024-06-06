#!/bin/env python3

"""
Author: Tche L., USTC, seistche@gmail.com
Created at: Sat 02 Jan 2021 09:24:16 PM CST
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

if 1:
  uyls = ['x', 'y', 'z']
  tyls = ['xz', 'yz', 'zz']
else:
  uyls = ['r', 't', 'z']
  tyls = ['rz', 'tz', 'zz']

if len(sys.argv) > 1:
  dfile = sys.argv[1]
else:
  dfile = 'out-model.gu'

dat = np.loadtxt(dfile, skiprows = 1)
t = dat[:, 0]
v = dat[:, 1:]

fig, axs = plt.subplots(3, 1, sharex = True)
if dfile[-1] == 'u':
  for i in range(3):
    axs[i].plot(t, v[:, i])
    axs[i].set_ylabel(r'$ u_{%s} $' % (uyls[i]))
else:
  for i in range(3):
    axs[i].plot(t, v[:, i])
    axs[i].set_ylabel(r'$ \tau_{%s} $' % (tyls[i]))

plt.show()

# vim:ft=python tw=80 ts=4 sw=2 et ai
