#!/bin/env python3

"""
Author: Tche L., USTC, seistche@gmail.com
Created at: Sat 02 Jan 2021 09:24:16 PM CST
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

toxyz = True
strain = False

if toxyz:
  uyls = ['x', 'y', 'z']
  tyls = ['xx', 'yy', 'zz', 'xy', 'xz', 'yz']
else:
  uyls = ['r', 't', 'z']
  tyls = ['rr', 'tt', 'zz', 'rt', 'rz', 'tz']

if len(sys.argv) > 1:
  dfile = sys.argv[1]
else:
  dfile = 'out-model.gu'

dat = np.loadtxt(dfile, skiprows = 1)
t = dat[:, 0]
v = dat[:, 1:]

Uyls = [    r'$ u_{%s} $' % (l) for l in uyls ]
if not strain:
  Tyls = [ r'$ \tau_{%s} $' % (l) for l in tyls ]
else:
  Tyls = [ r'$ \varepsilon_{%s} $' % (l) for l in tyls ]

if dfile[-1] == 'u':
  nv = 3
  yls = Uyls
else:
  nv = 6
  yls = Tyls

fig, axs = plt.subplots(nv, 1, sharex = True)
for i in range(nv):
  axs[i].plot(t, v[:, i])
  axs[i].set_ylabel(r'%s' % (yls[i]))
  axs[i].ticklabel_format(axis = 'y', style = 'sci', scilimits = (0, 0))
plt.show()

# vim:ft=python tw=80 ts=4 sw=2 et ai
