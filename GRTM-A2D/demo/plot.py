#!/bin/env python3

"""
Author: Tche L., USTC, seistche@gmail.com
Created at: Sat 02 Jan 2021 09:24:16 PM CST
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) > 1:
  dfile = sys.argv[1]
else:
  dfile = 'out.gu'

dat = np.loadtxt(dfile, skiprows = 1)
t = dat[:, 0]
v = dat[:, 1:]

fig, axs = plt.subplots(v.shape[1], 1)
if dfile[-1] == 'u':
  yls = ['x', 'z']
  for i in range(2):
    axs[i].plot(t, v[:, i])
    axs[i].set_ylabel(yls[i] + ' component')
else:
  axs.plot(t, v[:, 0])
  axs.set_title('Pressure')

plt.show()

# vim:ft=python tw=80 ts=4 sw=2 et ai
