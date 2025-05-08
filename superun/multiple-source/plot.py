#!/bin/env python

#===============================================================================
# Author: Tche L., USTC, seistche@gmail.com
# Created at: Mon 17 Jun 2024 06:51:56 PM CST
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import sys

ptype = 'u'
nfile = 7
strain = True
toxyz = True

if len(sys.argv) > 1:
  ptype = sys.argv[1]

t = np.loadtxt('out-0.g' + ptype, skiprows = 1)[:, 0]
nt = len(t)

dat = np.zeros((nfile, nt))
if ptype == 'u':
  if toxyz:
    uyls = ['x', 'y', 'z']
  else:
    uyls = ['r', 't', 'z']
  vts = [ r'$ u_{%s} $' % (l) for l in uyls ]
elif ptype == 't':
  if toxyz:
    tyls = ['xx', 'yy', 'zz', 'xy', 'xz', 'yz']
  else:
    tyls = ['rr', 'tt', 'zz', 'rt', 'rz', 'tz']
  if strain:
    vts = [ r'$ \tau_{%s} $' % (l) for l in tyls ]
  else:
    vts = [ r'$ \varepsilon_{%s} $' % (l) for l in tyls ]
else:
  raise ValueError('Unrecongnized plot type <%s>' % (ptype))

nv = len(vts)
for iv in range(nv):
  # reading
  for ir in range(nfile):
    dfile = 'out-%d.g%s' % (ir, ptype)
    dat[ir, :] = np.loadtxt(dfile, skiprows = 1)[:, iv + 1]
  sha = np.amax(np.abs(dat)) / 1.0
  # plotting
  for ir in range(nfile):
    plt.plot(t, dat[ir, :] + ir * sha, label = 'src #%d' % (ir))
  plt.xlabel('time (s)')
  plt.title('waveform of component ' + vts[iv])
  plt.legend()
  plt.show()

# vim:ft=python tw=80 ts=4 sw=2 et ai
