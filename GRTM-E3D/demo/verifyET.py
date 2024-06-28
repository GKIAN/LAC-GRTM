#!/bin/env python

#===============================================================================
# Author: Tche L., USTC, seistche@gmail.com
# Created at: Fri 28 Jun 2024 10:56:20 AM CST
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as pl

irec = 1 # layer number of the receiver layer
fE = './out.gt'
fT = './out-model.gt'
mfile = 'model.dat'

if 1:
  yls = ['xx', 'yy', 'zz', 'xy', 'xz', 'yz']
else:
  yls = ['rr', 'tt', 'zz', 'rt', 'rz', 'tz']

datM = np.loadtxt(mfile, skiprows = 1)
rho, bet, alp = datM[irec - 1, 2:5]

mu = rho * bet * bet
la = rho * alp * alp - 2.0 * mu
ka = la + 2.0 * mu

t = np.loadtxt(fE, skiprows = 1)[:, 0]

datT = np.loadtxt(fT, skiprows = 1)[:, 1:]
datE = np.loadtxt(fE, skiprows = 1)[:, 1:]

yls = [ r'$ \varepsilon_{%s} $' % (l) for l in yls ]

ftxx, ftyy, ftzz, ftxy, ftxz, ftyz = datT.T
fexx, feyy, fezz, fexy, fexz, feyz = datE.T

etxx = ka * fexx + la * (feyy + fezz)
etyy = ka * feyy + la * (fexx + fezz)
etzz = ka * fezz + la * (fexx + feyy)
etxy = 2.0 * mu * fexy
etxz = 2.0 * mu * fexz
etyz = 2.0 * mu * feyz
covT = np.vstack((etxx, etyy, etzz, etxy, etxz, etyz)).T

fig, axs = pl.subplots(6, 1, sharex = True)
for i in range(6):
  axs[i].plot(t, datT[:, i], label = 'dat', linewidth = 3)
  axs[i].plot(t, covT[:, i], label = 'cov', linewidth = 1)
  axs[i].set_ylabel(r'%s' % (yls[i]))
  axs[i].ticklabel_format(axis = 'y', style = 'sci', scilimits = (0, 0))
  axs[i].legend()
pl.show()

# vim:ft=python tw=80 ts=4 sw=2 et ai
