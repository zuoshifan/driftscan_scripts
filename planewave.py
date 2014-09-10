#!/usr/bin/env python

"""compute the spherical harmonic coefficients of plane wave in equatorial coordinate.

:Authors: Shifan Zuo
:Date: 2014-09-09
:email: sfzuo@bao.ac.cn
:usage:
    python planewave.py [-h] [-c {1,2}] [-a] [-f FREQ] [-n NSIDE] [-o OUT_FILE]
"""

import argparse
import numpy as np
import healpy
import h5py

from cora.util import hputil


# Read in arguments
parser = argparse.ArgumentParser(description="compute the spherical harmonic coefficients of plane wave in equatorial coordinate.")
parser.add_argument('-c', '--case', type=int, choices=[1, 2], default=1, help='Which array configuration, 1 for 32+32+32 case, 2 for 31+32+33 case.')
parser.add_argument('-a', '--auto_corr', action='store_false', help='Whether use auto correlation.')
parser.add_argument('-f', '--freq', type=float, default=750.0, help='Observing frequency.')
parser.add_argument('-n', '--nside', type=int, default=256, help='Healpix nside.')
parser.add_argument('-o', '--out_file', help='Output name file name.')
args = parser.parse_args()


# transform matrix from local coordinate (`x` points toward East, `y` points toward North, `z` points toward the Zenith) to equatorial coordinate (`X` points toward the local meridian in the `yz` plane, `Z` points toward the North pole, `Y` is perpendicular to `X` and `Z` and they together construct a right hand coordinate)
R = np.array([[0.0, -np.sqrt(2)/2.0, np.sqrt(2)/2.0],
              [1.0, 0.0, 0.0],
              [np.sqrt(2)/2.0, 0.0, np.sqrt(2)/2.0]])

# feeds positions in local coordinate
ncyl = 3
cyl_w = 15.0 # m
local_pos = []
# Case 1: 32 + 32 + 32
if args.case == 1:
    print 'Case1'
    nfeeds = 32
    feeds_spacing = 0.4 # m
    for cyl in range(ncyl):
        for feed in range(nfeeds):
            local_pos.append(np.array([cyl * cyl_w, feed * feeds_spacing, 0.0]))
# Case 2: 31 + 32 + 33
elif args.case == 2:
    print 'Case2'
    nfeeds = [31, 32, 33]
    feeds_spacing = [12.4/30, 0.4, 12.4/32]
    for cyl in range(ncyl):
        for feed in range(nfeeds[cyl]):
            local_pos.append(np.array([cyl * cyl_w, feed * feeds_spacing[cyl], 0.0]))

# feeds positions in equatorial coordinate
eq_pos = [np.dot(R, pos) for pos in local_pos]
pos_len = len(eq_pos)
# baseline vectors in equatorial coordinate
auto_corr = args.auto_corr
if auto_corr:
    bls = [(eq_pos[pos1] - eq_pos[pos2]) for pos1 in range(pos_len) for pos2 in range(pos1, pos_len)]
else:
    bls = [(eq_pos[pos1] - eq_pos[pos2]) for pos1 in range(pos_len) for pos2 in range(pos1+1, pos_len)]

# observing frequency
freq = args.freq # MHz
c = 3.0e8 # m/s
wavelen = c / (freq * 1.0e6)
uvecs = [(bl / wavelen) for bl in bls]

# healpix
nside = args.nside
nx, ny, nz = healpy.pix2vec(nside, range(12*nside**2))
enus = np.zeros_like(nx, dtype=np.complex128)
for uvec in uvecs:
    # enus += np.exp(2*np.pi*1.0J*np.abs(uvec[0] * nx + uvec[1] * ny + uvec[2] * nz))
    enus += np.exp(2*np.pi*1.0J*(uvec[0] * nx + uvec[1] * ny + uvec[2] * nz))

# spherical harmonic transform
alm = hputil.sphtrans_complex(enus, centered=True)

# save data
outfile = args.out_file or 'palm_%d_%.1f_%s_case%d.hdf5' % (args.nside, args.freq, args.auto_corr, args.case)
with h5py.File(outfile, 'w') as f:
    f.create_dataset('alm', data=alm)
