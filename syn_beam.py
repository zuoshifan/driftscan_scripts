#!/usr/bin/env python

"""Synthesized beam computation.

:Authors: Shifan Zuo
:Date: 2014-12-07
:email: sfzuo@bao.ac.cn
"""


import argparse
import numpy as np
import healpy
import h5py

from drift.telescope import cylinder
from drift.telescope import exotic_cylinder
from drift.core import visibility
from cora.util import hputil


# Read in arguments
parser = argparse.ArgumentParser(description="compute the synthesized beam of an cylinder array and its spherical harmonic coefficients.")
parser.add_argument('-c', '--case', type=int, choices=[1, 2, 3, 4], default=1, help='Which array configuration, 1 for 32+32+32 case, 2 for 31+32+33 case.')
parser.add_argument('--lat', type=float, nargs='?', default=45, help='Telescope latitude.')
parser.add_argument('--lon', type=float, nargs='?', default=90, help='Telescope longitude.')
parser.add_argument('-a', '--auto_corr', action='store_false', help='Whether use auto correlation.')
parser.add_argument('-m', '--mask_horizon', action='store_false', help='Whether mask out below horizon.')
parser.add_argument('-p', '--primary_beam', action='store_true', help='Multiply primary beam if True.')
parser.add_argument('-f', '--freq', type=float, default=750.0, help='Observing frequency.')
parser.add_argument('-n', '--nside', type=int, default=256, help='Healpix nside.')
parser.add_argument('-o', '--outfile', help='Output name file name.')
args = parser.parse_args()


# Case 1: 32 + 32 + 32
if args.case == 1:
    cyl = cylinder.PolarisedCylinderTelescope(args.lat, args.lon)
elif args.case == 2:
    cyl = exotic_cylinder.UnequalFeedsCylinder(args.lat, args.lon)
elif args.case == 3:
    cyl = exotic_cylinder.ArbitraryPolarisedCylinder(args.lat, args.lon)
elif args.case == 4:
    cyl = exotic_cylinder.UnequalFeedsCylinder(args.lat, args.lon)
else:
    raise Exception('Unsupported case: %d' % args.case)

# set the location of the telescope
# cyl.zenith = [45.0, 0.0]

# Set the measured frequencies of the telescope
cyl.num_freq = 3
cyl.freq_lower = args.freq - 10.0
cyl.freq_upper = args.freq + 10.0


# Set the properties of the cylinders
cyl.num_cylinders = 3
# cyl.num_cylinders = 1
cyl.cylinder_width = 15.0
if args.case == 1:
    cyl.num_feeds = 32
    cyl.feed_spacing = 0.4
    # cyl.feed_spacing = 1.2
elif args.case == 2:
    cyl.num_feeds = [31, 32, 33]
    cyl.feed_spacing = [12.4/30, 0.4, 12.4/32]
elif args.case == 3:
    cyl.num_feeds = 32
    D1 = 0.8
    D2 = 1.0
    D3 = 1.1
    cyl1_sp = [D1] * 10 + [D3] * 5 + [D2] + [D3] * 15 # spacing between adjacent feeds
    cyl2_sp = [D3] * 10 + [D1] * 5 + [D2] + [D1] * 5 + [D3] * 10
    cyl3_sp = [D3] * 15 + [D2] + [D3] * 5 + [D1] * 10
    cyl.feed_spacing = [np.cumsum(np.insert(cyl1_sp, 0, 0)).tolist(),
                        np.cumsum(np.insert(cyl2_sp, 0, 0)).tolist(),
                        np.cumsum(np.insert(cyl3_sp, 0, 0)).tolist()]
elif args.case == 4:
    cyl.num_feeds = [31, 32, 33]
    cyl.feed_spacing = [31.0/30, 1.0, 31.0/32]

# Set the thermal noise (T_sys flat across spectrum)
cyl.tsys_flat = 50.0
# cyl.ndays = 733
# cyl.l_boost = 1.0

cyl.auto_correlations = args.auto_corr

#cyl.in_cylinder = False


nside = args.nside
cyl._init_trans(nside)
angpos = cyl._angpos
horizon = cyl._horizon
wavelength = cyl.wavelengths[1] # central wavelength, corresponding to args.freq

Nbl = len(cyl.baselines) # number of baselines
enus = 0.0
for (n, bi) in zip(cyl.redundancy, range(Nbl)):
    if args.primary_beam:
        enus += n * np.array(cyl._beam_map_single(bi, (cyl.num_freq - 1)/2))
    else:
        uv = cyl.baselines[bi] / wavelength
        fringe = visibility.fringe(angpos, cyl.zenith, uv)
        if args.mask_horizon:
            enus += n * fringe * horizon
        else:
            enus += n * fringe

if args.outfile is not None:
    outenus = 'enus_' + args.outfile
    outpalms = 'palms_' + args.out_file
else:
    if args.primary_beam:
        outenus = 'enus_%d_%.1f_%.1f_%.1f_%s_case%d_%s.hdf5' % (args.nside, args.freq, args.lat, args.lon, args.auto_corr, args.case, 'p')
        outpalms = 'palms_%d_%.1f_%.1f_%.1f_%s_case%d_%s.hdf5' % (args.nside, args.freq, args.lat, args.lon, args.auto_corr, args.case, 'p')
    else:
        outenus = 'enus_%d_%.1f_%.1f_%.1f_%s_case%d_%s.hdf5' % (args.nside, args.freq, args.lat, args.lon, args.auto_corr, args.case, args.mask_horizon)
        outpalms = 'palms_%d_%.1f_%.1f_%.1f_%s_case%d_%s.hdf5' % (args.nside, args.freq, args.lat, args.lon, args.auto_corr, args.case, args.mask_horizon)

# save enus (synthesized beam?)
with h5py.File(outenus, 'w') as f:
    f.create_dataset('map', data=enus)

# spherical harmonic transform
if args.primary_beam:
    alm = hputil.sphtrans_complex_pol([enus[0], enus[1], enus[2], enus[3]], centered=True)
else:
    alm = hputil.sphtrans_complex(enus, centered=True)

# save data
with h5py.File(outpalms, 'w') as f:
    f.create_dataset('alm', data=np.array(alm))
