#!/usr/bin/env python

"""Create all zero value healpy map in hdf5 format.

:Authors: Shifan Zuo
:Date: 2014-05-28
:email: sfzuo@bao.ac.cn
:usage:
    python all0_map.py [-h] [-o [OUT_FILE]] [-f [NFREQ]] [-p [{1,4}]] [-n [NSIDE]]
"""

import argparse


def create_map(args):
    import numpy as np
    import h5py

    value = np.zeros((args.nfreq, args.npol, 12 * args.nside**2))
    out_file = args.out_file or 'all_zeros_%d_%d_%d.hdf5' % (args.nside, args.nfreq, args.npol)
    with h5py.File(out_file, 'w') as f:
        f.create_dataset('map', data=value)


parser = argparse.ArgumentParser(description='Create all zero value healpy map in hdf5 format.')
parser.add_argument('-o', '--out_file', type=str, nargs='?', help='Output h5py file name.')
parser.add_argument('-f', '--nfreq', type=int, nargs='?', default=10, help='Number of frequencies channels.')
parser.add_argument('-p', '--npol', type=int, nargs='?', default=4, choices=[1, 4], help='Number of polarization components.')
parser.add_argument('-n', '--nside', type=int, nargs='?', default=64, help='Healpix NSIDE.')
parser.set_defaults(func=create_map)

args = parser.parse_args()
args.func(args)
