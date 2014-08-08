#!/usr/bin/env python

"""Add healpy maps to create a total map in h5py format.

:Authors: Shifan Zuo
:Date: 2014-05-29
:email: sfzuo@bao.ac.cn
:usage:
    python add_maps.py [-h] [-o [OUT_FILE]] hpmaps [hpmaps ...]
"""

import argparse


def create_map(args):
    import numpy as np
    import h5py

    map_sum = None
    for hpmap in args.hpmaps:
        with h5py.File(hpmap, 'r') as f:
            if map_sum is None:
                map_sum = f['map'][...]
            else:
                map_sum += f['map'][...]
    out_file = args.out_file or 'total_map.hdf5'
    with h5py.File(out_file, 'w') as f:
        f.create_dataset('map', data=map_sum)


parser = argparse.ArgumentParser(description='Add healpy maps to create a total map in h5py format.')
parser.add_argument('hpmaps', type=str, nargs='+', help='Input healpy maps in hdf5 format.')
parser.add_argument('-o', '--out_file', type=str, nargs='?', help='Output h5py file name.')
parser.set_defaults(func=create_map)

args = parser.parse_args()
args.func(args)

