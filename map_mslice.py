#!/usr/bin/env python

"""Generate sky maps from a healpix map with a single `m` slice.

:Authors: Shifan Zuo
:Date: 2014-06-26
:email: sfzuo@bao.ac.cn
:usage:
    python map_mslice.py [-h] [-o [OUT_FILE]] [-m [M_INDEX]] [in_map]
"""

import argparse


def generate_map(args):
    import numpy as np
    import h5py
    import healpy
    from cora.util import hputil

    with h5py.File(args.in_map, 'r') as f:
        in_map = f['map'][...]
    nside = healpy.pixelfunc.get_nside(in_map[0])
    in_alm = hputil.sphtrans_sky(in_map)
    # lmax = in_alm.shape[-2] - 1
    mmax = in_alm.shape[-1] - 1
    print 'mmax = %d' % mmax
    # l_cut = min(args.l_cut, lmax)
    # m_cut = min(args.m_cut, mmax)
    m_index = min(args.m_index, mmax)
    temp_alm = np.zeros_like(in_alm, dtype=in_alm.dtype)
    temp_alm[:, :, :, m_index] = in_alm[:, :, :, m_index]
    out_map = hputil.sphtrans_inv_sky(temp_alm, nside)
    out_file = args.out_file or ('m_slice_%d_' % m_index + args.in_map)
    with h5py.File(out_file, 'w') as f:
        f.create_dataset('map', data=out_map)


parser = argparse.ArgumentParser(description='Generate sky maps from a healpix map with a single `m` slice.')
parser.add_argument('in_map', type=str, nargs='?', help='Input hdf5 sky map. If more than one, they will be first combined, i.e. added together.')
parser.add_argument('-o', '--out_file', type=str, nargs='?', help='Name of the hdf5 map file to save.')
# parser.add_argument('-l', '--l_cut', type=int, nargs='?', default=0, help='zeroing alm for l < l_cut.')
parser.add_argument('-m', '--m_index', type=int, nargs='?', default=0, help='m index.')
parser.set_defaults(func=generate_map)

args = parser.parse_args()
args.func(args)