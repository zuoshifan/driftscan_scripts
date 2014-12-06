#!/usr/bin/env python

"""Generate sky maps from a cut alm array.

:Authors: Shifan Zuo
:Date: 2014-06-16
:email: sfzuo@bao.ac.cn
:usage:
    python
"""

import argparse


def generate_map(args):
    import os
    import numpy as np
    import h5py
    import healpy
    from cora.util import hputil

    with h5py.File(args.in_file, 'r') as f:
        in_alm = f['alm'][...]
    lmax = in_alm.shape[-2] - 1
    mmax = in_alm.shape[-1] - 1
    print 'lmax = %d\nmmax = %d' % (lmax, mmax)
    lmin_cut = max(args.lmin, 0)
    lmax_cut = min(args.lmax, lmax) if args.lmax is not None else lmax
    mmin_cut = max(args.mmin, 0)
    mmax_cut = min(args.mmax, mmax) if args.mmax is not None else mmax
    if args.include:
        print 'Include lmin_cut = %d, lmax_cut = %d, mmin_cut = %d, mmax_cut = %d' % (lmin_cut, lmax_cut, mmin_cut, mmax_cut)
        temp_alm = np.zeros_like(in_alm, dtype=in_alm.dtype)
        temp_alm[:, :, lmin_cut:(lmax_cut+1), mmin_cut:(mmax_cut+1)] = in_alm[:, :, lmin_cut:(lmax_cut+1), mmin_cut:(mmax_cut+1)]
    else:
        print 'Exclude lmin_cut = %d, lmax_cut = %d, mmin_cut = %d, mmax_cut = %d' % (lmin_cut, lmax_cut, mmin_cut, mmax_cut)
        temp_alm = in_alm
        temp_alm[:, :, lmin_cut:(lmax_cut+1), mmin_cut:(mmax_cut+1)] = 0.0
    out_map = hputil.sphtrans_inv_sky(temp_alm, args.nside)
    incl = 'incl' if args.include else 'excl'
    out_file = args.out_file or ('lm_cut_%d_%d_%d_%d_%s_' % (lmin_cut, lmax_cut, mmin_cut, mmax_cut, incl) + os.path.basename(args.in_file).replace('alm', 'map'))
    with h5py.File(out_file, 'w') as f:
        f.create_dataset('map', data=out_map)


parser = argparse.ArgumentParser(description='Generate sky maps from a healpix map with large `l` and `m` only.')
parser.add_argument('in_file', type=str, nargs='?', help='Input hdf5 alm file.')
parser.add_argument('-o', '--out_file', type=str, nargs='?', help='Name of the hdf5 map file to save.')
# parser.add_argument('--maxl', type=int, nargs='?', help='Max l in the spherical transform.')
parser.add_argument('--lmin', type=int, nargs='?', default=0, help='Min cut threshold of l.')
parser.add_argument('--lmax', type=int, nargs='?', help='Max cut threshold of l.')
parser.add_argument('--mmin', type=int, nargs='?', default=0, help='Min cut threshold of m.')
parser.add_argument('--mmax', type=int, nargs='?', help='Max cut threshold of m.')
parser.add_argument('-n', '--nside', type=int, nargs='?', default=256, help='Healpix NSIDE')
parser.add_argument('-i', '--include', action='store_false', help='Include alms in range lmin<=l<=lmax, mmin<=m<=mmax if True, else exclude.')
parser.set_defaults(func=generate_map)

args = parser.parse_args()
args.func(args)