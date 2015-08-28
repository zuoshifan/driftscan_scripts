#!/usr/bin/env python

"""Scatter plot the :math:`a_{lm}` distribute of a healpix map.

:Authors: Shifan Zuo
:Date: 2014-06-16
:email: sfzuo@bao.ac.cn
"""

import argparse


def plot_alm(args):
    import os
    import numpy as np
    import h5py
    from cora.util import hputil
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    with h5py.File(args.in_map, 'r') as f:
        in_map = f['map'][...]
    alm = hputil.sphtrans_sky(in_map, lmax=args.maxl)

    # plot alm
    plt.figure(figsize=(args.figlength, args.figwidth))
    alms = alm[args.ifreq, args.pol]
    plt.scatter(alms.reshape(-1).real, alms.reshape(-1).imag)
    out_file = args.out_file or ('alm_scatter_of_' + os.path.basename(args.in_map).replace('.hdf5', '_%d_%s.%s' % (args.ifreq, ('{%d}' % args.pol).format('T', 'Q', 'U', 'V'), args.figfmt)))
    plt.savefig(out_file)

    if args.radius:
        if args.within:
            alms = np.where(np.abs(alms)<=args.radius, alms, 0.0)
            ratio = 1.0*np.prod(np.where(np.abs(alms))[0].shape) / np.prod(alms.shape)
            print 'r <= %f is %f%%' % (args.radius, 100.0*ratio)
            out_map_file = out_file.replace('.%s' % args.figfmt, '_in_%f.hdf5' % args.radius)
        else:
            alms = np.where(np.abs(alms)>args.radius, alms, 0.0)
            ratio = 1.0*np.prod(np.where(np.abs(alms))[0].shape) / np.prod(alms.shape)
            print 'r > %f is %f%%' % (args.radius, 100.0*ratio)
            out_map_file = out_file.replace('.%s' % args.figfmt, '_out_%f.hdf5' % args.radius)
        alms = alms.reshape((1, 1) + alms.shape)
        out_map = hputil.sphtrans_inv_sky(alms, int(np.sqrt(in_map.shape[-1]/12)))

        with h5py.File(out_map_file, 'w') as f:
            f.create_dataset('map', data=out_map)
            f.attrs['non_zero_percent'] = 100.0*ratio



parser = argparse.ArgumentParser(description='Scatter plot the alm distribute of a healpix map.')
parser.add_argument('in_map', type=str, nargs='?', help='Input hdf5 sky map. If more than one, they will be first combined, i.e. added together.')
parser.add_argument('-o', '--out_file', type=str, nargs='?', help='Name of the healpix map file to save.')
parser.add_argument('-f', '--figfmt', default='png', help='Output image format.')
parser.add_argument('-i', '--ifreq', type=int, nargs='?', default=0, help='Frequency channel index.')
parser.add_argument('-p', '--pol', type=int, nargs='?', default=0, help='Polarization index.')
parser.add_argument('--maxl', type=int, nargs='?', help='Max l in the spherical transform.')
parser.add_argument('-r', '--radius', type=float, nargs='?', help='Scatter radius, i.e., the absolute value.')
parser.add_argument('--within', action='store_true', help='preserve the values whose absolute value is within `radius`.')
parser.add_argument('-l', '--figlength', type=float, default=8, help='Output figure length.')
parser.add_argument('-w', '--figwidth', type=float, default=6, help='Output figure width.')
parser.set_defaults(func=plot_alm)

args = parser.parse_args()
args.func(args)
